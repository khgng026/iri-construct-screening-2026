#!/usr/bin/env python3
"""Separate label-free span prediction from model inference and held-out scoring."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.session_store import UnitStore,atomic_json,canonical,digest,file_hash
from iri_constructs.session_benchmark import load_suite,cosine_ranking,evaluate_query,summarize
from iri_constructs.session_segmentation import (fixed_windows,alignment_guided_segments,SegmentationConfig,
    fuse_rankings,per_segment_union,interval_metrics,select_configuration)
from session_experiment import check_deadline,public_query,mean_residues,model_identity,peak_mib


def open_store(path,session):
    """Open one complete immutable source snapshot, never latest loose units."""
    path=Path(path).resolve()
    store=UnitStore(path,json.loads((path/'identity.json').read_text()),session/'job_ledger.jsonl')
    completions=[(json.loads(p.read_text())['completed_utc'],p) for p in path.glob('completion-*.json')]
    if not completions:raise ValueError('Source job is not complete')
    descriptor_path=max(completions)[1];descriptor=json.loads(descriptor_path.read_text())
    expected=expected_keys(store.identity)
    if (descriptor['identity_hash']!=store.identity_hash or set(descriptor['units'])!=expected
        or descriptor['expected_count']!=len(expected)):raise ValueError('Source completion coverage/identity mismatch')
    store.completed_paths={}
    for key,entry in descriptor['units'].items():
        unit=(path/entry['path']).resolve()
        if not unit.is_relative_to(path):raise ValueError('Source completion path escapes job')
        if file_hash(unit/'complete.json')!=entry['descriptor_sha256']:raise ValueError('Source unit descriptor changed')
        checked=store.verify(unit)
        if checked['key']!=key:raise ValueError('Source completed unit key mismatch')
        store.completed_paths[key]=unit
    store.completion_sha256=file_hash(descriptor_path)
    return store


def expected_keys(identity):
    ids=identity['query_ids']
    if not ids or len(ids)!=len(set(ids)):raise ValueError('Expected query list is empty or duplicated')
    mode=identity['mode']
    if mode in ('model','alignment'):
        expected={'query-'+qid for qid in ids}|{'summary'}
        if mode=='model':expected|={'catalogue-'+a for a in identity['catalogue_accessions']}
    elif mode=='automatic_spans_v1':
        if not identity['configs']:raise ValueError('No prediction configurations')
        expected={qid+'--'+key for qid in ids for key in identity['configs']}
    elif mode=='automatic_segment_retrieval_v1':
        configs=[c['configuration_id'] for c in identity['configs']]
        if not configs or len(configs)!=len(set(configs)):raise ValueError('Empty or duplicate retrieval configurations')
        expected={qid+'--'+key for qid in ids for key in configs}|{'summary'}
    else:raise ValueError('Unsupported completed source job mode')
    return expected


def completed_read(store,key,arrays=False):
    path=store.completed_paths.get(key)
    if path is None:return None
    record=json.loads((path/'record.json').read_text())
    if not arrays:return record
    import numpy as np
    descriptor=json.loads((path/'complete.json').read_text())
    if 'arrays.npz' not in descriptor['files']:raise ValueError('Source arrays lack completion checksum')
    with np.load(path/'arrays.npz',allow_pickle=False) as z:data={k:z[k].copy() for k in z.files}
    return record,data


def validate_source_suite(identity,suite,challenge_hash):
    if identity.get('challenge_sha256')!=challenge_hash:raise ValueError('Source challenge mismatch')
    if identity.get('catalogue_accessions')!=sorted(suite['sequences']):raise ValueError('Source catalogue mismatch')
    expected={a:digest(seq.encode()) for a,seq in suite['sequences'].items()}
    if identity.get('catalogue_sequence_sha256')!=expected:raise ValueError('Source reference sequence hashes differ')


def validate_dev_representation(method,design):
    expected=design['selection_representation']
    for key in ('model_id','revision','device','dtype'):
        if method.get(key)!=expected[key]:raise ValueError('Development model differs from frozen selection representation: '+key)
    if (expected['layer']!='last' or type(method.get('layer')) is not int or method['layer']<1
        or method.get('layer')!=method.get('num_layers')):
        raise ValueError('Development selection requires the frozen final layer')


def setup(args):
    s=args.session.resolve();suite=load_suite(ROOT,s/'inputs');design=json.loads(args.design.read_text())
    if file_hash(s/'inputs/challenge_manifest.json')!=design['source_challenge_sha256']:raise ValueError('Design input mismatch')
    for relative,expected in design['implementation_sha256'].items():
        if file_hash(ROOT/relative)!=expected:raise ValueError('Frozen segmentation implementation changed: '+relative)
    ids=design['selection_query_ids'] if args.partition=='dev' else design['test_query_ids']+design['clean_sanity_query_ids']
    selected=None
    if args.partition=='test':
        if not args.selection:raise ValueError('Held-out scoring requires frozen dev selection')
        selected=json.loads(args.selection.read_text())
        if selected['design_sha256']!=file_hash(args.design):raise ValueError('Frozen selection/design mismatch')
        development=open_store(Path(selected['development_job']),s)
        if (development.identity['mode']!='automatic_segment_retrieval_v1'
            or development.identity['partition']!='dev'
            or development.identity['design_sha256']!=file_hash(args.design)
            or development.identity_hash!=selected['development_identity']
            or development.identity['query_ids']!=design['selection_query_ids']
            or development.identity['configs']!=design['candidate_configurations']):
            raise ValueError('Frozen selection lacks matching complete development evidence')
        validate_dev_representation(development.identity['model'],design)
        if completed_read(development,'summary')['selection']!=selected['selection']:
            raise ValueError('Selection differs from completed development summary')
        if file_hash(development.completed_paths['summary']/'complete.json')!=selected['development_summary_descriptor_sha256']:
            raise ValueError('Development selection summary descriptor changed')
    elif args.selection:
        raise ValueError('Development scoring cannot consume a pre-existing selection')
    configs=design['candidate_configurations']
    if selected:
        winner_map=selected['selection']['selected_by_segmenter']
        segmenters={c['segmenter'] for c in configs}
        if set(winner_map)!=segmenters:raise ValueError('Missing or unexpected segmenter winners')
        by_id={c['configuration_id']:c for c in configs}
        if any(winner not in by_id or by_id[winner]['segmenter']!=family for family,winner in winner_map.items()):
            raise ValueError('Unknown or wrong-family winning configuration')
        winners=set(winner_map.values())
        configs=[c for c in configs if c['configuration_id'] in winners]
    queries={q['query_id']:q for q in suite['queries']}
    if not ids or len(ids)!=len(set(ids)) or not set(ids)<=set(queries):raise ValueError('Invalid frozen partition query IDs')
    combined=design['selection_query_ids'] if args.partition=='dev' else design['test_query_ids']
    if any(design['partition']['queries'][qid]!=args.partition or queries[qid]['component_count']<2 for qid in combined):
        raise ValueError('Combined queries cross frozen partition boundary')
    if any(queries[qid]['component_count']!=1 for qid in design['clean_sanity_query_ids']):
        raise ValueError('Clean sanity panel contains a combined query')
    return s,suite,design,[queries[i] for i in ids],configs


def prediction_config(config):
    return {k:v for k,v in config.items() if k not in ['configuration_id','fusion','rrf_constant']}


def validate_proposal(proposal,query,config):
    if proposal['query']!=public_query(query) or proposal['prediction_config']!=config:
        raise ValueError('Frozen proposal query/configuration mismatch')
    prediction=proposal['prediction'];segments=prediction['segments'];length=query['length_aa']
    if not segments or prediction['boundaries_supplied'] is not False or prediction['segment_count']!=len(segments):
        raise ValueError('Missing or boundary-informed span prediction')
    for span in segments:
        if not (type(span['start0']) is int and type(span['end0']) is int and 0<=span['start0']<span['end0']<=length):
            raise ValueError('Predicted span lies outside exact query')
    encoded=sum(span['end0']-span['start0'] for span in segments)
    scan=0 if config['segmenter']=='fixed_windows' else prediction['total_queried_residues']
    if (proposal['logical_encoded_segment_residues']!=encoded or proposal['logical_candidate_scan_residues']!=scan
        or proposal['logical_total_queried_residues']!=scan+encoded):raise ValueError('Prediction work accounting mismatch')
    if proposal['boundary_evaluation']!=interval_metrics(segments,query['spans']):
        raise ValueError('Saved boundary evaluation does not match frozen spans')


def union_results(segment_rankings,query):
    unions={}
    for k in (1,5,10):
        union=per_segment_union(segment_rankings,k);present=set(union['accessions'])
        hits=[a in present for a in query['components']]
        union['evaluation']={'components':len(hits),'component_hits':sum(hits),
            'component_recall':sum(hits)/len(hits),'any_components':int(any(hits)),'all_components':int(all(hits)),
            'extraneous_to_exact_component_labels':len(present-set(query['components'])),
            'extraneous_note':'Outside exact target labels, not necessarily a false functional match'}
        unions[str(k)]=union
    return unions


def validate_retrieval_record(record,query,config,proposal,accessions,equivalence):
    if record['query']!=public_query(query) or record['configuration']!=config:
        raise ValueError('Retrieval record query/configuration mismatch')
    if (record['prediction']!=proposal['prediction'] or record['boundary_evaluation']!=proposal['boundary_evaluation']
        or record['logical_total_queried_residues']!=proposal['logical_total_queried_residues']):
        raise ValueError('Retrieval record differs from frozen label-free proposal')
    if len(record['segment_rankings'])!=len(proposal['prediction']['segments']):raise ValueError('Missing segment rankings')
    for ranking in record['segment_rankings']:
        if len(ranking)!=len(accessions) or {r['accession'] for r in ranking}!=set(accessions):
            raise ValueError('Segment ranking lacks complete catalogue')
        # The common evaluator validates finite scores, order, unique accessions,
        # and rank positions. Its labels have no role in prediction or fusion.
        evaluate_query(ranking,query,equivalence)
    fused=fuse_rankings(record['segment_rankings'],config['fusion'],rrf_constant=config['rrf_constant'],catalogue_accessions=accessions)
    if record['ranking']!=fused['ranking']:raise ValueError('Saved fusion does not match segment rankings')
    if record['evaluation']!=evaluate_query(record['ranking'],query,equivalence):raise ValueError('Saved segmentation metrics mismatch')
    if record['larger_budget']!=union_results(record['segment_rankings'],query):raise ValueError('Per-segment review-volume/metrics mismatch')


def spans(args):
    s,suite,design,queries,configs=setup(args)
    source=open_store(args.alignment_job.resolve(),s)
    if source.identity['mode']!='alignment':raise ValueError('Require raw alignment cache')
    validate_source_suite(source.identity,suite,design['source_challenge_sha256'])
    from iri_constructs.session_alignment import rank_alignment,alignment_metadata
    if source.identity['method']!=alignment_metadata():raise ValueError('Cached and fresh alignment methods differ')
    if source.identity['code_hashes']['src/iri_constructs/session_alignment.py']!=file_hash(ROOT/'src/iri_constructs/session_alignment.py'):
        raise ValueError('Cached and fresh alignment implementations differ')
    unique={digest(canonical(prediction_config(c)))[:12]:prediction_config(c) for c in configs}
    identity={'mode':'automatic_spans_v1','design_sha256':file_hash(args.design),'partition':args.partition,
              'selection_sha256':file_hash(args.selection) if args.selection else None,
              'alignment_identity':source.identity_hash,'alignment_completion_sha256':source.completion_sha256,
              'challenge_sha256':design['source_challenge_sha256'],
              'catalogue_accessions':sorted(suite['sequences']),
              'catalogue_sequence_sha256':{a:digest(seq.encode()) for a,seq in suite['sequences'].items()},
              'configs':unique,'query_ids':[q['query_id'] for q in queries],
              'code':{name:file_hash(ROOT/name) for name in ['scripts/session_segment_experiment.py','scripts/session_experiment.py','src/iri_constructs/session_segmentation.py','src/iri_constructs/session_alignment.py','src/iri_constructs/session_benchmark.py','src/iri_constructs/session_store.py']},
              'prediction_inputs':'sequence and catalogue only; boundary labels used only in subsequent interval evaluation'}
    store=UnitStore(s/'jobs'/(args.name+'-'+digest(canonical(identity))[:12]),identity,s/'job_ledger.jsonl')
    existing_by_hash={q['joined_sequence_sha256']:q['query_id'] for q in suite['queries']}
    expected=[];start=time.monotonic()
    for i,q in enumerate(queries):
        check_deadline(s);fragment_cache={};hits=0;fresh=0
        def ranker(sequence,catalogue):
            nonlocal hits,fresh
            sha=digest(sequence.encode())
            if sha in fragment_cache:hits+=1;return fragment_cache[sha]
            candidate=existing_by_hash.get(sha)
            saved=completed_read(source,'query-'+candidate) if candidate else None
            if saved is not None:
                if saved['query']['joined_sequence_sha256']!=sha or saved['query']['length_aa']!=len(sequence):
                    raise ValueError('Cached alignment sequence identity mismatch')
                result={'ranking':saved['rankings']['local_alignment']};hits+=1
            else:result=rank_alignment(sequence,catalogue);fresh+=1
            fragment_cache[sha]=result
            return result
        for key,cfg in unique.items():
            unit=q['query_id']+'--'+key;expected.append(unit)
            if store.find(unit):continue
            t=time.monotonic()
            if cfg['segmenter']=='fixed_windows':
                prediction=fixed_windows(len(q['sequence']),cfg['width'],cfg['stride'])
                scan=0
            else:
                prediction=alignment_guided_segments(q['sequence'],suite['sequences'],SegmentationConfig(**cfg['parameters']),ranker=ranker)
                scan=prediction['total_queried_residues']
            encoded=sum(v['end0']-v['start0'] for v in prediction['segments'])
            record={'query':public_query(q),'prediction':prediction,'prediction_config':cfg,
                 'boundary_evaluation':interval_metrics(prediction['segments'],q['spans']),
                 'logical_candidate_scan_residues':scan,'logical_encoded_segment_residues':encoded,
                 'logical_total_queried_residues':scan+encoded,'seconds':time.monotonic()-t,
                 'query_cumulative_cache_hits':hits,'query_cumulative_fresh_alignment_calls':fresh}
            validate_proposal(record,q,cfg)
            store.publish(unit,record)
        if i%25==0:print(json.dumps({'stage':'spans','query':i+1,'total':len(queries),'seconds':time.monotonic()-start}),flush=True)
    if set(expected)!=expected_keys(identity):raise ValueError('Span expected-unit coverage mismatch')
    for q in queries:
        for key,cfg in unique.items():validate_proposal(store.read(q['query_id']+'--'+key),q,cfg)
    store.complete(expected);print(json.dumps({'status':'complete','job':str(store.root)}),flush=True)


def infer(args):
    from iri_constructs.session_models import load_encoder
    s,suite,design,queries,configs=setup(args)
    predictions=open_store(args.spans_job.resolve(),s);parent=open_store(args.model_job.resolve(),s)
    if parent.identity['mode']!='model':raise ValueError('Catalogue must come from a completed model job')
    validate_source_suite(parent.identity,suite,design['source_challenge_sha256'])
    if parent.identity['code_hashes']['scripts/session_experiment.py']!=file_hash(ROOT/'scripts/session_experiment.py'):
        raise ValueError('Cached catalogue pooling implementation differs from current segment pooling')
    validate_source_suite(predictions.identity,suite,design['source_challenge_sha256'])
    expected_predictions={digest(canonical(prediction_config(c)))[:12]:prediction_config(c) for c in configs}
    if (predictions.identity['mode']!='automatic_spans_v1'
        or predictions.identity['design_sha256']!=file_hash(args.design)
        or predictions.identity['partition']!=args.partition
        or predictions.identity['query_ids']!=[q['query_id'] for q in queries]
        or predictions.identity['configs']!=expected_predictions
        or predictions.identity['selection_sha256']!=(file_hash(args.selection) if args.selection else None)):
        raise ValueError('Prediction job differs from frozen partition/configuration/selection')
    if args.partition=='dev':validate_dev_representation(parent.identity['method'],design)
    check_deadline(s)
    encoder=load_encoder(parent.identity['config'])
    try:
        if model_identity(encoder)!=parent.identity['method']:raise ValueError('Current model identity differs from cached catalogue')
        return infer_loaded(args,s,suite,design,queries,configs,predictions,parent,encoder)
    finally:encoder.close()


def infer_loaded(args,s,suite,design,queries,configs,predictions,parent,encoder):
    import numpy as np
    accessions=sorted(suite['sequences']);refs=[];catalogue_by_hash={}
    for a in accessions:
        record,arrays=completed_read(parent,'catalogue-'+a,arrays=True)
        vector=arrays['mean']
        if (set(arrays)!= {'mean'} or record['accession']!=a
            or record['sequence_sha256']!=digest(suite['sequences'][a].encode())
            or record['length']!=len(suite['sequences'][a])
            or vector.shape!=(record['encode']['residue_shape'][1],) or vector.dtype!=np.float32
            or not np.isfinite(vector).all()):raise ValueError('Source catalogue vector metadata/shape mismatch')
        refs.append(vector)
        if record['sequence_sha256'] in catalogue_by_hash and not np.array_equal(catalogue_by_hash[record['sequence_sha256']],vector):
            raise ValueError('Identical catalogue sequences have different cached vectors')
        catalogue_by_hash[record['sequence_sha256']]=arrays['mean']
    refs=np.stack(refs)
    identity={'mode':'automatic_segment_retrieval_v1','design_sha256':file_hash(args.design),'partition':args.partition,
              'selection_sha256':file_hash(args.selection) if args.selection else None,
              'predictions_identity':predictions.identity_hash,'parent_model_identity':parent.identity_hash,
              'predictions_completion_sha256':predictions.completion_sha256,'parent_completion_sha256':parent.completion_sha256,
              'challenge_sha256':design['source_challenge_sha256'],'catalogue_accessions':accessions,
              'catalogue_sequence_sha256':{a:digest(seq.encode()) for a,seq in suite['sequences'].items()},
              'model':model_identity(encoder),'query_ids':[q['query_id'] for q in queries],'configs':configs,
              'code_sha256':file_hash(__file__),'segmentation_code_sha256':file_hash(ROOT/'src/iri_constructs/session_segmentation.py'),
              'pooling_code_sha256':file_hash(ROOT/'scripts/session_experiment.py'),
              'metrics_code_sha256':file_hash(ROOT/'src/iri_constructs/session_benchmark.py'),
              'pooling':'float64 residue mean stored float32; cosine float64'}
    store=UnitStore(s/'jobs'/(args.name+'-'+digest(canonical(identity))[:12]),identity,s/'job_ledger.jsonl')
    cache_id={'model':model_identity(encoder),'pooling':identity['pooling'],'code_sha256':file_hash(__file__),
              'pooling_code_sha256':identity['pooling_code_sha256']}
    cache=UnitStore(s/'caches'/('segments-'+digest(canonical(cache_id))[:16]),cache_id,s/'job_ledger.jsonl')
    expected=[];start=time.monotonic()
    try:
        for i,q in enumerate(queries):
            check_deadline(s)
            rankings_by_prediction={}
            for cfg in configs:
                unit=q['query_id']+'--'+cfg['configuration_id'];expected.append(unit)
                if store.find(unit):continue
                prediction_key=digest(canonical(prediction_config(cfg)))[:12]
                proposal=completed_read(predictions,q['query_id']+'--'+prediction_key)
                if proposal is None:raise ValueError('Missing frozen span prediction')
                validate_proposal(proposal,q,prediction_config(cfg))
                t=time.monotonic();segment_rankings=[];cache_hits=0;new_encodes=0
                reused=prediction_key in rankings_by_prediction
                if prediction_key in rankings_by_prediction:
                    segment_rankings=rankings_by_prediction[prediction_key]
                else:
                    for span in proposal['prediction']['segments']:
                        seq=q['sequence'][span['start0']:span['end0']];sha=digest(seq.encode())
                        if sha in catalogue_by_hash:
                            vector=catalogue_by_hash[sha];cache_hits+=1
                        else:
                            saved=cache.read(sha,arrays=True)
                            if saved:
                                meta,arrays=saved;vector=arrays['mean'];cache_hits+=1
                                if (set(arrays)!= {'mean'} or meta['sequence_sha256']!=sha or meta['length']!=len(seq)
                                    or vector.shape!=(refs.shape[1],) or vector.dtype!=np.float32
                                    or not np.isfinite(vector).all()):raise ValueError('Segment cache metadata/shape mismatch')
                            else:
                                states=encoder.encode(seq);vector=mean_residues(states);new_encodes+=1
                                cache.publish(sha,{'sequence_sha256':sha,'length':len(seq),'encode':encoder.last_encode_metadata,
                                                   'source_query':q['query_id'],'source_span':span}, {'mean':vector})
                                del states
                        segment_rankings.append(cosine_ranking(vector,refs,accessions,dtype=np.float64))
                    rankings_by_prediction[prediction_key]=segment_rankings
                fusion=fuse_rankings(segment_rankings,cfg['fusion'],rrf_constant=cfg['rrf_constant'],catalogue_accessions=accessions)
                if fusion['segment_count']!=len(proposal['prediction']['segments']) or not fusion['ranking']:
                    raise ValueError('Incomplete segment ranking/fusion output')
                evaluation=evaluate_query(fusion['ranking'],q,suite['equivalence'])
                unions=union_results(segment_rankings,q)
                record={'query':public_query(q),'configuration':cfg,'prediction_unit':q['query_id']+'--'+prediction_key,
                  'prediction':proposal['prediction'],'boundary_evaluation':proposal['boundary_evaluation'],
                  'ranking':fusion['ranking'],'segment_rankings':segment_rankings,'evaluation':evaluation,'larger_budget':unions,
                  'logical_total_queried_residues':proposal['logical_total_queried_residues'],'seconds':time.monotonic()-t,
                  'ranking_reused_within_query':reused,'cache_hits':cache_hits,'new_encodes':new_encodes,'peak_process_rss_mib':peak_mib()}
                validate_retrieval_record(record,q,cfg,proposal,accessions,suite['equivalence'])
                store.publish(unit,record)
            if i%25==0:print(json.dumps({'stage':'segment_inference','query':i+1,'total':len(queries),'seconds':time.monotonic()-start}),flush=True)
        # Release model memory before reading and validating all saved outputs.
        encoder.close()
        summaries={};summaries_by_suite={};candidate_results=[]
        for cfg in configs:
            evaluations=[];selection_rows=[]
            prediction_key=digest(canonical(prediction_config(cfg)))[:12]
            for q in queries:
                record=store.read(q['query_id']+'--'+cfg['configuration_id'])
                proposal=completed_read(predictions,q['query_id']+'--'+prediction_key)
                validate_retrieval_record(record,q,cfg,proposal,accessions,suite['equivalence'])
                evaluations.append(record['evaluation'])
                if args.partition=='dev':
                    selection_rows.append({'query_id':q['query_id'],
                        'component_recall_at_5':record['evaluation']['query_metrics']['exact']['component_recall_at_5'],
                        'all_components_hit_at_5':record['evaluation']['query_metrics']['exact']['all_components_hit_at_5'],
                        'total_queried_residues':record['logical_total_queried_residues']})
            summaries[cfg['configuration_id']]=summarize(evaluations)
            summaries_by_suite[cfg['configuration_id']]={suite_name:summarize([e for e in evaluations if e['suite']==suite_name])
                for suite_name in sorted({e['suite'] for e in evaluations})}
            if args.partition=='dev':candidate_results.append({'configuration_id':cfg['configuration_id'],'rows':selection_rows})
        selection=select_configuration(candidate_results,design) if args.partition=='dev' else None
        expected.append('summary')
        scientific={'summaries':summaries,'summaries_by_suite':summaries_by_suite,'selection':selection}
        if store.find('summary') is None:store.publish('summary',{**scientific,'seconds':time.monotonic()-start})
        elif any(store.read('summary').get(key)!=value for key,value in scientific.items()):
            raise ValueError('Recomputed segmentation summary/selection mismatch')
        if set(expected)!=expected_keys(identity):raise ValueError('Segmentation expected-unit coverage mismatch')
        store.complete(expected)
        if selection:
            value={'design_sha256':file_hash(args.design),'development_job':str(store.root),'development_identity':store.identity_hash,
                   'development_summary_descriptor_sha256':file_hash(store.find('summary')/'complete.json'),'selection':selection}
            target=s/'inputs/segmentation_selection_v1.json'
            if target.exists() and json.loads(target.read_text())!=value:raise ValueError('Existing frozen selection differs')
            if not target.exists():atomic_json(target,value)
        print(json.dumps({'status':'complete','job':str(store.root),'selection':selection}),flush=True)
    except BaseException as exc:
        store.event(event='job_failed',state='failed',exception=f'{type(exc).__name__}: {exc}')
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['spans','infer']);p.add_argument('--session',type=Path,required=True)
    p.add_argument('--design',type=Path,required=True);p.add_argument('--partition',choices=['dev','test'],required=True)
    p.add_argument('--selection',type=Path);p.add_argument('--name',required=True);p.add_argument('--alignment-job',type=Path)
    p.add_argument('--spans-job',type=Path);p.add_argument('--model-job',type=Path);a=p.parse_args()
    spans(a) if a.stage=='spans' else infer(a)
