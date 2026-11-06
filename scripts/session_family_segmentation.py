#!/usr/bin/env python3
"""Exploratory unchanged-parameter segmentation transfer after reference removal.

Run spans and infer as separate guarded jobs. Prediction sees sequence and the
remaining common catalogue only. Family labels enter subsequent evaluation.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.session_store import UnitStore,canonical,digest,file_hash
from iri_constructs.session_benchmark import (load_suite,cosine_ranking,filter_ranking,evaluate_query)
from iri_constructs.session_segmentation import (SegmentationConfig,fixed_windows,alignment_guided_segments,
                                                fuse_rankings,per_segment_union,interval_metrics)
from session_experiment import check_deadline,public_query,mean_residues,model_identity,peak_mib
from session_segment_experiment import (open_store,completed_read,setup as primary_setup,
                                        prediction_config,validate_proposal,validate_source_suite,
                                        validate_dev_representation)
from session_catalogue_analysis import (validate_family_freeze,SummaryAccumulator,subset_denominators,
                                        denominator_check,latest_completion)
from session_independent_audit import verify_unit,validate_ranking


def require(value,message):
    if not value:raise ValueError(message)


def read(path):return json.loads(Path(path).read_text())


def transfer_setup(args):
    session=args.session.resolve();check_deadline(session)
    design=read(args.design)
    require(design['version']=='family-segmentation-transfer-v1'
            and design['status']=='prospective_frozen_before_transfer_scoring','Unknown or unfrozen transfer design')
    primary=session/'inputs'/design['primary_design_filename']
    family=session/'inputs'/design['family_labels_filename']
    selection=session/'inputs'/design['primary_selection_filename']
    require(file_hash(primary)==design['primary_design_sha256'],'Primary segmentation design changed')
    require(file_hash(family)==design['family_labels_sha256'],'Frozen family-label bytes changed')
    require(selection.exists(),'Primary development selection is not complete; transfer remains blocked')
    # This verifies complete development evidence, representation, every candidate,
    # summary descriptor and each winning configuration's correct segmenter family.
    session,suite,primary_design,primary_queries,configs=primary_setup(SimpleNamespace(
        session=session,design=primary,partition='test',selection=selection))
    require(len(primary_design['test_query_ids'])==design['required_combined_test_queries']==950,
            'Transfer combined panel differs from prospective freeze')
    require(set(design['policies'])=={'exact_removed','near90_removed'},'Transfer mask policies changed')
    require(design['retuning_allowed'] is False,'Transfer does not permit retuning')
    ids=list(primary_design['test_query_ids'])
    if args.include_clean:ids+=primary_design['clean_sanity_query_ids']
    query_map={q['query_id']:q for q in primary_queries}
    require(len(ids)==len(set(ids)) and set(ids)<=set(query_map),'Transfer query coverage invalid')
    require(all(query_map[q]['component_count']>1 for q in primary_design['test_query_ids']),
            'Primary test panel contains a clean query')
    labels=read(family)
    audit=validate_family_freeze(labels,suite,family,session/'audit/independent_family_label_review.json')
    return session,suite,design,primary_design,labels,[query_map[q] for q in ids],configs,selection,audit


def unit_key(qid,policy,configuration_id):return qid+'--'+policy+'--'+configuration_id


def mask_for_query(query,suite,labels,policy):
    removed=labels['removal_policies'][policy]['common_removed_accessions_by_group'][query['group_id']]
    independent=sorted(set().union(*(set(labels['removal_policies'][policy]['removed_accessions_by_target'][a])
                                    for a in query['components'])))
    require(removed==independent,'Mask must remove relatives of all components together')
    included=sorted(set(suite['sequences'])-set(removed))
    require(included,'No remaining references for transfer')
    return {'policy':policy,'removed_accessions':removed,'included_accessions':included,
            'included_accessions_sha256':digest(canonical(included)),
            'remaining_catalogue_count':len(included)}


def expected_transfer_keys(identity):
    require(identity['query_ids'] and len(identity['query_ids'])==len(set(identity['query_ids'])),
            'Transfer identity query IDs empty/duplicate')
    ids=[c['configuration_id'] for c in identity['configs']]
    require(ids and len(ids)==len(set(ids)),'Transfer configuration IDs empty/duplicate')
    policies=identity['policies'];require(set(policies)=={'exact_removed','near90_removed'} and len(policies)==2,
                                         'Transfer policies missing/duplicate')
    keys={unit_key(q,p,c) for q in identity['query_ids'] for p in policies for c in ids}
    if identity['mode']=='family_automatic_retrieval_v1':keys.add('summary')
    else:require(identity['mode']=='family_automatic_spans_v1','Unknown transfer mode')
    return keys


def open_transfer_store(path,session):
    path=Path(path).resolve();identity=read(path/'identity.json')
    store=UnitStore(path,identity,session/'job_ledger.jsonl')
    completion_path=latest_completion(path);completion=read(completion_path)
    expected=expected_transfer_keys(identity)
    require(completion['identity_hash']==store.identity_hash and completion['expected_count']==len(expected)
            and set(completion['units'])==expected,'Transfer completion lacks exact expected unit coverage')
    store.completed_paths={key:verify_unit(path,key,index,store.identity_hash)
                           for key,index in completion['units'].items()}
    store.completion_sha256=file_hash(completion_path)
    return store


def shared_identity(args,suite,design,queries,configs,selection,audit):
    return {'transfer_design_sha256':file_hash(args.design),'primary_design_sha256':design['primary_design_sha256'],
            'selection_sha256':file_hash(selection),'family_labels_sha256':design['family_labels_sha256'],
            'independent_family_audit':audit,'challenge_sha256':suite['verification']['challenge_manifest_sha256'],
            'catalogue_accessions':sorted(suite['sequences']),
            'catalogue_sequence_sha256':{a:digest(seq.encode()) for a,seq in suite['sequences'].items()},
            'query_ids':[q['query_id'] for q in queries],'configs':configs,'policies':design['policies'],
            'include_clean_sanity':args.include_clean,'parameter_transfer':'unchanged winners from completed original development only',
            'prediction_mask_policy':'same all-components removal mask is applied during prediction and retrieval; no target-boundary reuse',
            'code_hashes':{name:file_hash(ROOT/name) for name in [
                'scripts/session_family_segmentation.py','scripts/session_segment_experiment.py',
                'scripts/session_experiment.py','scripts/session_catalogue_analysis.py','scripts/session_independent_audit.py',
                'src/iri_constructs/session_segmentation.py','src/iri_constructs/session_alignment.py',
                'src/iri_constructs/session_models.py','src/iri_constructs/session_benchmark.py','src/iri_constructs/session_store.py']}}


def masked_alignment_ranker(remaining_catalogue,full_catalogue,source_lookup,fresh_ranker):
    """Exact fragment-hash cache reuse, with masking BEFORE candidate prediction.

    ``source_lookup`` returns a checksum-verified complete full-catalogue record
    or None. No target accession, constructed span, or expected answer enters.
    """
    included=sorted(remaining_catalogue);full=sorted(full_catalogue);cache={}
    counters={'memory_cache_hits':0,'source_full_alignment_hits':0,'fresh_remaining_catalogue_alignments':0,
              'fresh_reference_comparisons':0}
    def ranker(sequence,catalogue):
        require(catalogue==remaining_catalogue,'Prediction attempted to use another catalogue mask')
        sequence_hash=digest(sequence.encode());key=(sequence_hash,digest(canonical(included)))
        if key in cache:counters['memory_cache_hits']+=1;return cache[key]
        saved=source_lookup(sequence_hash,len(sequence))
        if saved is not None:
            require(saved['sequence_sha256']==sequence_hash and saved['length']==len(sequence),
                    'Raw alignment cache sequence identity differs')
            validate_ranking(saved['ranking'],full)
            ranking=filter_ranking(saved['ranking'],included)
            counters['source_full_alignment_hits']+=1
        else:
            result=fresh_ranker(sequence,remaining_catalogue)
            ranking=result['ranking'] if isinstance(result,dict) else result
            counters['fresh_remaining_catalogue_alignments']+=1
            counters['fresh_reference_comparisons']+=len(included)
        validate_ranking(ranking,included)
        require(all(row['accession'] in remaining_catalogue for row in ranking),'Removed accession leaked into predictor')
        result={'ranking':ranking};cache[key]=result;return result
    return ranker,counters


def validate_family_proposal(proposal,query,config,mask):
    validate_proposal(proposal,query,prediction_config(config))
    require(proposal['common_mask']==mask,'Proposal catalogue mask differs from frozen all-component mask')
    require(proposal['prediction_catalogue_accessions']==mask['included_accessions'],
            'Prediction did not use remaining catalogue')
    for span in proposal['prediction']['segments']:
        hit=span.get('seed_hit')
        if hit is not None:require(hit['accession'] in mask['included_accessions'],'Removed-reference hit informed predicted boundary')


def spans(args):
    session,suite,design,primary,labels,queries,configs,selection,audit=transfer_setup(args)
    require(args.alignment_job is not None,'Spans stage requires complete raw alignment cache')
    source=open_store(args.alignment_job,session)
    require(source.identity['mode']=='alignment','Require primary raw alignment source')
    validate_source_suite(source.identity,suite,suite['verification']['challenge_manifest_sha256'])
    from iri_constructs.session_alignment import rank_alignment,alignment_metadata
    require(source.identity['method']==alignment_metadata(),'Raw/fresh alignment methods differ')
    require(source.identity['code_hashes']['src/iri_constructs/session_alignment.py']==file_hash(ROOT/'src/iri_constructs/session_alignment.py'),
            'Cached/fresh alignment implementation differs')
    identity={**shared_identity(args,suite,design,queries,configs,selection,audit),
              'mode':'family_automatic_spans_v1','raw_alignment_identity':source.identity_hash,
              'raw_alignment_completion_sha256':source.completion_sha256}
    store=UnitStore(session/'jobs'/(args.name+'-'+digest(canonical(identity))[:12]),identity,session/'job_ledger.jsonl')
    store.event(event='job_started',state='running',pid=os.getpid(),stage='family_spans')
    existing={q['joined_sequence_sha256']:q['query_id'] for q in suite['queries']}
    expected=[];started=time.monotonic()
    def lookup(sha,length):
        qid=existing.get(sha);saved=completed_read(source,'query-'+qid) if qid else None
        if saved is None:return None
        return {'sequence_sha256':saved['query']['joined_sequence_sha256'],'length':saved['query']['length_aa'],
                'ranking':saved['rankings']['local_alignment']}
    try:
        for index,q in enumerate(queries):
            check_deadline(session);window_predictions={}
            for policy in design['policies']:
                mask=mask_for_query(q,suite,labels,policy)
                remaining={a:suite['sequences'][a] for a in mask['included_accessions']}
                ranker,counters=masked_alignment_ranker(remaining,suite['sequences'],lookup,rank_alignment)
                for config in configs:
                    key=unit_key(q['query_id'],policy,config['configuration_id']);expected.append(key)
                    old=store.read(key)
                    if old is not None:validate_family_proposal(old,q,config,mask);continue
                    start=time.monotonic();cfg=prediction_config(config)
                    if cfg['segmenter']=='fixed_windows':
                        geometry=digest(canonical(cfg))
                        if geometry not in window_predictions:
                            window_predictions[geometry]=fixed_windows(len(q['sequence']),cfg['width'],cfg['stride'])
                        prediction=window_predictions[geometry];scan=0
                    else:
                        prediction=alignment_guided_segments(q['sequence'],remaining,SegmentationConfig(**cfg['parameters']),ranker=ranker)
                        scan=prediction['total_queried_residues']
                    encoded=sum(v['end0']-v['start0'] for v in prediction['segments'])
                    record={'query':public_query(q),'prediction_config':cfg,'prediction':prediction,
                            'common_mask':mask,'prediction_catalogue_accessions':mask['included_accessions'],
                            'boundary_evaluation':interval_metrics(prediction['segments'],q['spans']),
                            'logical_candidate_scan_residues':scan,'logical_encoded_segment_residues':encoded,
                            'logical_total_queried_residues':scan+encoded,'mask_cumulative_cache_counts':dict(counters),
                            'seconds':time.monotonic()-start,'peak_process_rss_mib':peak_mib(),
                            'boundary_labels_used_for_prediction':False}
                    validate_family_proposal(record,q,config,mask);store.publish(key,record)
            if (index+1)%25==0:print(json.dumps({'stage':'family_spans','queries':index+1,'total':len(queries)}),flush=True)
        require(set(expected)==expected_transfer_keys(identity),'Family span completion omitted units')
        for q in queries:
            for policy in design['policies']:
                for cfg in configs:validate_family_proposal(store.read(unit_key(q['query_id'],policy,cfg['configuration_id'])),
                                                            q,cfg,mask_for_query(q,suite,labels,policy))
        store.complete(expected,{'elapsed_seconds':time.monotonic()-started})
        print(json.dumps({'status':'complete','job':str(store.root)}),flush=True)
    except BaseException as exc:
        store.event(event='job_failed',state='failed',stage='family_spans',exception=f'{type(exc).__name__}: {exc}');raise


def family_union_results(segment_rankings,query,acceptable_matches,remaining_accessions):
    """Evaluate the explicitly larger per-segment returned-hit union."""
    remaining=set(remaining_accessions);result={}
    for k in (1,5,10):
        union=per_segment_union(segment_rankings,k);present=set(union['accessions']);components=[]
        for accession in query['components']:
            accepted=set(acceptable_matches[accession])&remaining
            components.append({'accession':accession,'scorable':bool(accepted),
                               'accepted_surviving':sorted(accepted),'hit':int(bool(accepted&present)) if accepted else None,
                               'exact_reference_present':accession in remaining,
                               'exact_hit':int(accession in present) if accession in remaining else None})
        n=sum(c['scorable'] for c in components);hits=sum(c['hit'] for c in components if c['scorable'])
        all_scorable=n==len(components)
        union['evaluation']={'components':components,'original_components':len(components),'scorable_components':n,
                             'unscorable_components':len(components)-n,'component_hits':hits,
                             'component_recall':hits/n if n else None,'all_components_scorable':all_scorable,
                             'any_component_hit':int(hits>0) if all_scorable else None,
                             'all_components_hit':int(hits==len(components)) if all_scorable else None,
                             'unlabelled_returned_hits':len(present-set().union(*(set(c['accepted_surviving']) for c in components))),
                             'unlabelled_note':'Outside accepted source-family sets, not independently verified functional negatives.'}
        result[str(k)]=union
    return result


class UnionAccumulator:
    def __init__(self):
        self.ids=set();self.components=self.scorable=self.hits=self.full_queries=self.any_hits=self.all_hits=0
        self.occurrences=self.unique_hits=self.duplicates=0;self.targets={}
    def add(self,qid,value):
        require(qid not in self.ids,'Duplicate union evaluation');self.ids.add(qid);e=value['evaluation']
        self.components+=e['original_components'];self.scorable+=e['scorable_components'];self.hits+=e['component_hits']
        if e['all_components_scorable']:
            self.full_queries+=1;self.any_hits+=e['any_component_hit'];self.all_hits+=e['all_components_hit']
        self.occurrences+=value['returned_hit_occurrences'];self.unique_hits+=value['unique_returned_hits'];self.duplicates+=value['duplicate_returned_hits']
        for c in e['components']:
            if c['scorable']:
                target=self.targets.setdefault(c['accession'],[0,0]);target[0]+=1;target[1]+=c['hit']
    def result(self):
        n=len(self.ids)
        return {'original_queries':n,'original_components':self.components,'scorable_components':self.scorable,
                'unscorable_components':self.components-self.scorable,'component_hits':self.hits,
                'component_recall':self.hits/self.scorable if self.scorable else None,
                'per_target_macro_component_recall':math.fsum(h/n for n,h in self.targets.values())/len(self.targets) if self.targets else None,
                'all_components_scorable_queries':self.full_queries,'excluded_all_component_queries':n-self.full_queries,
                'any_component_hit':self.any_hits/self.full_queries if self.full_queries else None,
                'all_components_hit':self.all_hits/self.full_queries if self.full_queries else None,
                'returned_hit_occurrences':self.occurrences,'unique_hit_occurrences':self.unique_hits,
                'duplicate_hit_occurrences':self.duplicates,'mean_returned_hit_occurrences':self.occurrences/n if n else None,
                'retrieval_budget':'k per predicted segment'}


def masked_fusion(full_segment_rankings,included,config):
    """Critical order: restrict/rerank EACH segment BEFORE RRF or max-score."""
    masked=[filter_ranking(rows,included) for rows in full_segment_rankings]
    fused=fuse_rankings(masked,config['fusion'],rrf_constant=config['rrf_constant'],catalogue_accessions=included)
    require(masked and fused['ranking'],'Segment fallback must preserve a nonempty retrieval query')
    return masked,fused


def validate_transfer_record(record,query,config,proposal,mask,labels,equivalence,segment_rankings):
    validate_family_proposal(proposal,query,config,mask)
    require(record['query']==public_query(query) and record['configuration']==config and record['common_mask']==mask,
            'Family retrieval query/configuration/mask mismatch')
    require(record['prediction']==proposal['prediction'] and len(segment_rankings)==proposal['prediction']['segment_count'],
            'Family retrieval changed masked span predictions')
    require(record['boundary_evaluation']==proposal['boundary_evaluation']
            and all(record[name]==proposal[name] for name in ('logical_candidate_scan_residues',
                    'logical_encoded_segment_residues','logical_total_queried_residues')),
            'Family retrieval changed boundary evaluation or work accounting')
    require(record['prediction_unit']==unit_key(query['query_id'],mask['policy'],config['configuration_id']),
            'Family retrieval refers to a different prediction unit')
    masked,fused=masked_fusion(segment_rankings,mask['included_accessions'],config)
    require(record['ranking']==fused['ranking'],'Masked fusion does not regenerate')
    expected=evaluate_query(fused['ranking'],query,equivalence,acceptable_matches=labels['acceptable_matches'])
    require(record['evaluation']==expected,'Family retrieval metrics do not regenerate')
    require(expected['query_metrics']['exact']['scorable_components']==0
            and expected['query_metrics']['duplicate_equivalent']['scorable_components']==0,
            'Removed exact/duplicate references must remain unscorable')
    require(record['larger_budget']==family_union_results(masked,query,labels['acceptable_matches'],mask['included_accessions']),
            'Larger-budget family evaluation mismatch')
    frozen=next(x for x in labels['removal_policies'][mask['policy']]['query_scoring_status'] if x['query_id']==query['query_id'])
    require(expected['query_metrics']['labelled']['scorable_components']==frozen['scorable_components']
            and expected['query_metrics']['labelled']['all_components_scorable']==frozen['all_components_scorable'],
            'Transferred family denominators differ from frozen annotation panel')


def validate_parent_development_method(parent_method,selection_path):
    """Require the exact representation used by the already-verified dev winner.

    transfer_setup/primary_setup verify complete development/selection evidence.
    Recheck its immutable identity here before comparing full metadata: nominal
    model/revision/device/layer equality cannot detect backend or library drift.
    """
    selection=read(selection_path)
    development=read(Path(selection['development_job'])/'identity.json')
    require(digest(canonical(development))==selection['development_identity'],
            'Selected development identity changed before transfer inference')
    require(development['mode']=='automatic_segment_retrieval_v1' and development['partition']=='dev',
            'Selection does not refer to original development inference')
    require(parent_method==development['model'],
            'Parent representation differs from exact selected development method')


def infer(args):
    session,suite,design,primary,labels,queries,configs,selection,audit=transfer_setup(args)
    require(args.spans_job is not None and args.model_job is not None,'Inference needs completed masked spans and original8M model')
    predictions=open_transfer_store(args.spans_job,session);parent=open_store(args.model_job,session)
    require(parent.identity['mode']=='model','Require completed parent model catalogue')
    validate_source_suite(parent.identity,suite,suite['verification']['challenge_manifest_sha256'])
    validate_dev_representation(parent.identity['method'],primary)
    validate_parent_development_method(parent.identity['method'],selection)
    require(parent.identity['code_hashes']['scripts/session_experiment.py']==file_hash(ROOT/'scripts/session_experiment.py'),
            'Parent catalogue pooling code differs')
    base=shared_identity(args,suite,design,queries,configs,selection,audit)
    require(predictions.identity['mode']=='family_automatic_spans_v1','Require masked family transfer spans')
    for key,value in base.items():require(predictions.identity.get(key)==value,'Prediction transfer identity mismatch: '+key)
    from iri_constructs.session_models import load_encoder
    check_deadline(session);encoder=load_encoder(parent.identity['config'])
    try:
        require(model_identity(encoder)==parent.identity['method'],'Current model differs from parent catalogue')
        infer_loaded(args,session,suite,labels,queries,configs,predictions,parent,encoder,base)
    finally:encoder.close()


def infer_loaded(args,session,suite,labels,queries,configs,predictions,parent,encoder,base):
    import numpy as np
    accessions=sorted(suite['sequences']);vectors=[];catalogue_by_hash={}
    for accession in accessions:
        meta,arrays=completed_read(parent,'catalogue-'+accession,arrays=True);vector=arrays['mean'];sha=digest(suite['sequences'][accession].encode())
        require(set(arrays)=={'mean'} and meta['accession']==accession and meta['sequence_sha256']==sha
                and meta['length']==len(suite['sequences'][accession]) and vector.shape==(meta['encode']['residue_shape'][1],)
                and vector.dtype==np.float32 and np.isfinite(vector).all(),'Parent catalogue cache invalid')
        if sha in catalogue_by_hash:require(np.array_equal(vector,catalogue_by_hash[sha]),'Duplicate source vectors disagree')
        vectors.append(vector);catalogue_by_hash[sha]=vector
    refs=np.stack(vectors)
    identity={**base,'mode':'family_automatic_retrieval_v1','model':model_identity(encoder),
              'parent_model_identity':parent.identity_hash,'parent_completion_sha256':parent.completion_sha256,
              'prediction_identity':predictions.identity_hash,'prediction_completion_sha256':predictions.completion_sha256,
              'pooling':'float64 residue mean stored float32; cosine float64',
              'fusion_order':'filter and rerank each segment using common all-components mask BEFORE fusion'}
    store=UnitStore(session/'jobs'/(args.name+'-'+digest(canonical(identity))[:12]),identity,session/'job_ledger.jsonl')
    cache_identity={'mode':'family_transfer_segment_cache_v1','model':model_identity(encoder),
                    'pooling':identity['pooling'],'catalogue_sequence_sha256':base['catalogue_sequence_sha256'],
                    'code_hashes':base['code_hashes']}
    cache=UnitStore(session/'caches'/('family-segments-'+digest(canonical(cache_identity))[:16]),cache_identity,session/'job_ledger.jsonl')
    store.event(event='job_started',state='running',pid=os.getpid(),stage='family_infer')
    expected=[];started=time.monotonic()

    def get_segment(sequence,qid,span):
        sha=digest(sequence.encode());saved=cache.read(sha,arrays=True);source='segment_cache'
        if saved is None:
            if sha in catalogue_by_hash:vector=catalogue_by_hash[sha];source='parent_catalogue'
            else:
                states=encoder.encode(sequence);vector=mean_residues(states);del states;source='new_encode'
            ranking=cosine_ranking(vector,refs,accessions,dtype=np.float64)
            cache.publish(sha,{'sequence_sha256':sha,'length':len(sequence),'full_catalogue_ranking':ranking,
                               'source_query_id':qid,'source_span':span,'representation_identity_hash':digest(canonical(model_identity(encoder)))},
                          {'mean':vector})
            saved=cache.read(sha,arrays=True)
        meta,arrays=saved;vector=arrays['mean']
        require(set(arrays)=={'mean'} and meta['sequence_sha256']==sha and meta['length']==len(sequence)
                and meta['representation_identity_hash']==digest(canonical(model_identity(encoder)))
                and vector.shape==(refs.shape[1],) and vector.dtype==np.float32 and np.isfinite(vector).all(),
                'Segment cache representation/sequence/array mismatch')
        validate_ranking(meta['full_catalogue_ranking'],accessions)
        path=cache.find(sha)
        return meta['full_catalogue_ranking'],{'sequence_sha256':sha,'length':len(sequence),'cache_identity_hash':cache.identity_hash,
                 'unit_path':str(path),'descriptor_sha256':file_hash(path/'complete.json'),
                 'ranking_pointer':'/full_catalogue_ranking','source':source}

    def cached_rankings(record,query):
        require(len(record['segment_ranking_references'])==record['prediction']['segment_count'],'Missing segment ranking references')
        rows=[]
        for evidence,span in zip(record['segment_ranking_references'],record['prediction']['segments']):
            sequence=query['sequence'][span['start0']:span['end0']];sha=digest(sequence.encode())
            require(evidence['sequence_sha256']==sha and evidence['length']==len(sequence)
                    and evidence['cache_identity_hash']==cache.identity_hash,'Saved segment cache identity differs')
            path=Path(evidence['unit_path']).resolve()
            require(path.is_relative_to(cache.root.resolve()),'Segment cache evidence escaped cache')
            require(file_hash(path/'complete.json')==evidence['descriptor_sha256'],'Segment cache descriptor changed')
            descriptor=cache.verify(path);require(descriptor['key']==sha,'Segment cache key mismatch')
            meta=read(path/'record.json');require(meta['sequence_sha256']==sha,'Segment ranking sequence mismatch')
            ranking=meta['full_catalogue_ranking'];validate_ranking(ranking,accessions);rows.append(ranking)
        return rows

    try:
        for index,q in enumerate(queries):
            check_deadline(session)
            for policy in base['policies']:
                mask=mask_for_query(q,suite,labels,policy)
                for cfg in configs:
                    key=unit_key(q['query_id'],policy,cfg['configuration_id']);expected.append(key)
                    proposal=completed_read(predictions,key);require(proposal is not None,'Missing complete masked proposal')
                    validate_family_proposal(proposal,q,cfg,mask)
                    old=store.read(key)
                    if old is not None:
                        validate_transfer_record(old,q,cfg,proposal,mask,labels,suite['equivalence'],cached_rankings(old,q));continue
                    start=time.monotonic();rankings=[];references=[];counts=defaultdict(int)
                    for span in proposal['prediction']['segments']:
                        sequence=q['sequence'][span['start0']:span['end0']]
                        ranking,evidence=get_segment(sequence,q['query_id'],span)
                        rankings.append(ranking);references.append(evidence);counts[evidence['source']]+=1
                    masked,fused=masked_fusion(rankings,mask['included_accessions'],cfg)
                    evaluation=evaluate_query(fused['ranking'],q,suite['equivalence'],acceptable_matches=labels['acceptable_matches'])
                    record={'query':public_query(q),'configuration':cfg,'common_mask':mask,'prediction':proposal['prediction'],
                            'boundary_evaluation':proposal['boundary_evaluation'],'prediction_unit':key,
                            'prediction_unit_descriptor_sha256':file_hash(predictions.completed_paths[key]/'complete.json'),
                            'segment_ranking_references':references,'ranking':fused['ranking'],'evaluation':evaluation,
                            'larger_budget':family_union_results(masked,q,labels['acceptable_matches'],mask['included_accessions']),
                            'ranking_recipe':'Read each checksum-pinned full segment ranking; restrict/rerank to common included_accessions; THEN apply unchanged selected fusion. Full fused ranking is stored.',
                            'logical_candidate_scan_residues':proposal['logical_candidate_scan_residues'],
                            'logical_encoded_segment_residues':proposal['logical_encoded_segment_residues'],
                            'logical_total_queried_residues':proposal['logical_total_queried_residues'],
                            'cache_counts':dict(counts),'seconds':time.monotonic()-start,'peak_process_rss_mib':peak_mib()}
                    validate_transfer_record(record,q,cfg,proposal,mask,labels,suite['equivalence'],rankings);store.publish(key,record)
            if (index+1)%25==0:print(json.dumps({'stage':'family_infer','queries':index+1,'total':len(queries)}),flush=True)
        encoder.close()
        summaries={};counts={p:subset_denominators(labels,p,{q['query_id'] for q in queries}) for p in base['policies']}
        for policy in base['policies']:
            for cfg in configs:
                scopes={'all':SummaryAccumulator(),'combined':SummaryAccumulator(),'clean':SummaryAccumulator()};by_suite={}
                union_scopes={scope:{str(k):UnionAccumulator() for k in (1,5,10)} for scope in scopes}
                resources={'cache_counts':defaultdict(int),'logical_candidate_scan_residues':0,
                           'logical_encoded_segment_residues':0,'cap_bound_queries':0,'whole_query_fallback_queries':0}
                for q in queries:
                    key=unit_key(q['query_id'],policy,cfg['configuration_id']);record=store.read(key)
                    proposal=completed_read(predictions,key);mask=mask_for_query(q,suite,labels,policy)
                    validate_transfer_record(record,q,cfg,proposal,mask,labels,suite['equivalence'],cached_rankings(record,q))
                    scope='clean' if q['component_count']==1 else 'combined'
                    for target in ('all',scope):
                        scopes[target].add(record['evaluation'])
                        for k,value in record['larger_budget'].items():union_scopes[target][k].add(q['query_id'],value)
                    by_suite.setdefault(q['suite'],SummaryAccumulator()).add(record['evaluation'])
                    for name,value in record['cache_counts'].items():resources['cache_counts'][name]+=value
                    for name in ('logical_candidate_scan_residues','logical_encoded_segment_residues'):resources[name]+=record[name]
                    resources['cap_bound_queries']+=record['prediction']['cap_bound']
                    resources['whole_query_fallback_queries']+=record['prediction'].get('whole_query_fallback_used',False)
                output={scope:acc.result() for scope,acc in scopes.items()}
                for scope in ('all','combined'):denominator_check(output[scope],counts[policy][scope])
                suite_summaries={name:acc.result() for name,acc in by_suite.items()}
                for name,value in suite_summaries.items():denominator_check(value,counts[policy]['by_suite'][name])
                unions={scope:{k:acc.result() for k,acc in values.items()} for scope,values in union_scopes.items()}
                for scope in ('all','combined'):
                    expected_counts=counts[policy][scope]
                    for value in unions[scope].values():
                        require(value['original_queries']==expected_counts['original_queries']
                                and value['scorable_components']==expected_counts['scorable_component_observations']
                                and value['all_components_scorable_queries']==expected_counts['all_components_scorable_queries'],
                                'Larger union budget denominator drift')
                summaries.setdefault(policy,{})[cfg['configuration_id']]={**output,'by_suite':suite_summaries,
                    'larger_budget_summaries':unions,'resources':{**resources,'cache_counts':dict(resources['cache_counts'])}}
        report={'summaries':summaries,'frozen_family_denominators_for_panel':counts,'query_count':len(queries),
                'transfer_only_no_retuning':True,'common_mask_applied_during_prediction_and_retrieval':True,
                'interpretation':'Exploratory source-family retrieval on new combinations of same proteins; no exact-reference boundary reuse, functional equivalence, calibrated detection or pretrained holdout.'}
        old=store.read('summary')
        if old is None:store.publish('summary',report)
        else:require(old==report,'Recomputed family-transfer summary differs')
        expected.append('summary');require(set(expected)==expected_transfer_keys(identity),'Family retrieval completion omitted units')
        store.complete(expected,{'elapsed_seconds':time.monotonic()-started,'no_family_test_tuning':True})
        print(json.dumps({'status':'complete','job':str(store.root)}),flush=True)
    except BaseException as exc:
        store.event(event='job_failed',state='failed',stage='family_infer',exception=f'{type(exc).__name__}: {exc}');raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['spans','infer'])
    p.add_argument('--session',type=Path,required=True);p.add_argument('--design',type=Path,required=True)
    p.add_argument('--name',required=True);p.add_argument('--alignment-job',type=Path);p.add_argument('--spans-job',type=Path)
    p.add_argument('--model-job',type=Path);p.add_argument('--include-clean',action='store_true')
    args=p.parse_args();spans(args) if args.stage=='spans' else infer(args)
