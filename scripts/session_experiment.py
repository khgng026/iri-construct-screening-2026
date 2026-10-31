#!/usr/bin/env python3
"""Resumable shared-suite inference and alignment, one owned worker at a time."""
from __future__ import annotations
import argparse
import datetime as dt
import importlib.metadata
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
import traceback

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('HF_HUB_DISABLE_TELEMETRY','1')
os.environ.setdefault('TOKENIZERS_PARALLELISM','false')
from iri_constructs.session_store import UnitStore, atomic_json, append_event, canonical, digest, file_hash, utcnow


def peak_mib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(2**20 if sys.platform=='darwin' else 1024)


def code_identity():
    paths=[ROOT/'src/iri_constructs'/name for name in
           ['session_store.py','session_models.py','session_benchmark.py','session_alignment.py']]+[Path(__file__)]
    return {str(p.relative_to(ROOT)):file_hash(p) for p in sorted(paths)}


def check_deadline(session):
    end=dt.datetime.fromisoformat(json.loads((session/'launch.json').read_text())['hard_stop_utc'])
    if dt.datetime.now(dt.timezone.utc)>=end:
        raise TimeoutError('Session hard stop reached')


def public_query(query):
    return {k:v for k,v in query.items() if k!='sequence'}


def mean_residues(states):
    import numpy as np
    return np.mean(states,axis=0,dtype=np.float64).astype(np.float32)


def model_identity(encoder):
    return {k:v for k,v in encoder.metadata.items() if k!='load_seconds'}


def run(args):
    import numpy as np
    from iri_constructs.session_benchmark import load_suite, cosine_ranking, evaluate_query, evaluate_oracle_query, summarize
    from iri_constructs.session_alignment import rank_alignment, alignment_metadata
    session=args.session.resolve(); check_deadline(session)
    config=json.loads(args.config.read_text())
    suite=load_suite(ROOT,session/'inputs')
    queries=suite['queries']
    if args.subset=='repeat':
        wanted=set(suite['manifest']['repeat_query_ids']); queries=[q for q in queries if q['query_id'] in wanted]
    queries=sorted(queries,key=lambda q:({'clean':0,'pairs':1,'triples':2,'quadruples':3}[q['suite']],q['query_id']))
    if args.limit: queries=queries[:args.limit]
    accessions=sorted(suite['sequences'])
    encoder=None
    started=time.monotonic()
    if args.mode=='model':
        from iri_constructs.session_models import load_encoder
        encoder=load_encoder(config)
        method=model_identity(encoder)
    else:
        method=alignment_metadata()
    identity={'version':1,'mode':args.mode,'config':config,'method':method,
              'challenge_sha256':file_hash(session/'inputs/challenge_manifest.json'),
              'code_hashes':code_identity(),'code_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'catalogue_accessions':accessions,'catalogue_sequence_sha256':{a:digest(suite['sequences'][a].encode()) for a in accessions},
              'query_ids':[q['query_id'] for q in queries],
              'numerical_policy': 'residue reductions float64 then stored float32; cosine float64; raw means before aggregation',
              'environment':{p:importlib.metadata.version(p) for p in ['numpy','biopython']},
              'oracle_policy':'separately encoded complete component vectors cached by exact sequence; clean queries independently re-encoded'}
    job_id=args.name+'-'+digest(canonical(identity))[:12]
    store=UnitStore(session/'jobs'/job_id,identity,session/'job_ledger.jsonl')
    store.event(event='job_started',job_id=job_id,state='running',pid=os.getpid(),query_count=len(queries))
    atomic_json(store.root/'invocation-'+str(os.getpid())+'.json',{'argv':sys.argv,'started_utc':utcnow(),'pid':os.getpid()})
    expected=[]
    try:
        refs=None
        if encoder:
            means=[]
            for accession in accessions:
                check_deadline(session)
                key='catalogue-'+accession; expected.append(key)
                saved=store.read(key,arrays=True)
                if saved is None:
                    t=time.monotonic(); states=encoder.encode(suite['sequences'][accession]); v=mean_residues(states)
                    store.publish(key,{'accession':accession,'sequence_sha256':digest(suite['sequences'][accession].encode()),
                                      'length':len(states),'encode':encoder.last_encode_metadata,'seconds':time.monotonic()-t,
                                      'peak_process_rss_mib':peak_mib()}, {'mean':v})
                    del states
                else:
                    old,arrays=saved;v=arrays['mean']
                    if (v.ndim!=1 or not np.isfinite(v).all() or v.shape!=(old['encode']['residue_shape'][1],)
                        or old['accession']!=accession or old['sequence_sha256']!=digest(suite['sequences'][accession].encode())
                        or old['length']!=len(suite['sequences'][accession])):
                        raise ValueError('Invalid cached catalogue metadata/shape/value')
                means.append(v)
            refs=np.stack(means)
            reference_map=dict(zip(accessions,means))
            isolated_rankings={a:cosine_ranking(reference_map[a],refs,accessions,dtype=np.float64) for a in suite['manifest']['targets']}
        for index,query in enumerate(queries):
            check_deadline(session)
            key='query-'+query['query_id'];expected.append(key)
            if store.find(key): continue
            t=time.monotonic()
            record={'query':public_query(query),'rankings':{},'evaluations':{},'controls':{}}
            arrays={}
            if encoder:
                states=encoder.encode(query['sequence'])
                whole=mean_residues(states)
                contextual={s['accession']:mean_residues(states[s['start0']:s['end0']]) for s in query['spans']}
                lengths={s['accession']:s['end0']-s['start0'] for s in query['spans']}
                # Canonical summation order makes isolated mixture invariant to query order.
                component_ids=sorted(query['components']); total=query['length_aa']
                mixture=sum(lengths[a]*reference_map[a].astype(np.float64) for a in component_ids)/total
                reconstructed=sum(lengths[a]*contextual[a].astype(np.float64) for a in component_ids)/total
                error=float(np.max(np.abs(whole.astype(np.float64)-reconstructed)))
                if error>1e-6: raise ValueError(f'Weighted contextual identity failed: {error}')
                reverse=sum(lengths[a]*reference_map[a].astype(np.float64) for a in reversed(component_ids))/total
                order_error=float(np.max(np.abs(mixture-reverse)))
                if order_error>1e-6: raise ValueError('Isolated mixture order invariance failed')
                record['controls']={'whole_weighted_context_max_abs':error,'mixture_reverse_max_abs':order_error,
                                    'tolerance':1e-6,'context_vs_isolated':{}}
                arrays={'whole':whole,'isolated_mixture':mixture.astype(np.float32)}
                for a in component_ids:
                    arrays['contextual_'+a]=contextual[a];arrays['isolated_'+a]=reference_map[a]
                    va,vb=contextual[a].astype(np.float64),reference_map[a].astype(np.float64)
                    record['controls']['context_vs_isolated'][a]={'cosine':float(va@vb/(np.linalg.norm(va)*np.linalg.norm(vb))),
                              'relative_l2_change':float(np.linalg.norm(va-vb)/np.linalg.norm(vb)),
                              'contextual_norm':float(np.linalg.norm(va)),'isolated_norm':float(np.linalg.norm(vb))}
                record['rankings']={'whole':cosine_ranking(whole,refs,accessions,dtype=np.float64),
                   'isolated_mixture':cosine_ranking(arrays['isolated_mixture'],refs,accessions,dtype=np.float64),
                   'contextual_component':{a:cosine_ranking(contextual[a],refs,accessions,dtype=np.float64) for a in component_ids},
                   'isolated_component':{a:isolated_rankings[a] for a in component_ids}}
                record['encode']=encoder.last_encode_metadata
                if len(component_ids)==1:
                    diff=float(np.max(np.abs(whole-reference_map[component_ids[0]])))
                    record['controls']['independent_clean_reference_max_abs']=diff
                    if diff>1e-6: raise ValueError('Independent clean/reference encoding differs')
                del states
            else:
                alignment=rank_alignment(query['sequence'],suite['sequences'])
                record['rankings']={'local_alignment':alignment['ranking']}
            for arm,ranking in record['rankings'].items():
                evaluator=evaluate_oracle_query if isinstance(ranking,dict) else evaluate_query
                record['evaluations'][arm]=evaluator(ranking,query,suite['equivalence'])
            record.update(seconds=time.monotonic()-t,peak_process_rss_mib=peak_mib())
            store.publish(key,record,arrays if arrays else None)
            if index%25==0 or index+1==len(queries):
                print(json.dumps({'job_id':job_id,'completed_queries':index+1,'expected_queries':len(queries),
                                  'elapsed_seconds':time.monotonic()-started,'peak_rss_mib':peak_mib()}),flush=True)
        evaluations={}
        for query in queries:
            record=store.read('query-'+query['query_id'])
            for arm,value in record['evaluations'].items():evaluations.setdefault(arm,[]).append(value)
        summary={arm:{'all':summarize(values),
                       'by_suite':{s:summarize([v for v in values if v['suite']==s]) for s in sorted({v['suite'] for v in values})}}
                 for arm,values in evaluations.items()}
        # Final summary is itself checksum-validated and immutable.
        expected.append('summary')
        if store.find('summary') is None:
            store.publish('summary',{'summary':summary,'elapsed_seconds':time.monotonic()-started,
                                     'peak_process_rss_mib':peak_mib(),'query_count':len(queries),
                                     'model_load_seconds':encoder.load_seconds if encoder else None})
        elif store.read('summary')['summary'] != summary:
            raise ValueError('Recomputed scientific summary differs from saved summary; preserve job and investigate before completion')
        complete=store.complete(expected,{'query_count':len(queries),'peak_process_rss_mib':peak_mib()})
        print(json.dumps({'status':'complete','job_id':job_id,'path':str(store.root),'units':complete['expected_count']}),flush=True)
    except BaseException as exc:
        store.event(event='job_failed',job_id=job_id,state='failed',exception=f'{type(exc).__name__}: {exc}',traceback=traceback.format_exc())
        raise
    finally:
        if encoder: encoder.close()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--session',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
    p.add_argument('--name',required=True);p.add_argument('--mode',choices=['model','alignment'],default='model')
    p.add_argument('--subset',choices=['full','repeat'],default='full');p.add_argument('--limit',type=int)
    run(p.parse_args())
