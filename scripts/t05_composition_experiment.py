#!/usr/bin/env python3
"""Run the frozen ESM-2 composition experiment, offline; sequences stay local."""
import argparse
import csv
from datetime import datetime, timezone
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import resource
import subprocess
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('HF_HUB_DISABLE_TELEMETRY','1')
from iri_constructs.t05_data import sha256, build_query
from iri_constructs.t05_composition import (verify_freeze, load_reviewed, make_plan,
    equivalence_map, component_metrics, query_recall, fraction_bin, cache_read,
    cache_write, canonical_hash)

BASE = ROOT/'results/t05-expanded'
MANIFEST = ROOT/'data/manifests/t05_expanded_catalogue.json'
PLAN = ROOT/'config/t05_expanded_queries.json'
MODEL = ROOT/'config/t05_model.json'
FREEZE = ROOT/'config/t05_expanded_freeze.json'


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def csv_write(path, rows):
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w',newline='') as f:
        w = csv.DictWriter(f,fieldnames=fields)
        w.writeheader(); w.writerows(rows)


def safe_id(name):
    if not re.fullmatch('[A-Za-z0-9_-]+',name):
        raise ValueError('Invalid run ID')
    return name


def memory_mib():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss/(2**20 if sys.platform=='darwin' else 1024)


def save_query(spec, sequences):
    """Save local sequence references and known boundaries; no inferred annotation."""
    qid = spec['query_id']
    sequence = ''.join(sequences[a] for a in spec['components'])
    path = ROOT/f'data/processed/t05-expanded/queries/{qid}.fasta'
    data = f'>{qid}\n{sequence}\n'.encode()
    if path.exists() and path.read_bytes()!=data:
        raise ValueError('Local query changed; refusing overwrite')
    path.parent.mkdir(parents=True,exist_ok=True); path.write_bytes(data)
    features, start = [],0
    for a in spec['components']:
        features.append({'source_record_id':a,'start0':start,'end0':start+len(sequences[a]),
                         'coordinate_system':'protein_0_based_half_open','strand':None,
                         'boundary_source':'known construction, not segmentation',
                         'component_sequence_ref':{'path':f'data/processed/t05-expanded/{a}.fasta',
                         'record_id':a,'format':'fasta','length':len(sequences[a]),
                         'sha256':sha256(sequences[a].encode()),'alphabet':'protein'}})
        start += len(sequences[a])
    return {'query_id':qid,'artificial':spec['artificial'],'features':features,
            'sequence_ref':{'path':str(path.relative_to(ROOT)),'record_id':qid,'format':'fasta',
                            'length':len(sequence),'sha256':sha256(sequence.encode()),'alphabet':'protein'}}


def run(identifier, subset=False, reencode=False):
    out = BASE/safe_id(identifier); out.mkdir(parents=True,exist_ok=False)
    started = perf_counter()
    result = {'run_id':identifier,'started_utc':datetime.now(timezone.utc).isoformat(),
              'status':'running','failures':[],'evaluations':[],'component_rows':[],
              'timings_seconds':{},'subset':subset,'network':'offline local model and source files only'}
    stage = 'validate_inputs'
    try:
        freeze, manifest, plan, model = read(FREEZE),read(MANIFEST),read(PLAN),read(MODEL)
        verify_freeze(ROOT,freeze)
        sequences = load_reviewed(ROOT,manifest)
        rebuilt = make_plan(plan['targets'],sequences)
        if any(plan[k]!=v for k,v in rebuilt.items()):
            raise ValueError('Frozen queries differ from prespecified construction rule')
        ids = sorted(sequences); eq = equivalence_map(sequences)
        specs = [q for q in plan['queries'] if not subset or q['query_id'] in plan['repeat_query_ids']]
        result.update(model=model,freeze_sha256=sha256(FREEZE.read_bytes()),
            frozen_input_hashes=freeze['file_sha256'],catalogue_size=len(ids),
            unique_catalogue_sequences=len(set(sequences.values())),
            duplicate_sequence_groups=manifest['duplicate_sequence_groups'],
            excluded_records=manifest['excluded_records'],excluded_pairs=plan['excluded_pairs'],
            counts={'distinct_targets':len(plan['targets']),
                'distinct_pairs':len({q['pair_id'] for q in specs if q['artificial']}),
                'clean_queries':sum(not q['artificial'] for q in specs),
                'ordered_combined_queries':sum(q['artificial'] for q in specs),
                'distinct_query_strings':len({''.join(sequences[a] for a in q['components']) for q in specs})},
            code_hashes={str(p.relative_to(ROOT)):sha256(p.read_bytes()) for p in
                [Path(__file__),*sorted((ROOT/'src/iri_constructs').glob('*.py'))]},
            git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            git_dirty=bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()))
        stage = 'imports_and_model_load'
        import_start=perf_counter()
        import numpy as np
        import torch
        from huggingface_hub import hf_hub_download
        from iri_constructs.t05_embeddings import Encoder, rank_cosine
        result['timings_seconds']['numerical_imports']=perf_counter()-import_start
        result['environment']={'python':platform.python_version(),'platform':platform.platform(),
            'machine':platform.machine(),'packages':dict(sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions())),
            'cpu':subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip() if sys.platform=='darwin' else platform.processor(),
            'physical_memory_gib':int(subprocess.check_output(['sysctl','-n','hw.memsize'],text=True))/2**30 if sys.platform=='darwin' else None}
        t=perf_counter(); encoder=Encoder(model,offline=True)
        result['timings_seconds']['model_setup_and_load']=perf_counter()-t
        result['timings_seconds']['checkpoint_load']=encoder.load_seconds
        model_hashes = {name:sha256(Path(hf_hub_download(model['model_id'],name,revision=model['revision'],local_files_only=True)).read_bytes())
                        for name in ['model.safetensors','config.json','tokenizer_config.json','special_tokens_map.json','vocab.txt']}
        result['model_file_sha256']=model_hashes
        result['model_details']={'layers':encoder.model.config.num_hidden_layers,
            'embedding_dimensions':encoder.model.config.hidden_size,
            'loaded_parameters':sum(p.numel() for p in encoder.model.parameters()),
            'config_max_position_embeddings':encoder.model.config.max_position_embeddings,
            'operational_residue_cap':1022,'no_truncation':True,'device':'cpu','inference_mode':True}
        result['memory']={'after_load_peak_rss_mib':memory_mib(),'note':'Whole-process resident high-water mark; includes libraries; CPU only.'}
        key_data={'model':model,'model_files':model_hashes,
            'encoder_code':{k:v for k,v in result['code_hashes'].items() if k.endswith(('t05_embeddings.py','t05_data.py'))},
            'runtime':{k:result['environment'][k] for k in ['python','machine','platform']},
            'packages':{p:importlib.metadata.version(p) for p in ['torch','transformers','numpy','tokenizers','safetensors']},
            'sequences':{a:sha256(sequences[a].encode()) for a in ids}}
        cache_dir=BASE/'cache'; t=perf_counter()
        refs=None if reencode else cache_read(cache_dir,key_data,(len(ids),320))
        result['timings_seconds']['cache_lookup']=perf_counter()-t
        result['cache']={'key':canonical_hash(key_data),'hit':refs is not None,'forced_reencode':reencode,'key_data':key_data}
        stage='catalogue_encode'; t=perf_counter(); catalogue_info=[]
        if refs is None:
            vectors=[]
            for a in ids:
                pooled,_,info=encoder.encode([sequences[a]])
                vectors.append(pooled[0].numpy()); catalogue_info.append({'accession':a,**info})
            refs=np.stack(vectors)
            if not reencode:
                cache_write(cache_dir,key_data,refs)
        result['timings_seconds']['catalogue_encode_and_cache_write']=perf_counter()-t
        result['timings_seconds']['catalogue_forward']=sum(x['forward_seconds'] for x in catalogue_info)
        result['catalogue_encoding']=catalogue_info
        references=torch.from_numpy(refs)
        arrays={'catalogue':refs}
        stage='padding_check'; t=perf_counter()
        shortest=min(ids,key=lambda a:(len(sequences[a]),a)); longest=max(ids,key=lambda a:(len(sequences[a]),a))
        padded,_,_=encoder.encode([sequences[shortest],sequences[longest]])
        dif=float((padded[0]-references[ids.index(shortest)]).abs().max())
        result['padding_check']={'shortest':shortest,'longest':longest,'max_abs_difference':dif,'atol':1e-5,'rtol':1e-5,
                                'passed':bool(torch.allclose(padded[0],references[ids.index(shortest)],atol=1e-5,rtol=1e-5))}
        if not result['padding_check']['passed']: raise ValueError('Padding numerical check failed')
        result['timings_seconds']['padding_check']=perf_counter()-t
        clean={}; numerical_max=0.; clean_failures=[]; query_encode=0.; ranking_time=0.; construction_time=0.
        rank_rows=[]; mappings=[]
        refs64=refs.astype(np.float64)
        for number,spec in enumerate(specs):
            qid=spec['query_id']; stage='query_'+qid; t=perf_counter()
            mapping=save_query(spec,sequences); mappings.append(mapping)
            combined=''.join(sequences[a] for a in spec['components'])
            arms=[('whole' if spec['artificial'] else 'clean',combined,spec['components'],None)]
            if spec['artificial']:
                for i,a in enumerate(spec['components']):
                    _,oracle,interval=build_query(spec['components'],a,sequences)
                    arms.append(('oracle_'+str(i+1),oracle,[a],interval))
            construction_time+=perf_counter()-t
            for arm,sequence,expected,interval in arms:
                evaluation_id=qid+'--'+arm
                t=perf_counter(); pooled,_,info=encoder.encode([sequence]); query_encode+=perf_counter()-t
                vector=pooled[0]; arrays[evaluation_id]=vector.numpy()
                t=perf_counter(); ranking=rank_cosine(vector,references,ids); ranking_time+=perf_counter()-t
                # Independent float64 cosine check; ties may differ under higher precision.
                q64=vector.numpy().astype(np.float64)
                independent=refs64@q64/(np.linalg.norm(refs64,axis=1)*np.linalg.norm(q64))
                observed={r['accession']:r['cosine_similarity'] for r in ranking}
                numerical_max=max(numerical_max,max(abs(observed[a]-float(independent[i])) for i,a in enumerate(ids)))
                rows=[]
                for a in expected:
                    cm=component_metrics(ranking,a,eq,clean.get(a))
                    if arm=='clean':
                        clean[a]=cm
                        if cm['rank']!=1 or abs(cm['cosine_similarity']-1)>1e-5:
                            clean_failures.append({'accession':a,**cm})
                    position='alone' if not spec['artificial'] else ('N' if spec['components'].index(a)==0 else 'C')
                    fraction=len(sequences[a])/len(combined)
                    row={'evaluation_id':evaluation_id,'query_id':qid,'pair_id':spec['pair_id'],
                        'arm':arm,'accession':a,'gene':plan['target_genes'][a],
                        'partner':next((b for b in spec['components'] if b!=a),'') if spec['artificial'] else '',
                        'position_in_parent':position,'component_length_aa':len(sequences[a]),
                        'parent_length_aa':len(combined),'fraction_in_parent':fraction,'fraction_bin':fraction_bin(fraction),
                        'catalogue_size':len(ids),**cm}
                    rows.append(row); result['component_rows'].append(row)
                result['evaluations'].append({'evaluation_id':evaluation_id,'query_id':qid,'arm':arm,
                    'pair_id':spec['pair_id'],'parent_components':spec['components'],'expected':expected,
                    'sequence_sha256':sha256(sequence.encode()),'length_aa':len(sequence),
                    'oracle_interval_aa':interval,'encoding':info,'ranking':ranking,
                    'retrieval':query_recall(ranking,expected,eq)})
                rank_rows.extend({'evaluation_id':evaluation_id,'query_id':qid,'arm':arm,'catalogue_size':len(ids),**r} for r in ranking)
            if (number+1)%50==0: print(f'Encoded {number+1}/{len(specs)} fixed queries',flush=True)
        result['clean_query_failures']=clean_failures
        result['timings_seconds'].update(query_construction=construction_time,query_encoding=query_encode,ranking=ranking_time)
        result['numerical_cosine_check']={'method':'Independent NumPy float64 from saved float32 vectors; all reference/query scores',
            'scores_checked':len(rank_rows),'max_abs_difference':numerical_max,'atol':1e-6,'passed':numerical_max<=1e-6}
        result['counts'].update(evaluations=len(result['evaluations']),component_instances=len(result['component_rows']),
            combined_component_instances=sum(r['arm']=='whole' for r in result['component_rows']),complete_ranking_rows=len(rank_rows))
        np.savez_compressed(out/'embeddings.npz',**arrays)
        csv_write(out/'rankings.csv',rank_rows); csv_write(out/'components.csv',result['component_rows'])
        write(out/'sequence_references.json',mappings)
        if numerical_max>1e-6: raise ValueError('Independent cosine calculation differs beyond tolerance')
        result['status']='passed'
    except Exception as exc:
        result['status']='failed'; result['failures'].append({'stage':stage,'error_type':type(exc).__name__})
        raise
    finally:
        result['timings_seconds']['total_process_work']=perf_counter()-started
        result.setdefault('memory',{})['final_peak_rss_mib']=memory_mib()
        write(out/'results.json',result)
    print(json.dumps({'status':result['status'],'counts':result['counts'],'timings':result['timings_seconds'],
                      'peak_rss_mib':result['memory']['final_peak_rss_mib'],'clean_failures':len(clean_failures)},indent=2))


def compare(first, second):
    import numpy as np
    ad,bd=BASE/safe_id(first),BASE/safe_id(second)
    a,b=read(ad/'results.json'),read(bd/'results.json')
    for key in ['freeze_sha256','code_hashes','model_file_sha256','environment','catalogue_size']:
        if a[key]!=b[key]: raise ValueError('Repeatability provenance mismatch')
    if a['status']!='passed' or b['status']!='passed': raise ValueError('Failed runs cannot pass repeatability')
    lookup={e['evaluation_id']:e for e in a['evaluations']}
    max_score=0.; count=0
    for e in b['evaluations']:
        old=lookup[e['evaluation_id']]
        if e['sequence_sha256']!=old['sequence_sha256']: raise ValueError('Input mismatch')
        if [r['accession'] for r in e['ranking']]!=[r['accession'] for r in old['ranking']]:
            raise ValueError('Ranks changed in repeat process')
        for x,y in zip(e['ranking'],old['ranking'],strict=True):
            max_score=max(max_score,abs(x['cosine_similarity']-y['cosine_similarity'])); count+=1
    aa,bb=np.load(ad/'embeddings.npz',allow_pickle=False),np.load(bd/'embeddings.npz',allow_pickle=False)
    max_vector=0.
    for key in bb.files:
        if aa[key].shape!=bb[key].shape: raise ValueError('Array shapes changed')
        max_vector=max(max_vector,float(np.abs(aa[key]-bb[key]).max()))
    if not count or max(max_vector,max_score)>1e-6: raise ValueError('Repeatability tolerance exceeded')
    r={'passed':True,'first':first,'second':second,'evaluations':len(b['evaluations']),
       'scores':count,'arrays':len(bb.files),'ranks_identical':True,
       'max_score_abs_difference':max_score,'max_vector_abs_difference':max_vector,'atol':1e-6,
       'second_reencoded_catalogue':b['cache']['forced_reencode']}
    write(BASE/f'comparison-{first}-{second}.json',r); print(json.dumps(r,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='command',required=True)
    r=sub.add_parser('run'); r.add_argument('--run-id',required=True); r.add_argument('--subset',action='store_true'); r.add_argument('--reencode-catalogue',action='store_true')
    c=sub.add_parser('compare'); c.add_argument('first'); c.add_argument('second')
    args=p.parse_args()
    if args.command=='run': run(args.run_id,args.subset,args.reencode_catalogue)
    else: compare(args.first,args.second)
