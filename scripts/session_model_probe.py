#!/usr/bin/env python3
"""Bounded full-sequence model admission with persistent per-input evidence."""
import argparse
import json
import os
from pathlib import Path
import resource
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('HF_HUB_DISABLE_TELEMETRY','1')
from iri_constructs.session_store import UnitStore,canonical,digest,file_hash,atomic_json


def main(args):
    import numpy as np
    from iri_constructs.session_benchmark import load_suite
    from iri_constructs.session_models import load_encoder
    s=args.session.resolve();suite=load_suite(ROOT,s/'inputs');config=json.loads(args.config.read_text())
    clean=sorted([q for q in suite['queries'] if q['suite']=='clean'],key=lambda q:(q['length_aa'],q['query_id']))
    longest=max(suite['queries'],key=lambda q:(q['length_aa'],q['query_id']))
    panel=[clean[0],clean[len(clean)//2],clean[-1],longest]
    encoder=load_encoder(config)
    try:
        identity={'mode':'admission_probe_v1','config':config,
                  'model':{k:v for k,v in encoder.metadata.items() if k!='load_seconds'},
                  'probe_code_sha256':file_hash(__file__),'query_ids':[q['query_id'] for q in panel],
                  'query_hashes':{q['query_id']:q['joined_sequence_sha256'] for q in panel},
                  'repeats':3,'whole_residue_inputs_only':True,'warmup':'first forward retained but excluded from subsequent median'}
        job=UnitStore(s/'jobs'/(args.name+'-'+digest(canonical(identity))[:12]),identity,s/'job_ledger.jsonl')
        expected=[]
        for q in panel:
            key=q['query_id'];expected.append(key)
            if job.find(key):continue
            times=[];reference=None;max_difference=0
            for i in range(3):
                start=time.monotonic();states=encoder.encode(q['sequence']);times.append(time.monotonic()-start)
                if reference is None:reference=states.copy()
                else:max_difference=max(max_difference,float(np.max(np.abs(states-reference))))
            if max_difference>1e-6:raise ValueError('Within-configuration repeat exceeded CPU FP32 tolerance')
            info={'query_id':key,'length':len(q['sequence']),'times_seconds':times,
                  'warm_median_seconds':float(np.median(times[1:])),
                  'within_config_max_abs':max_difference,'encode':encoder.last_encode_metadata,
                  'peak_process_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,
                  'load_seconds':encoder.load_seconds}
            if config['family']=='esm3' and q==panel[0]:
                info['official_api_equivalence']=encoder.verify_official_equivalence(q['sequence'])
            job.publish(key,info,{'residues':reference,'mean':reference.mean(axis=0,dtype=np.float64).astype(np.float32)})
            print(json.dumps(info),flush=True)
        job.complete(expected)
        print(json.dumps({'status':'passed','job':str(job.root)}),flush=True)
    finally:encoder.close()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--session',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--name',required=True)
    main(p.parse_args())
