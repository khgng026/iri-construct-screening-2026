#!/usr/bin/env python3
"""Acquire only frozen official allowlists within the session's byte/disk limits."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('HF_HUB_DISABLE_TELEMETRY','1')
os.environ.setdefault('HF_HUB_DOWNLOAD_TIMEOUT','60')
from iri_constructs.session_store import append_event,atomic_json,file_hash


def main(args):
    from huggingface_hub import hf_hub_download,try_to_load_from_cache
    session=args.session.resolve()
    admission=json.loads((session/'audit/model_admission.json').read_text())
    model=next(m for m in admission['models'] if m['id']==args.model)
    if model['gated'] or model['private'] or model['disabled']:
        append_event(session/'job_ledger.jsonl',{'event':'download','state':'blocked','model':args.model,'reason':'official anonymous access unavailable'})
        return 77
    prior=session/'download_ledger.jsonl'
    events=[json.loads(x) for x in prior.read_text().splitlines()] if prior.exists() else []
    # Failed full-transfer reservations are conservatively counted; successful
    # files count only their actual newly cached payload bytes.
    consumed=sum(e.get('budget_bytes',0) for e in events)
    files={f['rfilename']:f for f in model['files']}
    new=[]
    for name in model['download_allowlist']:
        cached=try_to_load_from_cache(args.model,name,revision=model['sha'])
        if not isinstance(cached,str):new.append(name)
    proposed=sum(files[name]['size'] for name in new)
    if consumed+proposed>12000000000:raise RuntimeError('Authorized aggregate download limit would be exceeded')
    free=shutil.disk_usage(session).free
    if free-2*proposed<25*2**30:raise RuntimeError('Insufficient disk reserve plus conservative temporary headroom')
    records={}
    for name in model['download_allowlist']:
        item=files[name];cached=try_to_load_from_cache(args.model,name,revision=model['sha'])
        already=isinstance(cached,str)
        path=None
        for attempt in range(3):
            if attempt and consumed+item['size']>12000000000:raise RuntimeError('Retry would exceed conservative download budget')
            try:
                path=Path(hf_hub_download(args.model,name,revision=model['sha'],token=False))
                if path.stat().st_size!=item['size']:raise ValueError('Official file size mismatch')
                actual=file_hash(path);expected=item.get('lfs',{}).get('sha256')
                if expected and actual!=expected:raise ValueError('Official LFS SHA256 mismatch')
                charged=0 if already else path.stat().st_size
                append_event(prior,{'event':'file_complete','model':args.model,'revision':model['sha'],'file':name,
                                   'cached_before':already,'actual_new_payload_bytes':charged,'budget_bytes':charged,
                                   'sha256':actual,'attempt':attempt+1})
                consumed+=charged
                records[name]={'path':str(path),'bytes':path.stat().st_size,'sha256':actual,'new_payload_bytes':charged}
                break
            except Exception as exc:
                # A failed transfer's precise wire bytes are unavailable from Hub;
                # reserve the whole payload, rather than undercounting retries.
                charged=0 if already else item['size'];consumed+=charged
                append_event(prior,{'event':'file_failed','model':args.model,'revision':model['sha'],'file':name,
                                   'budget_bytes':charged,'actual_new_payload_bytes':None,'exception':f'{type(exc).__name__}: {exc}',
                                   'attempt':attempt+1,'accounting':'full payload reserved conservatively'})
                status=getattr(getattr(exc,'response',None),'status_code',None)
                if status in (401,403) or isinstance(exc,ValueError) or attempt==2:raise
                time.sleep(2*(attempt+1))
        if path is None:raise RuntimeError('No verified downloaded path')
    label=args.model.replace('/','--')
    atomic_json(session/'inputs/model_files'/f'{label}.json',{'model_id':args.model,'revision':model['sha'],'files':records,
            'aggregate_conservative_budget_bytes':consumed})
    print(json.dumps({'model_id':args.model,'status':'downloaded_verified','aggregate_budget_bytes':consumed}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--session',type=Path,required=True);p.add_argument('--model',required=True)
    sys.exit(main(p.parse_args()) or 0)
