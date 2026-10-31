#!/usr/bin/env python3
"""Actual owned-child interruption and restart acceptance test (numeric toy data)."""
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.session_store import UnitStore, atomic_json, file_hash
from iri_constructs.session_jobs import attempt_job, BlockedJob


def child(directory, pause):
    import numpy as np
    store=UnitStore(directory,{'test':'deliberate-process-interruption-v1'})
    for i in range(4):
        key=f'q{i}'
        if store.find(key): continue
        if pause and i == 1:
            stage=store.root/'units/q1/.staging-owned-interruption'
            stage.mkdir(parents=True,exist_ok=True)
            atomic_json(stage/'record.json',{'incomplete':True})
            atomic_json(store.root/'interrupt_ready.json',{'pid':os.getpid()})
            time.sleep(30)
        store.publish(key,{'index':i},{'vector':np.array([i,i+1],dtype=np.float32)})
    store.complete([f'q{i}' for i in range(4)])


def main(directory):
    store=UnitStore(directory,{'test':'deliberate-process-interruption-v1'})
    command=[sys.executable,__file__,'--directory',str(directory),'--child']
    p=subprocess.Popen(command+['--pause'])
    limit=time.monotonic()+20
    while not (store.root/'interrupt_ready.json').exists():
        if p.poll() is not None or time.monotonic()>limit:
            if p.poll() is None: p.terminate(); p.wait(timeout=5)
            raise RuntimeError('Interruption readiness failed')
        time.sleep(.05)
    before=store.find('q0'); before_hash=file_hash(before/'complete.json')
    p.terminate(); p.wait(timeout=5)
    assert store.find('q1') is None
    subprocess.run(command,check=True,timeout=20)
    assert store.find('q0')==before and file_hash(before/'complete.json')==before_hash
    paths={f'q{i}':store.find(f'q{i}') for i in range(4)}
    assert all(path is not None for path in paths.values())
    corrupt=store.find('q2'); (corrupt/'arrays.npz').write_bytes(b'deliberate corruption')
    assert store.find('q2') is None
    subprocess.run(command,check=True,timeout=20)
    assert store.find('q2')!=corrupt and corrupt.exists()
    assert store.find('q0')==before
    def fail(): raise RuntimeError('deliberate optional failure')
    def deny(): raise BlockedJob('simulated access denied; missing credentials; no prompt required')
    queue=[attempt_job(n,f,store.root/'queue.jsonl') for n,f in
           [('optional-failure',fail),('access-denied',deny),('next-job',lambda:store.complete([f'q{i}' for i in range(4)]))]]
    assert [x['state'] for x in queue]==['failed','blocked','complete']
    result={'status':'passed','terminated_owned_pid':p.pid,'interrupted_returncode':p.returncode,
            'completed_record_preserved':True,'incomplete_record_rerun':True,'corrupt_cache_detected_and_retained':True,
            'unique_logical_records':4,'missing_records':0,'queue_states':[x['state'] for x in queue],
            'paths_after_recovery':{f'q{i}':str(store.find(f'q{i}')) for i in range(4)}}
    atomic_json(store.root/'acceptance.json',result)
    print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--child',action='store_true');p.add_argument('--pause',action='store_true');a=p.parse_args()
    child(a.directory,a.pause) if a.child else main(a.directory)
