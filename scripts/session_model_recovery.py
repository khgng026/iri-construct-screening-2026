#!/usr/bin/env python3
"""Interrupt a real representative model run and resume with identical identity."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.session_store import UnitStore,atomic_json,file_hash,append_event
from iri_constructs.session_jobs import run_guarded


def main(args):
    s=args.session.resolve();name=args.name
    command=[sys.executable,str(ROOT/'scripts/session_experiment.py'),'--session',str(s),'--config',str(args.config.resolve()),'--name',name,'--subset','repeat']
    with (s/'logs'/(name+'-initial.log')).open('x') as log:
        p=subprocess.Popen([sys.executable,str(ROOT/'scripts/session_guard_entry.py'),'--timeout','300','--',*command],stdout=log,stderr=subprocess.STDOUT)
        awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-s','-t','300','-w',str(p.pid)])
        before={};job=None
        try:
            deadline=time.monotonic()+60
            while time.monotonic()<deadline:
                dirs=sorted((s/'jobs').glob(name+'-*'),key=lambda x:x.stat().st_mtime,reverse=True)
                if dirs:
                    job=dirs[0];paths=list((job/'units').glob('query-*/*/complete.json'))
                    if paths:
                        before={str(q):file_hash(q) for q in paths};break
                if p.poll() is not None:raise RuntimeError('Model trial exited before interruption; inspect log')
                time.sleep(.02)
            if not before:raise RuntimeError('No completed query before interruption deadline')
        finally:
            if p.poll() is None:p.terminate();p.wait(timeout=20)
            if awake.poll() is None:awake.terminate();awake.wait(timeout=5)
        atomic_json(s/'audit'/(name+'-before.json'),{'job':str(job),'guard_exit':p.returncode,'completed_descriptors':before})
        append_event(s/'job_ledger.jsonl',{'event':'deliberate_model_interruption','state':'failed','job':str(job),'guard_exit':p.returncode,'completed_query_records':len(before)})
        print('Interrupted owned model trial after',len(before),'completed queries',flush=True)
    result=run_guarded(command,interpreter=sys.executable,guard=ROOT/'scripts/session_guard_entry.py',timeout=300,
            deadline_utc=json.loads((s/'launch.json').read_text())['hard_stop_utc'],log_path=s/'logs'/(name+'-resumed.log'))
    store=UnitStore(job,json.loads((job/'identity.json').read_text()),s/'job_ledger.jsonl')
    assert all(file_hash(path)==value for path,value in before.items())
    expected=['query-'+q for q in store.identity['query_ids']]
    assert len(expected)==32 and all(store.find(k) for k in expected)
    result.update(status='passed',completed_query_records_preserved=len(before),expected_query_count=32,
                  missing_records=0,job_path=str(job),interrupted_guard_returncode=p.returncode)
    atomic_json(s/'audit'/(name+'-recovery.json'),result);print(json.dumps(result))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--session',type=Path,required=True);p.add_argument('--config',type=Path,required=True);p.add_argument('--name',required=True)
    main(p.parse_args())
