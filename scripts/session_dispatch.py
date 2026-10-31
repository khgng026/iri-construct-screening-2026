#!/usr/bin/env python3
"""Sequential guarded job queue; optional failure never erases other work."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time
import signal

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.session_store import UnitStore,atomic_json,append_event,file_hash,canonical,digest
from iri_constructs.session_jobs import run_guarded,attempt_job,BlockedJob,remaining_session_seconds


def stability(session,job_id):
    spec=importlib.util.spec_from_file_location('resource_guard',ROOT.parent/'machine-readiness-2026-10-31/resource_guard.py')
    guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
    start=time.monotonic();samples=[]
    while True:
        sample=guard.snapshot();samples.append(sample)
        if guard.reasons(sample) or sample['disk_free_gib']<25:
            atomic_json(session/'logs'/(job_id+'-stability.json'),samples)
            raise BlockedJob('Large-model 60-second stability gate failed; resources unsuitable')
        if time.monotonic()-start>=60:break
        time.sleep(2)
    atomic_json(session/'logs'/(job_id+'-stability.json'),samples)


def run_queue(args):
    s=args.session.resolve();queue=json.loads(args.queue.read_text());launch=json.loads((s/'launch.json').read_text())
    results=[]
    for job in queue['jobs']:
        job_id=job['id'];command=job['command']
        append_event(s/'job_ledger.jsonl',{'event':'queue_job_planned','job_id':job_id,'state':'planned','command':command,'timeout':job['timeout']})
        def action():
            for required in job.get('requires',[]):
                prior=[r for r in results if r['job_id']==required]
                if not prior or prior[-1]['state']!='complete':
                    raise BlockedJob('Required earlier job did not complete: '+required)
            if job.get('stable_before'):stability(s,job_id)
            output=run_guarded(command,interpreter=sys.executable,guard=ROOT/'scripts/session_guard_entry.py',
                    timeout=job['timeout'],deadline_utc=launch['hard_stop_utc'],log_path=s/'logs'/(job_id+'-guard.log'))
            if 'completion_glob' in job:
                matches=list((s/'jobs').glob(job['completion_glob']))
                valid=[]
                for path in matches:
                    if not list(path.glob('completion-*.json')):continue
                    store=UnitStore(path,json.loads((path/'identity.json').read_text()),s/'job_ledger.jsonl')
                    completions=sorted(path.glob('completion-*.json'),key=lambda p:p.stat().st_mtime)
                    completion=json.loads(completions[-1].read_text())
                    if completion['identity_hash']!=store.identity_hash:raise ValueError('Completion identity mismatch')
                    for key,entry in completion['units'].items():
                        unit=path/entry['path'];store.verify(unit)
                        if file_hash(unit/'complete.json')!=entry['descriptor_sha256']:raise ValueError('Completion descriptor mismatch')
                    if len(completion['units'])!=completion['expected_count']:raise ValueError('Completion count mismatch')
                    valid.append(str(path))
                if len(valid)!=1:raise ValueError(f'Expected one validated completed output, found {len(valid)}')
                output['verified_job']=valid[0]
                if job.get('audit'):
                    audit_path=s/'audit'/('independent_'+job_id+'.json')
                    subprocess.run([sys.executable,str(ROOT/'scripts/session_independent_audit.py'),valid[0],'--output',str(audit_path)],check=True,timeout=120)
                    output['independent_audit']=str(audit_path)
            if 'required_file' in job:
                file=s/job['required_file']
                if not file.exists():raise ValueError('Required verified acquisition record absent')
                output['required_file_sha256']=file_hash(file)
            return output
        result=attempt_job(job_id,action,s/'job_ledger.jsonl',optional=True);results.append(result)
        atomic_json(s/'logs'/(args.queue.stem+'-status.json'),{'queue':str(args.queue),'results':results})
        print(json.dumps(result),flush=True)
    return 0


def main(args):
    session=args.session.resolve();launch=json.loads((session/'launch.json').read_text())
    remaining=remaining_session_seconds(launch['hard_stop_utc'],session/'inputs/deadline_anchor.json')
    if remaining<=0:raise BlockedJob('Session hard deadline reached')
    awake=subprocess.Popen(['/usr/bin/caffeinate','-i','-s','-t',str(max(1,int(remaining))),'-w',str(__import__('os').getpid())])
    def terminate(signum,frame):
        signal.signal(signal.SIGTERM,signal.SIG_IGN)
        signal.signal(signal.SIGINT,signal.SIG_IGN)
        raise SystemExit(128+signum)
    signal.signal(signal.SIGTERM,terminate);signal.signal(signal.SIGINT,terminate)
    try:
        return run_queue(args)
    finally:
        if awake.poll() is None:awake.terminate();awake.wait(timeout=5)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--session',type=Path,required=True);p.add_argument('--queue',type=Path,required=True)
    sys.exit(main(p.parse_args()))
