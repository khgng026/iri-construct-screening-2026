#!/usr/bin/env python3
"""Signal-safe entry to the unchanged prepared hardware guard.

Turning external termination into SystemExit lets the guard's existing finally
terminate ONLY its owned child process group before this wrapper exits.
"""
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.session_jobs import cleanup_owned_group
from iri_constructs.session_store import atomic_json,utcnow

GUARD=Path('/Users/Rudra2.0/Desktop/IRI/machine-readiness-2026-10-31/resource_guard.py')
SPAWNING=False
PENDING_SIGNAL=None


def terminate(signum, frame):
    global PENDING_SIGNAL
    # Ignore subsequent termination while the guard performs bounded cleanup.
    signal.signal(signal.SIGTERM,signal.SIG_IGN)
    signal.signal(signal.SIGINT,signal.SIG_IGN)
    if SPAWNING:
        PENDING_SIGNAL=signum
        return
    raise SystemExit(128+signum)


def main():
    global SPAWNING
    signal.signal(signal.SIGTERM,terminate)
    signal.signal(signal.SIGINT,terminate)
    ownership_path=os.environ.get('IRI_GUARD_OWNERSHIP_PATH')
    record={'wrapper_pid':os.getpid(),'token':os.environ.get('IRI_GUARD_OWNERSHIP_TOKEN'),
            'guard_source_sha256':hashlib.sha256(GUARD.read_bytes()).hexdigest(),
            'child_groups':[],'cleanup':[],'started_utc':utcnow()}
    def publish():
        if ownership_path:atomic_json(Path(ownership_path),record)
    class GuardSubprocess:
        """Observe only direct Popen calls made by the unchanged guard module."""
        def __getattr__(self,name):return getattr(subprocess,name)
        def Popen(self,*args,**kwargs):
            global SPAWNING
            if kwargs.get('start_new_session') is not True:
                raise RuntimeError('Prepared guard worker must own a new process group')
            SPAWNING=True
            try:
                child=subprocess.Popen(*args,**kwargs)
                record['child_groups'].append({'pid':child.pid,'pgid':child.pid,
                    'start_new_session':True,'spawned_utc':utcnow(),
                    'command_sha256':hashlib.sha256(json.dumps(args[0],sort_keys=True).encode()).hexdigest()})
                publish()
                return child
            finally:
                SPAWNING=False
                if PENDING_SIGNAL is not None:raise SystemExit(128+PENDING_SIGNAL)
    publish()
    spec=importlib.util.spec_from_file_location('prepared_resource_guard',GUARD)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.subprocess=GuardSubprocess()
    try:
        return module.main()
    finally:
        try:
            for entry in record['child_groups']:
                record['cleanup'].append(cleanup_owned_group(entry['pgid']))
        except BaseException as exc:
            record['cleanup_error']=f'{type(exc).__name__}: {exc}'
            raise
        finally:
            record['finished_utc']=utcnow()
            publish()
            print(json.dumps({'owned_group_cleanup':record['cleanup']}),flush=True)


if __name__=='__main__':
    sys.exit(main())
