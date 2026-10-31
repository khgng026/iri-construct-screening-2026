#!/usr/bin/env python3
"""Signal-safe entry to the unchanged prepared hardware guard.

Turning external termination into SystemExit lets the guard's existing finally
terminate ONLY its owned child process group before this wrapper exits.
"""
import importlib.util
from pathlib import Path
import signal
import sys

GUARD=Path('/Users/Rudra2.0/Desktop/IRI/machine-readiness-2026-10-31/resource_guard.py')


def terminate(signum, frame):
    # Ignore subsequent termination while the guard performs bounded cleanup.
    signal.signal(signal.SIGTERM,signal.SIG_IGN)
    signal.signal(signal.SIGINT,signal.SIG_IGN)
    raise SystemExit(128+signum)


if __name__=='__main__':
    signal.signal(signal.SIGTERM,terminate)
    signal.signal(signal.SIGINT,terminate)
    spec=importlib.util.spec_from_file_location('prepared_resource_guard',GUARD)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    sys.exit(module.main())
