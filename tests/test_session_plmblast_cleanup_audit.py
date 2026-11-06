"""Bounded owned toy subprocess regressions; no upstream, JIT or model imports."""
import io
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch, Mock

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
import session_plmblast_experiment as runner
from iri_constructs.session_store import atomic_json, digest, file_hash
from iri_constructs.session_jobs import OwnedCleanupFailure


class OwnedChildCleanupAuditTests(unittest.TestCase):
    def setUp(self):
        evidence_root = os.environ.get('IRI_PLM_CLEANUP_EVIDENCE_DIR')
        self.temporary = None
        if evidence_root:
            Path(evidence_root).mkdir(parents=True, exist_ok=True)
            self.directory = Path(tempfile.mkdtemp(prefix=self._testMethodName + '-', dir=evidence_root))
        else:
            self.temporary = tempfile.TemporaryDirectory()
            self.directory = Path(self.temporary.name)
        self.handlers = {s: signal.getsignal(s) for s in (signal.SIGTERM, signal.SIGINT)}
        self.old_globals = runner.SPAWNING, runner.CLEANING, runner.PENDING_SIGNAL
        runner.SPAWNING = runner.CLEANING = False
        runner.PENDING_SIGNAL = None

    def tearDown(self):
        for sig, handler in self.handlers.items(): signal.signal(sig, handler)
        runner.SPAWNING, runner.CLEANING, runner.PENDING_SIGNAL = self.old_globals
        if self.temporary: self.temporary.cleanup()

    def test_leader_exit_does_not_leave_resistant_descendant(self):
        grandchild_code = ('import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); '
                           'print("ready",flush=True); time.sleep(20)')
        leader_code = ('import subprocess,sys,time,json; '
                       'child=subprocess.Popen([sys.executable,"-c",' + repr(grandchild_code) + '],stdout=subprocess.PIPE,text=True); '
                       'assert child.stdout.readline().strip()=="ready"; '
                       'print(json.dumps({"grandchild":child.pid}),flush=True); time.sleep(25)')
        process = runner.spawn_owned_child([sys.executable, '-c', leader_code], self.directory, dict(os.environ))
        try:
            ready = runner.child_message(process, 5)
            result = runner.stop_child(process)
            self.assertNotIn(process.pid, runner.OWNED_CHILDREN)
            self.assertTrue(result['group']['kill_sent'])
            self.assertTrue(result['group']['quiescent'])
            self.assertEqual(result['status'], 'passed')
            self.assertIsNotNone(process.returncode)
            self.assertTrue(process.stdin.closed and process.stdout.closed)
            state = subprocess.run(['/bin/ps', '-p', str(ready['grandchild']), '-o', 'stat='],
                                   capture_output=True, text=True, timeout=5).stdout.strip()
            self.assertTrue(not state or state[0] in 'ZX', state)
            atomic_json(self.directory / 'fixture-result.json', {'status': 'passed', 'grandchild_pid': ready['grandchild'],
                        'grandchild_final_status': state, 'cleanup': result, 'finite_toy_seconds': 25})
        finally:
            runner.stop_child(process)

    def test_pending_signal_after_spawn_is_delivered_after_ownership_and_cleanup(self):
        signal.signal(signal.SIGTERM, runner.interrupted)
        signal.signal(signal.SIGINT, runner.interrupted)
        original_popen = subprocess.Popen
        children = []
        def interrupt_before_return(*args, **kwargs):
            child = original_popen(*args, **kwargs)
            if not children:
                children.append(child)
                os.kill(os.getpid(), signal.SIGTERM)
                self.assertTrue(runner.SPAWNING)
                self.assertEqual(runner.PENDING_SIGNAL, signal.SIGTERM)
            return child
        with patch.object(runner.subprocess, 'Popen', side_effect=interrupt_before_return):
            with self.assertRaises(SystemExit) as caught:
                runner.spawn_owned_child([sys.executable, '-c', 'import time; time.sleep(20)'],
                                         self.directory, dict(os.environ))
        self.assertEqual(caught.exception.code, 143)
        self.assertEqual(len(children), 1)
        child = children[0]
        try:
            self.assertIsNotNone(child.returncode)
            self.assertNotIn(child.pid, runner.OWNED_CHILDREN)
            self.assertTrue(child.stdin.closed and child.stdout.closed)
            ownership = json.loads((self.directory / 'ownership.json').read_text())
            self.assertEqual(ownership['pid'], child.pid)
            cleanup = [json.loads(p.read_text()) for p in self.directory.glob('cleanup-*.json')]
            self.assertEqual(len(cleanup), 1)
            self.assertEqual(cleanup[0]['status'], 'passed')
            self.assertTrue(cleanup[0]['group']['quiescent'])
            atomic_json(self.directory / 'fixture-result.json', {'status': 'passed', 'signal_exit_code': caught.exception.code,
                        'ownership_recorded_before_unwind': True, 'cleanup': cleanup[0], 'finite_toy_seconds': 20})
        finally:
            runner.stop_child(child)

    def test_cleanup_failure_preserves_evidence_and_closes_streams(self):
        ownership = {'pid': 987654, 'pgid': 987654, 'start_new_session': True,
                     'directory': str(self.directory), 'token': 'test-only-no-real-process'}
        atomic_json(self.directory / 'ownership.json', ownership)
        process = SimpleNamespace(pid=987654, returncode=0, stdin=io.StringIO(), stdout=io.StringIO(),
                                  wait=Mock(return_value=0), _iri_ownership=ownership)
        with patch.object(runner.os, 'killpg'), patch.object(runner, 'cleanup_owned_group', side_effect=OwnedCleanupFailure('injected cleanup failure')):
            with self.assertRaises(OwnedCleanupFailure): runner.stop_child(process)
        self.assertTrue(process.stdin.closed and process.stdout.closed)
        evidence = [json.loads(p.read_text()) for p in self.directory.glob('cleanup-*.json')]
        self.assertEqual(len(evidence), 1)
        self.assertEqual(evidence[0]['status'], 'failed')
        self.assertEqual(evidence[0]['exception_type'], 'OwnedCleanupFailure')
        atomic_json(self.directory / 'fixture-result.json', {'status': 'passed', 'injected_failure_preserved': True,
                    'streams_closed': True, 'real_process_started': False})

    def test_registry_recovers_child_without_caller_assignment(self):
        ownership = {'pid': 987653, 'pgid': 987653, 'start_new_session': True,
                     'directory': str(self.directory), 'token': 'test-only-no-real-process'}
        atomic_json(self.directory / 'ownership.json', ownership)
        child = SimpleNamespace(pid=987653, returncode=0, stdin=io.StringIO(), stdout=io.StringIO(),
                                wait=Mock(return_value=0), _iri_ownership=ownership)
        runner.OWNED_CHILDREN[child.pid] = child
        try:
            with patch.object(runner.os, 'killpg'), patch.object(runner, 'cleanup_owned_group', return_value={'quiescent': True}) as cleanup:
                runner.stop_owned_children()
            cleanup.assert_called_once_with(987653, grace_seconds=2)
            self.assertNotIn(child.pid, runner.OWNED_CHILDREN)
            self.assertTrue(child.stdin.closed and child.stdout.closed)
        finally:
            runner.OWNED_CHILDREN.pop(child.pid, None)

    def test_saved_state_identity_rejects_same_length_substitution_or_failed_record(self):
        record = {'status': 'success', 'sequence_sha256': digest(b'0123'), 'length': 4}
        runner.validate_state_record(record, '0123')
        for changed, sequence in ((record, '0132'), ({**record, 'status': 'failed'}, '0123'),
                                  ({**record, 'length': 3}, '0123')):
            with self.subTest(record=changed, sequence=sequence), self.assertRaises(ValueError):
                runner.validate_state_record(changed, sequence)

    def test_existing_summary_preserves_and_revalidates_original_startup_snapshot(self):
        session = self.directory
        job = session / 'jobs/test'; job.mkdir(parents=True)
        store = SimpleNamespace(root=job, identity={'runtime': {'version': 'fixed'}})
        def startup(name, seconds):
            directory = session / 'tmp/plmblast-workers' / name; directory.mkdir(parents=True)
            ready = {'event': 'ready', 'runtime': store.identity['runtime'], 'fixtures': {'status': 'passed'},
                     'startup_seconds': seconds, 'peak_process_rss_mib': 10}
            path = directory / 'ready.json'; atomic_json(path, ready)
            record = {'ready': ready, 'ready_path': str(path), 'ready_sha256': file_hash(path)}
            atomic_json(job / ('startup-' + name + '.json'), record)
            return record
        first = startup('first', 1.)
        published = {'startup_evidence': [first]}
        second = startup('later-verification', 2.)
        self.assertEqual(runner.summary_startups(store, session, published, 3072), [first])
        self.assertEqual(runner.summary_startups(store, session, None, 3072), [first, second])
        Path(first['ready_path']).write_text('{}')
        with self.assertRaisesRegex(ValueError, 'changed'):
            runner.summary_startups(store, session, published, 3072)


if __name__ == '__main__': unittest.main()
