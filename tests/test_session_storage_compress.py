"""Native compression tests use only newly written, explicitly fictional files."""
from contextlib import redirect_stdout
import argparse
import copy
import io
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import session_storage_compress as storage

PREPARED_GUARD = ROOT / 'scripts/resource_guard.py'


class FixtureGuard:
    def __init__(self, *args, **kwargs):
        self.identity = {'path':str(PREPARED_GUARD), 'sha256':storage.sha(PREPARED_GUARD)}

    def check(self, **kwargs): pass


@unittest.skipUnless(sys.platform == 'darwin', 'Native macOS compression semantics required')
class StorageCompressionTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='iri-fictional-compression-')
        self.addCleanup(self.temporary.cleanup)
        self.session = Path(self.temporary.name).resolve()
        (self.session/'audit').mkdir(); (self.session/'jobs').mkdir()
        self.state = self.session/'audit'/'operational-state'
        storage.atomic(self.session/'launch.json', {'resource_guard':str(PREPARED_GUARD)})
        for name in storage.INITIAL_STORES:
            job = self.session/'jobs'/name; job.mkdir()
            identity = {'mode':'expanded_model_retrieval_v1','query_ids':['fiction1','fiction2','fiction3']}
            storage.atomic(job/'identity.json',identity)
            ih = storage.digest(storage.canonical(identity)); units = {}
            for key in sorted(storage.expected_units(identity)):
                unit = job/key; unit.mkdir()
                storage.atomic(unit/'record.json',{'fictional_only':True,'key':key,'text':'fictional example '*8192})
                os.chmod(unit/'record.json',0o640)
                os.utime(unit/'record.json',ns=(1700000000000000000,1700000000000000000))
                (unit/'arrays.npz').write_bytes(b'FICTIONAL opaque non-result bytes\n')
                files = {p.name:{'bytes':p.stat().st_size,'sha256':storage.sha(p)}
                         for p in (unit/'record.json',unit/'arrays.npz')}
                storage.atomic(unit/'complete.json',{'key':key,'identity_hash':ih,'files':files})
                units[key] = {'path':key,'descriptor_sha256':storage.sha(unit/'complete.json')}
            storage.atomic(job/'completion-fictional.json',{'identity_hash':ih,'expected_count':len(units),
                            'units':units,'completed_utc':'2026-01-01T00:00:00+00:00'})
        args = argparse.Namespace(session=self.session,state=self.state,timeout=60,approved_stores=None)
        with patch.object(storage,'Guard',FixtureGuard), patch.object(storage,'quiet_jobs'), redirect_stdout(io.StringIO()):
            storage.freeze(args)
        self.plan_sha = storage.sha(self.state/'plan.json')
        self.plan = storage.load_plan(self.state,self.plan_sha)
        self.entry = self.plan['files'][0]
        self.source = Path(self.entry['path'])
        self.guard = FixtureGuard()
        self.original_nonrecords = {str(p):storage.sha(p) for p in (self.session/'jobs').rglob('*')
                                   if p.is_file() and p.name != 'record.json'}

    def assert_nonrecords_unchanged(self):
        self.assertEqual(self.original_nonrecords,{name:storage.sha(name) for name in self.original_nonrecords})

    def tx(self):
        return storage.transactions(self.state,self.plan)[-1]

    def compress(self, **kwargs):
        return storage.compress_one(self.state,self.plan,self.entry,self.guard,**kwargs)

    def independent_report(self, *, succeeds=True):
        # A genuinely distinct process invokes the independent complete-unit
        # byte checker. Only machine resource sampling is replaced for fixtures.
        program = '''import sys
from pathlib import Path
sys.path.insert(0,sys.argv[1])
import session_storage_compress as s
class FixtureGuard:
    def check(self,**kwargs): pass
state=Path(sys.argv[2]); plan=s.load_plan(state,sys.argv[3])
s.verify(None,state,plan,FixtureGuard())
'''
        result = subprocess.run([sys.executable,'-c',program,str(ROOT/'scripts'),str(self.state),self.plan_sha],
                                capture_output=True,text=True,timeout=40)
        if not succeeds:
            self.assertNotEqual(result.returncode,0); return result
        self.assertEqual(result.returncode,0,result.stderr)
        value = json.loads(result.stdout.strip().splitlines()[-1])
        self.assertNotEqual(storage.read(value['report'])['pid'],os.getpid())
        return Path(value['report']),value['sha256']

    def finalize(self, report, checksum):
        args = argparse.Namespace(verification=report,verification_sha256=checksum)
        with patch.object(storage,'quiet_jobs'), redirect_stdout(io.StringIO()):
            storage.finalize(args,self.state,self.plan,self.guard)

    def test_native_copy_preserves_bytes_metadata_and_owned_backup_until_full_pass(self):
        self.assertEqual(self.compress(),'compressed')
        directory,intent,entry = self.tx()
        installed = storage.metadata(self.source)
        self.assertEqual(installed['flags'] & storage.UF_COMPRESSED,storage.UF_COMPRESSED)
        self.assertLess(installed['allocated_bytes'],entry['metadata']['allocated_bytes'])
        self.assertTrue((directory/'original.link').exists())
        storage.exact(self.source,entry); storage.stable(entry['metadata'],installed)
        self.assert_nonrecords_unchanged()
        with self.assertRaisesRegex(ValueError,'different process'):
            with patch.object(storage,'quiet_jobs'):
                storage.verify(None,self.state,self.plan,self.guard)
        report,checksum = self.independent_report()
        self.finalize(report,checksum)
        self.assertFalse((directory/'original.link').exists())
        self.assert_nonrecords_unchanged()
        self.assertEqual(self.compress(),'already_compressed')

    def test_fault_with_same_inode_backup_unlinks_backup_and_allows_retry(self):
        def fail(stage):
            if stage == 'backup_retained': raise RuntimeError('fictional stop before rename')
        with self.assertRaisesRegex(RuntimeError,'fictional stop'):
            self.compress(fault=fail)
        directory,intent,entry = self.tx()
        self.assertEqual(self.source.stat().st_ino,entry['metadata']['inode'])
        self.assertEqual(self.source.stat().st_nlink,1)
        self.assertFalse((directory/'original.link').exists())
        self.assertEqual(storage.terminal(directory),'rolled_back.json')
        self.assertEqual(self.compress(),'compressed')

    def test_fault_after_atomic_replace_restores_original_inode(self):
        def fail(stage):
            if stage == 'replaced': raise RuntimeError('fictional stop after rename')
        with self.assertRaises(RuntimeError): self.compress(fault=fail)
        directory,intent,entry = self.tx()
        self.assertEqual(self.source.stat().st_ino,entry['metadata']['inode'])
        self.assertEqual(self.source.stat().st_nlink,1)
        self.assertTrue((directory/'replaced-copy.link').exists())
        storage.exact(self.source,entry); self.assert_nonrecords_unchanged()

    def test_restart_recovers_an_unterminated_post_replace_intent(self):
        def fail(stage):
            if stage == 'replaced': raise RuntimeError('simulated process death')
        # Suppress exception cleanup only to construct the exact crash-window
        # disk state, then run ordinary restart reconciliation twice.
        with patch.object(storage,'rollback'), self.assertRaises(RuntimeError): self.compress(fault=fail)
        directory,intent,entry = self.tx()
        self.assertIsNone(storage.terminal(directory))
        with self.assertRaisesRegex(ValueError,'Unresolved intent'): storage.state_digest(self.state,self.plan)
        self.assertEqual(storage.recover_pending(self.state,self.plan),1)
        self.assertEqual(storage.recover_pending(self.state,self.plan),0)
        self.assertEqual(self.source.stat().st_ino,entry['metadata']['inode'])
        self.assertEqual(self.source.stat().st_nlink,1)
        storage.exact(self.source,entry)

    def test_restart_rejects_missing_backup_with_foreign_source_inode(self):
        def fail(stage):
            if stage == 'replaced': raise RuntimeError('simulated process death')
        with patch.object(storage,'rollback'), self.assertRaises(RuntimeError): self.compress(fault=fail)
        directory,intent,entry = self.tx()
        (directory/'original.link').unlink()
        with self.assertRaisesRegex(ValueError,'Missing rollback link'):
            storage.recover_pending(self.state,self.plan)
        self.assertIsNone(storage.terminal(directory)); self.assertTrue(self.source.exists())

    def test_interrupted_first_intent_is_retained_but_unknown_data_blocks(self):
        directory = self.state/'transactions'/(self.entry['id']+'-fictional'); directory.mkdir()
        partial = directory/'.intent.json.partial'; partial.write_bytes(b'{')
        self.assertEqual(storage.recover_pending(self.state,self.plan),0)
        self.assertTrue(partial.exists())
        (directory/'compressed.tmp').write_bytes(b'unknown')
        with self.assertRaisesRegex(ValueError,'unexpected data'):
            storage.recover_pending(self.state,self.plan)

    def test_external_hash_and_descriptor_path_binding_reject_modified_plans(self):
        with self.assertRaisesRegex(ValueError,'Root-reviewed plan'):
            storage.load_plan(self.state,'0'*64)
        changed = copy.deepcopy(self.plan)
        changed['files'][0]['path'] = str(self.source.with_name('arrays.npz'))
        storage.atomic(self.state/'plan.json',changed)
        with self.assertRaisesRegex(ValueError,'descriptor-bound record.json'):
            storage.load_plan(self.state,storage.sha(self.state/'plan.json'))

    def test_fixed_trial_and_allowlist_cannot_silently_expand(self):
        self.assertEqual(len(self.plan['trial_file_ids']),10)
        self.assertEqual([x['job_index'] for x in self.plan['files'][:3]],[0,1,2])
        changed = copy.deepcopy(self.plan); changed['trial_file_ids'].reverse()
        storage.atomic(self.state/'plan.json',changed)
        with self.assertRaisesRegex(ValueError,'First10'):
            storage.load_plan(self.state,storage.sha(self.state/'plan.json'))
        changed = copy.deepcopy(self.plan); changed['jobs'].append(copy.deepcopy(changed['jobs'][0]))
        storage.atomic(self.state/'plan.json',changed)
        with self.assertRaisesRegex(ValueError,'allowlist expanded'):
            storage.load_plan(self.state,storage.sha(self.state/'plan.json'))

    def test_admission_rejects_external_hardlinks_and_resourcefork(self):
        link = self.session/'external-link'; os.link(self.source,link)
        with self.assertRaisesRegex(ValueError,'Unexpected hard links'): self.compress()
        link.unlink()
        subprocess.run(['/usr/bin/xattr','-w','com.apple.ResourceFork','fictional',str(self.source)],check=True)
        with self.assertRaisesRegex(ValueError,'unexpected storage/resource-fork'): self.compress()
        self.assertEqual(storage.transactions(self.state,self.plan),[])

    def test_metadata_mismatch_is_preserved_as_failure_evidence(self):
        def changed_mode(source,target,guard,seconds):
            storage.native_copy(source,target,guard,seconds); os.chmod(target,0o600)
        with self.assertRaisesRegex(ValueError,'metadata changed'): self.compress(copier=changed_mode)
        directory,intent,entry = self.tx()
        self.assertTrue((directory/'compressed.tmp').exists())
        report,checksum = self.independent_report(); self.finalize(report,checksum)
        self.assertTrue((directory/'compressed.tmp').exists())
        self.assertIn('metadata_failure_evidence_retained',(self.state/'ledger.jsonl').read_text())

    def test_acl_changes_rejected_by_stable_metadata_check(self):
        original = storage.metadata(self.source); changed = copy.deepcopy(original)
        changed['acl'] = ['0: user:fictional deny delete']
        with self.assertRaisesRegex(ValueError,'metadata changed'): storage.stable(original,changed)

    def test_forced_low_aggregate_rss_cap_stops_before_transaction(self):
        guard = storage.Guard.__new__(storage.Guard)
        guard.started=time.monotonic(); guard.seconds=60; guard.state=self.state; guard.last=0
        with patch.object(storage,'MEMORY',1), self.assertRaisesRegex(ValueError,'aggregate worker/ditto memory'):
            storage.compress_one(self.state,self.plan,self.entry,guard)
        self.assertEqual(storage.transactions(self.state,self.plan),[])

    def test_metadata_cli_recover_starts_without_unsupported_data_rlimit(self):
        previous=signal.getsignal(signal.SIGTERM)
        self.addCleanup(signal.signal,signal.SIGTERM,previous)
        arguments=['session_storage_compress.py','recover','--state',str(self.state),'--plan-sha256',self.plan_sha]
        with patch.object(sys,'argv',arguments), patch.object(storage,'Guard',FixtureGuard), \
             patch.object(storage,'quiet_jobs'), redirect_stdout(io.StringIO()) as output:
            storage.main()
        self.assertEqual(json.loads(output.getvalue()),{'recovered':0})

    def test_separate_full_verifier_catches_untouched_array_corruption_and_retains_backup(self):
        self.compress(); directory,intent,entry = self.tx()
        arrays = self.source.with_name('arrays.npz'); original = arrays.read_bytes()
        arrays.write_bytes(b'changed fictional opaque file')
        result = self.independent_report(succeeds=False)
        self.assertIn('Unit bytes/hash mismatch',result.stderr)
        self.assertTrue((directory/'original.link').exists())
        arrays.write_bytes(original)
        report,checksum = self.independent_report(); self.finalize(report,checksum)
        self.assertFalse((directory/'original.link').exists())

    def test_wrong_verification_hash_or_partial_coverage_cannot_release_backups(self):
        self.compress(); directory,intent,entry = self.tx()
        report,checksum = self.independent_report()
        with self.assertRaisesRegex(ValueError,'Root-reviewed verification'):
            self.finalize(report,'0'*64)
        changed = storage.read(report); changed['verified_units'] -= 1
        altered = self.state/'fictional-incomplete-report.json'; storage.atomic(altered,changed)
        with self.assertRaisesRegex(ValueError,'omits store/files'):
            self.finalize(altered,storage.sha(altered))
        self.assertTrue((directory/'original.link').exists())

    def test_trial_full_gate_and_trial_selection(self):
        args = argparse.Namespace(quiescent=True,scope='full',per_file_seconds=10)
        with patch.object(storage,'quiet_jobs'), self.assertRaisesRegex(ValueError,'first10 trial'):
            storage.apply(args,self.state,self.plan,self.guard)
        args.scope='trial'; selected=[]
        def record_only(state,plan,entry,guard,**kwargs): selected.append(entry['id']); return 'fixture_only'
        with patch.object(storage,'quiet_jobs'), patch.object(storage,'compress_one',record_only), redirect_stdout(io.StringIO()):
            storage.apply(args,self.state,self.plan,self.guard)
        self.assertEqual(selected,self.plan['trial_file_ids'])
        self.assertEqual(len(selected),10)
        self.assertFalse((self.state/'trial_finalized.json').exists())


if __name__ == '__main__': unittest.main()
