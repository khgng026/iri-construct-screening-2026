import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from iri_constructs.session_store import UnitStore, snapshot_tree
from iri_constructs.session_jobs import attempt_job, BlockedJob


class StoreTests(unittest.TestCase):
    def test_completion_integrity_and_identity(self):
        with tempfile.TemporaryDirectory() as d:
            store = UnitStore(Path(d)/"job", {"model":"a", "revision":"fixed", "code":"x"})
            path = store.publish("q1", {"answer":1}, {"v": np.arange(3,dtype=np.float32)})
            with self.assertRaises(FileExistsError): store.publish("q1", {})
            with self.assertRaises(ValueError): store.complete(["q1","q2"])
            with self.assertRaises(ValueError): UnitStore(Path(d)/"job", {"model":"b"})
            self.assertEqual(store.complete(["q1"])["expected_count"],1)
            (path/"arrays.npz").write_bytes(b"corrupted")
            self.assertIsNone(store.find("q1"))
            new = store.publish("q1", {"answer":1}, {"v":np.arange(3)})
            self.assertNotEqual(new,path)
            self.assertTrue(path.exists())
            self.assertEqual(store.read("q1")["answer"],1)

    def test_incomplete_directory_and_snapshot(self):
        with tempfile.TemporaryDirectory() as d:
            store = UnitStore(Path(d)/"job", {"input":"a"})
            (store.root/"units/q/.staging-crash").mkdir(parents=True)
            self.assertIsNone(store.find("q"))
            store.publish("q", {"ok":True})
            index = snapshot_tree(store.root, Path(d)/"snapshot")
            self.assertTrue(index)
            self.assertFalse((Path(d)/"snapshot/units/q/.staging-crash").exists())

    def test_optional_failure_and_block_continue(self):
        with tempfile.TemporaryDirectory() as d:
            ledger=Path(d)/"ledger.jsonl"
            def failed(): raise ValueError("deliberate optional failure")
            def blocked(): raise BlockedJob("simulated denied access; credentials absent")
            results=[attempt_job(name,call,ledger) for name,call in
                     [("failure",failed),("denied",blocked),("next",lambda:{"validated":True})]]
            self.assertEqual([r["state"] for r in results],["failed","blocked","complete"])
            self.assertEqual(len(ledger.read_text().splitlines()),6)

    def test_supervisor_interruption_is_recorded_and_propagated(self):
        with tempfile.TemporaryDirectory() as d:
            ledger=Path(d)/'ledger.jsonl'
            def stop():raise SystemExit(143)
            with self.assertRaises(SystemExit):attempt_job('owned-stop',stop,ledger)
            events=[json.loads(line) for line in ledger.read_text().splitlines()]
            self.assertEqual(events[-1]['event'],'job_attempt_interrupted')
            self.assertEqual(events[-1]['state'],'failed')


if __name__ == "__main__": unittest.main()
