from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_dispatch import resolve_command, BlockedJob


class QueueDependencyTests(unittest.TestCase):
    def test_dependency_uses_verified_result_and_preserves_literal_arguments(self):
        job = {'command': ['python', '--source', {'verified_job_of': 'encode'}, 'a b;$x'], 'requires': ['encode']}
        prior = [{'job_id': 'encode', 'state': 'complete', 'outputs': {'verified_job': '/owned/jobs/exact-hash'}}]
        self.assertEqual(resolve_command(job, prior), ['python', '--source', '/owned/jobs/exact-hash', 'a b;$x'])

    def test_failed_or_missing_dependency_cannot_resolve(self):
        job = {'command': ['python', {'verified_job_of': 'encode'}], 'requires': ['encode']}
        for prior in ([], [{'job_id': 'encode', 'state': 'failed'}]):
            with self.assertRaises(BlockedJob):
                resolve_command(job, prior)

    def test_unverified_output_and_undeclared_dependency_rejected(self):
        job = {'command': ['python', {'verified_job_of': 'encode'}], 'requires': ['encode']}
        for outputs in ({}, {'verified_job': 'relative/path'}):
            with self.assertRaises(ValueError):
                resolve_command(job, [{'job_id': 'encode', 'state': 'complete', 'outputs': outputs}])
        with self.assertRaises(ValueError):
            resolve_command({**job, 'requires': []}, [])

    def test_malformed_objects_and_empty_command_rejected(self):
        for command in ([], [''], ['python', {'unexpected': 'encode'}], ['python', None]):
            with self.assertRaises(ValueError):
                resolve_command({'command': command}, [])


if __name__ == '__main__':
    unittest.main()
