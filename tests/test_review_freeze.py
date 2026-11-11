"""Evidence integrity checks; no fictional record is a scientific clearance."""
import copy
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.review_freeze import build_review_freeze, verify_review_freeze
from iri_constructs.runtime import InputError


class ReviewFreezeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / 'inputs.json').write_text('{"role": "proposed"}\n')
        (self.root / 'protocol.json').write_text('{"k": 5}\n')

    def test_order_independent_binding_never_authorizes_science(self):
        a = build_review_freeze(self.root, ['inputs.json', 'protocol.json'])
        b = build_review_freeze(self.root, ['protocol.json', 'inputs.json'])
        self.assertEqual(a, b)
        result = verify_review_freeze(self.root, a)
        self.assertEqual(result['artifact_count'], 2)
        self.assertFalse(result['scientific_execution_authorized'])

    def test_label_or_protocol_edit_invalidates_review_binding(self):
        for name in ['inputs.json', 'protocol.json']:
            with self.subTest(name=name):
                freeze = build_review_freeze(self.root, ['inputs.json', 'protocol.json'])
                path = self.root / name
                before = path.read_bytes()
                path.write_text('revised evidence')
                with self.assertRaisesRegex(InputError, 'evidence bytes changed'):
                    verify_review_freeze(self.root, freeze)
                path.write_bytes(before)

    def test_missing_source_refuses_rerun(self):
        freeze = build_review_freeze(self.root, ['inputs.json'])
        (self.root / 'inputs.json').unlink()
        with self.assertRaisesRegex(InputError, 'absent'):
            verify_review_freeze(self.root, freeze)

    def test_tampered_freeze_cannot_rename_evidence_or_claim_clearance(self):
        freeze = build_review_freeze(self.root, ['inputs.json'])
        changed = copy.deepcopy(freeze)
        changed['artifacts']['inputs.json']['sha256'] = '0' * 64
        with self.assertRaisesRegex(InputError, 'freeze content changed'):
            verify_review_freeze(self.root, changed)
        changed = copy.deepcopy(freeze)
        changed['scientific_clearance'] = 'CLEAR_FOR_E1'
        with self.assertRaisesRegex(InputError, 'cannot grant'):
            verify_review_freeze(self.root, changed)

    def test_empty_duplicates_and_escape_rejected(self):
        for paths in [[], ['inputs.json', './inputs.json'], ['../outside.json']]:
            with self.subTest(paths=paths), self.assertRaises(InputError):
                build_review_freeze(self.root, paths)


if __name__ == '__main__':
    unittest.main()
