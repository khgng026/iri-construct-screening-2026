"""Review gates use fictional metadata and temporary generated parser records only."""
import csv
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.bench_inputs import load_proteins, read_source, reviewed_sources
from iri_constructs.runtime import InputError, digest_file, local_output, ROOT
from iri_constructs.validate_manifest import REQUIRED_COLUMNS


class ReviewGateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.source = self.directory / 'fictional.fasta'
        self.source.write_text('>fictional-test\n' + 'A' * 4 + '\n', encoding='utf-8')
        self.row = dict(record_id='fictional-test', source='fictional-unit-test',
                        source_url='https://example.invalid/fictional-test', accession='FICTIONAL',
                        accession_version='1', archetype='single_cds', sequence_length_nt='12', feature_count='0',
                        annotation_status='approved', review_status='approved', licensing_status='approved',
                        safety_review='approved', retrieved_on='2026-11-05', checksum_sha256=digest_file(self.source),
                        notes='Fictional fixture, not an approved biological record')
        self.config = dict(schema_version='t05-inputs-1', review_manifest='review.csv',
                           records=[dict(input_id='p', record_id='fictional-test', path=self.source.name,
                                         format='fasta', molecule_type='protein', sequence_id='fictional-test')])

    def write(self):
        with (self.directory / 'review.csv').open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS)
            writer.writeheader()
            writer.writerow(self.row)
        path = self.directory / 'inputs.json'
        path.write_text(json.dumps(self.config), encoding='utf-8')
        return path

    def test_unreviewed_and_uncertain_records_rejected(self):
        for field in ('review_status', 'licensing_status', 'safety_review'):
            for status in ('pending','uncertain','rejected'):
                self.row[field] = status
                with self.subTest(field=field,status=status), self.assertRaises(InputError):
                    reviewed_sources(self.write())
            self.row[field] = 'approved'

    def test_changed_file_checksum_is_rejected(self):
        self.source.write_text('>fictional-test\n' + 'C' * 4 + '\n', encoding='utf-8')
        with self.assertRaisesRegex(InputError, 'checksum'):
            reviewed_sources(self.write())

    def test_duplicate_input_id_and_blank_audit_notes_rejected(self):
        self.config['records'] *= 2
        with self.assertRaisesRegex(InputError, 'Duplicate'):
            reviewed_sources(self.write())
        self.config['records'] = self.config['records'][:1]
        self.row['notes'] = ''
        with self.assertRaisesRegex(InputError, 'notes'):
            reviewed_sources(self.write())

    def test_output_confinement_checks_resolved_path(self):
        self.assertEqual(local_output(ROOT / 'results' / 'test'), ROOT / 'results' / 'test')
        for path in (ROOT / 'results', ROOT / 'data' / 'examples', ROOT / 'results' / '..' / 'not-ignored'):
            with self.assertRaises(InputError):
                local_output(path)

    @unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython')
    def test_fasta_id_binding_and_no_concatenation(self):
        self.source.write_text('>fictional-test\n' + 'A'*4 + '\n>other\n' + 'C'*4 + '\n', encoding='utf-8')
        self.row['checksum_sha256'] = digest_file(self.source)
        config, proteins, provenance, excluded = load_proteins(self.write())
        self.assertEqual(len(proteins['p'].sequence), 4)
        self.assertEqual(proteins['p'].construct_id, self.row['record_id'])
        self.assertEqual(provenance['sources'][0]['accession'], 'FICTIONAL')
        self.config['records'][0]['sequence_id'] = 'missing'
        with self.assertRaisesRegex(InputError, 'exactly one'):
            load_proteins(self.write())

    @unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython')
    def test_annotation_review_not_implicitly_approved(self):
        self.row['annotation_status'] = 'pending'
        _, proteins, _, excluded = load_proteins(self.write())
        self.assertEqual(proteins, {})
        self.assertEqual(len(excluded), 1)

    @unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython')
    def test_duplicate_fasta_ids_rejected(self):
        self.source.write_text(('>fictional-test\n' + 'A'*4 + '\n') * 2, encoding='utf-8')
        self.row['checksum_sha256'] = digest_file(self.source)
        _, sources, _ = reviewed_sources(self.write())
        with self.assertRaisesRegex(InputError, 'Duplicate sequence'):
            read_source(sources[0])
