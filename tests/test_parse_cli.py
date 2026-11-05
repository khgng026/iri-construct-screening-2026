"""End-to-end parser demonstration on generated, explicitly fictional test records."""
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.runtime import digest_file
from iri_constructs.schema_check import check_schema_example
from iri_constructs.validate_manifest import REQUIRED_COLUMNS


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython parser CLI integration')
class ParseCliTests(unittest.TestCase):
    def test_fasta_and_genbank_from_unrelated_working_directory(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        from Bio.SeqRecord import SeqRecord
        with tempfile.TemporaryDirectory(dir=ROOT / 'results') as directory:
            base = Path(directory)
            record = SeqRecord(Seq('ATG' + 'GCT'*2 + 'TAA'), id='FICTIONAL.1', description='Fictional unit test only')
            record.annotations = dict(molecule_type='DNA')
            record.features = [SeqFeature(SimpleLocation(0, len(record), strand=1), type='CDS',
                qualifiers={'translation': ['M' + 'A'*2], 'product': ['fictional fixture']})]
            records, rows = [], []
            for fmt in ('fasta','genbank'):
                path = base / ('fictional.' + ('fasta' if fmt == 'fasta' else 'gb'))
                SeqIO.write(record, path, fmt)
                rows.append(dict(record_id=fmt, source='fictional-test', source_url='https://example.invalid/fixture',
                    accession='FICTIONAL', accession_version='1', archetype='single_cds', sequence_length_nt=str(len(record)),
                    feature_count='1' if fmt == 'genbank' else '0', annotation_status='approved',
                    review_status='approved', licensing_status='approved', safety_review='approved', retrieved_on='2026-11-05',
                    checksum_sha256=digest_file(path), notes='Fictional test gate values; no real source approved'))
                records.append(dict(input_id=fmt, record_id=fmt, path=path.name, format=fmt,
                                    molecule_type='dna', sequence_id=record.id))
            with (base / 'review.csv').open('w',newline='',encoding='utf-8') as handle:
                writer = csv.DictWriter(handle, fieldnames=REQUIRED_COLUMNS)
                writer.writeheader()
                writer.writerows(rows)
            config = base / 'inputs.json'
            config.write_text(json.dumps(dict(schema_version='t05-inputs-1', review_manifest='review.csv', records=records)), encoding='utf-8')
            output = base / 'run'
            proc = subprocess.run([sys.executable, str(ROOT/'scripts'/'parse_annotations.py'),
                '--inputs', str(config), '--output', str(output)], cwd=base, capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            report = json.loads((output/'report.json').read_text())
            self.assertEqual(report['constructs'], 2)
            self.assertEqual(report['extracted_annotated_proteins'], 1)
            docs = json.loads((output/'normalized.json').read_text())
            self.assertEqual(docs[0]['features'], [])
            self.assertEqual(docs[0]['construct']['annotation_status'], 'unknown')
            for index, doc in enumerate(docs):
                check = base / ('check' + str(index) + '.json')
                check.write_text(json.dumps(doc), encoding='utf-8')
                self.assertEqual(check_schema_example(check), [])
            audit = json.loads((output/'annotation-audit.json').read_text())
            self.assertEqual(len(audit['genbank']['features']), 1)
            # Metadata-only reports/stdout must not contain raw sequences.
            self.assertNotIn(str(record.seq), (output/'report.json').read_text() + proc.stdout + proc.stderr)
