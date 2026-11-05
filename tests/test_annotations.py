"""Fictional records built in memory; never evidence for biological annotation accuracy."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.annotations import normalize_record
from iri_constructs.runtime import InputError
from iri_constructs.schema_check import check_schema_example

REVIEW = dict(record_id='fictional-test', source='fictional-test', source_url='https://example.invalid',
              accession='FICTIONAL', accession_version='1', archetype='single_cds', notes='Fictional unit-test metadata only')


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython parsing tests; install requirements-parsing.txt')
class AnnotationTests(unittest.TestCase):
    def setUp(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        self.Seq, self.Record, self.Feature, self.Location = Seq, SeqRecord, SeqFeature, SimpleLocation
        self.coding = 'ATG' + 'GCT' * 2 + 'TAA'
        self.protein = 'M' + 'A' * 2

    def record(self, sequence=None, location=None, qualifiers=None):
        record = self.Record(self.Seq(sequence or self.coding), id='FICTIONAL.1')
        record.annotations = {'molecule_type': 'DNA', 'topology': 'linear'}
        record.features = [self.Feature(location or self.Location(0, len(record), strand=1), type='CDS',
            qualifiers=qualifiers if qualifiers is not None else {'translation': [self.protein], 'product': ['fictional test label']})]
        return record

    def normalize(self, record):
        doc, audit, extracted = normalize_record(record, REVIEW, 'f' * 64)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'doc.json'
            path.write_text(json.dumps(doc), encoding='utf-8')
            self.assertEqual(check_schema_example(path), [])
        return doc, audit, extracted

    def test_exact_coordinates_translation_and_qualifiers(self):
        record = self.record()
        doc, audit, extracted = self.normalize(record)
        self.assertEqual(list(extracted.values()), [self.protein])
        fid = doc['features'][0]['feature_id']
        self.assertEqual(audit['features'][fid]['qualifiers'], record.features[0].qualifiers)
        self.assertEqual(doc['features'][0]['location'], dict(start0=0, end0=len(record), strand=1))
        self.assertEqual(doc['units'][0]['model_input_status'], 'pending_model_selection')

    def test_reverse_strand(self):
        reverse = str(self.Seq(self.coding).reverse_complement())
        record = self.record(reverse, self.Location(0, len(reverse), strand=-1))
        doc, _, extracted = self.normalize(record)
        self.assertEqual(list(extracted.values()), [self.protein])
        self.assertEqual(doc['features'][0]['location']['strand'], -1)

    def test_circular_join_and_reverse_join(self):
        from Bio.SeqFeature import CompoundLocation
        for reverse in (False, True):
            # Physical sequence is rotated; supplied segment order restores the annotated CDS.
            rotated = self.coding[6:] + self.coding[:6]
            spans = [self.Location(6,12,strand=1), self.Location(0,6,strand=1)]
            if reverse:
                rotated = str(self.Seq(rotated).reverse_complement())
                spans = [self.Location(0,6,strand=-1), self.Location(6,12,strand=-1)]
            record = self.record(rotated, CompoundLocation(spans, operator='join'))
            record.annotations['topology'] = 'circular'
            doc, _, extracted = self.normalize(record)
            self.assertEqual(list(extracted.values()), [self.protein])
            self.assertEqual(doc['features'][0]['location']['start0'], None)
            self.assertEqual(len(doc['units'][0]['source_segments']), 2)

    def test_fuzzy_unknown_and_remote_locations_withhold(self):
        from Bio.SeqFeature import BeforePosition
        for location in [self.Location(BeforePosition(0), 12, strand=1),
                         self.Location(0, 12), self.Location(0,12,strand=1,ref='REMOTE')]:
            doc, audit, extracted = self.normalize(self.record(location=location))
            self.assertEqual(extracted, {})
            self.assertEqual(doc['features'][0]['annotation_status'], 'ambiguous')
            self.assertTrue(audit['features'])
        record = self.record()
        record.features[0].location = None
        self.assertEqual(self.normalize(record)[2], {})

    def test_missing_conflicting_and_exceptional_annotations(self):
        for qualifiers, status in [({}, 'unknown'), ({'translation': ['M' + 'C' * 2]}, 'ambiguous'),
                                  ({'translation': [self.protein], 'pseudo': ['']}, 'ambiguous'),
                                  ({'translation': [self.protein], 'transl_except': ['fictional']}, 'ambiguous'),
                                  ({'translation': [self.protein, self.protein]}, 'ambiguous')]:
            doc, _, extracted = self.normalize(self.record(qualifiers=qualifiers))
            self.assertEqual(extracted, {})
            self.assertEqual(doc['features'][0]['annotation_status'], status)

    def test_out_of_bounds_rejected(self):
        with self.assertRaises(InputError):
            self.normalize(self.record(location=self.Location(0,99,strand=1)))

    def test_codon_start_and_translation_table_are_honored(self):
        record = self.record(sequence='C' + self.coding,
                             qualifiers={'translation': [self.protein], 'codon_start': ['2']})
        self.assertEqual(list(self.normalize(record)[2].values()), [self.protein])
        record = self.record(sequence='ATG' + 'TGA' + 'TAA',
                             qualifiers={'translation': ['M' + 'W'], 'transl_table': ['4']})
        self.assertEqual(list(self.normalize(record)[2].values()), ['M' + 'W'])
        record.features[0].qualifiers['transl_table'] = ['1']
        self.assertEqual(self.normalize(record)[2], {})

    def test_fasta_has_no_invented_annotations(self):
        record = self.record()
        record.features = []
        doc, _, extracted = self.normalize(record)
        self.assertEqual(doc['construct']['annotation_status'], 'unknown')
        self.assertEqual(doc['features'], [])
        self.assertEqual(extracted, {})

    def test_stable_ids_across_feature_reordering(self):
        record = self.record()
        record.features.append(self.Feature(self.Location(0,3,strand=1), type='promoter', qualifiers={'note': ['fictional']}))
        first = self.normalize(record)[0]
        record.features.reverse()
        self.assertEqual(first, self.normalize(record)[0])

    def test_duplicate_features_require_review(self):
        record = self.record()
        record.features.append(record.features[0])
        with self.assertRaises(InputError):
            self.normalize(record)
