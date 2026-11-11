"""Fictional parser cases; scientific admission is never granted by these tests."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.bench_inputs import detect_record_format, parse_sequence_record
from iri_constructs.runtime import InputError, ROOT, digest_file
from iri_constructs.week2_pipeline import load_preparatory_record, run_preparatory_qa


@unittest.skipUnless(importlib.util.find_spec('Bio'), 'Optional Biopython parsing dependency')
class Week2PipelineTests(unittest.TestCase):
    def setUp(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        from Bio.SeqRecord import SeqRecord
        self.SeqIO, self.Seq, self.Feature, self.Location, self.Record = SeqIO, Seq, SeqFeature, SimpleLocation, SeqRecord
        (ROOT / 'results').mkdir(exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix='fictional-week2-', dir=ROOT / 'results')
        self.addCleanup(self.temporary.cleanup)
        self.base = Path(self.temporary.name)
        self.output_number = 0
        self.coding = 'ATGGCTGCTTAA'

    def genbank(self, name='fixture', sequence=None):
        r = self.Record(self.Seq(sequence or self.coding), id='FICTIONAL.1', name='FICTIONAL')
        r.annotations = {'molecule_type': 'DNA', 'topology': 'linear'}
        r.features = [self.Feature(self.Location(0, len(r), strand=1), type='CDS',
                                  qualifiers={'translation': ['MAA'], 'product': ['fictional test product']})]
        path = self.base / (name + '.gb')
        self.SeqIO.write(r, path, 'genbank')
        return path, r

    def entry(self, path, source_id='fixture', **extra):
        return {'source_id': source_id, 'path': path.name, 'format': 'auto',
                'sequence_id': 'FICTIONAL.1', 'checksum_sha256': digest_file(path),
                'archetype': 'single_cds', 'source_url': 'https://example.invalid/fictional',
                'source': 'fictional-test-only', **extra}

    def run_records(self, entries):
        path = self.base / 'inputs.json'
        path.write_text(json.dumps({'schema_version': 'week2-source-inputs-1', 'records': entries}))
        self.output_number += 1
        output = self.base / f'qa-{self.output_number}'
        return run_preparatory_qa(path, output), output

    def test_supported_genbank_provenance_and_annotation_only_export(self):
        path, _ = self.genbank()
        report, output = self.run_records([self.entry(path, accession='FICTIONAL', accession_version='1')])
        self.assertEqual(report['prepared_records'], 1)
        self.assertFalse(report['scientific_admission'])
        self.assertFalse(report['scientific_scores_produced'])
        self.assertFalse(report['inference_executed'])
        self.assertEqual(report['exported_annotated_proteins'], 1)
        units = json.loads((output / 'unit-provenance.json').read_text())['units']
        self.assertEqual(units[0]['sequence_sha256'], hashlib.sha256(b'MAA').hexdigest())
        self.assertFalse(units[0]['scientific_admission'])
        self.assertEqual(units[0]['source_segments'], [{'start0': 0, 'end0': 12, 'strand': 1}])
        self.assertNotIn(self.coding, (output / 'report.json').read_text())
        self.assertNotIn('"sequence":', (output / 'report.json').read_text())

    def test_fasta_is_sequence_only_without_invented_features_or_proteins(self):
        path = self.base / 'fasta-with-wrong-extension.gb'
        path.write_text('>FICTIONAL.1\n' + self.coding + '\n')
        self.assertEqual(detect_record_format(path), 'fasta')
        report, output = self.run_records([self.entry(path)])
        self.assertEqual(report['prepared_records'], 1)
        self.assertEqual(report['exported_annotated_proteins'], 0)
        self.assertEqual(report['records'][0]['source_feature_count'], 0)
        self.assertEqual((output / 'units.fasta').read_text(), '')
        self.assertEqual(json.loads((output / 'normalized.json').read_text())[0]['construct']['annotation_status'], 'unknown')

    def test_content_format_disagreement_is_controlled(self):
        path = self.base / 'fixture.fa'
        path.write_text('>FICTIONAL.1\n' + self.coding + '\n')
        with self.assertRaisesRegex(InputError, 'Declared format'):
            parse_sequence_record(path, 'genbank', 'FICTIONAL.1', 'dna', digest_file(path))

    def test_integrity_and_exact_id_are_required(self):
        path, _ = self.genbank()
        for patch, reason in [({'checksum_sha256': 'f' * 64}, 'source_checksum_mismatch'),
                              ({'sequence_id': 'OTHER.1'}, 'identity_mismatch'),
                              ({'accession': 'FICTIONAL', 'accession_version': '2'}, 'identity_mismatch')]:
            report, _ = self.run_records([self.entry(path, **patch)])
            self.assertEqual(report['blocked_records'], 1)
            self.assertEqual(report['records'][0]['failure_class'], reason)

    def test_malformed_and_binary_inputs_do_not_leak_source_content(self):
        path = self.base / 'broken.gb'
        path.write_text('LOCUS this is a deliberately malformed local fixture\n' + self.coding)
        report, output = self.run_records([self.entry(path)])
        self.assertEqual(report['blocked_records'], 1)
        self.assertNotIn(self.coding, (output / 'report.json').read_text())
        path.write_bytes(b'\xff\xfe\x00')
        with self.assertRaisesRegex(InputError, 'UTF-8'):
            detect_record_format(path)

    def test_duplicate_source_id_rejected_and_aliases_visible(self):
        path, _ = self.genbank()
        with self.assertRaisesRegex(InputError, 'Duplicate preparation'):
            self.run_records([self.entry(path), self.entry(path)])
        report, _ = self.run_records([self.entry(path), self.entry(path, source_id='other-alias')])
        self.assertEqual(report['prepared_records'], 1)
        self.assertEqual(report['records'][1]['failure_class'], 'duplicate_identity')

    def test_identical_parent_sequence_is_a_link_not_independent_support(self):
        path, r = self.genbank()
        other = self.base / 'other.gb'
        r.id = 'OTHER.1'
        self.SeqIO.write(r, other, 'genbank')
        report, _ = self.run_records([self.entry(path), self.entry(other, source_id='other', sequence_id='OTHER.1')])
        self.assertEqual(report['prepared_records'], 2)
        self.assertEqual(report['duplicate_parent_sequence_groups'][0]['source_ids'], ['fixture', 'other'])
        self.assertEqual(len(report['duplicate_protein_sequence_groups'][0]['input_ids']), 2)

    def test_proposed_lineages_survive_preparation_without_claiming_independence(self):
        path, r = self.genbank()
        other = self.base / 'other.gb'
        r.id = 'OTHER.1'
        self.SeqIO.write(r, other, 'genbank')
        report, _ = self.run_records([
            self.entry(path, related_lineage_proposals=['shared-vector', 'shared-vector']),
            self.entry(other, source_id='other', sequence_id='OTHER.1',
                       related_lineage_proposals=['shared-vector'], lineage_group='legacy-link')])
        self.assertEqual(report['declared_lineage_groups'],
                         {'legacy-link': ['other'], 'shared-vector': ['fixture', 'other']})
        self.assertEqual(report['records'][0]['related_lineage_proposals'], ['shared-vector'])
        self.assertFalse(report['scientific_admission'])

    def test_deterministic_outputs_and_stable_ids(self):
        path, _ = self.genbank()
        report1, output1 = self.run_records([self.entry(path)])
        report2, output2 = self.run_records([self.entry(path)])
        self.assertEqual(report1['content_sha256'], report2['content_sha256'])
        for name in ('normalized.json', 'annotation-audit.json', 'unit-provenance.json', 'units.fasta'):
            self.assertEqual((output1 / name).read_bytes(), (output2 / name).read_bytes())
        with self.assertRaises(FileExistsError):
            run_preparatory_qa(self.base / 'inputs.json', output1)

    def test_circular_join_order_and_strand_preserved(self):
        from Bio.SeqFeature import CompoundLocation
        path, r = self.genbank()
        r.seq = self.Seq(self.coding[6:] + self.coding[:6])
        r.annotations['topology'] = 'circular'
        r.features[0].location = CompoundLocation([self.Location(6,12,strand=1), self.Location(0,6,strand=1)])
        self.SeqIO.write(r, path, 'genbank')
        report, output = self.run_records([self.entry(path)])
        self.assertEqual(report['exported_annotated_proteins'], 1)
        doc = json.loads((output / 'normalized.json').read_text())[0]
        self.assertEqual(doc['features'][0]['location']['start0'], None)
        self.assertEqual(doc['units'][0]['source_segments'], [{'start0': 6, 'end0': 12, 'strand': 1}, {'start0': 0, 'end0': 6, 'strand': 1}])

    def test_missing_translation_withheld_not_synthesized(self):
        path, r = self.genbank()
        del r.features[0].qualifiers['translation']
        self.SeqIO.write(r, path, 'genbank')
        report, _ = self.run_records([self.entry(path)])
        self.assertEqual(report['prepared_records'], 1)
        self.assertEqual(report['exported_annotated_proteins'], 0)
        self.assertEqual(report['records'][0]['annotation_states'], {'unknown': 1})

    def test_valid_translation_over_frozen_model_cap_is_not_model_eligible(self):
        path, r = self.genbank(sequence='ATG' + 'GCT' * 1022 + 'TAA')
        r.features[0].qualifiers['translation'] = ['M' + 'A' * 1022]
        self.SeqIO.write(r, path, 'genbank')
        report, output = self.run_records([self.entry(path)])
        self.assertEqual(report['exported_annotated_proteins'], 1)
        self.assertEqual(report['protein_interface_coverage'], {'blocked_length_cap': 1})
        units = json.loads((output / 'unit-provenance.json').read_text())
        self.assertEqual(units['model_interface']['max_residues'], 1022)
        self.assertEqual(units['units'][0]['selected_interface_eligibility'], 'blocked_length_cap')

    def test_igem_unknown_coordinates_and_missing_translation_preserved(self):
        path = self.base / 'part.json'
        path.write_text(json.dumps({'name': 'FICTIONAL.1', 'status': 'published', 'sequence': self.coding,
                                    'topology': {'label': 'Linear'}}))
        feature = {'uuid': 'fictional-feature-id', 'role': {'label': 'CDS', 'accession': 'SO:0000316'},
                   'label': 'fictional CDS', 'locations': [{'start': 1, 'end': 12}], 'strand': 'forward'}
        annotations = self.base / 'features.json'
        annotations.write_text(json.dumps({'data': [feature], 'page': 1, 'total': 1}))
        entry = self.entry(path, format='igem-json', annotations_path=annotations.name,
                           annotations_checksum_sha256=digest_file(annotations))
        report, output = self.run_records([entry])
        self.assertEqual(report['exported_annotated_proteins'], 0)
        doc = json.loads((output / 'normalized.json').read_text())[0]
        self.assertEqual(doc['features'][0]['annotation_status'], 'ambiguous')
        self.assertEqual(doc['features'][0]['segments'], [])
        audit = json.loads((output / 'annotation-audit.json').read_text())
        source_feature = next(iter(audit['fixture']['features'].values()))
        self.assertEqual(source_feature['source_feature_id'], 'fictional-feature-id')
        self.assertEqual(json.loads(source_feature['qualifiers']['registry_original_feature'][0]), feature)
        annotations.write_text(json.dumps({'data': [feature], 'page': 1, 'total': 2}))
        entry['annotations_checksum_sha256'] = digest_file(annotations)
        report, _ = self.run_records([entry])
        self.assertEqual(report['records'][0]['failure_class'], 'incomplete_annotations')
        for page, total in [(2, 1), (True, 1), (1, True)]:
            annotations.write_text(json.dumps({'data': [feature], 'page': page, 'total': total}))
            entry['annotations_checksum_sha256'] = digest_file(annotations)
            report, _ = self.run_records([entry])
            self.assertEqual(report['records'][0]['failure_class'], 'incomplete_annotations')

    def test_igem_verified_single_location_still_requires_source_translation(self):
        path = self.base / 'part.json'
        path.write_text(json.dumps({'name': 'FICTIONAL.1', 'status': 'published', 'sequence': self.coding}))
        annotations = self.base / 'features.json'
        annotations.write_text(json.dumps({'data': [{'uuid': 'f', 'role': {'label': 'CDS'},
                                                     'locations': [{'start': 1, 'end': 12}], 'strand': 'reverse'}], 'page': 1, 'total': 1}))
        entry = self.entry(path, format='igem-json', annotations_path=annotations.name,
                           annotations_checksum_sha256=digest_file(annotations),
                           coordinate_convention='one_based_inclusive', coordinate_evidence='https://example.invalid/fixture-only')
        record, _ = load_preparatory_record(entry, self.base)
        self.assertEqual(int(record.features[0].location.start), 0)
        self.assertEqual(record.features[0].location.strand, -1)
        report, _ = self.run_records([entry])
        self.assertEqual(report['exported_annotated_proteins'], 0)
        self.assertEqual(report['records'][0]['annotation_states'], {'unknown': 1})


if __name__ == '__main__':
    unittest.main()
