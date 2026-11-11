"""Fictional source and transport fixtures; never biological admission evidence."""
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import HTTPError
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.dataset_curation import (LICENSE_UUID, SourceClient, lineage_groups,
    parse_genbank, parse_igem, retry_after_seconds, safe_error, sha256, source_scope_hold, publish_review)
from iri_constructs.runtime import InputError


class Response(io.BytesIO):
    headers = {'Content-Type': 'application/json'}

    def geturl(self):
        return 'https://example.org/source'


class SourceClientTests(unittest.TestCase):
    def test_verified_cache_does_not_refetch(self):
        with tempfile.TemporaryDirectory() as tmp:
            calls = []
            def opener(*args, **kwargs):
                calls.append(1)
                return Response(b'original response')
            c = SourceClient(Path(tmp), opener=opener)
            data, metadata = c.get('https://example.org/source', 'fixture.json')
            self.assertEqual(sha256(data), metadata['sha256'])
            self.assertFalse(metadata['cache_hit'])
            cached, reused = c.get('https://example.org/source', 'fixture.json')
            self.assertEqual(data, cached)
            self.assertTrue(reused['cache_hit'])
            self.assertEqual(len(calls), 1)

    def test_corrupt_cache_is_not_silently_repaired(self):
        with tempfile.TemporaryDirectory() as tmp:
            c = SourceClient(Path(tmp), opener=lambda *a, **k: Response(b'original'))
            c.get('https://example.org/source', 'fixture.json')
            (Path(tmp) / 'fixture.json').write_bytes(b'tampered')
            with self.assertRaisesRegex(InputError, 'digest mismatch'):
                c.get('https://example.org/source', 'fixture.json')
            self.assertEqual((Path(tmp) / 'fixture.json').read_bytes(), b'tampered')

    def test_url_binding_and_incomplete_pairs(self):
        with tempfile.TemporaryDirectory() as tmp:
            c = SourceClient(Path(tmp), opener=lambda *a, **k: Response(b'original'))
            c.get('https://example.org/source', 'fixture.json')
            with self.assertRaises(InputError):
                c.get('https://example.org/changed', 'fixture.json')
            (Path(tmp) / 'fixture.json.retrieval.json').unlink()
            with self.assertRaisesRegex(InputError, 'Incomplete cache pair'):
                c.get('https://example.org/source', 'fixture.json')

    def test_offline_never_connects(self):
        with tempfile.TemporaryDirectory() as tmp:
            c = SourceClient(Path(tmp), offline=True,
                             opener=lambda *a, **k: self.fail('Offline mode attempted network'))
            with self.assertRaisesRegex(InputError, 'not present'):
                c.get('https://example.org/source', 'fixture.json')

    def test_retry_after_date_and_long_hold(self):
        now = datetime(2026, 11, 11, tzinfo=timezone.utc)
        self.assertEqual(retry_after_seconds('120'), 120)
        self.assertEqual(retry_after_seconds('Wed, 11 Nov 2026 00:02:00 GMT', now), 120)
        self.assertIsNone(retry_after_seconds('unknown'))
        with tempfile.TemporaryDirectory() as tmp:
            def opener(*args, **kwargs):
                raise HTTPError('https://example.org', 429, 'private response', {'Retry-After': '120'}, None)
            c = SourceClient(Path(tmp), opener=opener,
                             sleeper=lambda _: self.fail('Long provider wait should become resumable HOLD'))
            with self.assertRaisesRegex(InputError, '120 seconds'):
                c.get('https://example.org/source', 'fixture.json')

    def test_retry_honors_short_provider_wait(self):
        with tempfile.TemporaryDirectory() as tmp:
            calls, delays = [], []
            def opener(*args, **kwargs):
                calls.append(1)
                if len(calls) == 1:
                    raise HTTPError('https://example.org', 429, '', {'Retry-After': '2'}, None)
                return Response(b'ok')
            c = SourceClient(Path(tmp), opener=opener, sleeper=delays.append)
            self.assertEqual(c.get('https://example.org/source', 'fixture.json')[0], b'ok')
            self.assertIn(2, delays)

    def test_error_payload_and_path_not_exposed(self):
        self.assertEqual(safe_error(ValueError('SECRET_SEQUENCE')), 'ValueError')
        with tempfile.TemporaryDirectory() as tmp:
            c = SourceClient(Path(tmp))
            with self.assertRaises(InputError):
                c.get('https://example.org', '../escape')

    def test_published_source_free_manifest_resolves_same_original(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            output = root / 'results/curation'
            output.mkdir(parents=True)
            original = output / 'raw_sources/fixture.json'
            original.parent.mkdir()
            original.write_bytes(b'private original bytes')
            (output / 'source_review.json').write_text(json.dumps({'scientific_clearance_issued': False}))
            (output / 'inputs.preparatory.json').write_text(json.dumps({'records': [
                {'source_id': 'fiction', 'path': 'raw_sources/fixture.json',
                 'annotations_path': 'raw_sources/fixture.annotations.json'}]}))
            destination = root / 'data/manifests/week2'
            with patch('iri_constructs.dataset_curation.ROOT', root), patch('iri_constructs.dataset_curation.local_output', side_effect=lambda p: p):
                publish_review(output, destination)
            inputs = json.loads((destination / 'inputs.preparatory.json').read_text())
            self.assertEqual((destination / inputs['records'][0]['path']).resolve(), original)
            self.assertEqual({p.name for p in destination.iterdir()}, {'source_review.json', 'inputs.preparatory.json'})


class OriginalSourceTests(unittest.TestCase):
    def genbank(self):
        try:
            from Bio import SeqIO
            from Bio.Seq import Seq
            from Bio.SeqFeature import SeqFeature, SimpleLocation
            from Bio.SeqRecord import SeqRecord
        except ImportError:
            self.skipTest('Optional Biopython unavailable')
        record = SeqRecord(Seq('ATGGCTGCTTAA'), id='FICT00001.1', name='fictional',
                           description='Fictional software fixture')
        record.annotations['molecule_type'] = 'DNA'
        record.features = [SeqFeature(SimpleLocation(0, 12, strand=1), type='CDS',
                            qualifiers={'translation': ['MAA'], 'product': ['fictional fixture']})]
        stream = io.StringIO()
        SeqIO.write(record, stream, 'genbank')
        return stream.getvalue().encode()

    def test_exact_identity_and_translation_evidence(self):
        _, evidence = parse_genbank(self.genbank(), 'FICT00001.1')
        self.assertEqual(evidence['validated_source_translation_count'], 1)
        self.assertNotIn('translation', evidence['cds_evidence'][0])
        self.assertEqual(evidence['cds_evidence'][0]['translation_sha256'], sha256(b'MAA'))
        with self.assertRaisesRegex(InputError, 'accession.version mismatch'):
            parse_genbank(self.genbank(), 'FICT00001.2')

    def test_incomplete_original_rejected(self):
        with self.assertRaisesRegex(InputError, 'Incomplete'):
            parse_genbank(b'LOCUS incomplete\n', 'FICT00001.1')

    def part(self):
        return {'name': 'BBa_FICTION', 'status': 'published', 'sequence': 'ATGGCTGCTTAA',
                'licenseUUID': LICENSE_UUID, 'topology': {'label': 'Linear'},
                'source': 'Fictional source for software testing'}

    def features(self):
        return {'page': 1, 'total': 1, 'data': [{'uuid': 'fictional-feature',
                'role': {'accession': 'SO:0000316', 'label': 'CDS'},
                'label': 'fictional', 'locations': [{'start': 1, 'end': 12}], 'strand': 'forward'}]}

    def parse_part(self, part=None, features=None, license_data=None):
        return parse_igem(json.dumps(part or self.part()).encode(),
            json.dumps(features or self.features()).encode(), 'BBa_FICTION',
            license_data or {'uuid': LICENSE_UUID, 'spdxID': 'cc-by-4.0'})

    def test_source_cds_without_translation_not_synthesized(self):
        _, evidence = self.parse_part()
        self.assertEqual(evidence['validated_source_translation_count'], 0)
        self.assertEqual(evidence['source_cds_feature_count'], 1)
        self.assertEqual(evidence['cds_evidence'][0]['translation_state'], 'unknown')
        self.assertNotIn('sequence', evidence)

    def test_igem_pagination_and_license_not_assumed(self):
        features = self.features()
        features['total'] = 2
        with self.assertRaisesRegex(InputError, 'pagination incomplete'):
            self.parse_part(features=features)
        with self.assertRaisesRegex(InputError, 'license'):
            self.parse_part(license_data={'uuid': 'different', 'spdxID': 'cc-by-4.0'})
        for page, total in [(2, 1), (True, 1), (1, True)]:
            features = self.features()
            features.update(page=page, total=total)
            with self.assertRaisesRegex(InputError, 'pagination incomplete'):
                self.parse_part(features=features)

    def test_igem_identity_and_coordinates(self):
        part = self.part()
        part['name'] = 'BBa_OTHER'
        with self.assertRaisesRegex(InputError, 'identity'):
            self.parse_part(part=part)
        features = self.features()
        features['data'][0]['locations'][0]['end'] = 13
        with self.assertRaisesRegex(InputError, 'coordinates'):
            self.parse_part(features=features)

    def test_lineage_ids_are_not_independent_sample_counts(self):
        self.assertEqual(lineage_groups('BBa_25UIP3DW'), lineage_groups('BBa_25NPQIRU'))
        self.assertEqual(lineage_groups('unrelated-fictional'), [])

    def test_source_toxin_annotation_excluded_despite_vector_title(self):
        reason, basis = source_scope_hold('FICTIONAL.1', {'record_title': 'ordinary cloning vector',
            'cds_evidence': [{'product': 'ccdB', 'gene': None}]})
        self.assertEqual(reason, 'T02_source_annotated_toxin_component_excluded')
        self.assertIn('do not include', basis)
        self.assertEqual(source_scope_hold('FICTIONAL.1', {'cds_evidence': [{'product': 'reporter'}]}), (None, None))

    def test_publication_context_hold_does_not_claim_vector_hazard(self):
        for accession in ('KJ170897.1', 'PP098727.1'):
            reason, basis = source_scope_hold(accession, {'cds_evidence': []})
            self.assertEqual(reason, 'HOLD_context_review_source_publication_association')
            self.assertIn('review', basis)


if __name__ == '__main__':
    unittest.main()
