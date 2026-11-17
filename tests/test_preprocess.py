"""Public toolkit API/CLI behavior on fictional fixtures (software checks, not biology).

Runs against the checkout by default. With ``IRI_TEST_INSTALLED=1`` it refuses to
touch ``src/`` and asserts the imported package is an installed copy, so a
packaging failure cannot hide behind the source tree.
"""
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
INSTALLED = os.environ.get('IRI_TEST_INSTALLED') == '1'
if not INSTALLED:
    sys.path.insert(0, str(ROOT / 'src'))
FIXTURES = Path(os.environ.get('IRI_FIXTURES', ROOT / 'examples' / 'synthetic'))
HAVE_BIO = importlib.util.find_spec('Bio') is not None


def tsv(path):
    with Path(path).open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t', quoting=csv.QUOTE_NONE))


def cli(*args, cwd=None):
    env = dict(os.environ)
    if not INSTALLED:
        env['PYTHONPATH'] = str(ROOT / 'src')
    return subprocess.run([sys.executable, '-m', 'iri_constructs.cli', *map(str, args)],
                          capture_output=True, text=True, cwd=cwd, env=env)


@unittest.skipUnless(HAVE_BIO, 'Toolkit preprocessing needs Biopython (declared package dependency)')
class PreprocessTests(unittest.TestCase):
    def setUp(self):
        import iri_constructs
        self.api = iri_constructs
        if INSTALLED:
            self.assertNotIn(str(ROOT / 'src'), str(Path(iri_constructs.__file__).resolve()))
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = Path(self._tmp.name)
        self.addCleanup(self._tmp.cleanup)
        self.expected = json.loads((FIXTURES / 'EXPECTED.json').read_text())

    def run_api(self, *names, out='out', **options):
        return self.api.preprocess_files([FIXTURES / n for n in names], self.tmp / out,
                                         options=self.api.PreprocessOptions(**options))

    # -- the bundle and its traceability -------------------------------------------------
    def test_single_file_bundle_is_complete_and_traceable(self):
        result = self.run_api('single_cds.gb')
        self.assertEqual(result.exit_code, 0)
        out = result.output_dir
        self.assertEqual(sorted(p.name for p in out.iterdir()),
                         sorted(['normalized.json', 'annotations.json', 'proteins.fasta', 'features.tsv',
                                 'units.tsv', 'report.json', 'run.json']))
        normalized = json.loads((out / 'normalized.json').read_text())
        construct = normalized['constructs'][0]
        self.assertEqual(construct['construct_id'], 'FIXSINGLE.1')
        self.assertEqual(construct['source_input']['sha256'],
                         hashlib.sha256((FIXTURES / 'single_cds.gb').read_bytes()).hexdigest())
        self.assertNotEqual(construct['sequence_sha256'], construct['source_input']['sha256'])
        self.assertEqual(construct['provenance'], {'supplied': None, 'missing_fields': ['source', 'source_url', 'retrieved_on', 'terms_note']})
        annotations = json.loads((out / 'annotations.json').read_text())['constructs']['FIXSINGLE.1']
        features, units = tsv(out / 'features.tsv'), tsv(out / 'units.tsv')
        self.assertEqual(len(features), len(units))
        self.assertEqual([f['source_feature_type'] for f in features], ['source', 'promoter', 'CDS', 'terminator'])
        self.assertEqual(features[0]['feature_label'], 'other')        # normalized label, original type kept beside it
        for feature in features:                                        # every row resolves to its original annotation
            original = annotations['features'][feature['feature_id']]
            self.assertEqual(original['type'], feature['source_feature_type'])
        cds = next(u for u in units if u['molecule_type'] == 'protein')
        header, *lines = (out / 'proteins.fasta').read_text().splitlines()
        self.assertEqual(header.split()[0], '>' + cds['unit_id'])
        protein = ''.join(lines)
        self.assertEqual(hashlib.sha256(protein.encode()).hexdigest(), cds['protein_sha256'])
        original = annotations['features'][cds['feature_id']]['qualifiers']
        self.assertEqual(original['translation'], [protein])            # exported protein is the source translation
        self.assertEqual(original['gene'], ['fixA'])
        self.assertEqual(result.proteins, {cds['unit_id']: protein})

    def test_edge_cases_get_declared_reason_codes_and_ordered_segments(self):
        result = self.run_api('edge_cases.gb')
        self.assertEqual(result.exit_code, 0)            # withholding alone never fails a batch
        features = {f['gene']: f for f in tsv(result.output_dir / 'features.tsv') if f['gene']}
        units = {u['feature_id']: u for u in tsv(result.output_dir / 'units.tsv')}
        observed = {gene: units[f['feature_id']]['reason_code'] for gene, f in features.items()}
        self.assertEqual(observed, self.expected['FIXEDGE.1'])
        state = {gene: units[f['feature_id']]['extraction_state'] for gene, f in features.items()}
        self.assertEqual({g for g, s in state.items() if s == 'exported'}, {'fwd', 'rev', 'alt_table', 'revjoin', 'wrap'})
        self.assertEqual({g for g, s in state.items() if s == 'withheld'},
                         {'notrans', 'mismatch', 'nostop', 'pseudo', 'noncanon', 'partial5'})
        self.assertEqual(features['rev']['strand'], '-')
        # Origin-spanning and reverse joins keep their source order; no bounding interval is invented.
        self.assertRegex(features['wrap']['segments'], r'^\[\d+:1563\]\(\+\);\[0:\d+\]\(\+\)$')
        first, second = features['revjoin']['segments'].split(';')
        self.assertGreater(int(first[1:].split(':')[0]), int(second[1:].split(':')[0]))
        # An uncertain boundary is kept as uncertain text, never as an exact segment.
        self.assertEqual(features['partial5']['segments'], '')
        self.assertIn('<', features['partial5']['source_location'])
        normalized = json.loads((result.output_dir / 'normalized.json').read_text())['constructs'][0]
        detail = {u['unit_id']: u['detail'] for u in normalized['units']}
        self.assertEqual(detail[features['nostop']['unit_id']],
                         {'cds_defects': ['no_terminal_stop_codon'], 'in_frame_translation_agrees': True})
        # Identical duplicate features both survive with distinct IDs.
        dup_ids = [f['feature_id'] for f in tsv(result.output_dir / 'features.tsv') if f['gene'] == 'dup']
        self.assertEqual(len(dup_ids), 2)
        self.assertEqual(len(set(dup_ids)), 2)
        self.assertTrue(dup_ids[1].endswith('-dup2'))
        report = json.loads((result.output_dir / 'report.json').read_text())
        self.assertEqual(report['summary']['proteins_exported'], 5)
        self.assertEqual(report['summary']['cds_withheld'], 6)
        self.assertIn('duplicate_feature', [w['code'] for w in report['records'][0]['warnings']])

    def test_long_protein_is_exported_and_only_flagged_for_the_model(self):
        result = self.run_api('long_protein.gb', model_max_residues=1022)
        unit = next(u for u in result.units if u.molecule_type == 'protein')
        self.assertEqual((unit.extraction_state, unit.length_aa, unit.model_compatibility),
                         ('exported', 1100, 'exceeds_max_residues'))
        self.assertEqual(len(''.join((result.output_dir / 'proteins.fasta').read_text().splitlines()[1:])), 1100)
        summary = json.loads((result.output_dir / 'report.json').read_text())['summary']['model_compatibility']
        self.assertEqual((summary['compatible'], summary['exceeds_max_residues']), (0, 1))
        unflagged = self.run_api('long_protein.gb', out='plain')
        self.assertEqual(next(u for u in unflagged.units if u.sequence).model_compatibility, 'not_evaluated')

    # -- batches, failures and identity -------------------------------------------------
    def test_partial_batch_keeps_successes_and_reports_failures(self):
        result = self.run_api('multi_record_one_bad.gb', 'not_a_sequence.gb', 'single_cds.gb', 'sequence_only.fasta')
        self.assertEqual(result.exit_code, 2)
        self.assertEqual([(r.input_name, r.record_index, r.status) for r in result.records],
                         [('multi_record_one_bad.gb', 1, 'ok'), ('multi_record_one_bad.gb', 2, 'failed'),
                          ('multi_record_one_bad.gb', 3, 'ok'), ('not_a_sequence.gb', 0, 'failed'),
                          ('single_cds.gb', 1, 'ok'), ('sequence_only.fasta', 1, 'ok')])
        self.assertEqual([r.failure.code for r in result.failed_records], ['malformed_record', 'unrecognized_format'])
        self.assertEqual(result.failed_records[0].source_record_id, 'FIXBROKEN')
        units = tsv(result.output_dir / 'units.tsv')
        self.assertEqual({u['construct_id'] for u in units}, {'FIXMULTIA.1', 'FIXMULTIC.1', 'FIXSINGLE.1'})
        self.assertEqual(sum(u['extraction_state'] == 'exported' for u in units), 4)
        fasta = json.loads((result.output_dir / 'normalized.json').read_text())['constructs'][-1]
        self.assertEqual((fasta['construct_id'], fasta['features'], fasta['units']), ('FIXFASTA', [], []))
        # A parser message must not carry annotation payload into the report.
        report = (result.output_dir / 'report.json').read_text()
        self.assertNotIn('garbage(', report)
        self.assertNotIn('\\n', json.dumps([r.failure.message for r in result.failed_records]))
        self.assertEqual(json.loads(report)['status'], 'partial')

    def test_no_usable_record_is_exit_one_with_a_written_report(self):
        result = self.run_api('not_a_sequence.gb')
        self.assertEqual(result.exit_code, 1)
        report = json.loads((result.output_dir / 'report.json').read_text())
        self.assertEqual((report['status'], report['summary']['failures_by_code']), ('no_usable_records', {'unrecognized_format': 1}))

    def test_record_level_failures_are_isolated_and_coded(self):
        good = (FIXTURES / 'single_cds.gb').read_text()
        cases = {
            'truncated.gb': (good.replace(good.split('ORIGIN')[1], '\n        1 acgtacgtac\n//\n'), 'declared_length_mismatch'),
            'no_sequence.gb': (good.split('ORIGIN')[0] + '//\n', 'malformed_record'),
            'protein.gb': (good.replace(' bp    DNA ', ' aa            ', 1), None),
            'empty.gb': ('', 'empty_input'),
        }
        for name, (text, _) in cases.items():
            (self.tmp / name).write_text(text)
        (self.tmp / 'latin1.gb').write_bytes(good.replace('fictional promoter', 'fictional promotér').encode('latin-1'))
        result = self.api.preprocess_files([self.tmp / n for n in [*cases, 'latin1.gb']], self.tmp / 'out')
        by_name = {r.input_name: r for r in result.records}
        for name, (_, code) in cases.items():
            if code:
                self.assertEqual((by_name[name].status, by_name[name].failure.code), ('failed', code), name)
        self.assertEqual(by_name['latin1.gb'].status, 'ok')
        self.assertIn('non_utf8_input', [w.code for w in by_name['latin1.gb'].warnings])
        self.assertEqual(result.exit_code, 2)

    def test_duplicate_and_colliding_records_never_overwrite(self):
        a = self.tmp / 'a'
        a.mkdir()
        shutil.copy(FIXTURES / 'single_cds.gb', a / 'copy.gb')
        text = (FIXTURES / 'single_cds.gb').read_text()
        changed = self.tmp / 'changed.gb'
        changed.write_text(text.replace('fictional promoter', 'edited note'))
        self.assertNotEqual(text, changed.read_text())
        result = self.api.preprocess_files([FIXTURES / 'single_cds.gb', a / 'copy.gb', changed], self.tmp / 'out')
        self.assertEqual([r.status for r in result.records], ['ok', 'skipped_duplicate', 'ok'])
        self.assertEqual(result.exit_code, 0)
        first, duplicate, other = result.records
        self.assertEqual(duplicate.duplicate_of, first.construct_id)
        self.assertNotEqual(first.construct_id, other.construct_id)
        for record in (first, other):                      # same identifier, different content: both disambiguated
            self.assertRegex(record.construct_id, r'^FIXSINGLE\.1~[0-9a-f]{10}$')
            self.assertIn('shared_record_identifier', [w.code for w in record.warnings])
        unit_ids = [u['unit_id'] for u in tsv(result.output_dir / 'units.tsv')]
        self.assertEqual(len(unit_ids), len(set(unit_ids)))
        self.assertEqual(len(unit_ids), 8)
        export = self.run_api('no_accession_export.gb', out='export')
        self.assertRegex(export.records[0].construct_id, r'^no_accession_export~[0-9a-f]{10}$')

    def test_outputs_are_deterministic_apart_from_run_metadata(self):
        names = ['edge_cases.gb', 'multi_record_one_bad.gb', 'long_protein.gb']
        one = self.run_api(*names, out='one', model_max_residues=1022)
        two = self.run_api(*names, out='two', model_max_residues=1022)
        for name in ('normalized.json', 'annotations.json', 'proteins.fasta', 'features.tsv', 'units.tsv', 'report.json'):
            self.assertEqual((one.output_dir / name).read_bytes(), (two.output_dir / name).read_bytes(), name)
        self.assertIn('started_at_utc', json.loads((one.output_dir / 'run.json').read_text()))
        self.assertNotIn(str(self.tmp), (one.output_dir / 'report.json').read_text())
        # A record's IDs do not depend on which other files share the batch.
        alone = self.run_api('edge_cases.gb', out='alone')
        self.assertEqual({u.unit_id for u in alone.units}, {u.unit_id for u in one.units if u.construct_id == 'FIXEDGE.1'})

    def test_one_normalizer_toolkit_agrees_with_historical_interface(self):
        from Bio import SeqIO
        from iri_constructs.annotations import normalize_record
        record = SeqIO.read(FIXTURES / 'single_cds.gb', 'genbank')
        review = dict(record_id='FIXSINGLE.1', source='fictional', source_url='https://example.invalid',
                      accession='FIXSINGLE', accession_version='1', archetype='single_cds', notes='fictional')
        document, _, extracted = normalize_record(record, review, '0' * 64)
        result = self.run_api('single_cds.gb')
        self.assertEqual(extracted, result.proteins)
        self.assertEqual({f['feature_id'] for f in document['features']}, {f.feature_id for f in result.records[0].features})

    # -- safety of the output location ---------------------------------------------------
    def test_output_location_is_protected(self):
        api = self.api
        self.run_api('single_cds.gb')
        before = (self.tmp / 'out' / 'report.json').read_bytes()
        with self.assertRaises(api.OutputError):
            self.run_api('single_cds.gb')
        self.assertEqual((self.tmp / 'out' / 'report.json').read_bytes(), before)
        replaced = self.run_api('long_protein.gb', overwrite=True)
        self.assertEqual(replaced.records[0].construct_id, 'FIXLONG.1')
        other = self.tmp / 'other'
        other.mkdir()
        (other / 'notes.txt').write_text('keep me')
        with self.assertRaises(api.OutputError):          # --overwrite never replaces a non-bundle directory
            api.preprocess_files([FIXTURES / 'single_cds.gb'], other, options=api.PreprocessOptions(overwrite=True))
        self.assertEqual((other / 'notes.txt').read_text(), 'keep me')
        work = self.tmp / 'work'
        work.mkdir()
        source = shutil.copy(FIXTURES / 'single_cds.gb', work / 'single_cds.gb')
        original = Path(source).read_bytes()
        with self.assertRaises(api.OutputError):          # an input inside the output directory is refused
            api.preprocess_files([source], work, options=api.PreprocessOptions(overwrite=True))
        self.assertEqual(Path(source).read_bytes(), original)
        with self.assertRaises(api.PreprocessError):
            api.preprocess_files([self.tmp / 'missing.gb'], self.tmp / 'new')
        self.assertFalse((self.tmp / 'new').exists())
        self.assertEqual([p.name for p in self.tmp.iterdir() if p.name.startswith('.')], [])   # no staging leftovers

    def test_tsv_escaping_and_supplied_provenance(self):
        from Bio import SeqIO
        import io
        record = SeqIO.read(io.StringIO((FIXTURES / 'single_cds.gb').read_text()), 'genbank')
        next(f for f in record.features if f.type == 'CDS').qualifiers['product'] = ['has\ttab and \\ backslash']
        path = self.tmp / 'odd.gb'
        SeqIO.write(record, path, 'genbank')
        result = self.api.preprocess_files([path], self.tmp / 'out', options=self.api.PreprocessOptions(
            provenance={'odd.gb': {'source': 'fictional', 'source_url': 'https://example.invalid/x'}}))
        raw = (result.output_dir / 'features.tsv').read_text().splitlines()
        self.assertTrue(all(line.count('\t') == raw[0].count('\t') for line in raw))
        construct = json.loads((result.output_dir / 'normalized.json').read_text())['constructs'][0]
        self.assertEqual(construct['provenance']['supplied']['source'], 'fictional')
        self.assertEqual(construct['provenance']['missing_fields'], ['retrieved_on', 'terms_note'])

    # -- command line -------------------------------------------------------------------
    def test_cli_exit_codes_help_and_isolation(self):
        outside = self.tmp / 'elsewhere'
        outside.mkdir()
        ok = cli('preprocess', FIXTURES / 'single_cds.gb', '--out', self.tmp / 'c0', cwd=outside)
        self.assertEqual(ok.returncode, 0, ok.stderr)
        self.assertIn('proteins exported: 1', ok.stdout)
        partial = cli('preprocess', FIXTURES, '--out', self.tmp / 'c2', '--max-model-residues', '1022', cwd=outside)
        self.assertEqual(partial.returncode, 2, partial.stderr)
        self.assertIn('malformed_record', partial.stderr)
        self.assertTrue((self.tmp / 'c2' / 'proteins.fasta').stat().st_size > 0)
        again = cli('preprocess', FIXTURES / 'single_cds.gb', '--out', self.tmp / 'c0', cwd=outside)
        self.assertEqual(again.returncode, 1)
        self.assertIn('not empty', again.stderr)
        for bad in (['preprocess', FIXTURES / 'single_cds.gb'], ['preprocess', '--out', self.tmp / 'x'], ['nonsense'],
                    ['preprocess', self.tmp / 'absent.gb', '--out', self.tmp / 'c1'],
                    ['preprocess', FIXTURES / 'not_a_sequence.gb', '--out', self.tmp / 'c3']):
            self.assertEqual(cli(*bad, cwd=outside).returncode, 1, bad)
        helped = cli('preprocess', '--help', cwd=outside)
        self.assertEqual(helped.returncode, 0)
        self.assertIn('--out', helped.stdout)

    def test_core_path_loads_no_model_or_experiment_code_and_opens_no_socket(self):
        script = (
            'import json, socket, sys\n'
            'def refuse(*a, **k): raise RuntimeError("network use in core preprocessing")\n'
            'socket.socket.connect = refuse; socket.create_connection = refuse\n'
            'import iri_constructs\n'
            'r = iri_constructs.preprocess_files([sys.argv[1]], sys.argv[2])\n'
            # numpy is a Biopython dependency, not a model library, so it is not listed here.
            'heavy = sorted(m for m in sys.modules if m.split(".")[0] in ("torch", "transformers", "esm", "huggingface_hub", "psutil"))\n'
            'legacy = sorted(m for m in sys.modules if m.startswith("iri_constructs.") and m.split(".")[1] not in '
            '("preprocess", "annotations", "runtime", "schema_check", "sequences"))\n'
            'print(json.dumps([r.exit_code, heavy, legacy, iri_constructs.__file__]))\n')
        env = dict(os.environ)
        if not INSTALLED:
            env['PYTHONPATH'] = str(ROOT / 'src')
        done = subprocess.run([sys.executable, '-c', script, str(FIXTURES / 'edge_cases.gb'), str(self.tmp / 'iso')],
                              capture_output=True, text=True, cwd=self.tmp, env=env)
        self.assertEqual(done.returncode, 0, done.stderr)
        code, heavy, legacy, location = json.loads(done.stdout)
        self.assertEqual((code, heavy, legacy), (0, [], []))
        if INSTALLED:
            self.assertNotIn(str(ROOT), location)


if __name__ == '__main__':
    unittest.main()
