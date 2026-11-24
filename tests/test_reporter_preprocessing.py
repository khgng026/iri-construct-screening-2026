"""Check the new reporter inputs without model weights or external packages."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from experiments.paired_preprocessing.reporter_pilot import (
    MAX_RESIDUES, _gene_caller_candidates, _rank_by_start, run,
)


class ReporterPreprocessingTests(unittest.TestCase):
    def test_start_evidence_reorders_but_preserves_every_candidate(self):
        units = [
            {'id': 'long_without_start', 'protein': 'A' * 80},
            {'id': 'slightly_longer', 'protein': 'A' * 55},
            {'id': 'early_start', 'protein': 'M' + 'A' * 43},
        ]
        ranked = _rank_by_start(units)
        self.assertEqual([u['id'] for u in ranked], ['long_without_start', 'early_start', 'slightly_longer'])
        self.assertCountEqual(ranked, units)

    def test_gene_caller_spans_follow_translated_chunks_on_reverse_strand(self):
        # The caller reports a translated product plus a stop in its source span.
        protein = 'A' * (MAX_RESIDUES + 30)
        end = len(protein) * 3 + 3
        gene = SimpleNamespace(begin=1, end=end, strand=-1,
                               translate=lambda *, include_stop: protein + ('*' if include_stop else ''))
        finder = SimpleNamespace(find_genes=lambda sequence: [gene])
        first, second = _gene_caller_candidates('N' * end, finder)
        self.assertEqual(first['spans'], [[93, end]])
        self.assertEqual(second['spans'], [[3, 93]])
        self.assertEqual(first['protein'] + second['protein'], protein)
        self.assertEqual((first['strand'], second['strand']), (-1, -1))

    def test_gene_caller_rejects_a_span_shorter_than_its_translation(self):
        gene = SimpleNamespace(begin=1, end=86, strand=1, translate=lambda **kwargs: 'A' * 30)
        finder = SimpleNamespace(find_genes=lambda sequence: [gene])
        with self.assertRaisesRegex(ValueError, 'source coordinates disagree'):
            _gene_caller_candidates('N' * 90, finder)

    def test_runner_rejects_changed_inputs_before_loading_model(self):
        with TemporaryDirectory() as folder:
            path = Path(folder) / 'constructs.fasta'
            path.write_text('different input', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'Sequence input differs'):
                run('start_ranked', path, path, {'sequence_input_sha256': 'different'}, path)


if __name__ == '__main__':
    unittest.main()
