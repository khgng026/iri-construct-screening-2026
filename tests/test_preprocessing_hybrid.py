"""Regression checks for the sequence-free hybrid replay."""
from unittest import TestCase

from experiments.preprocessing_exploration.hybrid import (
    Candidate,
    family_predictions,
    same_frame_overlap,
    select_order,
    selection_orders,
)


class HybridCandidateTest(TestCase):
    def setUp(self):
        self.gene_a = Candidate("gene_caller", 0, 1, 0, 300, 100)
        self.gene_b = Candidate("gene_caller", 1, 1, 600, 900, 100)
        self.duplicate = Candidate("start_ranked", 0, 1, 0, 300, 100)
        self.outside = Candidate("start_ranked", 1, 1, 330, 480, 50)
        self.pools = {"gene_caller": [self.gene_a, self.gene_b],
                      "minimal": [], "start_ranked": [self.duplicate, self.outside]}

    def test_distinct_supplements_use_coordinate_and_frame(self):
        other_frame = Candidate("start_ranked", 2, 1, 1, 301, 100)
        self.assertEqual(same_frame_overlap(self.gene_a, self.duplicate), 1)
        self.assertEqual(same_frame_overlap(self.gene_a, other_frame), 0)
        genes_then_distinct = selection_orders(self.pools, 0.8, 1)["gene_first_distinct"]
        self.assertEqual(genes_then_distinct, [self.gene_a, self.gene_b, self.outside])

    def test_selector_enforces_both_budgets_and_skips_exact_duplicate(self):
        order = [self.gene_a, self.duplicate, self.gene_b, self.outside]
        selected = select_order(order, 2, 150, unique_loci=True)
        self.assertEqual(selected, (self.gene_a, self.outside))
        selected = select_order(order, 2, 200, unique_loci=True)
        self.assertEqual(selected, (self.gene_a, self.gene_b))

    def test_reserved_slot_displaces_only_the_last_primary_candidate(self):
        order = selection_orders(self.pools, 0.8, 1)["reserve_one_distinct"]
        self.assertEqual(select_order(order, 2, None, unique_loci=True),
                         (self.gene_a, self.outside))

    def test_output_cap_cycles_through_candidates(self):
        candidates = [self.gene_a, self.outside]
        ranks = {self.gene_a: ["a", "b"], self.outside: ["q", "s"]}
        self.assertEqual(
            family_predictions(candidates, ranks, 3),
            [{"query_id": self.gene_a.query_id, "family": "a"},
             {"query_id": self.outside.query_id, "family": "q"},
             {"query_id": self.gene_a.query_id, "family": "b"}],
        )
