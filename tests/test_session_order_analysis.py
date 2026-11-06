"""Order-analysis fixtures; no model, alignment or production saved outcomes."""
from collections import Counter
from copy import deepcopy
import itertools
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import session_order_analysis as order


def queries(names=("A", "B")):
    lengths = {a: i + 2 for i, a in enumerate(sorted(names))}
    result = []
    for arrangement in sorted(order.allowed_orders(names)):
        total, offset, spans = sum(lengths.values()), 0, []
        for a in arrangement:
            spans.append({"accession": a, "start0": offset, "end0": offset + lengths[a], "fraction": lengths[a] / total})
            offset += lengths[a]
        result.append({"query_id": "q_" + "_".join(arrangement), "components": list(arrangement),
                       "component_count": len(names), "length_aa": total, "group_id": "__".join(sorted(names)),
                       "suite": {1: "clean", 2: "pairs", 3: "triples", 4: "quadruples"}[len(names)], "spans": spans})
    return result


def ranking(scores):
    values = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    ties = Counter(scores.values())
    return [{"accession": a, "score": value, "cosine_similarity": value, "rank": i,
             "exact_score_tie_count": ties[value]} for i, (a, value) in enumerate(values, 1)]


def observations(group, accession="A", ranks=(5, 6), scores=(.6, .4)):
    return [{"query_id": q["query_id"], "order": q["components"], "position1": q["components"].index(accession) + 1,
             "rank": rank, "score": score, "hit_at_5": int(rank <= 5)} for q, rank, score in zip(group, ranks, scores)]


class FrozenPanelTests(unittest.TestCase):
    def test_quadruples_are_eight_orders_not_all_twenty_four(self):
        planned = order.allowed_orders(("A", "B", "C", "D"))
        self.assertEqual(len(planned), 8)
        self.assertNotEqual(planned, set(itertools.permutations("ABCD")))
        for target in "ABCD":
            self.assertEqual(Counter(arrangement.index(target) for arrangement in planned), {0: 2, 1: 2, 2: 2, 3: 2})

    def test_group_orders_and_fixed_fraction_are_checked(self):
        group = queries(("A", "B", "C"))
        self.assertEqual(len(order.group_queries(group)["A__B__C"]), 6)
        bad = deepcopy(group); bad[1]["spans"][0]["fraction"] += .01
        with self.assertRaisesRegex(ValueError, "fraction/geometry"):
            order.group_queries(bad)
        with self.assertRaisesRegex(ValueError, "Missing, duplicated"):
            order.group_queries(group[:-1])
        with self.assertRaisesRegex(ValueError, "Duplicate query"):
            order.group_queries(group + group[:1])

    def test_changed_group_composition_rejected(self):
        group = queries()
        group[1]["components"] = ["A", "C"]
        with self.assertRaisesRegex(ValueError, "composition"):
            order.group_queries(group)

    def test_full_query_coverage_rejects_partial_and_extra_sources(self):
        ids = ["q" + str(i) for i in range(1146)]
        accessions = ["a" + str(i) for i in range(105)]
        source = {"mode": "model", "query_ids": ids, "catalogue_accessions": sorted(accessions)}
        order.require_full_panel(source, ids, accessions)
        for altered in (ids[:-1], ids + ["unexpected"], ids[:-1] + [ids[0]]):
            with self.assertRaisesRegex(ValueError, "intersections are forbidden"):
                order.require_full_panel({**source, "query_ids": altered}, ids, accessions)


class DescriptiveCellTests(unittest.TestCase):
    def test_rank5_boundary_flip_is_retained_with_position_support(self):
        group = queries()
        cell = order.describe_cell(group, "A", "whole", observations(group))
        self.assertEqual((cell["rank_min"], cell["rank_max"], cell["rank_range"]), (5, 6, 1))
        self.assertTrue(cell["r5_flip"])
        self.assertEqual(cell["r5_hit_fraction"], .5)
        self.assertEqual(cell["position_means"][0]["within_cell_centered_means"]["rank"], -.5)
        self.assertEqual(cell["position_means"][1]["within_cell_centered_means"]["rank"], .5)
        self.assertTrue(cell["position_support_balanced"])

    def test_missing_position_support_is_null(self):
        obs = [{"query_id": "q", "position1": 1, "rank": 5, "score": .5, "hit_at_5": 1}]
        values = order.position_means(obs, 2)
        self.assertEqual(values[1]["observation_count"], 0)
        self.assertIsNone(values[1]["means"]["rank"])
        self.assertIsNone(values[1]["within_cell_centered_means"]["hit_at_5"])
        absent = order.aggregate_positions([])
        self.assertTrue(all(r["group_target_cells"] == 0 and r["equal_cell_position_means"]["score"] is None for r in absent))

    def test_group_macro_does_not_weight_more_targets_more_heavily(self):
        group = queries()
        cell = order.describe_cell(group, "A", "whole", observations(group, ranks=(1, 11)))
        cells = [cell]
        for accession in ("B", "C", "D"):
            clone = deepcopy(cell); clone.update(group_id="other", accession=accession, rank_range=0)
            cells.append(clone)
        summary = order.aggregate_cells(cells)
        self.assertEqual(summary["cell_weighted"]["rank_range"]["mean"], 2.5)
        self.assertEqual(summary["group_weighted_means"]["rank_range"], 5)
        self.assertEqual(summary["group_target_cells"], 4)

    def test_clean_single_order_is_excluded_from_permutation_summaries(self):
        group = queries(("A",))
        cell = order.describe_cell(group, "A", "whole", observations(group, ranks=(1,), scores=(1.,)))
        self.assertTrue(cell["single_order_control"])
        with self.assertRaisesRegex(ValueError, "Clean/single-order"):
            order.aggregate_cells([cell])

    def test_cell_cannot_drop_an_order(self):
        group = queries()
        with self.assertRaisesRegex(ValueError, "lost an observed order"):
            order.describe_cell(group, "A", "whole", observations(group)[:1])


class MixtureSanityTests(unittest.TestCase):
    def test_near_tie_inversion_is_labelled_and_never_forced_equal(self):
        anchor = ranking({"A": .5, "B": .5 - 1e-7, "C": .1})
        observed = ranking({"A": .5 - 1e-7, "B": .5, "C": .1})
        result = order.mixture_changes(anchor, observed, ["A", "B", "C"])
        self.assertEqual(len(result["changed_candidate_ranks"]), 2)
        self.assertEqual(len(result["inverted_candidate_pairs"]), 1)
        self.assertTrue(result["all_rank_changes_within_near_ties"])
        self.assertFalse(result["scores_exactly_equal"])
        self.assertEqual(result["changed_candidate_ranks"][0]["observed_rank"], 2)

    def test_large_rank_inversion_fails_sanity_instead_of_reranking(self):
        result = order.mixture_changes(ranking({"A": .8, "B": .2}), ranking({"A": .2, "B": .8}), ["A", "B"])
        self.assertFalse(result["all_rank_changes_within_near_ties"])
        self.assertFalse(result["scores_within_tolerance"])
        self.assertEqual(result["changed_candidate_ranks"][0]["observed_rank"], 2)

    def test_corrupt_tie_metadata_or_nonfinite_scores_are_rejected(self):
        rows = ranking({"A": .5, "B": .5}); rows[0]["exact_score_tie_count"] = 1
        with self.assertRaises(ValueError): order.rank_map(rows, ["A", "B"])
        rows = ranking({"A": .5, "B": .2}); rows[0]["score"] = float("nan")
        with self.assertRaises(ValueError): order.rank_map(rows, ["A", "B"])

    def toy_saved_group(self, perturb=False):
        import numpy as np
        group = queries()
        refs = {"A": np.array([1, 0], dtype=np.float32), "B": np.array([0, 1], dtype=np.float32),
                "X": np.array([-1, -1], dtype=np.float32)}
        records, arrays = {}, {}
        for index, query in enumerate(group):
            mixture = np.array([.4, .6], dtype=np.float32)
            if perturb and index == 1: mixture[0] += .01
            vector = {"whole": mixture.copy(), "isolated_mixture": mixture,
                      "isolated_A": refs["A"], "isolated_B": refs["B"],
                      "contextual_A": refs["A"], "contextual_B": refs["B"]}
            def scores(value): return ranking({a: order.independent_cosine(value, ref) for a, ref in refs.items()})
            key = "query-" + query["query_id"]
            arrays[key] = vector
            records[key] = {"query": query, "rankings": {"whole": scores(mixture), "isolated_mixture": scores(mixture),
                            "contextual_component": {a: scores(refs[a]) for a in ("A", "B")},
                            "isolated_component": {a: scores(refs[a]) for a in ("A", "B")}}}
        with patch.object(order, "read_arrays", side_effect=lambda snapshot, key: arrays[key]):
            return order.summarize_group(SimpleNamespace(record=lambda key: records[key]), group, refs, SimpleNamespace(check=lambda: None))

    def test_saved_vectors_and_all_candidate_scores_agree_with_canonical_mixture(self):
        result = self.toy_saved_group()
        self.assertEqual(result["mixture_sanity"]["status"], "passed")
        self.assertTrue(all(c["vector_exactly_equal"] and c["scores_exactly_equal"] for c in result["mixture_sanity"]["checks"]))
        self.assertTrue(result["no_rank_or_vector_equality_forced"])

    def test_vector_perturbation_is_persisted_as_failed_with_observed_scores(self):
        result = self.toy_saved_group(perturb=True)
        self.assertEqual(result["mixture_sanity"]["status"], "failed")
        cell = next(c for c in result["cells"] if c["arm"] == "isolated_mixture" and c["accession"] == "A")
        self.assertGreater(cell["score_range"], 1e-6)
        self.assertGreater(result["mixture_sanity"]["checks"][1]["maximum_absolute_vector_difference"], 1e-6)


if __name__ == "__main__":
    unittest.main()
