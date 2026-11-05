"""Independent audit edge cases; no model, alignment or production metric import."""
from collections import Counter
from copy import deepcopy
import datetime as dt
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import session_extension_audit as audit


def query(names, qid="q"):
    return {"query_id": qid, "group_id": "g-" + qid, "suite": "quadruples" if len(names) == 4 else "pairs",
            "components": names, "component_count": len(names),
            "spans": [{"accession": a, "start0": i * 10, "end0": (i + 1) * 10, "fraction": 1 / len(names)}
                      for i, a in enumerate(names)]}


def ranking(pairs):
    pairs = sorted(pairs, key=lambda item: (-item[1], item[0]))
    counts = Counter(value for _, value in pairs)
    return [{"accession": a, "rank": i, "score": value, "cosine_similarity": value, "exact_score_tie_count": counts[value]}
            for i, (a, value) in enumerate(pairs, 1)]


def toy_context():
    return {"hashes": {a: a for a in "ABXY"}, "equivalence": {a: [a] for a in "ABXY"},
            "families": {"A": ["A", "X"], "B": ["B", "Y"]},
            "catalogues": {"masks": {"accession105": list("ABXY")}, "aliases": {}},
            "removals": {"exact_removed": {"A": ["A"], "B": ["B"]},
                         "near90_removed": {"A": ["A", "X"], "B": ["B", "Y"]}}}


def toy_family_record(ctx, q, raw):
    recipes = audit.expected_recipes(ctx, q)
    values = {}
    for name, recipe in recipes.items():
        rows = audit.filtered(raw, recipe["included_accessions"])
        values[name] = {"whole": audit.expected_evaluation(q, rows, ctx["equivalence"],
                       ctx["families"] if name.startswith("family_") else None)}
    return {"query": q, "ranking_recipes": recipes, "source_ranking_pointers": {"whole": "/rankings/whole"},
            "evaluations": values, "exact_reference_removal_diagnostics": {"whole": audit.removal_diagnostics(
                raw, audit.filtered(raw, recipes["family_exact_removed"]["included_accessions"]), q)}}


class MissingLabelTests(unittest.TestCase):
    def test_absent_exact_is_null_surviving_family_is_scorable(self):
        q = query(["A", "B"])
        result = audit.expected_evaluation(q, ranking([("X", .9), ("Y", .5)]),
                                           {"A": ["A"], "B": ["B"]}, {"A": ["A", "X"], "B": ["B"]})
        self.assertFalse(result["components"][0]["exact"]["scorable"])
        self.assertIsNone(result["components"][0]["exact"]["hit_at_1"])
        self.assertEqual(result["components"][0]["labelled"]["rank"], 1)
        metrics = result["query_metrics"]["labelled"]
        self.assertEqual(metrics["scorable_components"], 1)
        self.assertEqual(metrics["component_recall_at_1"], 1)
        self.assertIsNone(metrics["all_components_hit_at_1"])
        self.assertIsNone(metrics["any_component_hit_at_1"])
        self.assertEqual(metrics["any_scorable_component_hit_at_1"], 1)

    def test_zero_all_scorable_quadruple_denominator_is_null_not_zero(self):
        q = query(list("ABCD"))
        value = audit.expected_evaluation(q, ranking([("X", .9)]), {a: [a] for a in "ABCD"}, {a: ["X"] if a == "A" else [a] for a in "ABCD"})
        accumulator = audit.Rollup()
        accumulator.add(value)
        result = accumulator.result()["kinds"]["labelled"]
        self.assertEqual(result["original_components"], 4)
        self.assertEqual(result["scorable_components"], 1)
        self.assertEqual(result["unscorable_components"], 3)
        qm = result["by_component_count"]["4"]["query_metrics"]
        self.assertEqual(qm["all_components_scorable_queries"], 0)
        self.assertEqual(qm["excluded_unscorable_queries"], 1)
        self.assertIsNone(qm["all_components_hit_at_10"])
        bad = deepcopy(result)
        bad["by_component_count"]["4"]["query_metrics"]["all_components_hit_at_10"] = 0
        with self.assertRaises(ValueError):
            audit.same(bad, result)

    def test_all_unscorable_summary_preserves_targets_support_and_nulls(self):
        q = query(["A", "B"])
        acc = audit.Rollup()
        acc.add(audit.expected_evaluation(q, [], {"A": ["A"], "B": ["B"]}))
        data = acc.result()["kinds"]["exact"]
        self.assertEqual(data["original_components"], 2)
        self.assertEqual(data["scorable_targets"], 0)
        self.assertEqual(data["parent_queries"], 1)
        self.assertEqual(set(data["by_target"]), {"A", "B"})
        for field in ("mean_rank", "mrr", "recall_at_1", "per_target_macro_recall_at_1"):
            self.assertIsNone(data[field])

    def test_summary_strata_and_target_macro_have_distinct_denominators(self):
        acc = audit.Rollup()
        equivalence = {a: [a] for a in "ABX"}
        for qid, pairs in [("q1", [("A", .9), ("B", .8), ("X", .1)]),
                           ("q2", [("B", .9), ("X", .8), ("A", .1)])]:
            acc.add(audit.expected_evaluation(query(["A", "B"], qid), ranking(pairs), equivalence))
        data = acc.result()["kinds"]["exact"]
        self.assertEqual(data["mean_rank"], 1.75)
        self.assertAlmostEqual(data["mrr"], (1 + .5 + 1 / 3 + 1) / 4)
        self.assertEqual(data["hits_at_1"], 2)
        self.assertEqual(data["per_target_macro_recall_at_1"], .5)
        self.assertEqual(data["query_metrics"]["any_component_hit_at_1"], 1)
        self.assertEqual(data["query_metrics"]["all_components_hit_at_1"], 0)
        self.assertEqual(data["by_position"]["C"]["mean_rank"], 1.5)
        self.assertEqual(data["by_fraction"]["0.40_to_lt_0.60"]["original_components"], 4)
        bad = deepcopy(data)
        bad["by_position"]["C"]["scorable_components"] = 1
        with self.assertRaises(ValueError):
            audit.same(bad, data)


class MaskAndFusionTests(unittest.TestCase):
    def test_common_mask_unions_every_component(self):
        ctx, q = toy_context(), query(["A", "B"])
        masks = audit.family_masks(q, ctx["hashes"], ctx["removals"])
        self.assertEqual(masks["family_exact_removed"]["included_accessions"], ["X", "Y"])
        self.assertEqual(masks["family_near90_removed"]["included_accessions"], [])
        rows = audit.denominator_rows([q], ctx["hashes"], ctx["families"], ctx["removals"], "near90_removed")
        totals = audit.denominator_totals(rows)
        self.assertEqual(totals["all"]["unscorable_component_observations"], 2)
        self.assertEqual(totals["combined"]["all_components_scorable_queries"], 0)

    def test_tampered_recipe_rejected_even_if_hash_and_evaluation_are_forged(self):
        ctx, q = toy_context(), query(["A", "B"])
        raw = ranking([("A", 1.), ("B", .9), ("X", .8), ("Y", .7)])
        record = toy_family_record(ctx, q, raw)
        audit.check_family_query(record, {"whole": (raw, "/rankings/whole")}, q, ctx)
        mask = record["ranking_recipes"]["family_exact_removed"]
        mask["included_accessions"] = ["B", "X", "Y"]
        mask["removed_accessions"] = ["A"]
        mask["included_accessions_sha256"] = audit.canonical_hash(mask["included_accessions"])
        record["evaluations"]["family_exact_removed"]["whole"] = audit.expected_evaluation(
            q, audit.filtered(raw, mask["included_accessions"]), ctx["equivalence"], ctx["families"])
        with self.assertRaisesRegex(ValueError, "ranking masks/recipes"):
            audit.check_family_query(record, {"whole": (raw, "/rankings/whole")}, q, ctx)

    def test_filtered_ranking_renumbers_before_component_metrics(self):
        rows = audit.filtered(ranking([("A", 1.), ("B", .9), ("X", .8)]), ["X"])
        self.assertEqual(audit.match(rows, ["X"])["rank"], 1)
        with self.assertRaises(ValueError):
            audit.filtered(rows, ["Y"])

    def test_totalK_fusion_is_not_per_component_oracle_budget(self):
        q = query(["A", "B"])
        lists = [ranking([("A", 1.), ("B", .2)]), ranking([("B", 1.), ("A", .2)])]
        equivalence = {a: [a] for a in "AB"}
        oracle = audit.expected_evaluation(q, dict(zip(q["components"], lists)), equivalence)
        self.assertEqual(oracle["query_metrics"]["exact"]["all_components_hit_at_1"], 1)
        for method in ("max_score", "rrf"):
            fused = audit.fused_ranking(lists, method, list("AB"))
            self.assertEqual([row["accession"] for row in fused], ["A", "B"])
            self.assertEqual(fused[0]["exact_score_tie_count"], 2)
            total = audit.expected_evaluation(q, fused, equivalence)
            self.assertEqual(total["query_metrics"]["exact"]["all_components_hit_at_1"], 0)
            self.assertEqual(total["returned_hits_ceiling"]["1"], 1)
            self.assertEqual(total["retrieval_budget"], "k total per original query")
            self.assertEqual(total["exact_recall_ceiling"]["1"], .5)
        self.assertAlmostEqual(audit.fused_ranking(lists, "rrf", list("AB"))[0]["score"], 1 / 61 + 1 / 62)

    def test_fusion_rejects_incomplete_or_nonfinite_original_rankings(self):
        rows = ranking([("A", 1.), ("B", .1)])
        with self.assertRaises(ValueError):
            audit.fused_ranking([rows[:1]], "rrf", list("AB"))
        rows[0]["score"] = float("nan")
        with self.assertRaises(ValueError):
            audit.fused_ranking([rows], "max_score", list("AB"))

    def test_frozen_unique_catalogues_rebuilt_and_tampering_rejected(self):
        session = ROOT / "results/iri-session-20261031-032506"
        manifest = audit.read(ROOT / "config/reproduction/packet/challenge_manifest.json")
        catalogue = audit.read(ROOT / "data/manifests/t05_expanded_catalogue.json")
        hashes = {r["accession"]: r["sequence_sha256"] for r in catalogue["records"]}
        masks = audit.rebuild_catalogues(manifest, hashes)
        self.assertEqual(len(masks["masks"]), 8)
        self.assertEqual(len(masks["masks"]["unique25_rep0"]), 25)
        self.assertEqual(len(masks["masks"]["unique50_rep2"]), 50)
        self.assertEqual(masks["aliases"]["unique101_rep2"], "unique101_rep0")
        manifest["catalogue_subsets"][0]["entries"][0]["representative_accession"] = "forged"
        with self.assertRaises(ValueError):
            audit.rebuild_catalogues(manifest, hashes)


class SnapshotTests(unittest.TestCase):
    def build(self, root):
        session = root / "session"
        identity = {"mode": "alignment", "query_ids": ["q"]}
        digest = audit.canonical_hash(identity)
        job = session / "jobs" / ("source-" + digest[:12])
        job.mkdir(parents=True)
        future = dt.datetime.now(dt.timezone.utc) + dt.timedelta(hours=1)
        (session / "launch.json").write_text(json.dumps({"hard_stop_utc": future.isoformat()}))
        (job / "identity.json").write_text(json.dumps(identity))
        units = {}
        for key in ("query-q", "summary"):
            unit = job / "units" / key
            unit.mkdir(parents=True)
            record = unit / "record.json"
            record.write_text(json.dumps({"saved": key}))
            descriptor = {"key": key, "identity_hash": digest, "files": {"record.json": {"bytes": record.stat().st_size, "sha256": audit.sha(record)}}}
            (unit / "complete.json").write_text(json.dumps(descriptor))
            units[key] = {"path": str(unit.relative_to(job)), "descriptor_sha256": audit.sha(unit / "complete.json")}
        completion = job / "completion-fixture.json"
        completion.write_text(json.dumps({"identity_hash": digest, "expected_count": 2, "units": units, "completed_utc": future.isoformat()}))
        return session, job, completion

    def test_complete_snapshot_resolves_alias_and_rejects_changed_record(self):
        with tempfile.TemporaryDirectory() as temporary:
            session, job, completion = self.build(Path(temporary))
            snap = audit.Snapshot(job, session, audit.Deadline(60, session), audit.sha(completion))
            self.assertEqual(snap.record("query-q"), {"saved": "query-q"})
            (snap.units["query-q"] / "record.json").write_text('{"saved":"tampered"}')
            with self.assertRaises(ValueError):
                snap.record("query-q")
            with self.assertRaises(ValueError):
                audit.Snapshot(job, session, audit.Deadline(60, session))

    def test_missing_unit_cannot_be_hidden_by_consistent_expected_count(self):
        with tempfile.TemporaryDirectory() as temporary:
            session, job, completion = self.build(Path(temporary))
            value = audit.read(completion)
            del value["units"]["query-q"]
            value["expected_count"] = 1
            completion.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError, "unit coverage"):
                audit.Snapshot(job, session, audit.Deadline(60, session))

    def test_external_job_and_wrong_completion_hash_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            session, job, completion = self.build(Path(temporary))
            with self.assertRaisesRegex(ValueError, "escaped session"):
                audit.Snapshot(job, session / "elsewhere", audit.Deadline(60, session))
            with self.assertRaises(ValueError):
                audit.Snapshot(job, session, audit.Deadline(60, session), "0" * 64, completion)

    def test_expired_deadline_and_atomic_report(self):
        with tempfile.TemporaryDirectory() as temporary:
            session, _, _ = self.build(Path(temporary))
            (session / "launch.json").write_text(json.dumps({"hard_stop_utc": "2000-01-01T00:00:00+00:00"}))
            with self.assertRaises(TimeoutError):
                audit.Deadline(60, session).check()
            output = session / "audit/report.json"
            audit.atomic_report(output, {"status": "failed", "missing": None})
            self.assertEqual(audit.read(output), {"status": "failed", "missing": None})
            self.assertEqual(list(output.parent.iterdir()), [output])


if __name__ == "__main__":
    unittest.main()
