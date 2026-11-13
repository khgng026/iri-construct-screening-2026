# IRI Fall 2026: annotated construct decomposition

The cleared E1–E4 descriptive pilot has run. Read the [current results, limitations,
workflow status and handoff](docs/EXPERIMENT_STATUS.md). [PR #11](https://github.com/khgng026/iri-construct-screening-2026/pull/11) was merged on November 13, 2026;
the study has not established a general decomposition benefit.

The research question is whether converting documented benign DNA constructs
into traceable protein components improves protein-family retrieval. Existing
protein models receive validated annotated proteins, never whole-plasmid DNA,
inferred ORFs or concatenated separate proteins. No model is trained or tuned.

The original E2 endpoint gives longest eligible CDS 2/2 supported parent families
and all eligible CDSs 1/2. Component attribution exposes off-family baseline
credit and competition for the shared five-hit budget. E3 alignment and pinned
ESM-2 each recover 4/4 query families. E4 separates ESM-2 loading from inference
and alignment work. Four queries, three parents and two dependence clusters cannot
establish a winner, equivalence or a population effect. Function prediction and
autonomous safety decisions are outside this evidence.

## Install and test

Python 3.10+ supports dependency-free checks. Parsing uses the separately pinned
optional requirements; model execution uses its recorded environment and exact
cached checkpoint. Installation requires ordinary package-provider access.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-parsing.txt
.venv/bin/python scripts/smoke_test.py
.venv/bin/python -m unittest discover -s tests -v
```

The smoke command is a fictional metadata check. Tests do not establish new
biological evidence or clear new panels. Model setup is documented in
[the experiment interface](experiments/README.md), with exact run commands and
resource measurements in the current handoff.

## Source-free reproducible results demo

From a clean checkout, the standard library is sufficient to reconstruct the
cleared E4 table and derived E3 archetype report from committed approved evidence:

```bash
python3 scripts/analyze_cleared_pilot.py --stage E4
```

Outputs are versioned under `docs/results/analysis-2026-11-12/`. Exact approved E3
and E4 protocol hashes are checked before calculation. No downloads, private
sequences or model inference occur. This demonstrates reporting, not an independent
biological rerun.

The existing annotation-to-unit demo uses the separately admitted historical
pBR322 snapshot and produces ignored preparatory outputs:

```bash
.venv/bin/python scripts/validate_week2.py \
  --inputs data/manifests/week2/demo.preparatory.json --output results/demo-new
```

For E1 parsing and E2/E3 retrieval/replay commands, use the
[reproduction instructions](docs/EXPERIMENT_STATUS.md#reproduce-and-maintain).
Restore the exact hash-identified private source/reference/model caches first.
Missing or modified scientific inputs block execution. Use a fresh run directory.
Cached-score reuse verifies the immutable committed experimental fingerprints;
it performs zero additional model calls and is not a new independent observation.

## Evidence and governance

- [Original four-week workbook](planning/4-week-workflow.xlsx) defines task criteria.
- [T02](docs/week1/T02_data_governance.md) and [T07](docs/week1/T07_benchmark_protocol.md),
  their versioned addenda and actual [authorization](docs/oversight/EXPERIMENT_RUN_AUTHORIZATION.json)
  govern source admission, frozen panels, labels and scored execution.
- [E4 clearance](docs/oversight/E4_RUN_AUTHORIZATION.json) records the actual narrow
  external delegated decision; exploratory diagnostics do not amend primary criteria.
- [Original results](docs/results/cleared-2026-11-11/) remain unchanged.
- [Prospective panel v2](data/manifests/prospective-2026-11-12/panel.v2.json) is prepared,
  unscored and pending actual delegated clearance. Renilla, ESM-C scored retrieval
  and stricter relatedness tiers remain held.

Source sequences, full qualifiers, embeddings, weights and full logs remain in
ignored private storage. Public accessibility does not itself grant redistribution
rights. Source-free manifests retain original identities and checksums. The
historical review freeze and removed duplicate diagnostics are recoverable at
reviewed commit `3fb4262`. The launcher preserves the complete original authority
digest and checks a separately pinned dependency inventory for the selected stage.
Replay requires its source/reference/cache evidence, without model weights;
inference additionally requires the pinned ESM-2 weights. Held ESM-C and iGEM
originals are not required. Earlier dated reports are historical snapshots; current results
are linked above.
