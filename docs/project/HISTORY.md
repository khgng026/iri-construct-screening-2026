# Historical records and superseded instructions

Active research planning and authority live in [SCOPE.md](SCOPE.md), [AUTONOMY.md](AUTONOMY.md), [PLAN.md](PLAN.md), [STATUS.md](STATUS.md) and root [AGENTS.md](../../AGENTS.md). The [toolkit contract](TOOLKIT_CONTRACT.md) describes the existing annotation-aware package.

The following are retained for provenance and reproduction, not active release requirements or instructions for new toolkit work:

- `planning/4-week-workflow.xlsx` and `docs/week2/WORKFLOW_ACCEPTANCE.json`: original workbook and its extraction; T01–T28 ordering, 150–250-construct target and benchmark-before-packaging dependencies are superseded.
- `docs/week1/STATUS.md`, `docs/coordination/CODEX_STATUS.md`, `DECISIONS_NEEDED.md`, the execution-readiness/approval packet, `docs/week2/COMPLETION_REPORT.md` and `TASK_STATUS.json`: dated snapshots. Statements about zero admitted files, pending E1–E4 or external review describe their observation time and must not be treated as the current state.
- `docs/week1/T01_scope_and_acceptance.md`, T02/T07 documents, `docs/oversight/` and their amendments/mission records: historical approvals. Active authority replaces per-record/per-panel human or external-ChatGPT gates for new in-scope work. Preserved approvals and frozen conditions are not rewritten.
- [EXPERIMENT_STATUS.md](../EXPERIMENT_STATUS.md), `docs/runs/`, `docs/results/cleared-2026-11-11/`, `docs/results/analysis-2026-11-12/` and [panel-v2 evidence](../results/panel-v2-2026-11-13/README.md): scientific/technical observations with their original inputs and limitations.
- Dated `config/execution_*`, `config/week2/`, pilot/extension freezes and associated manifests: historical experiment definitions. For a new run, create a new input/configuration manifest; do not edit old hashes or approval fields to make the old launcher accept a new panel.
- `scripts/execute_cleared_retrieval.py`, `scripts/panel_v2_*_in_memory.py`, old A–D/E1–E4 runners and `experiments/`: reproduction/research entry points, not the new toolkit's required interface. Retain useful shared implementation and tests.

Files referenced by hashes remain at their original paths and byte content. Historical code/doc snapshots are recoverable from their Git commits; current README and active project documentation may evolve independently. No new toolkit path should require these historical Git objects or private experimental caches.

## Added November 17, 2026

`docs/week1`, `docs/week2`, `docs/coordination`, `docs/oversight` and `docs/runs` each carry a `README.md` notice saying their contents are historical. The three dated status pages that are not bound by any recorded checksum (`docs/week1/STATUS.md`, `docs/coordination/CODEX_STATUS.md`, `docs/coordination/DECISIONS_NEEDED.md`) also carry a one-line banner. Checksum-bound files (for example `docs/week1/T02_data_governance.md`, `T07_benchmark_protocol.md`, `docs/week2/COMPLETION_REPORT.md`, `TASK_STATUS.json`) and the original scientific report `docs/EXPERIMENT_STATUS.md` were left byte-for-byte unchanged; the directory notices and this index are their redirect.

New work lives elsewhere: the toolkit in `src/iri_constructs/preprocess.py` and `cli.py`; the separately versioned [component-coverage study v1](../results/component-study-v1/README.md) in `experiments/component_study_v1/` and `docs/results/component-study-v1/`. That study uses its own specification, labels, panel and endpoint. It does not re-score or amend E2/E3, whose historical parent-family endpoint and "longest eligible" baseline remain as recorded.

## Research scope reset, November 22, 2026

The active research scope supersedes the toolkit-first research plan and its restrictions on sequence-only decomposition and model comparisons. The former scope, plan, status, and agent instructions remain available in Git history; completed studies, release evidence, and frozen records retain their original content and meaning.
