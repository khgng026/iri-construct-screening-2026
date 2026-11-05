# T05 implementation handoff — 2026-11-05

T01–T04 remain approved. T05 remains **in progress**. The implementation is ready
for code review; reviewed-data pretrained experiments and a feasibility decision
remain outstanding. Implementation-stage verification preceded Git publication;
publication is recorded in the repository's Git history.

## Changes and file inventory

| Files | Change |
| --- | --- |
| `src/iri_constructs/runtime.py` | New standard-library run metadata, code/input hashes, controlled errors, ignored-output confinement and optional resource sampling |
| `src/iri_constructs/bench_inputs.py` | New shared T02 review/checksum gate, source binding, canonical sequence validation and traceable protein inputs |
| `src/iri_constructs/model_adapters.py` | New common residue interface; fixed ESM-2, ESM-C 300M and ESM-3 Small adapters; immutable local snapshots, strict pretrained loads and token checks |
| `src/iri_constructs/metrics.py` | New within-model cosine/pooling and independent-family retrieval, margins, coverage and abstentions |
| `src/iri_constructs/fusion.py`, `src/iri_constructs/fusion_legacy.py` | New isolated/full-fusion/contextual-target comparison with reviewed coordinates and per-condition failures/resources; upstream validation/numerical APIs retained for compatibility |
| `src/iri_constructs/annotations.py` | New Biopython T03 normalization and separate qualifier/location audit contract |
| `src/iri_constructs/experiment_cli.py` | New shared CLI, sequential execution, per-model results and exact common query/reference panels |
| `src/iri_constructs/config.py`, `src/iri_constructs/validate_manifest.py` | Improved configuration errors/BOM support and malformed CSV/duplicate-header validation; accepted schemas retained |
| `experiments/t05_esm2_fusion/compare_fusion.py`, `experiments/t05_model_comparison/compare_models.py` | New path-independent experiment entry points |
| `scripts/t05_model_probe.py` | Existing hard-coded probe replaced with the reviewed-input CLI |
| `scripts/parse_annotations.py`, `scripts/t05_environment.py` | New parsing demonstration and dependency/cache/device/offline-load preflight |
| `requirements-parsing.txt`, `experiments/requirements-esm2.txt`, `experiments/requirements-esm-sdk.txt` | New separate parser/model dependency instructions with direct pins |
| `tests/test_annotations.py`, `tests/test_bench_inputs.py`, `tests/test_experiment_cli.py`, `tests/test_fusion.py`, `tests/test_model_adapters.py`, `tests/test_parse_cli.py`, `tests/test_retrieval.py` | New meaningful parsing, coordinate, governance, adapter and retrieval regression tests |
| `tests/test_config.py`, `tests/test_validate_manifest.py` | Extended existing standard-library metadata tests |
| `.github/workflows/ci.yml` | Existing dependency-free job retained; separate optional parser job added |
| `.gitignore` | Additional local model environment, weight/embedding and legacy-output exclusions |
| `README.md`, `docs/week1/T01_scope_and_acceptance.md`, `docs/week1/T03_archetypes_and_schema.md` | Updated stale status text to the user's current approvals; four-week scope and T03 v0.1 retained |
| `docs/week1/STATUS.md`, `docs/week1/T05_annotation_audit_extension.md`, this report | New current status, separate audit/input contracts and handoff evidence |
| `experiments/README.md`, both experiment `README.md` files | New setup, CLI commands, protocols, provenance and limitations |
| `planning/4-week-workflow.xlsx` | Final workbook starts from the newer upstream workbook, preserves all approvals/probe notes, and updates only T05 status/note cells G9/L9 |
| `iri_t05_model_comparison/compare_models.py`, `README.md`, `requirements.txt` | Pre-existing untracked prototype retained for review; executable redirects to the shared reviewed-input runner before optional imports; requirements/docs point to the new workflow |

The T02 policy text, T03 example JSON files and baseline schema checker are
unchanged; upstream approval-record additions are preserved. No raw source sequence or embedding artifact was added to Git.
Prototype inspection informed the adapter/pooling and leave-construct-out
distinct-sequence retrieval policies; automatic candidate downloading, heuristic
label approval and winner selection are excluded from current executable use.

## Verification actually executed

| Check | Observed result |
| --- | --- |
| Initial local suite before changes | 11 passed; the subsequently fetched upstream branch added seven fusion tests, all preserved |
| Combined `python -m unittest discover -s tests` | 58 passed, including all upstream fusion regression tests; no skips in the local installed environment |
| `python -S -m unittest discover -s tests` | 40 passed; 18 optional Biopython/PyTorch tests skipped |
| `python -S scripts/smoke_test.py` | Passed; dependency-free Week 1 checks preserved |
| Six entry points with `python -S <absolute script> --help`, unrelated working directory | All returned 0 without model/parser packages |
| Parsing CLI integration | FASTA and GenBank successfully parsed from generated fictional test records; T03 validation and qualifier audit verified |
| Numerical fusion/retrieval/SDK-contract tests | Passed; these are hand-computable software fixtures, not pretrained-model measurements |
| Offline cached ESM-2 load preflight | Encoder/tokenizer loaded at `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`; no missing encoder weights, no added pooler and no sequence inference |
| Missing reviewed-input CLI attempts | Fusion, comparison and parser returned blocked reports/nonzero exit, before sequence processing/inference |
| Workbook integrity | ZIP valid; only Weekly To-Dos XML differs from upstream, with its T05 status/note updated and earlier evidence retained |
| Git whitespace/schema/data checks | `git diff --check` passed; accepted examples/schema/T02 unchanged |

Biopython 1.86 and psutil 7.1.3 were installed locally to run parser/resource tests.
The existing environment is Python 3.14.7 with PyTorch 2.15.0+cpu; CUDA is
unavailable and the ESM SDK is absent. The separate Python 3.11 direct model pins
have not been installed and clean-tested here. Following publication, both GitHub
Actions jobs (`test` and `parsing`) passed for implementation commit `d8980cb`:
[verified CI run](https://github.com/khgng026/iri-construct-screening-2026/actions/runs/42073187942).

## Experiments and blockers

Successfully exercised: actual local record parsing on fictional software
fixtures, metadata-only smoke checks, environment/resource inspection, and an
offline cached ESM-2 **load check**. No pretrained sequence representations,
fusion similarities, family retrieval scores or scientific model comparison
were produced. Numerical unit-test scores must not be presented as experimental
results. The earlier researcher-reported GFP probe is retained as attributed
prior evidence in STATUS.md; it was not independently reproduced in this session.

The requested scientific experiments are blocked by the absence of a reviewed
local input configuration, approved source manifests/checksums, reviewed fusion
boundaries and an independent annotation-family panel. The existing local FASTA
does not include review evidence. ESM-C/ESM-3 additionally require the supported
SDK environment and accessible checkpoints. No CUDA hardware is available here;
CPU feasibility for those larger models has not been measured.

Local machine-readable evidence is ignored by Git:

- `results/t05_environment/report.json`: dependency/device/cache and actual ESM-2 load status/resources, with no sequence inference.
- `results/t05_fusion_blocked/20261105/report.json` and `results/t05_comparison_blocked/20261105/report.json`: actual blocked attempts against the absent reviewed-input config.
- `results/t05_parsing_blocked/20261105/report.json`: corresponding missing-input parser report.
- `results/workbook_before_t05.xlsx`: initial local pre-edit workbook copy for review.
- `results/workbook_remote_before_publish.xlsx`: newer upstream workbook retained locally before reconciliation.

## Review and next execution

Review the implementation and [separate annotation/input contract](T05_annotation_audit_extension.md).
Then supply the T02-reviewed local records with source permissions, notes and
checksums; approve supported family labels/independence groups and documented
fusion target coordinates. Install and freeze a compatible model environment,
validate inputs, and execute the fixed-checkpoint commands in
[experiments/README.md](../../experiments/README.md). Review actual coverage,
failures, costs and annotation agreement before deciding whether T05 has enough
evidence to close. No later benchmark task is marked complete by this handoff.

## GitHub integration

Before publication, `origin/main` was six commits ahead of the initial local
checkout. The implementation was committed locally, then merged with upstream
history without force pushing. Original T01/T03/T04 approval dates (2026-10-29),
T02 reaffirmation, earlier researcher-reported GFP evidence, pair-template metadata
and all seven upstream fusion tests were retained. The earlier pure input-validation
and numerical APIs are available through `fusion_legacy.py` and the fusion entry
point; current scientific inference uses the source-checksum reviewed-input CLI.
The old pair-manifest CLI migration is documented in the fusion README.
