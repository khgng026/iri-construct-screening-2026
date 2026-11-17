> **Historical snapshot — superseded.** This page describes the project on the date it was written. Current status and authority: [docs/project/STATUS.md](../project/STATUS.md). Nothing below is an active requirement.

# Project status - 2026-11-06 UTC

## T07 v1.1 update (2026-11-06 UTC)

The exact ten candidate identities are selected, with zero biological source files
admitted. T07 v1.0 is approved in principle and remains unchanged. E1-E4 now have
tested software implementations in a separate v1.1 addendum: annotation quality,
longest/all-CDS parent comparison, reviewed reference tiers and strict measured
cost reporting. Real experiments remain blocked. See the
[implementation report](T07_v1_1_implementation_report.md),
[milestone handoff](../coordination/CODEX_STATUS.md) and
[pending evidence](../coordination/DECISIONS_NEEDED.md).

PR #2 was merged upstream by the researcher at 2026-11-06T00:40:46Z;
main is now `284a7008f0420da4557eccf434cabef64696d787`, whose tree is identical
to PR #2's `2c9a5bf`. The v1.1 PR therefore targets main directly. This
implementation has made no merge, force push or change to main or PR #1/#2 heads.

## Earlier implementation and approval history

T01-T04 approvals and their original 2026-10-28/29 dates are preserved from
the newer GitHub history. Approval of scope/schema/governance remains separate
from approval of any source record or label. The implementation was reconciled
with the upstream fusion experiment before publication.

| Task | Status | Evidence / remaining work |
| --- | --- | --- |
| T01 | Approved | Existing MVP question and four-week scope retained |
| T02 | Approved | Policy v1.0 retained; record-level reviews remain pending |
| T03 | Approved | Baseline v0.1 and all three existing examples retained |
| T04 | Approved / implemented | Dependency-free smoke/tests retained; optional parsing is separate |
| T05 | In progress | Fusion, probe, annotation parsing and modular model comparison implemented; reviewed-data pretrained experiments remain blocked |
| T06 | Inventory implemented; identities selected | 20 discovery candidates; exact ten-identity v1.1 shortlist; zero source files admitted |
| T07 | Design/scoring implemented, scientifically unvalidated | A priori A-D protocol/configurations frozen before pilot scoring; reference labels/panel await review |
| T08 | Pipeline implemented; real-data admission blocked | Ten-record pending queue, annotation-preservation audit and regression tests; zero real records parsed/scored |
| T09–T28 | Pending | No dataset expansion, full paired benchmark, final evaluation or research delivery claimed |

The initial local checkout lacked the fusion experiment and status file;
publication now includes the newer upstream versions and their validation tests.
The existing untracked `iri_t05_model_comparison/` prototype was
reviewed and retained; its executable now delegates to the reviewed local-input
runner. The original hard-coded probe has been replaced with the same validated
CLI. No biological input files were changed and no source records/labels were
approved. Git publication is recorded separately in repository history.

Implementation and verification details, commands, dependency pins, input review
requirements and measurement limitations are in [experiments/README.md](../../experiments/README.md).
The parser's separate audit contract is documented in
[T05_annotation_audit_extension.md](T05_annotation_audit_extension.md); T03's
accepted schema was not extended silently.

Local verification: the original 11 tests passed before changes. The expanded
suite passes 58 tests, including all seven upstream fusion regression tests and
a real FASTA/GenBank parsing CLI run on generated, explicitly
fictional software fixtures, plus numerical adapter/retrieval tests. These are
software checks, not source approval or model-performance evidence. Optional
Biopython 1.86 and psutil 7.1.3 were installed in the local environment.

The dependency-free run passes 40 tests and skips 18 optional parsing/tensor tests;
the standard-library smoke check passes. All six entry points provide dependency-free
CLI help from another working directory. After publication, both GitHub CI jobs
(`test` and `parsing`) passed for implementation commit `d8980cb`:
[verified Actions run](https://github.com/khgng026/iri-construct-screening-2026/actions/runs/42073187942).

Actual environment preflight: Python 3.14.7, PyTorch 2.15.0+cpu, CUDA unavailable
(zero devices), ESM SDK absent. The local environment differs from the separate
Python 3.11 model pins, whose full clean installation has not been validated here.
The pinned cached ESM-2 8M encoder/tokenizer loaded successfully offline, with no
missing encoder parameters and no added random pooler; no sequence was passed to it.
No reviewed fusion/panel configuration was supplied. The lone local protein
FASTA has no accompanying record-level review evidence. No pretrained sequence
inference or scientific model comparison was executed; no numerical performance
claim or T05 completion is justified. The local preflight report lives under
ignored `results/t05_environment/`.

Fusion/comparison/parser CLI attempts against the missing reviewed-input config
produced blocked JSON reports under ignored `results/t05_*_blocked/20261105/`.
The planning workbook preserves upstream approval dates and probe notes;
T05 remains In progress, with its current implementation note appended.

Review remaining: implementation and separate audit contract; T02 source/licensing/
safety approvals with checksums; annotation family and independence policy;
documented fusion target boundaries; a compatible model environment and model
access; then actual fixed-checkpoint runs before a T05 feasibility decision.

## Preserved earlier evidence

The 2026-10-29 upstream status records a researcher-reported ESM-2
`facebook/esm2_t6_8M_UR50D` probe on public GFP UniProt `P42212`: 238 amino
acids and output shape `[1, 240, 320]`. That prior report was not independently
reproduced in this implementation session. It remains preliminary evidence,
with no annotation-family benchmark or T05 completion implied; the local source
still needs its record-specific T02 review evidence. Upstream T01/T03/T04
approvals are dated 2026-10-29; T02 approval is dated 2026-10-28 and reaffirmed
2026-10-29.

## T06-T08 feature-branch progress

The [candidate inventory](T06_candidate_inventory.md) and
[ten-record queue](../../data/manifests/t08_pilot_review_queue.json) are discovery
metadata only. Source versions, permissions, checksums, full annotations and
separate label/independence decisions remain pending or uncertain. The
[T07 protocol](T07_benchmark_protocol.md) predefines budgets, valid protein
baselines, claims, leakage exclusions and paired uncertainty. The
[T08 pipeline](T08_pilot_pipeline.md) preserves T03 v0.1 and original qualifiers
without source approval or unsupported feature inference.

Actual feature-branch checks: inventory/queue metadata validation succeeded;
preparation and A returned `no_inputs`, and B-D returned
`blocked_no_reviewed_inputs` with no source parsing/inference. No real pilot score
or plot was produced. Software fixtures do not replace those experiments. Tests
and limitations are recorded in the [report](T06_T08_implementation_report.md).

New contributions are confined to `codex/real-construct-pilot`, based on
`origin/main` at `f39666d`. The main worktree and draft PR #1 were not modified.
PR #1 was inspected read-only at `7b15f02`; its artificial-combination results
inform controls, not present-day real-construct claims. No task approval or
dataset admission is inferred from code completion.
