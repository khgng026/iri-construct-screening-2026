# Project status — 2026-11-05

Current task approvals follow the user's implementation request. Earlier README,
scope/schema headings and workbook notes were stale and have been aligned with
that instruction. Approval of the scope/schema/governance is separate from
approval of any source record or label.

| Task | Status | Evidence / remaining work |
| --- | --- | --- |
| T01 | Approved | Existing MVP question and four-week scope retained |
| T02 | Approved | Policy v1.0 retained; record-level reviews remain pending |
| T03 | Approved | Baseline v0.1 and all three existing examples retained |
| T04 | Approved / implemented | Dependency-free smoke/tests retained; optional parsing is separate |
| T05 | In progress | Fusion, probe, annotation parsing and modular model comparison implemented; reviewed-data pretrained experiments remain blocked |
| T06–T28 | Not completed by this work | Helpers support later tasks without asserting dataset, evaluation, clean-model rerun or delivery acceptance |

The requested `experiments/t05_esm2_fusion/` and this status file were absent at
inspection. The existing untracked `iri_t05_model_comparison/` prototype was
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
suite passes 51 tests, including a real FASTA/GenBank parsing CLI run on generated, explicitly
fictional software fixtures, plus numerical adapter/retrieval tests. These are
software checks, not source approval or model-performance evidence. Optional
Biopython 1.86 and psutil 7.1.3 were installed in the local environment.

The dependency-free run passes 33 tests and skips 18 optional parsing/tensor tests;
the standard-library smoke check passes. All six entry points provide dependency-free
CLI help from another working directory. The optional CI parsing job is configured
but has not been run on GitHub in this session.

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
The planning workbook now records T01/T03 as Done (approved) and T05 as In progress;
only six status/note cells in the Weekly To-Dos worksheet changed.

Review remaining: implementation and separate audit contract; T02 source/licensing/
safety approvals with checksums; annotation family and independence policy;
documented fusion target boundaries; a compatible model environment and model
access; then actual fixed-checkpoint runs before a T05 feasibility decision.
