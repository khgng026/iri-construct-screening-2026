# E1–E4 execution assessment — November 10, 2026

**Outcome: new E1–E4 experiments are blocked.** The repository does not contain
the completed approved readiness packet or panel-bound scientific clearance
referenced in the execution request. Existing implementation passed software
checks, but there is no newly cleared scientific experiment to execute.

This assessment records the executor's findings, not a scientific approval. Dates
in this title use America/New_York; machine timestamps use UTC.

## Authorization and evidence

Inspected main: `4b032feda0672dc2903246ee897a034d665653fc`.
Inspected proposed governance: `cb31db46bc56cf0eab600b5058ef82dfc400b4c8`.
The exact source-document digests and missing-path checks are in
[assessment.json](assessment.json). The remote branch inventory, PR listing and
PR #9 discussion were also inspected for a newer packet or clearance.

- The researcher's [recorded scope authorization](https://github.com/khgng026/iri-construct-screening-2026/blob/cb31db46bc56cf0eab600b5058ef82dfc400b4c8/docs/oversight/2026-11-10_RESEARCHER_EXPERIMENT_AUTHORIZATION.md)
  states that PR #9 must merge and a dated `CLEAR_FOR_E1/E2/E3/E4` decision must
  identify the exact panel before new scoring. It explicitly says that no new
  experiment has yet been cleared. The current execution request also directs
  the executor to report experiments lacking scientific clearance as blocked.
- [PR #9](https://github.com/khgng026/iri-construct-screening-2026/pull/9)
  was still **draft/open, unmerged** at the saved observation time. Its latest
  researcher-authorization comment distinguishes permission from scientific
  clearance. See [the saved GitHub snapshot](governance-pr-snapshot.json).
- The latest [readiness assessment on that branch](https://github.com/khgng026/iri-construct-screening-2026/blob/cb31db46bc56cf0eab600b5058ef82dfc400b4c8/docs/coordination/EXPERIMENT_EXECUTION_READINESS_2026_11_10.md)
  requires source admission, family/reference evidence, a prospective panel and
  model/resource validation. It is a list of unmet requirements, not an approved
  experimental packet. Neither inspected tree contains
  `docs/coordination/EXPERIMENT_APPROVAL_PACKET.md` or the proposed
  `docs/oversight/EXPERIMENT_RUN_AUTHORIZATION.md` clearance record.
- Active [T02 v1.0/v1.1](../../week1/T02_data_governance.md),
  [T07 v1.0](../../week1/T07_benchmark_protocol.md) and
  [T07 v1.1](../../week1/T07_v1_1_evaluation_addendum.md) retain their source,
  independent-reference and pre-score freeze conditions. The proposed v1.2
  workflow has not become active through this assessment.
- The [11-record GenBank audit](../../archive-review-2026-11-06/README.md)
  records acquisition without admission. The
  [ten selected iGEM records](../../igem-review-2026-11-10/README.md) likewise
  lack complete per-source admission and reviewed independent references.
  Previously admitted pBR322 is a separate historical source-only E1 release.

## Experiment outcomes

| Stage | Current request | Required evidence still absent |
| --- | --- | --- |
| E1 | **BLOCKED**, no new experiment | Panel-specific scientific clearance and exact-source T02 admission/freeze for the new subset |
| E2 | **BLOCKED**, no inference | Clearance, reviewed parent-family labels, independent positive and competing references, frozen exclusions and panel |
| E3 | **BLOCKED**, no alignment or model scores | Clearance, reviewed reference tiers, exact model snapshots and available model runtime; strict tier 3 also lacks its reviewed relatedness policy |
| E4 | **BLOCKED**, no cost/quality comparison | Clearance and actual comparable E2/E3 performance/resource measurements |

For all four stages, new performance metrics and experimental compute costs are
null. There were zero biological model calls. A blocked stage is not a measured
model failure, zero accuracy, or a scientific result. No model ranking follows.

The prior [pBR322 E1 report](../pilot-2026-11-06/E1/report.json) remains unchanged:
48/48 observable features and qualifiers preserved, 3/3 source translations
validated, zero erroneous inventions. Its original measured wall time was
0.4881994 seconds and sampled process RSS was 53,567,488 bytes on its original
Windows environment. These are historical annotation-audit measurements, not a
new run or model performance/cost measurements in this cloud environment. The
source has 50 total features; two source-bookkeeping features are outside the
48-feature observable denominator. One backbone-heavy record does not establish
multi-archetype performance, biological function or retrieval generalization.

## Technical validation and costs

The new [validation harness](validate.py) ran existing checks without downloading
models or scoring biological inputs. Production parser/scorer code and frozen
protocols needed no change based on these checks. Further experimental
implementation depends on a cleared panel and configuration; no unsupported
source, label, reference, model snapshot or threshold was supplied to fill a gap.

| Check | Observed outcome |
| --- | --- |
| Smoke | Passed, fictional metadata |
| Full regression suite | 187 run, 182 passed, five optional skips, no failures/errors |
| Historical release verification | Passed: exact original pBR322 bytes/translations and published artifact hashes verified; historical HOLD boundaries retained |
| `pip check` | Passed, no broken installed requirements |
| Runtime | Python 3.12.14; Biopython 1.86, NumPy 2.5.4, psutil 7.1.3 |
| ESM-2 offline load | **Blocked**, missing model dependencies |
| ESM-C | Unavailable: ESM SDK absent; no validated checkpoint or execution |
| Model stack | PyTorch, Transformers and Hugging Face Hub absent; GPU availability not established |

The five command invocations took **1.9973 seconds wall time** in total and
**2.4216 child CPU seconds**. The OS-reported child-process peak RSS was
**40,112,128 bytes**, a high-water mark across checks, not an additive or
per-model measurement. These costs exclude repository inspection, tool overhead,
and evidence preparation. GPU time and monetary cost were not measured and
remain null. The offline-load failure itself took 0.00583 seconds with a sampled
16,883,712-byte process RSS; that is failed preflight work, not inference cost.

Raw logs, command arguments, source-code hashes, package inventory, resources and
timestamps are saved under [technical-validation](technical-validation/validation.json).
The environment command exits zero after recording its findings even when
model loading is blocked; its
[nested load status](technical-validation/environment/report.json) controls the
availability conclusion. Passing software checks does not establish scientific
readiness. Optional model/plot tests were skipped; no model environment was
installed or checkpoint substituted in the absence of a cleared experiment.

## Reproduction and remaining conditions

From the repository root, use a new output directory:

```bash
.venv/bin/python docs/runs/execution-assessment-2026-11-10/validate.py \
  --output results/execution-assessment-2026-11-10/rerun01
```

The harness uses the selected interpreter and standard POSIX resource accounting.
It reproduces technical checks, not the external scientific review or live PR
state. The committed validation artifacts are byte-preserving copies of the
original ignored outputs, indexed by their original `checksums.json`.
[artifact_checksums.json](artifact_checksums.json) covers this complete assessment.
Wall time/RSS vary across machines; re-check live governance state before any
future execution.

Execution requires an active delegated-review policy, a completed source/panel
evidence packet and a separately recorded dated scientific clearance bound to
that exact panel. E2/E3 additionally require reviewed independent families and
references; E3 requires pinned available models and supported reference tiers;
E4 requires real compatible measurements. These remain unresolved evidence
conditions. This assessment neither merges PR #9 nor creates an approval record,
and does not request another routine approval from the researcher.
