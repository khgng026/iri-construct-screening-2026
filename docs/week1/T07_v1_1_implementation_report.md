# T07 v1.1 implementation and validation report

Evidence date: 2026-11-06 UTC. Feature branch: `codex/t07-evaluation-extensions`.
Base implementation: PR #2 at `2c9a5bf6e157b05326199a6de0c345598cb22035`.
The researcher merged PR #2 upstream at 2026-11-06T00:40:46Z as `284a700`;
its merged tree is identical to that base. The new PR therefore targets main.
This work made no merge, force push or change to main, PR #1/#2 heads or oversight.

## Delivered work and ownership

Three actual delegated agents saved their implementations/tests in separate
files. All three subsequently hit the workspace spend cap before final handoff.
The lead reviewed their saved files, completed integration and corrections,
ran all reported checks, and published the milestones. No simulated agent
participation or completed agent handoff is claimed.

| Owner | Delivered contribution | Lead integration |
| --- | --- | --- |
| Agent A | E1 annotation-quality strata, source-preservation/abstention metrics, selection/queue validation and future expansion plan | Explicit evidence classification, fixture nonpromotion, cap propagation and full regression review |
| Agent B | E2 annotated protein-unit loading, longest/all-CDS parent panels, deterministic max merge, matched-unit controls and phase measurements | Parent-cluster uncertainty, saturated-panel exclusions and E4 record identity/resource alignment |
| Agent C | E3 reviewed reference tiers and E4 strict cost/quality reporting, optional plotting and focused tests | Alignment/model execution on common tiers, frozen model identities, exact scored support, ambiguous-label exclusions and common CLI |
| Lead | Versioned configs/addendum, orchestration, cross-module tests, Git milestones and delivery | Frozen-contract checks, documentation/workbook, source/artifact review and publication |

E1 describes agreement with observable source annotations, with per-record,
archetype and quality-stratum denominators. Missing/unsupported evidence causes
abstention; a source-only record cannot become perfect annotation agreement.
No inferred CDS, ORF completion or biological function is introduced.

E2 selects the longest eligible source translation with stable ID tie breaking;
the other condition represents every eligible reviewed CDS independently.
Maximum score per reference precedes one deduplicated five-hit parent ranking.
Full documented expectations and independently supported family subsets remain
separate. Matched-unit controls prevent additional component coverage from being
misrepresented as improved per-protein accuracy. Parent-cluster paired intervals
require eight independent groups and exclude catalogue-saturated observations.

E3 retains every original T07 independence safeguard in all tiers. Additional
linkage and stricter sequence-relatedness exclusions require predeclared reviewed
pair evidence bound to exact fingerprints. Missing strict policy blocks tier 3.
Raw/reference-normalized local alignment, ESM-2 and ESM-C use identical per-tier
panels; missing model support cannot silently become a successful comparison.

E4 retains actual/null measurements, source/panel/protocol/runtime identities,
cache effects and exact scored support. Incompatible records remain descriptive.
Shared execution IDs disclose reused costs; there is no automatic model winner.
Plotting requires actual admitted comparable measurements and an optional plotting
dependency. No scientific plot was produced in this implementation.

## Preserved contracts and candidate selection

T02, T03 v0.1/schema/examples, the entire `config/pilot/` freeze, the written
T07 v1.0 protocol and original twenty-candidate inventory are unchanged.
The original queue is archived in `t08_pilot_review_queue_v1_0.json`; the current
queue follows the exact ten approved identities, including pGreen2x35SLUC.
pUC19c is only an optional annotation-limited control after its own admission.
Candidate selection admits zero source files, annotations or labels.

`config/extensions/design_freeze.json` independently hashes the v1.1 addendum
and E1-E4 designs using UTF-8/LF text normalization. Biological checksums remain
raw-file checksums. Every actual runner freezes its chosen protocol before source
admission and its panel/model identity before scoring/forward execution.

The workflow workbook changed only `Weekly To-Dos!L10`, `L11`, `L12` notes.
All other cells and ZIP members were verified unchanged; task status cells,
approval dates and prior notes were retained. T05/T06/T07 remain in progress
where real-data assessment is incomplete; T08 remains blocked.

## Validation actually executed

Local environment: Python 3.14.7, Biopython 1.86, CPU-only PyTorch; no CUDA,
ESM SDK or Matplotlib. Model-like outputs in tests are fictional numerical
adapters, never pretrained research measurements. No weights were downloaded or
passed real inputs by these checks.

| Check | Observed result |
| --- | --- |
| Fresh unchanged PR #2 baseline | 116 tests passed |
| Final full suite with installed optional parser packages | 187 tests run: 186 passed, one skipped |
| Dependency-free suite (`python -S -m unittest discover -s tests`) | 187 tests run: 143 passed, 44 optional-dependency skips |
| Dependency-free smoke check (`python -S scripts/smoke_test.py`) | Passed; three fictional schema examples |
| Metadata inventory plus explicit v1.1 shortlist/queue validation | Exit 0; metadata only, zero admitted sources |
| E1 pending CLI | Exit 1, `no_inputs` |
| E2 pending CLI | Exit 1, `blocked_no_reviewed_inputs` |
| E3 pending CLI | Exit 1, `blocked_no_reviewed_inputs` |
| E4 pending CLI | Exit 1, `blocked_no_measurements` |

The full-suite skip is optional Matplotlib fixture rendering. Standard-library
execution also skips tests needing optional parsing/numerical packages. Tests
cover source gates/checksums, annotation uncertainty/invention, missing/duplicate
families, deterministic longest/max/tie behavior, empty parents, strict five-hit
budgets, score/panel tampering, same-panel reference tiers, strict policy binding,
unavailable methods, matched support, null costs and fixture nonpromotion.

Actual pending reports and sequence-free verification summaries are local ignored
artifacts under `results/extensions/verification-20261106/`. All four report
`inference_executed: false` and `scientific_scores_produced: false`.
There are zero admitted real records, zero real scored parents and zero measured
scientific improvements. Test fixture values are not biological findings.

The E1/E2 and E3/E4 milestone pushes both passed hosted metadata and parsing CI.
The final PR/check links and publication evidence are maintained in
[CODEX_STATUS.md](../coordination/CODEX_STATUS.md).

## Scientific limitations and first real pilot actions

Annotation preservation does not establish biological truth. Family retrieval
does not establish function, unseen-protein generalization or operational
screening utility. Parent coverage gains do not imply per-unit accuracy gains.
Small or saturated panels cannot establish method advantage. Strict relatedness
reviews and model environments have not been validated on a real admitted panel.
The alignment helper conservatively withholds a query when a base-panel reference
fails; this attrition remains visible even if a stricter tier removes that reference.
Cold initialization is distinct from cold filesystem caches, and shared inference
phases are not independent timing replicates. Missing resource fields stay null.

Before the first real scored pilot:

1. Supply permitted exact local source files and a complete T02 reviewed manifest:
   immutable version/URL, reviewer/date, licensing/reuse scope, safety rationale,
   annotation review and original-byte SHA256. Keep pending/uncertain sources held.
2. Review original CDS annotations/caveats, independently supported family labels,
   query/reference roles and relatedness groups. Retain unsupported/unknown units
   and select enough independent positives and competitors without increasing K.
3. Supply an independently justified, predeclared tier-3 relatedness policy bound
   to those exact fingerprints. Without it, that tier remains blocked.
4. Create a separate local reviewed input configuration; explicitly declare
   `admitted_research` only for admitted research inputs. Run E1/source audits and
   manually spot-check admitted records before any protein-model scoring.
5. Validate compatible model dependencies/access/memory, freeze actual checkpoint
   identities, and run the fixed E2/E3 methods sequentially, retaining failures,
   coverage, paired support and measured cache/resource context.
6. Feed recorded E2/E3 measurement batches into E4. Compare or plot only matched
   contexts; seek research review before interpreting scientific gains.

Commands and policy contracts are in the
[addendum](T07_v1_1_evaluation_addendum.md), [E1](E1_annotation_quality.md),
[E2](E2_multi_cds.md) and [E3/E4 guide](E3_E4_generalization_cost.md).
Specific human evidence is tracked in
[DECISIONS_NEEDED.md](../coordination/DECISIONS_NEEDED.md).
The [150-250-record expansion plan](FUTURE_DATASET_PLAN.md) is prospective and
conditional on assessing the pilot; no bulk acquisition was performed.
