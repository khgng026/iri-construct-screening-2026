# T06-T08 pilot handoff — 2026-11-05

This feature contribution moves candidate discovery, review admission and
experimental design forward without fabricating an approved dataset. There are
**20 public metadata candidates, ten proposed pilot records, zero admitted real
records and no real biological/model scores**. T01-T04 stay approved; T05 remains
in progress. T06 inventory and T07 design are implemented but await researcher
review; T08's pipeline is tested while scientific admission remains blocked.

## Code and artifact inventory

| Task | Files / behavior |
| --- | --- |
| T06 | `data/manifests/t06_candidate_inventory.json`, `candidate_inventory.py`, `scripts/inventory_candidates.py`, `tests/test_candidate_inventory.py`, `T06_candidate_inventory.md`: 20 source-linked candidates, explicit unknown versions/counts/checksums, source terms and exclusions; no approval-capable discovery schema |
| T08 review | `data/manifests/t08_pilot_review_queue.json`: ten prioritized records with blockers/checklist, canonical inventory hash and zero admissions; `config/pilot/inputs.pending.json` intentionally contains no scientific inputs |
| T08 / A | `pilot.py`, `scripts/prepare_pilot.py`, `tests/test_pilot.py`, `T08_pilot_pipeline.md`: reviewed DNA preparation, unchanged T03 v0.1 normalization, independent qualifier/location/strand/translation/traceability audit and per-archetype counts/rates |
| T07 / A-D | `pilot_scoring.py`, `tests/test_pilot_scoring.py`, `T07_benchmark_protocol.md`, four `config/pilot/experiment_*.json` designs and `design_freeze.json`: exact frozen panel identities, fixed total hit budgets, independent references, abstentions and paired cluster uncertainty |
| B-D execution | `pilot_experiments.py`, `scripts/run_pilot_experiment.py`, `tests/test_pilot_experiments.py`: shared design/panel/model freezes, existing fusion adapter reuse, matched retrieval, sequential model calls and truthful failure/admission statuses |
| C alignment | `pilot_alignment.py`, `tests/test_pilot_alignment.py`: local BLOSUM62 raw alignment and predeclared reference-self-normalized sensitivity, identical panel/cap/budget, recorded matrix identity and resource/call counts |
| Shared repairs | `bench_inputs.py`, `tests/test_bench_inputs.py`: reject source-record aliases and verify raw source SHA256 immediately before/after parsing; `fusion.py`: optionally expose pooled values to an in-memory sink without serializing embeddings |
| Integration | README and STATUS link the candidate/protocol/pipeline; workbook updates only T06-T08 status/note cells G10-G12 and L10-L12. Approved specs/examples remain unchanged. |

Module paths in this table are under `src/iri_constructs/`; task documents are
under `docs/week1/`. All source qualifiers, normalized biological objects, source
files, protein exports, model weights and execution outputs remain ignored local
artifacts. Only public discovery metadata, code, fictional software fixtures,
designs and documentation are proposed for GitHub.

## Source provenance and manual review

The [inventory document](T06_candidate_inventory.md) links each primary source and
records its observation limits. Coverage is four backbone-heavy, six principal
reporters, three multi-reporter and seven tagged/fusion candidates. These are
proposed labels from public descriptions, not verified feature-table assignments.
The ten-record queue covers all four categories (2/3/2/3 respectively).

Addgene material IDs do not supply an admitted immutable full-record sequence
version. Source-specific sequence reuse permissions and exact file variants need
researcher review; no login or sequence acquisition was attempted. Direct NCBI
pages presented browser challenges, so those rows explicitly identify weaker
indexed/official supporting metadata. Source counts, original CDS coordinates,
translations and file checksums remain unknown rather than fabricated.

Every candidate needs researcher decisions under the complete T02 checklist.
The discovery schema rejects `approved` values; it cannot become an input dataset
by editing a counter. Source approval, annotation review, family labels, fusion
bounds and independence groups are separate decisions. A real source can enter
the scored pipeline only through a complete reviewed T02 manifest, matching local
bytes and the existing T05 input contract. All twenty are currently excluded.

## Reuse and collaboration provenance

The feature worktree was created from `origin/main` at
`f39666dca62380dab16847034efcfd5525f92300`. Existing main supplies T02/T03
validators, normalization, qualifier audit semantics, resource metadata, model
adapters, fusion context extraction and reviewed-label/candidate eligibility.
New code extends those modules rather than importing a second research framework.

Draft [PR #1](https://github.com/khgng026/iri-construct-screening-2026/pull/1) was
read-only at `7b15f02ac929cd06a8dee76ccec0541641844ccf`: START_HERE, REPRODUCE,
the October 31 report and selected alignment/benchmark interfaces. It informed
the alignment comparator, total-budget controls, independent-family endpoint,
oracle-boundary distinctions and limits of artificial combination experiments.
No PR #1 codebase, historical results, model artifacts or input archive was copied,
and no operation altered or merged its branch. Its results are not present-day
real-construct evidence.

Three explicitly requested agents owned candidate/provenance, parsing/audits and
protocol/scoring files. Root owned shared execution, integration, regression
repairs and publication. Cross-track read-only reviews identified and fixed:
full-fusion exact-copy retrieval leakage; incorrect isolated-parent summaries;
false independence among targets sharing a fusion; no-annotation A success; and
budget-saturated panels that cannot distinguish methods. These fixes have
regression tests. Heavy pretrained inference was never launched concurrently—or
at all—because no real panel passed admission.

## Verification actually executed

| Check | Observed result |
| --- | --- |
| Existing main suite before changes | 58 passed |
| Integrated feature suite | 116 passed, no skips in the installed environment |
| Same suite with `python -S` | 86 passed; 30 optional Biopython/tensor checks skipped |
| Dependency-free smoke check | Passed |
| Three new CLI `--help` checks with `-S`, from another directory | Passed without parser/model packages |
| Fictional parsing/audit fixtures | Real Biopython parsing exercised; reverse/multipart/fuzzy/conflicting/missing annotations and traceability regressions passed |
| Numerical scoring/adapters/alignment fixtures | Fixed-panel, budget, leakage, dependency-failure, grouping and insufficient-panel tests passed; no scientific measurement claimed |
| Inventory and checked-in queue verification | `metadata_validated_pending_review`, 20 candidates, ten proposed, zero admitted; successful metadata exit code only |
| Real-input preparation and shared A pending attempts | `no_inputs`, nonzero exit, no source reads/inference/scientific scores |
| Shared B/C/D pending attempts | `blocked_no_reviewed_inputs`, nonzero exits, no inference/scientific scores |
| Workbook integrity | Valid ZIP; only Weekly To-Dos XML changed, with six intended T06-T08 cells verified |

The environment preflight was actually run: Python 3.14.7, CPU PyTorch,
Biopython/psutil available, no ESM SDK and no CUDA devices. It loaded no sequences
or model weights. Existing separate Python 3.11 model pins were not installed or
clean-validated in this assignment. No new scientific model inference or plot was
produced. Source-agreement tests on fictional records are software checks, not
real-record biological validation or an actual pretrained-model experiment.

Local evidence is in `results/t06_inventory/validated/`,
`results/t08_pilot/pending-admission/`, `results/pilot_validation_20261105/`, and
`results/t05_environment/`. No real biological source/output was fabricated or
released. Discovery queue hashes use canonical JSON; checked-in design
hashes normalize UTF-8/LF for Git portability. Biological source digests continue
to identify exact original file bytes.

## Scientific limits and next actions

The new shared B runner tests documented fusion targets, not an entire DNA
construct or the prospective multi-CDS parent-score merge. It excludes target
and fusion source parents, exact target/full-fusion copies and linked groups
before freezing one common reference panel. Shared source relationships form
conservative uncertainty clusters, and per-construct/archetype summaries refer
to the documented fusion. Known boundaries remain an annotation/oracle control.

C and D have runnable common scoring infrastructure but no admitted panel.
Model-specific cosines are never mixed. Whole/reference budgets and actual model
calls/lengths are reported separately. Missing references abstain; fewer than eight
independent groups produces null confidence intervals. Returning all eligible
references is descriptive only, not method-advantage evidence. No advantage over
alignment, functional prediction, model winner or real-construct improvement is
asserted.

The researcher must review source permissions, benign rationale, exact records,
annotation quality, family/reference definitions, boundaries and independence.
Then acquire permitted local inputs, create a separate approved T02/configuration,
run A before B, and review actual coverage/costs/limitations before C/D or a T05
completion decision. T09 expansion and later scoring/benchmark tasks remain pending.
The four-week MVP is preserved; no training, generated biology or operational
screening decision system was added.

Git delivery is a new `codex/real-construct-pilot` branch and an unmerged pull
request to main. The original main worktree, main branch and PR #1 branch are not
modified. Publication commits and CI status are available in the feature PR.
