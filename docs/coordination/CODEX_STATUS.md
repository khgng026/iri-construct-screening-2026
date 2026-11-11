# Current Week 2 handoff — November 11, 2026

The autonomous mission has delivered actual source curation, deterministic
preparatory QA, parser and scoring corrections, and pinned offline ESM-2/ESM-C
software validation. PR #9 is merged; its earlier draft/open snapshot below is
historical. No new panel-specific scientific clearance has been issued.

Start with [the completion report](../week2/COMPLETION_REPORT.md),
[task status](../week2/TASK_STATUS.json),
[scientific-review packet](EXPERIMENT_APPROVAL_PACKET.md), and
[run evidence](../runs/week2-2026-11-11/README.md).
New scoring remains blocked pending the separate actual ChatGPT decision;
Codex has not signed its own scientific clearance. Source-specific HOLDs and
unsupported reference/strict-tier conditions remain visible in the packet.

---

The handoff entries below are preserved historical evidence.

# Codex implementation handoff

## Pilot execution — 2026-11-06

Branch: `codex/pilot-verified-subset-2026-11-06`. Base main:
`4352009703dcd1825e181d203edf6088c882182e`, including merged researcher sign-off
PR #5. **One source admitted for E1; nine exact-source HOLDs; zero independent
references or retrieval scores.** The researcher's approval of all ten remains
valid. The older handoff below is historical; its zero-source/no-parsing state
is superseded by this execution. The original inventory, frozen protocols,
discovery queue and researcher sign-off are preserved as historical evidence.

The [release manifest](../../data/manifests/pilot_execution_2026_11_06/release.json)
and ten per-record manifests are the current technical admission ledger.
Downloaded complete NCBI pBR322 J01749.1: 24,806 original bytes, 4,361 nt,
50 features, three documented CDSs. Original-file SHA-256:
`7660b0fe95a91d492b102192dc72924508a69fda762988aa0923b3d5cc2073f6`.
The accession version identifies the sequence revision; retrieval date/checksum
identify the complete annotated-file snapshot. All three provided translations
agree with their exact source coordinates, frame and table.

The [source review](../../data/manifests/pilot_execution_2026_11_06/records/ncbi-J01749.1.json)
documents all feature locations/types, original CDS evidence, benign cloning
context including antibiotic markers, reviewer/date and private-copy rights.
NCBI's molecular-data policy imposes no NCBI distribution restriction and no
conflicting notice was identified in this exact file. Attribution to NCBI,
Will Gilbert and the original references is retained, alongside NCBI's
submitter-IP caveat. This is a T02 v1.1 private-project hosting determination,
not a claim of CC0, patent clearance or permission for public publication.

All nine Addgene public catalogues and sequence indexes returned HTTP 200;
every sequence index explicitly requires login to download/copy sequences.
Each lists a full-file category, but none was downloaded or checksum/annotation
verified. Public metadata hashes are clearly labeled HTML response hashes,
never biological-source checksums. No login gate was bypassed, credentials used,
or generic backbone/insert accession/partial read substituted. General
noncommercial-use terms were inspected; exact-file notices and hosting scope
remain unverified. See the [access log](../../data/manifests/pilot_execution_2026_11_06/source_access.json).

| Execution | Actual outcome |
| --- | --- |
| E1, unchanged v1.1 implementation | `ok`; one well-annotated backbone-heavy source; 48/48 substantive features and qualifiers preserved, excluding two source bookkeeping features; 3/3 validated source translations; zero invention/parser failures |
| E2 | `blocked_insufficient_parent_panel`; three parsed CDSs, zero reviewed family-eligible units and zero independent references; no model calls/retrieval scores |
| E3 | `blocked_insufficient_panel`; no independent positive/competing-family panel; no alignment/model scores; no strict relatedness policy invented |
| E4 | `blocked_no_measurements`; consumed empty E3 measurements; no cost/quality plot or ranking |
| Main parsing/regression suite | 187 tests: 182 passed, five optional skips, zero failures/errors; software fixtures are not biological results |
| Offline ESM-2 preflight | Existing cached snapshot `c731040fcd8d73dceaa04b0a8e6329b345b0f5df` loaded on CPU; no forward pass; weight SHA-256 `24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189` |

E1 measures agreement with one source's annotations, not biological correctness,
protein function, model accuracy or generalization. Abstention/uncertainty
denominators are zero here, so those rates remain null. E1's report field
`scientific_scores_produced: true` means source-annotation agreement only, as
its `evidence_class` states; no pretrained inference occurred.

The isolated parser environment has Python 3.14.7, Biopython 1.86, NumPy 2.5.4
and psutil 7.1.3. Four optional PyTorch adapter tests and one Matplotlib fixture
were skipped. The separate existing global model environment has PyTorch
2.15.0+cpu and Transformers 5.22.0, no ESM SDK and no CUDA. ESM-C was not
validated; no new checkpoint downloaded. ESM-2 availability cannot supply
missing scientific labels/references.

[Run evidence](../runs/pilot-2026-11-06/README.md) includes original sequence-free
reports, freezes, resources, test log, artifact digests and an attempt ledger.
Parser/evaluation modules and protocols are unchanged from base main; original
reports retain source-code hashes and dirty working-tree status from the added
dataset/configuration. Generated T03 objects, qualifier sidecars and protein
exports remain in ignored local results.

The [independence review](../../data/manifests/pilot_execution_2026_11_06/independence_review.json)
places all three pBR322 proteins in one verified same-parent exclusion group.
Standalone protein accessions from this plasmid would be source copies, not
independent references. Product names are retained without guessing homology
labels. Addgene reporter/backbone linkages remain review flags, not validated
protein clusters. No relatedness cutoff was invented. Remaining decisions:
legitimate Addgene access and exact-file rights; evidence-supported family
definitions/independent positive and competing references; justified pre-score
strict E3 exclusions. See [DECISIONS_NEEDED.md](DECISIONS_NEEDED.md).

Only `data/approved_public/J01749.1.gb` is allowlisted, with byte-preserving Git
attributes. No raw directory, biological extension, weights, embeddings or
credentials are broadly enabled. PR #1 was inspected read-only at `7b15f02`;
its code/data/results were neither modified nor imported as fresh pilot evidence.
The new branch is for review; no merge is performed.

## Previous implementation handoff (historical)

## Final integration and publication

UTC timestamp: 2026-11-06T03:40:52Z.
Branch: `codex/t07-evaluation-extensions`.
Recorded implementation HEAD before this status commit:
`9a83686ee73538302155b124bae6119d0a259718`.
The containing status commit is identified by this file's Git history.
Current review: [PR #3](https://github.com/khgng026/iri-construct-screening-2026/pull/3),
open and unmerged, with no merge conflicts reported by GitHub.

E1-E4 are implemented, integrated and tested as a separately frozen T07 v1.1
addendum. The exact approved ten-identity selection is reflected in the current
pending queue, while source admissions remain zero. No real sequence parsing,
pretrained inference, scientific scores or scientific plots were produced.
The completed work is published for review; researcher evidence is still required
before the first real scored pilot.

### GitHub state and scope

PR #2 was merged **upstream by the researcher**, not by this implementation, at
2026-11-06T00:40:46Z. Its merge commit is
`284a7008f0420da4557eccf434cabef64696d787`, and its tree matches the original
base `2c9a5bf6e157b05326199a6de0c345598cb22035` exactly. Consequently PR #3 now
targets main directly, applying the requested post-merge retargeting rule.
The original main worktree is clean; main itself was never pushed or merged by
this agent. PR #1 stays open/unmerged at `7b15f02`; PR #2's head stays `2c9a5bf`;
the oversight reference stays `4fd8256`. There was no force push.

### Agent and lead status

| Owner | Final state | Contribution |
| --- | --- | --- |
| Agent A | Saved deliverables integrated; final handoff interrupted by workspace spend cap | E1 quality/abstention, provenance-aware selection/queue, regression tests and conditional 150-250-record plan |
| Agent B | Saved deliverables integrated; final handoff interrupted by workspace spend cap | E2 longest/all eligible CDS evaluation, maximum-score merge, fixed five-hit parent budget, coverage/controls and measured phase costs |
| Agent C | Saved deliverables integrated; final handoff interrupted by workspace spend cap | E3 independently reviewed reference tiers and E4 strict cost/quality tables/optional plots and tests |
| Lead | Complete | Shared configs/CLI and design freeze, execution orchestration, review/fixes, parent uncertainty, cross-module tests, documentation/workbook and GitHub publication |

All three were actual delegated agents. The lead finished from saved files after
the cap; no final agent handoff is claimed and no replacement agents were spawned.

### Validation and actual experiments

| Check | Actual outcome |
| --- | --- |
| Fresh unchanged PR #2 baseline | 116 passed |
| Final full suite | 187 run: 186 passed, one optional Matplotlib fixture skipped |
| Dependency-free suite (`python -S`) | 187 run: 143 passed, 44 optional-dependency skips |
| Dependency-free smoke | Passed; three fictional examples |
| Explicit v1.1 selection/queue validation | Exit 0, metadata only |
| E1 pending command | Exit 1, `no_inputs` |
| E2/E3 pending commands | Exit 1, `blocked_no_reviewed_inputs` |
| E4 pending command | Exit 1, `blocked_no_measurements` |

All four pending reports explicitly record no inference and no scientific scores.
Observed scientific denominator: zero admitted real records, zero real scored
parents. Fictional numerical adapter/annotation tests are software validation.
The local environment is Python 3.14.7/Biopython 1.86, CPU-only, without ESM SDK
or Matplotlib; no scientific runtime compatibility or model-access claim follows.

Hosted CI passed both E1/E2 and E3/E4 milestones, plus the
[final implementation push](https://github.com/khgng026/iri-construct-screening-2026/actions/runs/42157933453).
The status-only publication push has its own CI run visible on PR #3.
No code changed after the final local suite; subsequent edits record publication.

Diff review checked 40 changed files and found no raw sequence files, restricted
annotation exports, weights, embeddings or credentials. T02, T03/schema/examples,
the original twenty-candidate inventory, all `config/pilot/` and frozen v1.0 text
are unchanged. The original pending queue is archived byte-for-byte as a Git blob.
V1.1 has a separate five-file design hash freeze. Workbook verification confirms
only `Weekly To-Dos!L10`, `L11`, `L12` notes changed; task statuses, prior notes,
approval history and every other archive member were preserved.

Full file ownership, acceptance evidence, limitations and reproduction commands:
[implementation report](../week1/T07_v1_1_implementation_report.md),
[addendum](../week1/T07_v1_1_evaluation_addendum.md),
[E1](../week1/E1_annotation_quality.md), [E2](../week1/E2_multi_cds.md),
[E3/E4](../week1/E3_E4_generalization_cost.md).
Sequence-free local check reports are ignored under
`results/extensions/verification-20261106/`.

### Scientific decisions and exact remaining actions

Engineering choices under existing authority: deterministic longest-CDS ties,
maximum per-reference merging, all-component leakage exclusions, conservative
parent clusters, fixed total K=5, null uncertainty below eight independent groups,
matched scored support, no comparison of saturated panels, and measured shared
cost attribution. No source approval, function inference, tuned threshold or
scientific winner was introduced. Annotation agreement is not biological truth;
additional CDS coverage is not higher per-unit model accuracy.

1. Supply permitted exact local files and complete T02 admission rows with immutable
   source/version, reviewer/date, reuse scope, safety rationale, annotation approval
   and original-byte checksums. Hold every pending/uncertain record.
2. Review original CDS annotations/caveats, family labels, query/reference roles
   and independence groups; retain unknowns and ensure independent positives and
   competing families without enlarging the five-hit budget.
3. Provide an independently justified, predeclared tier-3 sequence-relatedness
   policy bound to those exact fingerprints; otherwise that tier stays blocked.
4. Create a separate local reviewed input config with explicit research evidence
   classification; run E1 and manual spot checks before protein-model scoring.
5. Validate model dependencies/access/memory, freeze actual model identities and
   run fixed E2/E3 checkpoints sequentially with all costs/failures/coverage saved.
6. Feed actual measurement batches to E4; interpret or plot only comparable panels
   and conditions, then request scientific review of the observed results.

[DECISIONS_NEEDED.md](DECISIONS_NEEDED.md) records those remaining evidence gates.
Expansion to 150-250 reviewed constructs is prospective, conditional on pilot
assessment, and has not started. ChatGPT does not automatically monitor this
readable handoff; the researcher can explicitly request its review.

---

The milestone entries below are historical snapshots at their stated timestamps.

## E3/E4 completion milestone

UTC timestamp: 2026-11-06T03:33:20Z. Branch: `codex/t07-evaluation-extensions`.
Recorded HEAD before this milestone commit: `c5ba9f39e66ddc2b3821b1a02e954e9de145d642`.
PR [#2](https://github.com/khgng026/iri-construct-screening-2026/pull/2)
remains the unmodified stacked base; PR #1 and main remain read-only.

E3 and E4 software is integrated. E3 freezes all reference tiers before any
scores and preserves every v1.0 independence safeguard; a missing strict policy
blocks tier 3. E4 retains null costs and rejects unmatched panels, measurement
conditions, scored support and ambiguous method configurations. Scientific
measurements and plots remain pending; there is no automatic winning model.

Agent A/B deliverables have been reviewed by the lead; Agent C's saved E3/E4
modules and tests have now been integrated by the lead. All three agents' final
handoffs were interrupted by the workspace spend cap; their saved files were
available, and no replacement agents were launched. Lead work is now final
regression, documentation/workbook updates, diff review and GitHub publication.

Focused evidence: eight reference-tier tests passed; 14 cost/quality tests ran
(13 passed, one skipped for optional Matplotlib); six cross-module tests passed.
These cover strict policy binding before inference, identical per-tier panels,
missing model common-support failure, fixture-to-E4 nonpromotion, and the common
CLI's four pending outcomes: E1 `no_inputs`, E2/E3 `blocked_no_reviewed_inputs`,
E4 `blocked_no_measurements`, all with non-success exit codes. Zero real records
were parsed or scored; all model-like scores in tests came from fictional adapters.

Files added: `reference_tiers.py`, `cost_quality.py`, `extension_experiments.py`,
the shared and standalone reporting CLIs, their focused/cross-module tests and
`E3_E4_generalization_cost.md`. E2-to-E4 normalized records are also checked.
Implemented under delegated authority: explicit pair-bound reviewed policies,
shared execution attribution, exact support checks and conservative comparison
exclusions. Source and family admission, tier-3 scientific justification and
actual comparable measurements remain the same pending researcher evidence;
no new approval or interpretation has been inferred.

## E1/E2 completion milestone

UTC timestamp: 2026-11-06T03:31:18Z. Branch: `codex/t07-evaluation-extensions`.
Recorded HEAD before this milestone commit: `fcfb882a98fbe198bc3133b796c7d6b59d3245b3`.
Related PRs remain [#1](https://github.com/khgng026/iri-construct-screening-2026/pull/1)
and the read-only base [#2](https://github.com/khgng026/iri-construct-screening-2026/pull/2).

E1 and E2 software is complete and tested. The current pending queue now follows
the ten approved identities, including pGreen2x35SLUC; the original inventory
and original queue are preserved. Selection approves zero source files.
No real source was acquired, admitted or scored, and no pretrained inference ran.

| Owner | State | Evidence |
| --- | --- | --- |
| Agent A | Files delivered; final handoff blocked by workspace spend cap | Annotation-quality module, source-preservation/abstention tests, exact shortlist/queue and conditional expansion plan |
| Agent B | Files delivered; final handoff blocked by workspace spend cap | Longest versus all-CDS execution, maximum per-reference merge, five total parent hits, matched-unit controls and measured stage costs |
| Agent C | Files delivered; final handoff blocked by workspace spend cap | E3/E4 modules and tests saved; lead review underway |
| Lead | Integration and testing | Added fixture evidence safeguards, paired parent-cluster uncertainty and saturated-panel exclusions; finishing integration without respawning capped agents |

Fresh unchanged PR #2 baseline: 116 tests passed. Initial integrated suite:
183 tests run, 182 passed, one skipped (optional Matplotlib fixture).
After the E1/E2 review fixes: 46 focused tests passed (16 E1, 18 E2, 12 inventory).
Tests use fictional records/numerical adapters; these are software checks,
not measured improvements on real constructs. No biological result denominator
exists: zero admitted records and zero real scored parents.

Files: `annotation_quality.py`, `multi_cds.py`, their CLI/tests and E1/E2 docs;
`candidate_inventory.py`, selection/queue manifests and metadata CLI/tests;
`FUTURE_DATASET_PLAN.md`; shared `extension_config.py`, four versioned extension
configs and the v1.1 addendum. E3/E4 integration/CLI publication follows next.
Original T02/T03/T07 v1.0 files and `config/pilot/` are unchanged.

Decisions under existing authority: stable longest-CDS tie breaking, maximum
per-reference parent aggregation, all-component leakage exclusions, conservative
reviewed parent clustering and explicit null intervals below eight clusters.
Source-specific admission, independently supported family/reference evidence and
the strict tier-3 relatedness policy remain pending in
[DECISIONS_NEEDED.md](DECISIONS_NEEDED.md). Future expansion stays conditional.

## Initial milestone (historical)

### Git state and timestamp

UTC timestamp: 2026-11-06T00:12:31.363348+00:00

Working branch: `codex/t07-evaluation-extensions`.
Recorded HEAD before this status commit: `2c9a5bf6e157b05326199a6de0c345598cb22035`
(the current PR #2 implementation). The commit containing each status revision
is available in this file's Git history; a file cannot contain its own commit SHA.
Related PRs: [PR #1, read-only](https://github.com/khgng026/iri-construct-screening-2026/pull/1)
and [PR #2, read-only base](https://github.com/khgng026/iri-construct-screening-2026/pull/2).
The new stacked PR will target `codex/real-construct-pilot`, then be retargeted
to `main` after PR #2 is merged by the researcher. No merge or force-push is planned.

## Executive summary

The researcher approved ten candidate identities and the T07 v1.0 design in
principle, with no biological file, permission or family-label admission implied.
Three real subagents are implementing the approved E1-E4 extensions in separate
files, with E1/E2 first. Existing T02, T03 v0.1 and frozen T07 v1.0 contracts
remain the baseline. Scientific results remain pending reviewed local sources
and independent reference evidence; implementation will proceed using explicitly
fictional software fixtures and transparent blocked outcomes.

## Agent status

| Owner | State | Deliverable and acceptance |
| --- | --- | --- |
| Agent A | Ongoing | E1 annotation-quality strata, denominators and appropriate abstention; approved identity shortlist/pending queue; conditional 150-250-record expansion plan; parser/selection regressions |
| Agent B | Ongoing | E2 longest-CDS versus every eligible CDS, maximum score per reference, five total parent hits, matched-unit controls, coverage/resource records and leakage/budget/empty-parent tests |
| Agent C | Ongoing | E3 frozen relatedness tiers with all existing independence guards; E4 strict comparable cost/quality records and optional plotting of actual measurements; null/blocked/mismatch tests |
| Lead | Integration and testing | Shared v1.1 configurations/addendum, common CLI, cross-module tests, frozen-v1 regression, milestone publication and stacked PR |

## Evidence at this milestone

- Inspected PR #2 at `2c9a5bf`, PR #1 at `7b15f02`, and oversight branch at
  `4fd825691f5441f88fd0419996689f8bb76a98d9` using read-only references.
- Read the oversight decisions/specification, PR #2 inventory, protocol, pipeline
  and implementation report; retained explicit separation of candidate selection
  and T02 source/annotation/label admission.
- PR #2's existing validation record is 116 tests passed; a fresh baseline rerun
  is being recorded before integration. No new experiment has been executed.
- No sequence acquisition, model training, pretrained inference, invented
  annotations, result figures or scientific improvement claim at this milestone.
- Source review, independent family references and a justified stricter
  relatedness policy remain unresolved scientific inputs, not implementation blockers.

## Scientific decisions and controls

Implemented under delegated authority: separate feature worktree; three-agent file
ownership; conservative fixed-panel safeguards; descriptive-only saturated panels;
no threshold tuning on test results; sequential model execution; no source approval.
The v1.1 addendum extends rather than replaces frozen v1.0. Shortlist selection
will follow the ten identities in the researcher instruction; pUC19c becomes an
optional annotation-limited control, still subject to T02.

Specific researcher actions before real scoring are tracked in
[DECISIONS_NEEDED.md](DECISIONS_NEEDED.md). This file is a readable handoff, not a
command channel; ChatGPT does not automatically monitor GitHub.
