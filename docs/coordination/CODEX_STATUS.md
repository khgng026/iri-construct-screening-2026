# Codex implementation handoff

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
