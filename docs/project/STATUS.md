# Current project status

Updated November 17, 2026. This is the only active status ledger. [Scope](SCOPE.md),
[autonomy](AUTONOMY.md), [contract](TOOLKIT_CONTRACT.md) and [plan](PLAN.md) are
authoritative; the researcher's most recent instructions take precedence over
them. Historical reports keep the dates and limits of their own observations.
All review recorded here is agent self-review.

## Where things stand

| Area | State |
| --- | --- |
| Toolkit U1–U5 | **Done.** Installed API and `iri-constructs preprocess` CLI, output bundle, fixtures, tests, clean-wheel gate. [Release assessment](RELEASE_ASSESSMENT.md): R1–R7 pass on Python 3.10–3.13 (Linux); Windows and macOS not run |
| Research U6 | **Done as a pre-registered study plus one pre-registered extension.** [Component-coverage study v1](../results/component-study-v1/README.md) and [extension 1](../results/component-study-v1/extension-1/README.md): specifications and labels frozen before scoring; ESM-2 8M and alignment scored; the extension replicates v1 on unseen unit sets |
| Documentation | README shows working commands; superseded directories carry notices; this ledger and HISTORY.md are current |
| Open pull request | [#13](https://github.com/khgng026/iri-construct-screening-2026/pull/13), draft, branch `claude/gifted-gauss-i9jbhz`. Not merged; merging needs the researcher |

## Evidence

- **Toolkit.** `python scripts/verify_wheel_install.py` passes: wheel built,
  installed outside the checkout, imported from `site-packages`, no model
  package installed, CLI exit codes 0 and 2 as specified, 14 toolkit tests rerun
  against the installed copy with sockets disabled (13 pass, 1 checkout-only skip).
  Repository suite: 271 tests, 266 pass, 5 optional historical skips. The 251
  tests that existed before this work pass unchanged;
  `scripts/verify_pilot_release.py` passes.
- **Real records.** Two disjoint seeded samples of NCBI vector records, 1,000 in
  all: 1,000 processed, 0 failed, no parser warning; 2,338 of 2,491 CDSs
  exported (93.9%); 153 withheld (138 partial boundaries, 12 exceptional
  annotations, 3 non-canonical translations); 47 exported proteins over the
  1022-residue model limit.
- **Study v1** (110 distinct multi-protein parents, 185 supported components,
  five hits per parent; components recovered through their own unit):

  | Condition | ESM-2 8M | Alignment |
  | --- | --- | --- |
  | Longest protein, label-blind | 55 | 59 |
  | Longest supported protein | 106 | 110 |
  | All proteins, maximum-score merge | 151 | 135 |
  | All proteins, round-robin merge | 178 | 184 |
  | Per-unit ceiling (larger budget) | 180 | 185 |

  The gain is coverage from asking more queries (3.2× model calls) and depends
  on the merge rule. Per-unit retrieval is identical across conditions. Alignment
  matches or exceeds ESM-2 8M per protein and under round-robin; under
  maximum-score merging it is worse, because raw alignment scores grow with length. Whole fusion proteins are where
  mean-pooled ESM-2 fails (3 of 8 family slots against 7 of 8). Six dependence
  clusters, so no population claim. Model revision `c731040f…5df`, weight
  SHA-256 `24c5fa47…4189`, identical to the historical runs.
- **Extension 1** (363 distinct unit sets not seen in v1, 614 components): the
  same five rows give 200, 352, 509, 588, 599 for ESM-2 and 209, 362, 474, 605,
  613 for alignment. With 20 dependence clusters the paired cluster bootstrap
  gives round-robin − longest supported = +0.12 [0.04, 0.20] and maximum merge −
  longest supported = +0.05 [−0.08, 0.17] (ESM-2). Fusion-class family slots
  over both samples: ESM-2 28 of 53, alignment 50 of 53; an exploratory oracle
  cut to the anchor-aligned segment recovers 53 of 53.

## What is not established

- That preprocessing changes how the model represents any individual protein.
  It cannot under this design; the study measures coverage and aggregation.
- Any population-level rate. Components recur across vectors; the 78 distinct
  supported sequences are variants of 30 anchor proteins.
- Behavior on plasmid-editor exports with missing or inconsistent translations
  (fixtures only).
- Anything about noncoding content or the 230 exported proteins without a
  family label.

## Actual blockers

None for in-scope work. Private inputs of the v1 run (records, UniProt
responses, score matrices, weights) exist only in the ignored cache of the
environment that ran it; a later environment must re-acquire them, and a changed
provider snapshot is a new input that the freeze will refuse rather than accept
silently. The historical E2/E3 private caches remain unavailable, as before.

## Next useful task

Decide between the two follow-ups the v1 results motivate (details and
recommendation in the study README and the pull request):

1. **Fusion proteins.** The oracle diagnostic shows the loss is entirely due to
   pooling over unrelated sequence and is fully recoverable in principle. The
   open question is whether boundaries can be found *without* the answer key
   (source domain annotations, or a residue-level/local comparison in the
   manner of EBA/PEbA) and whether that beats plain local alignment, which
   already recovers 50 of 53. This is domain-aware decomposition: a scope
   expansion that needs the researcher's decision.
2. **A harder, larger reference tier** so that alignment is no longer at ceiling
   and the model can be distinguished from it.

Smaller in-scope items: confirm the hosted CI result of the new wheel job;
obtain a modest public sample of plasmid-editor exports under suitable terms to
replace fixture-only coverage of missing translations.

## Reproduce

```bash
python -m pip install -e . && python -m unittest discover -s tests
python scripts/verify_wheel_install.py
```

Study commands are listed at the end of the
[study README](../results/component-study-v1/README.md). Environment of this
checkpoint: Linux x86-64; toolkit Python 3.10–3.13 with Biopython 1.86 (dev) and
1.88 (wheel gate); model environment Python 3.11, torch 2.3.0 CPU,
transformers 4.48.1, huggingface-hub 0.28.1, numpy 1.26.4.

## Commits of this stretch

`c27d261` normalization core split · `10b8fcf` API and CLI · `2814371` fixtures,
tests, wheel gate · `17ea54b` study specification and panel freeze (before
scoring) · `e8e4a9d` study results · `f12d8cf` documentation, CI, notices ·
`b2fe310` extension specification (before fetching) · later commits: extension
freeze and results.
