# IRI research handoff — November 5, 2026

## Run the core benchmark

Follow [REPRODUCE.md](REPRODUCE.md) for a clean environment, verified shared input bundle and bounded CPU runs. Optional historical extensions remain separately scoped.

## Read first

[October 31 research report](docs/research/2026-10-31/REPORT.md) contains the results, limits and links. [Shared evidence manifest](docs/research/2026-10-31/SHARED_MANIFEST.json) lists copied tables, figures and audit records. Code and tests are on this branch; main is preserved pending review.

## Main finding

Five protein models completed 1,146 artificial queries built from known benign proteins. Whole-query component recall@5 ranged from 29.26% (ESM-2 8M) to 44.94% (ESMC 600M). Raw alignment achieved 82.56%; exploratory reference-self-normalization reached 100% on this exact-reference panel. On matched held-out combinations, ESMC 300M windows improved recall from 42.27% to 68.21%, without exceeding alignment.

These are component retrieval results, not biological function predictions or real-construct screening performance. Queries reuse 16 proteins. Known boundaries are oracle information; held-out combinations share proteins with development data.

## Relationship to the project plan

T01–T04 approval is preserved from shared main. Local T05 work now extends to multi-model retrieval and controls. The reviewed real-construct pilot and complete annotation-to-units parser evaluation remain open. Suggested meeting decision: prioritize a small annotated benign pilot, agree a valid protein-model baseline, and choose one local embedding comparator.

## What is shared and what remains local

Shared: committed implementation and tests, T05 documentation, October 31 report, continuation/lessons, selected aggregate CSV tables and PNG figures, historical audit summaries and checksums.

Local only: raw biological sequences, downloaded source records, model weights, embeddings, large per-query ranking stores and recovery snapshots. Six unreviewed casebook/unit-normalization draft files remain outside this branch. Nishit's confidential presentation is not included. The original report contains some local-only links and execution paths; those are evidence locators on Rudra's Mac, not portable downloads. This package is a review handoff, not a self-contained reproduction bundle. Model downloads and admitted input acquisition are needed for fresh inference.

## Review order

1. Read the report's headline table and limitations.
2. Inspect the matching CSVs and figures under docs/research/2026-10-31.
3. Review T05 design/rerun documents under docs/week1.
4. Review code and tests for the selected next experiment. Optional DCT/pLM adapters and unrun extensions are not completed scientific findings; see CONTINUE.md.

Historical validation: 370 supported tests, 369 passed and one SciPy-dependent skip subsequently checked in the compatible environment. New handoff verification is recorded separately. No new inference results are claimed.

## Initial handoff validation before portability fixes

On November 5, clean-worktree discovery ran 367 tests: 352 passed, 14 errored and one skipped. Errors are retained in [handoff test log](docs/research/handoff-tests.txt); they expose dependencies on local-only session fixtures and resource-guard files. The historical 370-test closeout used the original local environment and is not a claim that this checkout passes all tests. These initial errors are preserved as historical evidence; see REPRODUCE.md and the subsequent validation record for portability fixes.

Shared Drive entry point: https://docs.google.com/document/d/1fAH0zEmMkPr_5EqWrB3ck0R-OEZcnDk42QoMw_IkhSw/edit
