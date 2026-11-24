# Results

This directory is an **ignored working area**: only this README is tracked
(see `.gitignore`). Preprocessing bundles, model caches, score matrices and other
run artifacts are written here locally and are not committed.

Committed, source-free results live under version control elsewhere:

- Toolkit release evidence: [docs/project/RELEASE_ASSESSMENT.md](../docs/project/RELEASE_ASSESSMENT.md).
- Component-coverage study and its extension:
  [docs/results/component-study-v1/](../docs/results/component-study-v1/README.md).
- Historical E1–E4 pilot and panel v2:
  [docs/EXPERIMENT_STATUS.md](../docs/EXPERIMENT_STATUS.md) and
  [docs/results/](../docs/results/).

Each committed report records the model version and revision, input/panel
identities and hashes, denominators, failures and the compute environment. Raw
biological inputs, embeddings and weights are never committed.

## Sequence-free run outputs (cohort v2)

With the project owner's approval (23 November 2026), the sequence-free outputs of these runs are committed here despite the ignore rule:

- `paired-cohort-v2/`: the cohort freeze, evaluator answer keys, group manifests, construct map and file hashes, the final development traces (`dev-run3`, `dev-fusion`), and every holdout trace and analysis report. Two superseded development runs from before the pre-lock fixes are not included. See [cohort study](../experiments/paired_preprocessing/COHORT_STUDY.md) and [results](../experiments/paired_preprocessing/COHORT_RESULTS.md).
- `paired-corpus-batch-v1/`: the PR #19 acquisition, audit and shortlist metadata, reproduced in this run.
- `paired-reporter-pilot-v0/prepared/evaluator/answer-key.json`, `paired-reporter-characterizers-v2-repro/` and `cohort-runner-regression/`: the reporter reproduction and runner regression traces.

Each file was scanned before committing and holds no nucleotide or protein sequence. Traces record loci, residue counts, predictions and family rankings. Files named `*.private.json` keep their original names so that the code paths still work. Construct FASTA inputs, reference sequences, annotation-derived and documented-region protein units, and GenBank and UniProt snapshots remain untracked. NCBI notes that submitters may hold rights in deposited records, so those sequence files are reacquirable from the pinned accessions and hashes instead. Because the holdout answer keys are now public, any confirmatory retest needs fresh held-out groups.

## Q3 follow-up

The four JSON files in [q3-followup](q3-followup/) record the aggregate
holdout feasibility audit, the matched-domain sensitivity replay, the Open
Reporters source audit, and the historical rank combination. See the
[methods, limitations and replay commands](../docs/project/Q3_FOLLOWUP.md).
The new source audit outputs contain identifiers, hashes, label metadata and
candidate reach summaries. Raw records and protein sequences are kept in
the source checkouts.
