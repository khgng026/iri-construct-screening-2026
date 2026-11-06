# Rerun the fixed-model composition pilot

Work from `/Users/Rudra2.0/Desktop/IRI/iri-construct-screening-2026`.
Local branch: `rudra/t05-composition-pilot`. The original experiment remains on
`rudra/t05-esm2-benign-pilot`, and its files/results remain in this working tree.
Nothing was pushed, merged or sent to anyone.

## Existing local environment

```bash
cd /Users/Rudra2.0/Desktop/IRI/iri-construct-screening-2026
../.venv-t5/bin/python -m pip check
../.venv-t5/bin/python -m unittest discover -s tests -v
../.venv-t5/bin/python scripts/smoke_test.py
```

For environment recovery, use Python 3.12 and install
`requirements-t05-expanded.lock.txt` in an isolated environment. It includes the
unchanged T05 numerical packages plus Matplotlib/Biopython. Do not upgrade the
model packages mid-comparison. No GPU or account token is needed for these runs.

## Run and verify

Use new names: the runner refuses to overwrite an existing run directory.

```bash
../.venv-t5/bin/python scripts/t05_composition_experiment.py run --run-id my-full-01
../.venv-t5/bin/python scripts/t05_composition_experiment.py run --run-id my-repeat-01 --subset --reencode-catalogue
../.venv-t5/bin/python scripts/t05_composition_experiment.py compare my-full-01 my-repeat-01
../.venv-t5/bin/python scripts/t05_composition_analyze.py my-full-01
```

The first command validates frozen hashes and local sources, loads the fixed model
offline, uses or creates a validated catalogue cache, and separately encodes every
clean, whole and oracle input. `--reencode-catalogue` bypasses cache reuse and leaves
existing cache files intact. Each command launches a separate Python process.
`--subset` uses the previously frozen 16 clean queries plus both orders of four
pairs: TRX1–RPS28A, ENO1–HSP60, SSA1–RPS28A and PRE1–TSA1. No score-based selection.

The comparison requires identical reference rankings and score/vector absolute
differences ≤1e-6, with matching model/input/code/environment provenance. New code
or a different platform can therefore require a separately documented validation;
do not hide a mismatch by increasing tolerance. Clean-query failures are recorded
without silently removing them. A failed run retains its stage/error-type record.

`analyze` writes descriptive grouped tables and three figures in PNG, PDF and SVG.
It recomputes the complete global sequence-identity audit and explicitly labels
the duplicate-collapse sensitivity as post hoc. Re-running analysis regenerates
derived analysis files for that run, without touching its recorded model outputs.

## Files and recovery

- Fixed selection/model/query rules: `config/t05_expanded_selection.json`,
  `config/t05_model.json`, `config/t05_expanded_queries.json`.
- Provenance and approvals: `data/manifests/t05_expanded_catalogue.json`.
  `config/t05_expanded_freeze.json` pins their bytes, design and dependency locks.
- Local original HTTP response pages, extracted source records and licence:
  `data/raw/t05-expanded/`. An extracted-record hash is distinct from its original
  HTTP page hash. Both are retained; release IDs unavailable from EBI remain null.
- Local canonical FASTA and artificial query FASTA: `data/processed/t05-expanded/`.
  Query-to-component sequence references and known coordinates are also saved in
  each run's `sequence_references.json`; no guessed biological annotation.
- Content-addressed catalogue cache: `results/t05-expanded/cache/`. Keys include
  model revision/file hashes, sequence hashes, encoder code, settings and numerical
  runtime versions. Array checksums catch corruption; invalid caches fail explicitly.
- Per run: `results.json` contains provenance, settings, timings, memory, complete
  rankings, components and failures; `rankings.csv` is the full ranking table;
  `components.csv` contains ranks, reciprocal ranks, scores, changes, ties and hits;
  `embeddings.npz` contains pooled vectors (no pickled objects).
- Observed evidence: `expanded-01`, `repeat-fresh-01`, `cached-02`, their comparison
  JSONs, and `preflight/` checks. `expanded-01/analysis/summary.json` contains the
  aggregate results; companion CSVs cover targets, bins, positions and reversals.

Source acquisition is deliberately separate from scoring. The acquisition script
`scripts/t05_expanded_acquire.py` resolves the prespecified genes through the public
official API; it does not grant review approval or update frozen manifests. Current
source annotations may differ from frozen bytes. Keep the original local raw files
for exact reruns; if missing, retrieve and verify the exact expected bytes before
scoring. Changed/missing sources must be reported, not silently replaced or accepted.
The original T05 preparation route can restore the pinned public model cache.

Do not commit raw/processed sequences, downloaded annotation pages, embeddings,
cache files or result directories. Git stores code, metadata, design and reports.
A Git commit saves a local checkpoint; a push uploads it; a pull request proposes
merging it elsewhere. This task uses only the local checkpoints.
