# Run T05 again

## This Mac: existing environment and downloaded inputs

Open Terminal and run:

```bash
cd /Users/Rudra2.0/Desktop/IRI/iri-construct-screening-2026
source ../.venv-t5/bin/activate
python scripts/t05_experiment.py prepare
python scripts/smoke_test.py
python -m unittest discover -s tests -v
HF_HUB_OFFLINE=1 python scripts/t05_experiment.py run --run-id my-run-01
HF_HUB_OFFLINE=1 python scripts/t05_experiment.py run --run-id my-run-02
python scripts/t05_experiment.py compare my-run-01 my-run-02
```

Use new run IDs each time: existing result folders are never overwritten. `prepare`
without `--download` checks local files and creates/checks the artificial inputs.
`run` always loads the model locally and verifies frozen source hashes. The extra
`HF_HUB_OFFLINE=1` makes the no-network intent explicit. `deactivate` leaves the
Python environment when finished.

## Recreate the environment if needed

Tested with Python **3.12.14**, macOS arm64. This Mac's ordinary `python3` is 3.14;
the experiment instead uses the bundled Python 3.12 interpreter. If the environment
has been removed, recreate it from the repository directory:

```bash
/Users/Rudra2.0/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 -m venv ../.venv-t5
source ../.venv-t5/bin/activate
python -m pip install -r requirements-t05.lock.txt
python -m pip check
```

On another machine with Python 3.12, use `python3.12 -m venv .venv` and activate
`.venv/bin/activate`. Other platforms/package builds are not validated by this run.
The environment is disposable; source files and results are separate from it.

For missing source files or model weights, run:

```bash
python scripts/t05_experiment.py prepare --download
```

This downloads only the pinned model files and the five approved source records.
Existing files are reused. Changed upstream JSON or licence bytes fail the frozen
checksum check; the command will not silently replace versions or change labels.
Original source bytes are retained locally. Fresh download after upstream changes
may need recovery of the archived inputs or a separately reviewed/versioned
manifest; exact recovery from the mutable public endpoints is not guaranteed.

## Where to look

- `config/t05_model.json`: immutable model revision, CPU settings, 1,022-residue cap.
- `config/t05_queries.json`: the fixed component recipes and expected target.
- `data/manifests/t05_benign_catalogue.json`: provenance, review, terms and checksums.
- `data/raw/t05/`: original API JSON and licence evidence, ignored by Git.
- `data/processed/t05/`: canonical FASTA, artificial inputs and normalized mappings, ignored.
- `results/t05/<run-id>/results.json`: shapes, full rankings, target ranks/scores,
  timings, environment, model/source/code hashes and failures.
- `results/t05/<run-id>/rankings.csv`: all 40 ranking rows (4 queries × 2 arms × 5 references).
- `results/t05/<run-id>/embeddings.npz`: catalogue means and query residue/mean arrays.
- `results/t05/comparison-<first>-<second>.json`: repeatability check.

`numpy.load(path)` reads the numeric `.npz` arrays without pickle. Results JSON
contains no sequences. Timing includes model setup, checks and output operations
inside `run`, but excludes downloads and the CLI's initial Python/Torch imports;
model load and forward timings are recorded separately.

## Git: what happened locally

Terminal could not authenticate a normal private-repository `git clone`. Safari
downloaded the source ZIP; its unchanged files became local baseline `imported-main`
(commit `2969374`). The upstream snapshot was
`685df3d451aada1b5bd13cd347c0144c0dd274c5`. No repository `AGENTS.md` was present.
The earlier `meeting-starter/` folder was left untouched.

- **Clone:** normally copies files plus Git history. This ZIP-based copy has local history only.
- **Branch:** `rudra/t05-esm2-benign-pilot` is our separate line of changes.
- **Commit:** a local checkpoint. `2562106` fixed the design before inference.
- **Pull request:** an online request to review/merge a branch. None was created.

`git status`, `git log --oneline`, and `git diff imported-main..HEAD` let you inspect
the changes. No remote is configured, and nothing was pushed, merged or sent.
Do not push this imported branch as if it shared upstream history. For a future PR,
first obtain an authenticated clone and apply the exported local patch
`results/t05/t05-changes.patch` on a branch based on upstream. Data/weights/results
are intentionally absent from that patch and need separate local preparation.

The original `scripts/t05_model_probe.py` is preserved as received. Use the new
`t05_experiment.py` entry point for pinned, checked and reproducible T05 runs.
