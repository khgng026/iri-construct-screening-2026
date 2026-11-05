# Running the shared core experiments

Use Python 3.11 or 3.12 on macOS or Linux. On Windows use WSL2/Linux; the historical runners use POSIX process groups and are not native-Windows compatible. CPU is the supported default. No CUDA configuration has been validated.

## 1. Get the research branch and test dependencies

```sh
git clone --branch rudra/research-handoff-2026-11-05 https://github.com/khgng026/iri-construct-screening-2026.git
cd iri-construct-screening-2026
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -e . -r requirements-checks.txt
python scripts/smoke_test.py
python -m unittest discover -s tests -v
```

The checks use NumPy, SciPy, Biopython, psutil and PyTorch (toy tensors only; no model weights). On Linux, install the CPU-only PyTorch wheel from the official PyTorch CPU index first if you do not need CUDA. Default tests deliberately omit the three private-data integration checks and the optional DCT upstream-source fixture when that third-party source is absent. macOS-only compression tests skip on Linux. Missing numerical dependencies are errors, not silently hidden skips.

## 2. Install the exact reviewed core data

Download [the 1.8 MB reviewed input archive](https://drive.google.com/file/d/1mDn2pI4NjCNAyw1wgMRaKWNO1naOQWwW/view) from the shared project Drive using your signed-in account. This contains only the approved original core inputs and source/license evidence, not model weights or confidential material. Attribution and release scope are in `config/reproduction/DATA_NOTICE.txt`; historical manifests are unchanged.

```sh
python scripts/reproduce.py install-inputs ~/Downloads/iri-reviewed-core-inputs-2026-10-31.zip
IRI_RUN_DATA_TESTS=1 python -m unittest discover -s tests -v
```

The importer checks the whole archive, exact membership, each file checksum and existing destinations before writing. It never overwrites conflicting data. The subsequent suite check verifies all 105 records, 101 unique sequences, 1,146 queries, frozen metadata and boundaries. Opting into integration tests with missing or corrupt data fails rather than skipping. No symlinks to Rudra's workspace are needed.

## 3. Create a fresh bounded session

Connect power and keep at least 20 GiB free. The guard monitors resource pressure, swap growth, power and a finite time limit. Linux uses psutil measurements; unavailable measurements fail rather than silently bypassing the guard. Do not reuse historical launch files or queues.

```sh
python scripts/reproduce.py init results/my-core-run --hours 2
python scripts/reproduce.py run results/my-core-run --model alignment --limit 17 --timeout 300
```

The first 17 sorted queries include clean controls and a joined query. `--limit 0` selects all 1,146 queries. Outputs and checksummed completion records go under your new session's `jobs/`; guard logs go under `scripts/logs/`. Existing output identities are validated before resume. A fresh session must use a new directory. Session deadline and individual timeout are both enforced.

## 4. Run the five protein models

Install `requirements-models.txt` in the environment (the historical model-interface versions are pinned there). This installs model libraries but does not download weights. Each model's immutable revision and checkpoint hash are in `config/reproduction/models/`. Download only the model you intend to run; the five historical weight payloads together were about 6.84 GB, and installations/caches require additional space. CPU performance depends on hardware. Do not expect bit-identical scores across platforms/library kernels; outputs record the configuration and environment.

```sh
python -m pip install -r requirements-models.txt
python scripts/reproduce.py download-model esm2_t6_8M_UR50D-cpu
python scripts/reproduce.py run results/my-core-run --model esm2_t6_8M_UR50D-cpu --limit 17 --timeout 600
```

Available configuration names:

- `esm2_t6_8M_UR50D-cpu`
- `esm2_t30_150M_UR50D-cpu`
- `esm2_t33_650M_UR50D-cpu`
- `ESMC-300M-cpu`
- `ESMC-600M-cpu`

Download each separately, then use its name with `run`. Increase the bounded session duration and timeout based on a small run before selecting `--limit 0`. Inference is offline after the pinned files are cached. No ESM3 experiment is claimed.

## Scope of portability

This entry point reproduces the core whole-query, oracle/control and raw-alignment experiment over the frozen core catalogue. Shared CSVs and figures cover the wider historical work. The 1,000-reference expansion, optional DCTdomain/pLM-BLAST installations and historical analysis queues are **not** turnkey cross-platform workflows in this entry point. Their source code is shared, but they need their own input/configuration admission and optional dependencies; do not mistake their availability for a completed validated method. Native macOS compression is an optional storage utility, never required to run the core benchmark.

Historical report paths and historical hashes are retained as evidence. New execution uses repository-relative metadata and fresh session paths. Original results are not changed by this portability work.
