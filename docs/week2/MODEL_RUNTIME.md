# Model interface validation and reproducible CPU runtime

Both fixed protein-model interfaces now load and execute **fully offline** in a
separate Python 3.11.17 environment. Each encoded the same clearly fictional
20-residue alphabet fixture twice, with finite residue-level output, verified
mean pooling, and exactly identical repeated representations. Empty, ambiguous
and over-cap inputs were rejected before a forward pass.

These are software-interface checks. No biological source protein was encoded,
no family retrieval was scored, and no new E1–E4 experiment or scientific
clearance follows from this validation. The CLI accepts no biological input file.

## Fixed identities

| Model | Exact snapshot | Weight SHA-256 | Output on fixture |
| --- | --- | --- | --- |
| ESM-2 8M | `c731040fcd8d73dceaa04b0a8e6329b345b0f5df` | `24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189` | 20 × 320 |
| ESM-C 300M legacy SDK checkpoint | `7f10b20ae75017b2dbc884070e03434515709a8d` | `323dff9fbf3fef297a74f4f18b6528e6f2e599b0bcf72b6927516804015becea` | 20 × 960 |

The original ESM-C repository request,
`EvolutionaryScale/esmc-300m-2024-12`, resolves to the provider's renamed
[`biohub/esmc-300m-2024-12`](https://huggingface.co/biohub/esmc-300m-2024-12).
Its saved model card identifies the retained GitHub/SDK-format checkpoint as a
legacy compatibility release. The newer Hugging Face-compatible checkpoint
recommended by that card was **not substituted**. Provider metadata, card/notice
digests, exact model revisions and measured checks are recorded in
[model-freeze.json](model-freeze.json). Model weights remain ignored local cache
artifacts and are not redistributed in this contribution. The captured
third-party notice describes dependency licenses; it is not a complete legal
determination of model-weight rights. Further use retains the provider's model
terms. No public model release or additional paid service is used.

The adapter now rejects a resolved immutable revision that differs from an
explicitly pinned requested revision. Regression tests cover a mismatched SHA
and a mutable/non-SHA snapshot directory. Historical protocol files and model
identities are preserved; a new scientific run must separately freeze its actual
configuration and obtain the required panel-specific review.

## Environment and commands

The parser environment remains separate. The model environment uses PyTorch
2.3.0+cpu, Transformers 4.48.1, Hugging Face Hub 0.28.1, ESM SDK 3.2.1,
NumPy 1.26.4, Biopython 1.86 and psutil 7.1.3. The full transitive package
inventory is [model-runtime.lock.txt](model-runtime.lock.txt).

One dependency defect was confirmed: `esm==3.2.1` imports `httpx` through its SDK
package at local-encoder import time without declaring it. A successful
`pip check` therefore did not establish that ESM-C could import. The additive
[requirements-model-runtime.txt](../../experiments/requirements-model-runtime.txt)
retains the historical direct pins and adds `httpx==0.28.1`; no remote inference
client is called. The model environment's final `pip check`, five runtime tests
and four adapter tests all pass.

For a fresh CPU environment with Python 3.11 available:

```bash
python3.11 -m venv .venv-models
.venv-models/bin/python -m pip install torch==2.3.0 torchvision==0.18.0 \
  --index-url https://download.pytorch.org/whl/cpu
.venv-models/bin/python -m pip install -r experiments/requirements-model-runtime.txt
.venv-models/bin/python -m pip check
```

For the first permitted weight acquisition, set `HF_HOME` to an ignored local
cache and run the fixture CLI with `--allow-download`, retaining the exact
revisions below. No scientific source data are downloaded or read by this CLI.
After caching, reproduce the validated offline command with a **new** output
directory:

```bash
HF_HOME="$PWD/results/week2/model-cache" HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 \
  .venv-models/bin/python scripts/validate_model_runtime.py \
  --models esm2,esmc \
  --esm2-revision c731040fcd8d73dceaa04b0a8e6329b345b0f5df \
  --esmc-revision 7f10b20ae75017b2dbc884070e03434515709a8d \
  --threads 2 --output results/week2/model-runtime/rerun01
```

The shared cap remains 1022 residues. Inputs are canonical proteins only,
float32 weights without autocast, deterministic algorithms, sequence-only SDK
tracks, final adapter residue states and residue-mean pooling. BOS/EOS/PAD are
excluded; token/residue counts, finite values and model shapes are checked.
DNA is not passed to a protein-only encoder. Model dimensions and embedding
spaces differ; their cosine values or coordinates are not compared directly.

## Measured fixture costs

The final offline run used CPU, two PyTorch intra-op threads and three inter-op
threads on the selected three-logical-CPU cloud environment. The models were
initialized sequentially in one process without a fixture warm-up.

| Model | Initialization wall s | First fixture wall s | Repeat fixture wall s | Highest sampled process RSS bytes |
| --- | ---: | ---: | ---: | ---: |
| ESM-2 | 3.12688 | 0.01831 | 0.00875 | 480,595,968 |
| ESM-C | 8.16608 | 0.33019 | 0.34409 | 3,115,356,160 |

Initialization includes cache lookup, weight hashing, package work and CPU
loading. The filesystem/OS cache was uncontrolled; these are not independent
cold-system benchmarks. ESM-C's initialization peak includes loading allocations;
its first-forward sampled RSS was 1,807,491,072 bytes. Absolute process RSS is
sampled every 20 ms and can miss brief peaks. Process CPU times for each phase
are retained in the machine-readable freeze/report. GPU execution time and
monetary cost were not measured and are not inferred from wall time.

Final original artifacts are in ignored
`results/week2/model-runtime/final-offline/report.json`, its environment lock,
`final-offline.log`, `final-runtime-tests.txt`, `final-adapter-tests.txt` and
`final-pip-check.txt`. The integrated run evidence copies selected source-free
artifacts for review. Fixture costs are not E4 biological cost/quality
measurements, nor evidence that either model is a better biological retriever.

T05 and T15/T16 now have validated fixed local interfaces and deterministic
software evidence. End-to-end benign-source feasibility, admitted biological
baseline/unit results, paired retrieval and E1–E4 remain subject to source,
reference and external scientific-clearance conditions. Small-fixture success
does not establish memory or latency for a complete biological panel.
