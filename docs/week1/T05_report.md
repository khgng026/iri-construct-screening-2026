# T05 observed report — 2026-10-29

**Outcome:** one existing pretrained model runs locally end to end on a fixed, reviewed benign catalogue. The scope is protein-interface feasibility and an artificial retrieval comparison; it is not the full construct benchmark. No model was trained.

## Method and provenance

Independently implemented from the supplied presentation’s conceptual method (Nishit Kalani, methods pp. 11–14): encode a whole query or a known component, mean-pool residue embeddings, then rank references by cosine similarity. His presentation uses ESM-2 650M and reports much larger evaluations. We have not inspected his code or reproduced his numerical results. The confidential PDF remains outside this repository.

Our model is `facebook/esm2_t6_8M_UR50D`, revision `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`: six layers, 320 hidden dimensions. `EsmModel` loads 7,409,081 parameters without the masked-language head/pooler. Checkpoint identity and file hashes are in each run JSON. Settings: CPU float32, four threads, seed 0, eval plus inference mode; no gradients or training.

Hardware: Apple M4 MacBook Air, 10 CPU cores, 16 GB RAM, macOS 26.6.2 arm64. MPS is available but was not used. Isolated Python 3.12.14 environment at `../.venv-t5`; PyTorch 2.15.0, Transformers 5.22.0, NumPy 2.5.4, huggingface-hub 2.4.0, safetensors 0.8.0. `requirements-t05.lock.txt` pins all installed packages including pip.

Five UniProtKB/Swiss-Prot references, selected before scoring: GFP P42212 (238 aa), yeast PGK1 P00560 (416), ENO1 P00924 (437), TPI1 P00942 (248), ADH1 P00330 (348). Source-specific entry/sequence versions, original-file and sequence SHA-256, public URLs, review rationale, reviewer identity and licence evidence are frozen in `data/manifests/t05_benign_catalogue.json`. No exclusions or accession substitutions. Canonical sequences remain unmodified. Reviewed status refers to UniProt curation; benign inclusion is a separate documented Codex review under the user’s explicit T05 instruction, not advisor sign-off.

UniProt data are attributed to the UniProt Consortium under CC BY 4.0. The licence’s other-rights disclaimer is retained. Four yeast records came through EMBL-EBI’s official Proteins API after UniProt REST errors. Their API release identifier was unavailable; exact entry/sequence versions and downloaded bytes are pinned instead. Source dates are retained verbatim. Raw annotations/sequences, derived FASTA, weights and numeric outputs remain outside Git. No hazardous or proprietary sequence material was used.

## Observed outputs

The first GFP forward pass took **0.0475 s** and returned `[1,240,320]`: 238 residue positions plus two special tokens. Saved residue array: `[238,320]`; mean vector: `[320]`. Pooling masks both padding and special tokens. The real GFP/PGK1 padded-batch check matched the single-input GFP mean exactly (maximum difference 0).

Each artificial query joins full proteins directly, without a linker or DNA translation. The oracle slices the known GFP interval and re-encodes it. These boundaries are supplied answers, not automatic segmentation. The clean input is a self-retrieval control.

| Query | Length aa | GFP fraction | Whole rank / 5 | Whole cosine | Oracle rank / 5 | Oracle cosine | Forward ms: whole / oracle |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| clean_gfp | 238 | 1.0000 | 1 | 1.000000 | 1 | 1.000000 | 9.05 / 8.68 |
| gfp_pgk | 654 | 0.3639 | 3 | 0.972122 | 1 | 1.000000 | 27.46 / 9.53 |
| pgk_gfp | 654 | 0.3639 | 3 | 0.967909 | 1 | 1.000000 | 29.22 / 9.76 |
| gfp_pgk_adh | 1002 | 0.2375 | 3 | 0.951354 | 1 | 1.000000 | 51.33 / 9.27 |

Final run totals: **2.100 s** and **1.432 s**; model loads: 0.127 s and 0.090 s. Initial pinned-model download: 2.116 s. Totals exclude downloads and initial Python/Torch imports; per-stage times are recorded. These are single-run timings, not a speed benchmark. No inference failures, truncations or dropped queries occurred.

## Checks and repeatability

All **20 unit tests pass**, including the original 11. Checks cover masked pooling, length limits, multi-record rejection, review gates, corrupted source checksums, path escapes, exact oracle boundaries, schema references and cosine ranking/ties. The metadata smoke check and package dependency check pass. An independent NumPy float64 recalculation confirms all 40 scores within 1e-6 and all eight rankings.

`run-final-01` and `run-final-02` ran in separate processes, offline, from clean implementation commit `9416138`. All **40 scores, eight rankings and 21 embedding arrays match exactly** (maximum absolute difference 0; acceptance tolerance 1e-6). Model, input, package and source-code hashes match. Earlier first-forward and development runs are retained separately.

## Limits, failures and unresolved work

- **Input contract:** canonical 20-letter protein strings only, 1–1,022 residues. Oversize/invalid inputs fail; no silent truncation. Longest observed query was 1,002 residues. The 1,022 cap is this pilot’s conservative policy, not a proven architectural ceiling. Model config says 1,026 positions with rotary embeddings; the tokenizer reports a very large sentinel rather than a useful length limit. DNA and noncoding regions are unsupported here.
- **Interpretation:** whole-query GFP ranks fell from 1 to 3 in this fixed catalogue, while the known GFP span self-retrieved at rank 1. High cosine values are not probabilities or functional/safety labels. Changing the query also changes residue context, so this does not isolate pooling as the only mechanism.
- **Generalization:** three correlated artificial examples, one target, five references, exact target inclusion and no family holdout cannot establish general performance. Oracle success is expected by construction. Model training-set overlap is unknown. Automatic segmentation, DNA parsing/translation, real construct evaluation and held-out protein/family testing remain unresolved.
- **Retrieval/setup failures:** terminal clone authentication failed; Safari source ZIP supplied the working files. UniProt REST returned HTTP 500/timeouts for the yeast accessions; the official EBI route resolved acquisition without changing examples. The model loader’s unused language-head notice is expected; a macOS distributed-log notice did not affect inference.
- **Git limitation:** local branch `rudra/t05-esm2-benign-pilot` is based on an imported source snapshot at upstream `685df3d451aada1b5bd13cd347c0144c0dd274c5`, not its 27-commit history. Original files are recoverable on `imported-main`; the original probe is unchanged. No push, merge, pull request or message occurred. Future upstream work needs an authenticated clone and the local patch.
- **Recovery:** local repeatability is verified. Fresh source downloads may fail frozen hashes when public annotations change; original bytes are retained locally. Cross-platform and GPU numerical reproducibility are untested.

## Files and rerun

Run commands: [T05 rerun guide](T05_rerun.md). Fixed design: [T05 design](T05_design.md). T3 bridge: [sequence references](T05_sequence_refs.md). Main entry point: `scripts/t05_experiment.py`. Final evidence: `results/t05/run-final-01/`, `run-final-02/`, `comparison-run-final-01-run-final-02.json` and `independent-numpy-check.json`. Local sequences and normalized mappings: `data/processed/t05/`.

## Official documentation consulted

- [Model card and MIT licence designation](https://huggingface.co/facebook/esm2_t6_8M_UR50D)
- [Transformers ESM model/tokenizer documentation](https://huggingface.co/docs/transformers/model_doc/esm)
- [Meta ESM repository: residue averaging excluding special tokens](https://github.com/facebookresearch/esm)
- [PyTorch installation](https://docs.pytorch.org/get-started/locally/), [MPS](https://docs.pytorch.org/docs/stable/notes/mps.html) and [inference mode](https://docs.pytorch.org/docs/stable/generated/torch.autograd.grad_mode.inference_mode.html)
- [Python venv](https://docs.python.org/3.12/library/venv.html)
- [UniProt licence evidence](https://rest.uniprot.org/help/license) and [EMBL-EBI Proteins API](https://www.ebi.ac.uk/proteins/api/doc/)
