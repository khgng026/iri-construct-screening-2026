# IRI T05: ESM-2 vs ESM-C vs ESM-3 (benign construct pilot)

**2026-11-05 update:** The executable now delegates to the reviewed local-input
workflow in [experiments/t05_model_comparison](../experiments/t05_model_comparison/README.md).
Use the commands there. `--prepare`, automatic accession downloads, heuristic
family approval and automated winner selection are disabled in executable use.
The original prototype functions and the historical description below are retained
for review; they are not the current CLI or evidence of an executed experiment.
The dependency file now points to the separate pinned model environment.

This executable experiment compares **representation-based functional-family retrieval** for protein coding sequences extracted from public, benign annotated GenBank plasmid/vector records. It **does not** test direct function prediction from intact DNA, nor does it test detection of hazardous proteins.

## What it actually measures

1. Parse annotated `CDS` translations (do not infer novel ORFs) from a convenience set of benign public GenBank plasmid/vector and reporter records (including pEGFP variants, pGL3 variants, pBR322, pUC19, pUC18, and independent fluorescent/luciferase reporter constructs).
2. Extract one L2-normalized protein embedding per unique sequence per model, averaging representations over amino-acid residues. Different model embedding spaces are never mixed.
3. Leave one plasmid/construct out and classify each annotated protein by the **closest reference protein** in another plasmid, using cosine similarity. Reference labels are heuristically derived from GenBank annotation text.
4. Report retrieval accuracy, macro F1 and construct-level **eligible** trait recall. The primary `distinct_sequence` condition excludes cross-plasmid exact duplicate protein sequences. The `allow_reused_sequence` condition shows how duplicate vector sequences can inflate results.
5. Compare actual median inference time and GPU allocated memory. Report insufficient evidence, download failures, and model failures explicitly.

**Important:** Real-world biochemical activity is not established. These are annotation-matching scores. The labels are heuristic and the initial plasmid set is small/nonrandom; inspect `data/parts.csv` before treating scores as evidence. The pipeline deliberately declines to select a winner when coverage is weak.

## Environment and usage

Use a recent Python 3.10 or 3.11 environment, preferably Linux with NVIDIA CUDA and **at least ~12 GB free GPU RAM**, with more preferable for ESM-3. This memory figure is planning guidance, not measured in this environment. A CPU run of ESM-3 may be extremely slow. ESM model weights are downloaded on first run (multiple GB), requiring internet access. Installations vary by platform; begin from a GPU-enabled PyTorch environment (https://pytorch.org/get-started/locally/) before installing the requirements below.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python compare_models.py --self-test
python compare_models.py --prepare --email you@example.edu
python compare_models.py --run --models esm2,esmc,esm3 --device auto
```

The optional `--email` is passed to NCBI E-utilities. It is not saved. You can run `--prepare` without it. For a lower-resource first test:

```bash
python compare_models.py --run --models esm2 --device auto
```

If the optional ESM-3 weights require Hugging Face authentication, use `hf auth login` after accepting any access conditions on the model page. Always follow model license and use conditions.

Outputs in `results/`:

- `model_metrics.csv`: metrics by model and condition.
- `predictions.csv`: audited nearest-reference matches and cosine scores.
- `model_run_info.json`: checkpoint, device, measured time and failure metadata.
- `single_cds_baseline.csv`: deterministic trait coverage bound of retaining only one annotated CDS, not a model score.
- `embeddings_esm2.npz`, `embeddings_esmc.npz`, `embeddings_esm3.npz`: pooled embedding vectors indexed by sequence SHA256.
- `decision.md`: written decision under explicit data sufficiency and 0.05 macro-F1 tie criteria.

Data downloaded on `--prepare` live under `data/genbank/` with CSVs for candidate parts and omitted features. Edit `ACCESSIONS` in the script to add additional **benign, manually verified** annotated vectors. Curation of labels and independent homolog test sets are needed before using this for a formal benchmark.

## Model-specific technical notes

- `esm2`: Meta ESM-2 150M using `transformers.AutoModel` and `AutoTokenizer` with special tokens excluded from pooling.
- `esmc`: ESM-C 300M using `ESMC.from_pretrained("esmc_300m")`, `logits(..., return_embeddings=True)` via the official `esm` Python package.
- `esm3`: ESM-3 small using `ESM3.from_pretrained("esm3-sm-open-v1")` with the same ESM SDK embedding workflow.
- The adapter fails loudly if the token count does not match the sequence; no silent truncation. All datasets and models score on the same leave-one-construct-out policy, and none are fine-tuned.

Official references:

- Meta ESM-2: https://huggingface.co/facebook/esm2_t30_150M_UR50D
- Biohub ESM repository: https://github.com/Biohub/esm
- ESM-3 small: https://huggingface.co/biohub/esm3-sm-open-v1
- NCBI E-utilities: https://www.ncbi.nlm.nih.gov/books/NBK25499/
- Example record: https://www.ncbi.nlm.nih.gov/nuccore/U55762.1

## Interpreting the result

`INCONCLUSIVE` is not a failed run. Shared vector coding genes and near-identical constructs can leave very few **distinct** test examples. If this happens, expand the public benign panel with independent proteins and have a human review the annotation mappings.

A reasonably scoped follow-on choice would be ESM-C 300M if it matches ESM-3 within 5 macro-F1 points and inference cost remains acceptable; choose ESM-2 instead if it scores comparably and is clearly lighter/faster. ESM-3 needs a measurable gain on the specific task to justify its larger footprint. The automated `decision.md` algorithm uses measured performance and does not hard-code a favorite.
