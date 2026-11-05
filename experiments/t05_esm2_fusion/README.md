# T05: ESM-2 isolated-protein versus fusion-protein embeddings

## Scope
Compare representations of benign, publicly documented, reviewed protein fusions. ESM-2 accepts proteins, **not intact DNA**. This is a feasibility study, not a biosecurity classifier or a validation of protein function.

## Required inputs
Create local `data/raw/proteins.fasta` containing reviewed amino-acid FASTA entries with unique IDs, and local `data/raw/pairs.csv` using `pairs.template.csv` as a header template.

Each row pairs an independently annotated isolated protein with a longer fusion containing the **identical isolated sequence**. Coordinates are 1-based inclusive positions in the fusion protein. The script rejects mismatched subsequences and unapproved data. Required accession, source URL, and boundary-evidence fields must be populated; the researcher must verify the source, licensing and benign status manually under T02. Do not commit sequences or generated embeddings.

Optionally provide two other FASTA IDs: `reference_id` (independent, annotated related benign protein, but not sequence-identical) and `decoy_id` (different benign family). Avoid reused duplicates across pairs.

## Run
From the repository root (Python 3.10+):

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r experiments/t05_esm2_fusion/requirements.txt
python experiments/t05_esm2_fusion/compare_fusion.py   --fasta data/raw/proteins.fasta   --pairs data/raw/pairs.csv --validate-only
python experiments/t05_esm2_fusion/compare_fusion.py   --fasta data/raw/proteins.fasta   --pairs data/raw/pairs.csv   --model facebook/esm2_t6_8M_UR50D   --outdir results/t05_esm2_fusion
```

On a machine with more resources, also try `facebook/esm2_t30_150M_UR50D`. Record the model version for each run. Inference requires internet access to download uncached pretrained weights; the validation mode does not.

## Measures
- **Isolated vs whole fusion cosine:** similarity of full-sequence residue-mean-pooled representations.
- **Isolated vs contextual target cosine:** similarity between the isolated protein and the target-region residue embeddings *from the fusion-model pass*.
- **Contextual-minus-whole:** the paired difference between the two measures.
- **Optional two-candidate reference margin:** related reference cosine minus unrelated decoy cosine, computed for all three representations, to check whether the benign related reference is ranked higher.

Pooling excludes BOS/EOS special tokens. There is **no truncation**; overlength proteins raise errors. Identical isolated and separately extracted subsequences inevitably yield identical embeddings when encoded alone, and therefore cannot be evidence of improvement. These descriptive cosine differences do not prove protein function; independent annotation retrieval on a broader panel is necessary for stronger claims.

## Outputs
All outputs are local and gitignored under `results/`: `pair_metrics.csv`, `embeddings.npz`, `embedding_index.csv`, and `run_metadata.json` (device, versions, runtime, input checksums, aggregates). For a useful pilot, use multiple independent reviewed benign pairs and note whether the fusion is N-terminal or C-terminal. T05 remains in progress until actual model inference is completed and reviewed.

## Tests
Run dependency-free checks with `python -m unittest discover -s tests -v`. No biological sequence fixtures are committed.
