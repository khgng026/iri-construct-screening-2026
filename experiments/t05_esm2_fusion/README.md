# ESM-2 fusion feasibility

This experiment uses one fixed ESM-2 instance for:

1. Mean-pooled isolated target protein residues.
2. Mean-pooled complete fusion protein residues.
3. Mean-pooled target-region residues from the full-fusion forward pass.

Condition 3 slices the reviewed amino-acid interval **after** full-context
inference. It does not re-embed an extracted sequence. All three conditions use
the same weights, final-layer output, residue pooling and exclusion of special
tokens. Only two forward passes are needed per case. There is no generated
fusion, alignment-based boundary inference, reconstructed fragment or training.

```powershell
python experiments/t05_esm2_fusion/compare_fusion.py --inputs data/processed/reviewed-inputs.json --validate-only --output results/fusion-check/run01
python experiments/t05_esm2_fusion/compare_fusion.py --inputs data/processed/reviewed-inputs.json --allow-download --device cpu --output results/fusion/run01
```

See the [shared setup and reviewed-input contract](../README.md). Coordinates
are 0-based half-open amino-acid bounds in the fusion. The isolated target must
exactly equal that reviewed span; both source records require T02 approval.

The JSON report contains three pairwise pooled cosines, pooled cosine distance
between isolated/contextual targets, mean per-residue cosine and RMS difference
of residue representations, plus condition-specific coverage, failures, timing
and resources. RMS is scale dependent and meaningful only within this model.
Contextual timing reuses the full-fusion inference; it is not charged as another
forward pass. If the fusion exceeds the cap, isolated inference can still succeed
and coverage is reported as partial. No optional reference-retrieval margin is
implemented here; independent-family retrieval lives in the comparison module.

These measurements quantify representation changes. Similarity or context
dependence does not establish protein function, biological activity, a model
improvement or a safety classification. No reviewed fusion was available during
implementation, so no pretrained fusion measurements are claimed.

## Earlier pair-manifest workflow

The upstream FASTA/pair-validation and numerical comparison helpers remain
importable from `compare_fusion.py`; its original regression tests are preserved.
`pairs.template.csv` remains a record of the earlier CSV interface. To migrate
from the former `--fasta/--pairs` CLI, create the shared reviewed-input JSON and
T02 manifest described above, with source checksums and review notes. Convert
1-based inclusive bounds explicitly: `start0 = start_aa - 1`, `end0 = end_aa`.
Pair-level approval flags alone do not replace source-file T02 reviews. Current
scientific execution uses the `--inputs` commands above; the compatibility
helpers do not load models or download records.
