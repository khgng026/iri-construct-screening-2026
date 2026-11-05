# ESM-2 / ESM-C 300M / ESM-3 Small feasibility

```powershell
python experiments/t05_model_comparison/compare_models.py --inputs data/processed/reviewed-inputs.json --validate-only --output results/comparison-check/run01
python experiments/t05_model_comparison/compare_models.py --inputs data/processed/reviewed-inputs.json --models esm2,esmc,esm3 --allow-download --device cpu --output results/comparison/run01
```

Install the [separate model dependencies and reviewed inputs](../README.md).
Start with `--models esm2` for a lower-resource feasibility run. Weights and
biological data are not downloaded by default. Inference runs sequentially.

The pre-existing untracked `iri_t05_model_comparison/` prototype was inspected.
This implementation reuses its local sequence-only adapter approach, residue
mean pooling, leave-construct-out nearest-neighbor retrieval, distinct-sequence
exclusion and explicit distinction between annotation retrieval and function
prediction. Its historical code is retained for review; its executable redirects
here. Automatic source downloads, heuristic label acceptance and automatic
winner selection are removed from executable use.

One shared set of reviewed input units and annotation-family labels is used by
all model adapters. References exclude the query's parent construct, exact
sequence duplicates, and reviewer-linked independence groups. At least one
independent same-family reference and one competing-family reference are
required. The nearest cosine reference supplies a family label; margins are
`best same-family cosine - best other-family cosine`. Ties use stable input-ID
order and are flagged. Missing references cause abstention rather than shrinking
that model's full reference panel.

The report saves each model's full-panel predictions, accuracy/macro-F1 over
scored annotation families, per-archetype results, embeddings coverage, retrieval
coverage/abstentions, individual failures, first-pass latency and process/GPU
memory measurements. It also re-scores an **exact common query and reference
panel** using the intersection of units embedded by every requested model. Both
sets of identities are saved. If any requested model cannot load, the common
panel is empty. Examine coverage and panel restrictions before comparing scores.

No cosine is computed between different models' vectors. No score is generated
from unreviewed annotations, and no overall model winner is selected. Family
retrieval is annotation agreement, not an intrinsic function assay. Shared
backbones, related proteins and possible pretraining overlap remain confounders
even after the explicit exclusions. T07/T08 must review family definitions,
independence groups and panel adequacy before scientific interpretation.

No approved comparison panel or SDK was available during implementation; no
pretrained comparison measurements are claimed. Numerical unit-test fixtures
only verify scoring, coverage and orchestration.
