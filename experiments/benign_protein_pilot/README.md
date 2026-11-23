# Low similarity pilot on housekeeping proteins

This exploratory pilot asks whether a pretrained model retrieves a documented enzyme label when its known protein input has distant reference matches. The comparison isolates characterization from construct preprocessing. It uses naturally occurring housekeeping proteins from the reviewed part of UniProtKB. Some competing labels belong to the same broad enzyme superfamily, so a similar fold alone may be insufficient.

## Panel and comparison

- Require the single documented EC label and the corresponding InterPro annotation listed in `study.py`. These define the target before retrieval. Pick at most one reviewed entry per taxon and label by accession, and remove exact duplicates. Training exposure and experimental validation of every source annotation are unknown.
- Use the two listed model bacteria for references and the listed eukaryotes for queries. Every method ranks the complete reference panel for every query. No closer positive is removed after scoring. A query with a close positive stays in its appropriate tier.
- The primary remote tier requires less than 40% identity to the closest same-label reference among alignments covering at least 70% of **each** protein. Identity is the fraction of identical paired residues. Intermediate and high tiers use 40–<60% and ≥60%. Report insufficient coverage separately. The tiers describe full-length comparisons, with local matches also included in retrieval.
- Compare mean-pooled ESM-2 8M plus cosine against raw local BLOSUM62 alignment. Also report reference-self-normalized alignment as a sensitivity analysis. Use the same accession tie-break, and report first correct label at ranks 1 and 5, paired differences at rank 1, and denominators. The sample supports descriptive counts.

Run acquisition and scoring separately. Acquisition records the source release, request URI, entry and sequence version, and hashes before scoring. Source responses and sequences are stored in ignored local research storage, with terms from [UniProt](https://www.uniprot.org/help/license). The [ESM-2 model](https://huggingface.co/facebook/esm2_t6_8M_UR50D) uses the immutable revision pinned in the existing adapter and verifies the safetensors hash.

The exact selection for the completed run is recorded without sequence payloads in [panel.freeze.json](panel.freeze.json). To identify the same inputs after reacquisition, compare `sources_by_label` and `entries`. The retrieval date can change independently. See the [pilot handoff](../../docs/project/HANDOFF.md) for results.

The [related enzyme family comparison](RELATED_FAMILIES.md) adds a fresh panel with more low-identity queries and a residue alignment arm.

The [structural family study](STRUCTURAL_FAMILIES.md) tests very distant domains against pairwise and profile sequence methods. The [experimentally documented reaction study](EXPERIMENTAL_REACTIONS.md) tests whether an advantage carries over to exact biochemical labels on known full proteins.

```sh
python experiments/benign_protein_pilot/study.py acquire results/benign-protein-pilot-v0
# Use a Python environment with the project installed and the pinned ESM-2 model cached.
python experiments/benign_protein_pilot/study.py score results/benign-protein-pilot-v0 \
  --out results/benign-protein-pilot-v0/report.json
```
