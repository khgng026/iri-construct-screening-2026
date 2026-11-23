# Native multi-enzyme protein comparison

The reporter fusions in the initial construct pilot had whole-product labels. This experiment tests simple segmentation on known proteins with two deposited enzyme regions. Use four reviewed natural proteins representing three architectures. The two UMP synthases share an architecture. The other examples are yeast Ade5,7 and the tryptophan synthesis protein TrpCF. All are ordinary housekeeping proteins. [UniProt terms](https://www.uniprot.org/help/license)

## Comparison fixed before scoring

- Freeze the four query records and their deposited EC region labels. Use the six fixed single-enzyme references and three unrelated housekeeping references listed in the runner. Query and reference records have separate accessions.
- Compare the whole protein with two halves split at the midpoint, using the same total number of input residues per query. Compare both against the two documented regions as an annotation-derived diagnostic. Run reference-normalized local BLOSUM62 and mean-pooled ESM-2 on each choice of inputs. The midpoint case uses two model calls and the whole sequence uses one.
- For ESM-2, also pool the two midpoint windows from the already encoded whole protein. This uses the same full sequence encoding as the whole-protein arm and two reference comparisons per label.
- In every arm, retain the best score per label across the input spans and return the top two distinct labels. Record label hits, and separately count hits whose best scoring input covers at least 80% of the documented region. A whole input covers both regions by construction. The region check identifies when a proposed split cannot supply the documented region as one input, even if its label is recovered.

This is a descriptive protein-level diagnostic. It tests a fixed midpoint and pooled residue representation on known native proteins. The original construct-level fusion comparison needs its own verified inputs and key.

## Reproduce

Use the cached ESM-2 environment from the [pilot handoff](../../docs/project/HANDOFF.md). The runner verifies record roles, species, and labels during acquisition. Exact replays require the recorded sequence hashes and region annotations. Source and sequence payloads stay in ignored results storage.

```sh
run=results/native-housekeeping-fusions-v1
python -m experiments.benign_protein_pilot.native_fusions acquire "$run"
python -m experiments.benign_protein_pilot.native_fusions score "$run" --out "$run/report.json"
```
