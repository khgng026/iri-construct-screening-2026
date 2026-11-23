# Native multi-enzyme protein comparison

The reporter fusions in the initial construct pilot had whole-product labels. This experiment tests simple segmentation on known proteins with two deposited enzyme regions. Use four reviewed natural proteins representing three architectures. The two UMP synthases share an architecture. The other examples are yeast Ade5,7 and the tryptophan synthesis protein TrpCF. All are ordinary housekeeping proteins. [UniProt terms](https://www.uniprot.org/help/license)

## Comparison fixed before scoring

- Freeze the four query records and their deposited EC region labels. Use the six fixed single-enzyme references and three unrelated housekeeping references listed in the runner. Query and reference records have separate accessions.
- Compare the whole protein with two halves split at the midpoint, using the same total number of input residues per query. Compare both against the two documented regions as an annotation-derived diagnostic. Run reference-normalized local BLOSUM62 and mean-pooled ESM-2 on each choice of inputs. The midpoint case uses two model calls and the whole sequence uses one.
- For ESM-2, also pool the two midpoint windows from the already encoded whole protein. This uses the same full sequence encoding as the whole-protein arm and two reference comparisons per label.
- In every arm, retain the best score per label across the input spans and return the top two distinct labels. Record label hits, and separately count hits whose best scoring input covers at least 80% of the documented region. A whole input covers both regions by construction. The region check identifies when a proposed split cannot supply the documented region as one input, even if its label is recovered.

This is a descriptive protein-level diagnostic. It tests a fixed midpoint and pooled residue representation on known native proteins. The original construct-level fusion comparison needs its own verified inputs and key.

## Results

The panel and protocol were committed at `f0f4786` before scoring. The table gives correct labels out of eight. In this run every correct label was also supported by an input covering the corresponding region. [Frozen panel](native_fusion_panel.freeze.json)

| Input or representation | Mean-pooled ESM-2 | Normalized BLOSUM62 |
| --- | ---: | ---: |
| Whole protein | 3/8 | 8/8 |
| Two halves encoded separately | 4/8 | 8/8 |
| Two documented regions | 4/8 | 8/8 |
| Pool two halves from the whole encoding | 3/8 | — |

The midpoint split changed mean-pooled ESM-2 by one additional recovered label at the same residue budget and with two model calls per query. Even when the inputs followed the documented regions, it recovered four labels. Local pooling from the same encoding recovered three, the same total as whole-protein pooling.

An exploratory replay applied the previously fixed residue alignment scorer to the intact proteins. It recovered **7 of 8** labels with one encoding per query. This used `_encode` and `_residue_local_score` from the related-family comparison at `f4d0bfa`, with its existing gap settings. Each intact query was scored against all nine frozen references, and the top two distinct labels were retained. This result concerns label recovery. All calls in the replay came from the whole input. [Residue alignment scorer](https://github.com/khgng026/iri-construct-screening-2026/pull/21)

For these four proteins, residue alignment was more effective than changing the input boundaries for mean pooling. Normalized BLOSUM62 remained the strongest comparator. The sample is small and represents three architectures, so this result guides method selection for the next paired study.

## Reproduce

Use the cached ESM-2 environment from the [pilot handoff](../../docs/project/HANDOFF.md). The runner verifies record roles, species, and labels during acquisition. Exact replays require the recorded sequence hashes and region annotations. Source and sequence payloads stay in ignored results storage.

```sh
run=results/native-housekeeping-fusions-v1
python -m experiments.benign_protein_pilot.native_fusions acquire "$run"
python -m experiments.benign_protein_pilot.native_fusions score "$run" --out "$run/report.json"
```

The exploratory replay report is `results/native-housekeeping-fusions-v1/exploratory-residue-local.json` in the run workspace. The local suite ran 278 tests with five optional skips. The smoke and historical evidence checks passed.
