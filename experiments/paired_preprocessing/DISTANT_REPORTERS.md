# Distant reporters from whole deposited DNA records

This pilot connects the known-protein experiments to the original preprocessing question. Given the whole DNA record, can a model identify a documented reporter's protein group where the available references are distant? The targets are familiar fluorescent proteins and firefly luciferase in public GenBank records. The source records are unchanged. [Earlier construct study](COHORT_RESULTS.md) · [Known-protein follow-up, PR #28](https://github.com/khgng026/iri-construct-screening-2026/pull/28)

## Protocol fixed before scoring

- Use the ten pinned records and 14 named targets in [the panel](distant_reporter_panel.json). These records are absent from the earlier 1,000-record audit and five-reporter pilot. Validate each deposited target's translation against the DNA. Evaluate the GFP-like barrel and the adenylate-forming enzyme superfamily. Reporter publications and source annotations supply the labels. The scoring references come from a separate source. Group related examples by parent reporter ancestry.
- The reference library has ten proteins in five groups. The GFP group uses jellyfish proteins from *Aequorea* and *Clytia*. The adenylate-forming group uses two natural *E. coli* K-12 enzymes with the same broad fold as firefly luciferase. The remaining reference groups are avidin, lipocalin and phosphoglycerate kinase. The result concerns broad protein groups in this reference collection. [Jellyfish GFP study](https://doi.org/10.1039/c0pp00023j) · [Firefly structure](https://www.rcsb.org/structure/1LCI)
- Freeze the primary low-identity tier below 30% global identity to the best same-group reference, counting gaps in the alignment length. For comparison with earlier studies, record whether the same group also has a local alignment covering at least 70% of both proteins and what that local identity is. Interpret the strength of the result in light of both measurements.
- Reuse the established minimal translation and Pyrodigal arms. Both receive the DNA alone in its deposited orientation. Score the same selected protein candidates for all methods at caps of one, two and four, with two primary. The methods return at most one group per candidate. Count each candidate and target at most once, with the correct label and at least 80% overlap of the documented coding bases on the same strand.
- Compare raw and reference-normalized BLOSUM62, HMMER profiles, and both mean and local scoring with the pinned ESM-2 8M and 150M checkpoints. Rank against the same reference proteins for all methods. Repeated identical protein inputs are cached across arms. Record candidate counts and residues separately from encoder cache reuse.
- Score each documented target directly as a diagnostic for known boundaries. At each budget record whether the candidate pool covers the locus, whether the selected candidates cover it, and whether any selected candidate exactly matches the documented translation. Locus coverage alone can include a different translation. Report component counts and equal-ancestry recall. Interpretation is descriptive because the panel has few independent reporter ancestries.

The fluorescent labels are supported by the deposited product names and published reporter protein studies, independently of the library used to score them. Only the target reporters enter the denominator. Supporting vector features are outside the evaluation. [mCherry](https://doi.org/10.1038/nbt1037) · [mNeonGreen](https://doi.org/10.1038/nmeth.2413) · [GFP family](https://www.uniprot.org/uniprotkb/P42212/entry) · [Firefly luciferase](https://www.uniprot.org/uniprotkb/P08659/entry)

## Reference-library sensitivity fixed before scoring it

The first audit found ten targets below 30% global identity, each with a same-group local match covering at least 70% of both proteins at over 30% local identity. To probe another set of natural relationships, repeat the same comparison with a [second prespecified library](distant_reporter_alternate_references.json). Use two structurally characterized copepod GFPs and two reviewed plant 4-coumarate-CoA ligases for the positive reference groups. Keep the three other groups fixed. This library is selected from published natural proteins by taxon and structural group. Compare the results for both libraries and give the below-25% global tier separately. The two experiments use the same reporter constructs. [Copepod structure](https://doi.org/10.1111/j.1365-2443.2009.01305.x) · [Copepod GFP superfamily](https://doi.org/10.1016/j.jmb.2006.04.002) · [Plant ligase](https://www.uniprot.org/uniprotkb/Q42524/entry)

The second audit found eight targets below 30% global identity. Four were below 25% and lacked an adequately covering local match. All four come from the same *Aequorea* GFP ancestry. This subset can illustrate an end-to-end result for one laboratory reporter lineage.

## Additional reporter ancestry fixed before scoring it

Use [three separately deposited Renilla reporter vectors](distant_renilla_panel.json) to extend the comparison to a different reporter ancestry and fold. The deposited Renilla products belong to the broad alpha/beta hydrolase fold, supported by an experimentally solved Renilla luciferase structure. The two same-fold references are a fungal and a plant lipase with the corresponding structural annotation. Keep the second library's other groups as competitors. Apply the same DNA-only arms, loci, caps and scorers. The endpoint is recognition of the fold group. The three vectors represent one reporter ancestry. [Renilla structure](https://www.rcsb.org/structure/2PSF) · [Fold annotation](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR029058/)

The source audit validated the three deposited translations. Their best same-group global identities range from 18.3% to 18.8%. None has a local match covering 70% of both the reporter and reference.

## Results on the ten original records

The first comparison was fixed at `dffb342` and the second library at `783f61d`. All original decisions reproduced after reacquiring the public data. For each arm and cap, the fingerprints of the selected inputs are identical across the libraries. These records contain nine distinct target translations. Seven of the 14 target occurrences have a translation also present in the earlier construct holdout. This is a test on different vector records, with repeated reporter proteins and shared ancestry.

Each row below is the primary cap of two protein inputs per construct. Counts are correct broad-group calls at a matching locus. Both DNA-only arms made 20 queries. The minimal arm used 7,660 residues and gene calling used 6,891.

| Reference library and tier | DNA preparation | Targets reached by selected candidates | Raw alignment | Profile | ESM-2 8M local | ESM-2 150M local | ESM-2 150M mean |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Initial, below 30% (10) | Minimal | 6 | 5 | 6 | 5 | 5 | 6 |
| Initial, below 30% (10) | Gene calling | 6 | 6 | 6 | 6 | 6 | 6 |
| Alternate, below 30% (8) | Minimal | 2 | 1 | 2 | 0 | 0 | 2 |
| Alternate, below 30% (8) | Gene calling | 4 | 4 | 4 | 4 | 4 | 4 |
| Alternate, below 25% (4) | Minimal | 2 | 1 | 2 | 0 | 0 | 2 |
| Alternate, below 25% (4) | Gene calling | 3 | 3 | 3 | 3 | 3 | 3 |

With the initial library and the primary cap, the profile recovered six targets in each arm. The successful targets differed. The minimal arm recovered all five firefly instances and the sole lancelet instance. Gene calling recovered the five firefly instances and one of four *Discosoma* instances. This gives the minimal arm 66.7% recall with each reporter ancestry weighted equally, versus 41.7% for gene calling. With raw alignment those numbers are 33.3% and 41.7%. The small sample makes these estimates sensitive to a single reporter.

With the alternate library and the primary cap, gene calling moved three of the four very distant *Aequorea* reporters into the selected inputs. Every scorer correctly labeled those three. The minimal arm selected covering candidates for two cases, neither with the exact documented translation. Mean ESM-2 150M and the profile recovered both. Raw alignment found one, and the two local embedding scorers found zero. When given the four documented proteins directly, every scorer found all four. Three of the four occurrences repeat protein translations seen in the earlier construct study.

With a cap of four, gene calling recovered all ten initial low-identity targets and all eight alternate low-identity targets with every scorer. It used 34 queries across the ten records, compared with 40 for the minimal arm. The minimal arm in the alternate below-25% subset recovered two with raw alignment, three with the profile and 150M mean, and one with 150M local. The strongest model representation varied across comparisons.

## Results on the three Renilla records

The Renilla protocol and sources were fixed at `2557b3a` before scoring. The three target translations are distinct and absent from the earlier construct holdout. Both preparation arms reached all three target loci with just one candidate per construct. Raw BLOSUM62, HMMER, and both local ESM-2 variants recovered all three at every cap. Reference-normalized BLOSUM62 recovered zero, including when given the documented proteins directly. Mean-pooled ESM-2 150M recovered zero with minimal preparation and one with gene calling at the primary cap. One of the three gene-caller inputs exactly matched its documented protein, and the other two received correct fold calls from the strongest scorers. An analysis using only the normalized baseline would have missed the parity with raw alignment and HMMER.

## What this answers

The two below-25% comparisons contain seven target occurrences from two reporter ancestries. Every target had a compatible candidate somewhere in each preparation arm's pool. At a budget of two, the gene caller selected compatible candidates for six, and raw alignment, HMMER and local ESM-2 each recovered those same six. Minimal preparation recovered five with HMMER, four with raw alignment, and three with either local ESM-2 scorer. The two comparisons use different reference libraries, so these totals describe the examples tested. The tracked file stores source and identity details under `targets` and outcomes under each query budget. [Paired outcomes](distant_reporter_results.json)

Low identity across the full protein and the absence of a broadly covering local alignment can coexist with useful local signal. Raw alignment and HMMER recognized all seven of these very distant reporters when given their documented proteins. On the prepared inputs, the dominant remaining miss was a candidate that fell outside the query cap. This pilot establishes parity for these two reporter ancestries after gene calling. It leaves open whether embeddings add recoveries for fresh construct components where both strong conventional comparators miss the documented protein group. A later study needs more independent ancestries, credible alternatives from similar folds, and a fixed reference library.

The result concerns recall of named reporters among the tested groups. Estimating false-positive rates would require verified labels for the candidate outputs from the other vector features.

## Reproduce

Use the [model environment](../../docs/project/HANDOFF.md), Pyrodigal 3.6.3, PyFAMSA 0.7.0, and PyHMMER 0.12.3. Cache both ESM-2 checkpoints at the immutable revisions in `cohort_scorers.py`. The runner verifies their weight hashes and loads the models offline. NCBI's batch and individual record responses differ in their terminal line breaks. This pilot canonicalizes only those line breaks before verifying each GenBank snapshot. [NCBI data policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) · [UniProt terms](https://www.uniprot.org/help/license)

```sh
study=experiments.paired_preprocessing.distant_reporter_study
for panel in initial alternate renilla; do
  run="results/distant-reporters-replay-$panel"
  python -m "$study" --panel "$panel" acquire "$run"
  python -m "$study" --panel "$panel" audit "$run"
  OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
    python -m "$study" --panel "$panel" score "$run" --out "$run/report.json"
done
```
