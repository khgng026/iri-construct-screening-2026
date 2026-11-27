# Conventional reporter controls, November 27, 2026

This supplement applies the three methods from [PR 36 at the scoring revision](https://github.com/khgng026/iri-construct-screening-2026/tree/342c49846cb33e3cf07e52c9f3d8855e8a57f6c6) to conventional plasmids with documented products. The [aggregate table](counts.tsv) records the unit and denominator for each measure, including earlier comparisons from [PR 34](https://github.com/khgng026/iri-construct-screening-2026/pull/34).

| Primary endpoint | Raw DNA + BLASTx | Minimal translation + ESM-2 | Gene calling + ESM-2 |
|---|---:|---:|---:|
| Optical plasmids with every documented product recovered | 5/11 | 7/11 | 10/11 |
| Optical plasmids with complete recovery and zero uncredited calls | 5/11 | 7/11 | 6/11 |
| Documented product occurrences recovered | 20/26 | 22/26 | 25/26 |
| Uncredited calls among all calls returned | 36/56 | 14/36 | 6/31 |
| External yeast stocks with the documented marker recovered | 3/4 | 2/4 | 4/4 |

Each optical series was evaluated once using thresholds fitted on the other series. The primary rule allowed one uncredited call per training plasmid on average. Returned call counts varied on the held-out stocks. Five of eleven optical plasmids belonged to one related series. When each of the four series received equal weight, complete recovery was 25%, 75% and 75%, respectively. Across 17 reporting rules, restricting attention to those with at most eleven uncredited calls across the eleven held-out stocks gave 19/26 recovered occurrences for BLASTx, 14–17/26 for minimal translation and 20–25/26 for gene calling. No rule yielded identical observed extra-call counts for all three methods.

The four yeast stocks share one series. Each had one documented marker and each marker had appeared in earlier evaluations. Thresholds were frozen on the optical stocks before this check. The marker-only endpoint leaves other yeast calls unadjudicated.

## Database comparisons

Six distinct supplied reporter proteins came from twelve stocks and four families. Both BLASTp and phmmer gave a qualifying top-five identification for 1/6 on the 134-reference panel and 2/6 on the 139-reference panel. Full-index BLASTp gave 6/6 on UniProtKB 2025_03 and phmmer gave 5/6 on release 2025_01. The separate November 2026 provider run has its own study cutoff and is excluded from these completed totals.

One physical plasmid overlapping the optical series was also searched as whole DNA. The small panel supported 2/2 documented products and the first 100 subjects returned by NCBI nr supported 1/2. Output depth, software and index content differed. These counts describe that finite report.

## Verification and scope

Run the [offline recount](portable-recount/README.md) with Python 3.10 or later to reproduce optical and yeast aggregate outcomes from anonymized scored calls. The [provenance record](portable-recount/provenance.json) identifies the scoring revision and archived input and result hashes. The supplied-protein and single-plasmid database counts are transcribed from separately verified provider outputs.

Khang supplied the scoring implementation and whole-construct design. Rudra directed the preprocessing, inventory and reference-coverage comparisons. This supplement used AI-assisted computation and reanalysis. The diverse independent whole-construct comparison and proposed new-insert experiment remain open.
