# Distant protein recognition with same-fold alternatives

This study addresses the characterization part of [Q3](../../docs/project/SCOPE.md) on experimentally resolved natural domains involved in ordinary cellular functions. The prior structural study showed different results across three folds. Here the structural labels come from additional folds, and the competing labels share the query's broad fold. [Previous comparison](STRUCTURAL_FAMILIES.md)

## Protocol fixed before scoring

- Take contiguous domains from curated SCOPe 2.08 structure assignments in the model-organism taxa and fine families listed in `fold_study.py`. Include domains from 55 to 620 residues, using the archived boundaries. Give priority to distinct proteins within each fine family, take at most six unique domains per fine family, and exclude exact sequences in the earlier protein panels. Both the source and the label are fixed before evaluating the scorers.
- Predict the structural superfamily. Hold out entire fine families with the same deterministic three-way assignment as the previous study. Compare only within one fold, retaining folds with at least two candidate superfamilies. Every method receives the same reference domains, including the same-label relatives and all competing labels.
- The primary group is below **25% global identity to every available reference**, counting gaps in the alignment length. Also report the 25% to below 30% tier and a check for a same-label local alignment covering at least 70% of both sequences. The two measures distinguish low whole-domain identity from absence of broadly covering local alignment.
- Use raw and reference-normalized BLOSUM62, HMMER profiles per fine family and per superfamily, and both mean and residue-local ESM-2 at 8M and 150M. Score the first-ranked label, with ties treated as missed identifications. The main paired comparisons use residue-local 150M against raw alignment and each profile baseline. Count model rescues among cases missed by all four conventional methods.
- Summarize success per domain, by fold, and as an average giving folds equal weight. Within each fold give superfamilies equal weight, then fine families within their superfamily. Report paired differences, a descriptive whole-fold bootstrap interval, and sensitivity when dropping one fold. Treat the intervals as panel summaries. Several groups share reference data and the domain sequences may have appeared in ESM-2 pretraining.

The measured label is a structural superfamily, which can contain enzymes with different functions. Evidence that a scorer can recover this label after receiving the domain is one part of Q3. The separate construct experiments address whether a useful input reaches the scorer.

## Reference-depth check fixed before running it

The first comparison holds out several fine families at once. As a check on reference availability, give every method all the other fine families within the query's fold. Keep the entire query family out. Use the same 132 queries, candidate labels and scoring rules. Recalculate identity against the larger reference set and compare decisions on both the original below-25% cohort and the cases still below 25%. Some superfamilies have only two fine families, so this adds same-label references for only some queries. This check reuses the original sample.

## Results on the five additional folds

The [source manifest](fold_panel.freeze.json) and original protocol were committed at `fce0b53` before scoring. The panel contains 132 unique domains from 55 fine families and 15 structural superfamilies. In the three-way split, 128 domains from 52 fine families qualify for the below-25% tier. Their highest global identities to any reference range from 11.8% to 24.6%. Eleven have a same-label local alignment covering at least 70% of both domains. [Aggregate results and domain-level decisions](fold_results.json)

| Method | Correct first label, below 25% | Accuracy with folds weighted equally |
| --- | ---: | ---: |
| Raw BLOSUM62 | 86/128 | 70.1% |
| Reference-normalized BLOSUM62 | 60/128 | 60.0% |
| HMMER fine-family profiles | 81/128 | 66.4% |
| HMMER superfamily profiles | 86/128 | 70.3% |
| ESM-2 8M mean | 76/128 | 56.8% |
| ESM-2 8M local | 82/128 | 61.1% |
| ESM-2 150M mean | 85/128 | 62.8% |
| ESM-2 150M local | **99/128** | **83.5%** |

The 150M local model ranked the correct label first for 8 of the 17 domains missed by all four conventional methods. Those eight recoveries represent four folds. Against the superfamily profile alone, the model had 25 unique correct calls and the profile had 12. The equal-fold gain was 13.1 percentage points, with a descriptive whole-fold bootstrap interval of −3.1 to +27.5 points. The equal-fold gain over raw BLOSUM62 was 13.4 points, with an interval of −3.3 to +34.2. The model still led in the average after dropping any single fold. With five folds, the intervals leave substantial uncertainty about an average on new folds.

| Fold | Domains below 25% | Raw alignment | Fine-family profile | Superfamily profile | ESM-2 150M local |
| --- | ---: | ---: | ---: | ---: | ---: |
| a.118 | 49 | 34 | 25 | 29 | **39** |
| b.34 | 31 | 18 | **24** | **24** | 20 |
| b.43 | 17 | **11** | 7 | 7 | 10 |
| d.58 | 18 | 15 | 14 | 15 | **17** |
| d.79 | 13 | 8 | 11 | 11 | **13** |

Giving more weight to distant cases leaves only four queries from 25% to below 30% and zero at 30% or above. The panel documents successes and misses within the low-identity range. It has too little data at higher identity to estimate a crossover. The two model sizes received identical inputs. The 150M encoder took 133 seconds and the 8M encoder took 9 seconds on two CPU threads. For the original split, each scorer considered 2,666 protein pairs or the HMMER profiles fitted from the same references. Conventional methods built 140 profiles across the two profile arms.

## Results with more references per query

The alternate reference rule was committed at `7732f6e` before scoring. It increased the total comparisons from 2,666 to 3,792 pairs. Forty-one of the 128 primary examples gained a closer reference and all remained below 25%. The count of queries with a broadly covering same-label local match increased from 11 to 15. [Reference-depth outcomes](fold_wider_reference_results.json)

| Reference policy | Raw alignment | Fine-family HMM | Superfamily HMM | Local ESM-2 150M | Local model correct among all conventional misses |
| --- | ---: | ---: | ---: | ---: | ---: |
| Three-way split | 86/128 | 81/128 | 86/128 | **99/128** | 8/17 |
| Hold out one fine family | 78/128 | 73/128 | 85/128 | **98/128** | 10/19 |

The aggregate gain over the strongest profile persisted. The equal-fold difference narrowed to +8.9 points, with a descriptive interval of −6.8 to +21.7 points. The same three folds favored the model and the other two favored a conventional method. Fourteen examples were missed by all four conventional methods under both reference policies. The model correctly identified six of those fourteen under both policies.

The small joint-miss sets make the fraction of model recoveries uncertain. A uniformly random choice among the two to four possible labels would average 5.9 correct on the 17 original misses and 6.7 on the 19 wider-reference misses. Those expectations are the sum of one divided by the number of labels for each case. The larger primary cohort gives the stronger evidence for a model contribution. The two reference settings reuse the same domains and constitute a sensitivity comparison.

## Reproduce

Use the [CPU model environment](../../docs/project/HANDOFF.md) with `pyfamsa==0.7.0` and `pyhmmer==0.12.3`. The run used Python 3.12, Biopython 1.86, PyTorch 2.6.0 CPU, and Transformers 4.48.1. The runner verifies both model checkpoints and both the FASTA and sequence hashes for each published domain. The final JSON files contain the methods' decisions and summaries. Their timing fields may differ on replay. The source snapshots remain in ignored storage. [SCOPe 2.08](https://scop.berkeley.edu/astral/ver=2.08) · [SCOPe terms](https://scop.berkeley.edu/about)

A fresh download matched all 132 source and sequence hashes. Replaying the three-way split reproduced both model sizes, both alignments, the superfamily profile, and the all-conventional-miss counts. One correct/incorrect decision changed in the fine-family profile, taking its below-25% total from 81 to 82. Eight further evaluations of that query block confirmed that the best two label scores can switch order. The tracked result retains the original run. The main comparison uses the stronger superfamily profile, whose decisions all matched.

```sh
study=experiments.benign_protein_pilot.fold_study
run=results/q3-folds-source-replay
python -m "$study" restore "$run"
for policy in three-way leave-one-family-out; do
  out="$run/$policy"
  mkdir -p "$out"
  for arm in conventional 8m 150m; do
    OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python -m "$study" score "$run" \
      --arm "$arm" --reference-policy "$policy" --out "$out/$arm.json"
  done
  python -m "$study" summarize --conventional "$out/conventional.json" --small "$out/8m.json" \
    --large "$out/150m.json" --out "$out/report.json"
done
```
