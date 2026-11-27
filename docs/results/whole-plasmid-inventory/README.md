# Whole-plasmid inventory recovery at a fixed budget

Archived predictions for 777 documented vectors were reanalyzed. Nine vectors lacked a corroborated component, leaving 768 evaluable plasmids and 1,221 corroborated component occurrences. At most two candidate proteins were submitted per plasmid and five calls were returned. A plasmid was counted as complete when every corroborated occurrence was recovered. Additional calls could remain uncredited. [Source cohort and evaluation definitions](https://github.com/khgng026/iri-construct-screening-2026/blob/43c5c1f4a0599700305de981bb78df09c656f754/experiments/paired_preprocessing/COHORT_RESULTS.md)

| Input preparation | Components recovered with alignment / 1,221 | Complete with alignment / 768 | Complete with ESM-2 / 768 |
|---|---:|---:|---:|
| Minimal preprocessing | 594 | 322 | 315 |
| Start ranking | 618 | 342 | 337 |
| Gene calling | 746 | 383 | 377 |
| Supplied coding annotations | 1,134 | 696 | 686 |

The supplied-annotation condition received additional information. In the planned comparison with minimal preprocessing, 61 more inventories were completed with gene calling and alignment. Combining alignment and ESM-2 by rank sum added three completions to gene calling with alignment alone.

In an exploratory paired comparison, 45 inventories were completed only with gene calling and four only with start ranking. On all 45 gene-calling wins, two proteins were submitted by each method and fewer amino acids were submitted by gene calling. The effect of protein choice cannot be separated from the effect of input length in this comparison.

Every deposited coding feature was represented in the answer key for 246 of 768 plasmids. Within those 246, 161 inventories were completed with gene calling and 145 with start ranking. The other 222 gene-calling completions were scored against partial keys. Previously acquired records were reused. Agreement with the deposited annotations does not establish expression or performance on an independent population.

## Reproduction

From this directory in a full repository clone, run the following with Git and Python 3.10 or later. The first command recomputes paired counts, input volumes and key coverage directly from the pinned repository revision. The second verifies the hashes of the frozen predictions and recounts the two output budgets across four preparations and six aggregation methods. Model inference and network access are unnecessary.

```sh
python3 audit.py
python3 sealed_outcomes.py
```

[Paired results](audit-results.json) · [Aggregation results](sealed-outcome-results.json) · [Prediction hashes](sources/prediction-seals.json)
