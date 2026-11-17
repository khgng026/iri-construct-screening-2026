# Component-coverage study v1 — extension 1 results

Run 2026-11-17 under [SPEC_EXTENSION.md](../../../../experiments/component_study_v1/SPEC_EXTENSION.md),
committed (`b2fe310`) before any extension record was fetched; the extension
labels were frozen (`panel.freeze.json`, reusing the v1 reference panel
unchanged) before scoring. Agent self-review only. [v1 results](../README.md)
are not changed by anything here.

## Inputs and work

A second seeded sample of 800 records, disjoint from v1 ([acquisition](acquisition.json)).
Same toolkit command, model (`esm2_t6_8M_UR50D`, same revision and weight hash),
alignment scorer, family table and 134-reference panel. No labeled unit was a
near copy of a reference. ESM-2: 832 encoder calls, 68 s, 585 MB peak RSS.
Alignment: 93,532 pair scores, 77 s, 62 MB. (With 698 queries instead of 248,
alignment was the slower of the two; v1's ordering was a property of its size.)

## 1. Toolkit behavior on 800 more real records

800 of 800 records processed, no record failure, no parser warning. 1,850 of
1,963 annotated CDSs exported (94.2%); withheld: 102 partial/fuzzy boundary, 8
exceptional annotation, 3 non-canonical residues in the source translation. 35
exported proteins exceed the 1022-residue model limit. Mean share of construct
length inside exported CDSs: 32%.

Together with v1: **1,000 records, 0 failures; 2,338 of 2,491 CDSs exported
(93.9%)**, 153 withheld (138 fuzzy, 12 exceptional, 3 non-canonical), 47 over
the model limit. Still only NCBI flat files.

## 2. Parent-level replication on unit sets not seen in v1

505 multi-unit parents; 114 were dropped because their set of unit sequences had
already occurred in v1, leaving **363 distinct, previously unseen unit sets with
614 supported components**. Components recovered through their own unit, five
hits per parent:

| Condition | ESM-2 8M | Alignment | v1, ESM-2 (of 185) |
| --- | --- | --- | --- |
| `single_longest` (label-blind) | 200 / 614 (33%) | 209 (34%) | 55 (30%) |
| `single_longest_supported` | 352 (57%) | 362 (59%) | 106 (57%) |
| `all_max` | 509 (83%) | 474 (77%) | 151 (82%) |
| `all_round_robin` | 588 (96%) | 605 (98.5%) | 178 (96%) |
| Per-unit ceiling (larger budget) | 599 (98%) | 613 (99.8%) | 180 (97%) |

The v1 pattern reappears with almost the same proportions. Query-side cost of
the all-unit conditions: 1,107 instead of 363 model calls (3.05×), 326,299
instead of 145,676 residues (2.24×). The cost of extra units also reappears:
among the 176 distinct extension parents with exactly one supported component
(counted before removing sets seen in v1) the label-aware single protein
recovers 175, `all_max` 159 and round-robin 172 (ESM-2). The 127 single-unit
control parents give identical results in every condition.

**Uncertainty.** These parents form 20 dependence clusters, above the
pre-declared minimum of eight, so the project's paired cluster bootstrap was
applied (each cluster one vote, 2,000 replicates, 95% percentile interval). One
cluster contains 332 of the 363 parents — almost every vector shares a common
marker — so this is a deliberately conservative view in which 19 small, unusual
clusters outweigh the bulk of the sample:

| Difference in attributed recall | ESM-2 8M | Alignment |
| --- | --- | --- |
| round-robin − longest supported | +0.12 [0.04, 0.20] | +0.13 [0.05, 0.22] |
| maximum merge − longest supported | +0.05 [−0.08, 0.17] | +0.00 [−0.13, 0.10] |
| round-robin − maximum merge | +0.07 [0.01, 0.16] | +0.13 [0.03, 0.25] |
| round-robin − longest label-blind | +0.57 [0.39, 0.75] | +0.58 [0.41, 0.75] |

Per parent (unweighted) the same differences are +0.27, +0.17, +0.11, +0.62 for
ESM-2. Read conservatively: keeping all units with round-robin merging beats the
best single protein; keeping all units with maximum-score merging is not
distinguishable from the best single protein at cluster level. The population is
the set of records matching the fixed NCBI query, nothing wider.

## 3. Per-unit retrieval for sequences new in the extension

| Distinct supported sequences not in v1 | ESM-2 8M hit@1 / hit@5 | Alignment hit@1 / hit@5 |
| --- | --- | --- |
| Single-family whole protein (`full`, 121) | 119 / 119 | 121 / 121 |
| Family plus extra sequence (`embedded`, 29) | 14 / 17 | 28 / 28 |
| Two families in one CDS (`multi_family`, 5 units, 10 slots) | 5 of 10 slots at 5 | 9 of 10 slots at 5 |

The condition for withdrawing the v1 fusion observation is not met: ESM-2 finds
the family for 59% of new tagged/extended proteins against 98% of ordinary ones.
Over v1 and the extension together, fusion-class family slots recovered within
five hits: **ESM-2 28 of 53, alignment 50 of 53.**

## 4. Exploratory oracle diagnostic (not a method)

For every fusion-class unit, only the segment aligned to its labeling anchor was
embedded ([exploratory.oracle-segment.json](exploratory.oracle-segment.json)).
The cut uses the label's own alignment — the answer key.

| Family slots | Whole protein, hit@1 / hit@5 | Oracle segment, hit@1 / hit@5 |
| --- | --- | --- |
| Fusion classes (53) | 20 / 28 | 53 / 53 |
| Single-family control (189) | 186 / 186 | 186 / 186 |

Every family lost by pooling over a whole fusion protein is recovered at rank 1
once the extra sequence is removed, and trimming changes nothing for ordinary
proteins. The loss is therefore caused by mean pooling over unrelated sequence,
not by the model failing to recognize the domain. This bounds what splitting
*could* recover; it does not show that boundaries can be found without labels,
and it is not evidence about any splitting method.

## What the extension adds to the answer

- The coverage gain of representing every validated protein, its dependence on
  the merge rule, and its cost all replicate on new data and survive a
  conservative cluster-level interval for round-robin merging.
- The historical maximum-score merge gives no cluster-level advantage over the
  best single protein. If one rule is carried forward, it should be turn-taking.
- A fusion kept as one protein is the clearest limit of the current approach
  for an embedding model: about half of its documented families are lost to
  pooling, all of them recoverable in principle, and local alignment already
  recovers nearly all of them without any splitting.

## Limits

All v1 limits apply. Additionally: labels and segment boundaries both come from
alignment to the same anchors, so the oracle result is close to a best case; the
53 fusion slots are dominated by fluorescent-protein fusions; the bootstrap's
equal cluster weights and the single dominant cluster make its intervals a
statement about cluster means, not about a typical vector.

## Reproduce

```bash
E=docs/results/component-study-v1/extension-1
python experiments/component_study_v1/fetch_ncbi_sample.py --seed 20261118 --size 800 \
  --exclude docs/results/component-study-v1/acquisition.json \
  --cache data/raw/component_study_v1_ext1/ncbi --manifest $E/acquisition.json
iri-constructs preprocess data/raw/component_study_v1_ext1/ncbi \
  --out results/component_study_v1_ext1/bundle --max-model-residues 1022
python experiments/component_study_v1/build_panel.py --bundle results/component_study_v1_ext1/bundle \
  --private results/component_study_v1_ext1/panel.private.json --freeze $E/panel.freeze.json \
  --reuse-freeze docs/results/component-study-v1/panel.freeze.json \
  --reuse-private results/component_study_v1/panel.private.json
X="--bundle results/component_study_v1_ext1/bundle --private results/component_study_v1_ext1/panel.private.json --freeze $E/panel.freeze.json --out results/component_study_v1_ext1 --record $E"
python experiments/component_study_v1/score.py --method alignment $X
python experiments/component_study_v1/score.py --method esm2 $X          # model environment
python experiments/component_study_v1/analyze.py --scores results/component_study_v1_ext1 --record $E --compact
(cd experiments/component_study_v1 && python analyze_extension.py && python oracle_segment_diagnostic.py)
```
