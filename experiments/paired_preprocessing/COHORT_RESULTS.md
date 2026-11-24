# Cohort v2 results: Pyrodigal preprocessing at a matched query budget

## Answer

**Yes, for this cohort. At two protein inputs per construct, Pyrodigal preprocessing raised documented-component recall by about nine percentage points over minimal six-frame conversion, with either scorer.**

- **ESM-2:** recall rose from 46.9% to 59.4% pooled (+12.4 points). The equal-lineage-group gain was **+9.2 points (95% whole-group bootstrap 6.6–11.8)**.
- **Normalized BLOSUM62:** **+9.0 points (6.5–11.6)**. The gain persists with conventional alignment.
- **Cohort:** 1,221 corroborated documented components in 777 held-out public benign constructs, spanning 422 lineage groups and 27 families.
- **Budget:** Pyrodigal used slightly fewer inputs than allowed and 7.6% fewer input residues than the baseline.
- **Locked rule:** both gains meet the criterion for a practically meaningful benefit. The lower bound exceeds five points, and the leave-one-group-out minimum is 9.0 points for ESM-2 and 8.8 for BLOSUM62.

**The gain comes from candidate selection, not characterization.** Both scorers gained on almost the same components. The difference in gain between scorers was 0.2 points (−0.6 to 1.0). Pyrodigal moved 156 of the baseline's below-cap documented genes into the two-input budget.

**On identical inputs, BLOSUM62 found slightly more components than ESM-2.** The difference was 2.3 points with Pyrodigal inputs (0.9–3.7) and 2.5 points with minimal inputs (1.2–3.9). The gap is concentrated in tagged or fused proteins and in the low-similarity tier.

Qualifications that bound the claim:

1. **The gain is concentrated in two families.** Of the 152 extra ESM-2 recoveries, 124 are class A β-lactamases and aminoglycoside phosphotransferases. These are bacterial resistance genes that a prokaryotic gene caller finds well. The pre-specified two-way bootstrap, which resamples lineage groups and families together, gives a pooled interval of 2.8–23.6 points. Without those two families the equal-group gain is still 7.4 points (4.0–11.0), but the pooled gain falls to 4.9.
2. **Under conservative shared-protein grouping** the holdout is effectively one group, holding 96% of components. That view is descriptive: ESM-2 +5.8 (−4.7 to 16.6) and BLOSUM62 +3.1 (−5.4 to 12.2) with equal group weights.
3. **Tagged or fused constructs are inconclusive.** ESM-2 gained +4.6 (−4.0 to 13.5) and BLOSUM62 +3.1 (−4.1 to 10.5) across 43 groups.
4. **The gain depends on the budget.** At one input per construct it was +4.0 points (2.1–6.0). At four inputs it was +9.1 with ESM-2 and +6.8 with BLOSUM62.
5. **Pyrodigal also loses some genes.** It has no compatible candidate for 45 documented components, compared with 11 for the minimal pool.

Clean documented boundaries still leave about 27 points of headroom. Annotation-derived inputs reached 91.4% with ESM-2 and 92.9% with BLOSUM62 at the same cap.

These are research results for this cohort of deposited cloning and expression vectors. They do not establish biological function, risk, or performance on other construct sources. Self-review only; no independent specialist review is implied.

## Cohort

Rules, partitions and hashes are in the [locked protocol](COHORT_STUDY.md), committed at `deaf88b` before holdout inference. The holdout ran once with that commit's code.

| Holdout | Count |
| --- | ---: |
| Constructs | 777 |
| Corroborated components (primary denominator) | 1,221 |
| Anchor-only components (exploratory) | 33 |
| Lineage groups with a corroborated component | 422 (largest 20% of components; 23.9 effective) |
| Conservative groups | 37 (largest 96%) |
| Families | 27 |
| Archetype records: multi-CDS, backbone-heavy, single label + other CDSs, single-CDS, tagged/fused | 355, 144, 128, 99, 51 |

## Primary paired table (cap of two inputs per construct)

Each recall cell gives pooled recall, followed by equal-lineage-group recall in parentheses. Inputs are actual / allowed. Residues are the same for both scorers because the selected inputs are identical: fingerprints matched for every arm.

| Arm | ESM-2 recall | BLOSUM62 recall | Inputs | Residues | Candidate pool |
| --- | ---: | ---: | ---: | ---: | ---: |
| Minimal (baseline) | 46.9% (56.5%) | 48.6% (59.0%) | 1,554 / 1,554 | 596,911 | 101,407 |
| Start ranking (control) | 49.0% (58.9%) | 50.6% (61.1%) | 1,554 / 1,554 | 596,161 | 101,407 |
| **Pyrodigal** | **59.4% (65.7%)** | **61.1% (68.0%)** | 1,537 / 1,554 | 551,358 | 4,194 |
| Annotation-derived (diagnostic) | 91.4% (93.4%) | 92.9% (95.0%) | 1,169 / 1,554 | 344,108 | 1,250 |

| Quantity at cap 2 | Equal-group gain (95% CI) | Pooled gain (95% CI) | Leave-one-group-out, equal-group | Groups up / down / tied |
| --- | ---: | ---: | ---: | ---: |
| 1. Pyrodigal − minimal, ESM-2 | **+9.2 (6.6, 11.8)** | +12.4 (6.7, 17.9) | 9.0 to 9.5 | 74 / 7 / 341 |
| 2. Pyrodigal − minimal, BLOSUM62 | **+9.0 (6.5, 11.6)** | +12.4 (6.8, 17.9) | 8.8 to 9.3 | 71 / 5 / 346 |
| Difference in gain (ESM-2 − BLOSUM62) | +0.2 (−0.6, 1.0) | 0.0 (−0.6, 0.5) | — | — |
| 3a. ESM-2 − BLOSUM62, minimal inputs | −2.5 (−3.9, −1.2) | −1.7 (−3.2, −0.8) | — | 3 / 20 / 399 |
| 3b. ESM-2 − BLOSUM62, Pyrodigal inputs | −2.3 (−3.7, −0.9) | −1.7 (−3.2, −0.7) | — | 5 / 21 / 396 |
| 3c. ESM-2 − BLOSUM62, annotation inputs | −1.6 (−2.7, −0.7) | −1.5 (−2.9, −0.5) | — | 4 / 17 / 401 |
| Start ranking − minimal, ESM-2 / BLOSUM62 | +2.4 (1.0, 4.0) / +2.2 (0.9, 3.6) | +2.0 / +2.0 | — | — |
| Annotation − Pyrodigal, ESM-2 / BLOSUM62 | +27.7 (24.3, 31.1) / +27.0 (23.5, 30.4) | +32.0 / +31.8 | — | — |

Under the locked rule, Pyrodigal shows a **practically meaningful benefit** with both scorers. Start ranking's gain is positive, but a **gain of target size is unsupported** (upper bound below five points). Omitting the largest lineage group (240 components) leaves a pooled Pyrodigal gain of 8.4 points.

## Budget curves

Pooled recall (equal-group gain of Pyrodigal over minimal, 95% CI). Inputs and residues are totals over 777 constructs.

| Cap | Minimal ESM-2 / BLOSUM62 | Pyrodigal ESM-2 / BLOSUM62 | Gain ESM-2 | Gain BLOSUM62 | Inputs (minimal / Pyrodigal) | Residues (minimal / Pyrodigal) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 22.9% / 23.8% | 26.5% / 27.5% | +4.0 (2.1, 6.0) | +4.0 (2.1, 6.0) | 777 / 777 | 345,865 / 334,550 |
| **2** | 46.9% / 48.6% | 59.4% / 61.1% | +9.2 (6.6, 11.8) | +9.0 (6.5, 11.6) | 1,554 / 1,537 | 596,911 / 551,358 |
| 4 | 83.0% / 86.2% | 90.0% / 91.5% | +9.1 (6.7, 11.6) | +6.8 (4.6, 9.1) | 3,108 / 2,772 | 971,259 / 789,199 |

Pyrodigal never used more input work than the baseline. Its pools are about 24 times smaller, but gene calling took 9.0 s over the holdout, against 3.0 s for six-frame translation. ESM-2 made 1,688 encoder calls for the minimal arm at cap 2 and 1,671 for Pyrodigal, both including 134 reference encodings. BLOSUM62 made 208,236 and 205,958 reference comparisons, with 27.6 and 25.5 billion alignment cells. Each ESM-2 stage took 51–199 s of wall time and each BLOSUM62 stage 34–194 s, on two threads and one thread respectively, with peak memory of about 0.5 GB.

## Failure breakdown at the primary cap

| Status of each documented component | Minimal | Start ranking | Pyrodigal | Annotation-derived |
| --- | ---: | ---: | ---: | ---: |
| Recovered (ESM-2 / BLOSUM62) | 573 / 594 | 598 / 618 | 725 / 746 | 1,116 / 1,134 |
| Compatible candidate below the cap | 600 | 577 | 423 | 78 |
| No compatible candidate in the pool | 11 | 11 | 45 | 4 |
| Covering input reached the scorer, label missed (ESM-2 / BLOSUM62) | 37 / 16 | 35 / 15 | 28 / 7 | 23 / 5 |

Selection dominates the misses. Of the minimal arm's 600 below-cap misses, Pyrodigal recovered 156, left 412 below the cap, and dropped 29 from its pool. Start ranking recovered 36 of them from the same pool. The annotation-derived arm recovered 528. Its remaining 78 below-cap misses come from constructs with more than two documented components, where two inputs cannot reach every locus. When a covering input reached the scorer and the label was missed, the correct family usually ranked 4th to 16th. The scorer, not the locus, was the problem in those cases. Clean annotation boundaries reduced ESM-2's count from 37 to 23.

## Where it helps

| Stratum (components, lineage groups) | Minimal → Pyrodigal ESM-2 | Minimal → Pyrodigal BLOSUM62 | Gain ESM-2 (CI) | Gain BLOSUM62 (CI) | ESM-2 − BLOSUM62 on Pyrodigal inputs |
| --- | ---: | ---: | ---: | ---: | ---: |
| Backbone-heavy (142, 107) | 52 → 74 | 52 → 73 | +14.4 (6.6, 21.8) | +13.4 (6.0, 20.9) | +0.9 (0.0, 2.8) |
| Multi-CDS (750, 183) | 302 → 409 | 302 → 409 | +9.0 (5.9, 12.2) | +9.2 (6.2, 12.5) | −0.4 (−1.6, 0.8) |
| Single label, other CDSs (125, 72) | 112 → 122 | 112 → 122 | +8.1 (2.8, 14.3) | +8.1 (2.8, 14.3) | 0.0 |
| Single-CDS (98, 70) | 69 → 78 | 69 → 78 | +7.6 (2.4, 13.6) | +7.6 (2.4, 13.6) | 0.0 |
| Tagged or fused (106, 43) | 38 → 42 | 59 → 64 | +4.6 (−4.0, 13.5) | +3.1 (−4.1, 10.5) | **−24.8 (−35.4, −14.8)** |
| Low similarity (250, 154) | 121 → 137 | 137 → 152 | +7.6 (2.9, 12.9) | +7.4 (3.0, 12.3) | **−8.0 (−12.5, −3.9)** |
| Intermediate similarity (426, 199) | 159 → 243 | 159 → 244 | +7.7 (3.7, 12.1) | +8.2 (4.0, 12.7) | −0.5 (−1.5, 0.0) |
| High similarity (508, 262) | 262 → 311 | 270 → 321 | +9.9 (6.7, 13.5) | +10.5 (7.2, 14.3) | −2.2 (−3.7, −1.0) |
| Insufficient same-family coverage (37, 30) | 31 → 34 | 28 → 29 | +8.3 (0.0, 18.3) | +3.3 (0.0, 10.0) | +15.0 (3.3, 28.3) |

**Family concentration** is an exploratory breakdown of the locked holdout data; see `exploratory.by_family.json`. These are ESM-2 recoveries at cap 2 (minimal → Pyrodigal → annotation-derived):

- class A β-lactamase: 240 → 287 → 348 of 357
- aminoglycoside phosphotransferase: 98 → 175 → 278 of 297
- GFP-like: 31 → 40 → 89 of 128
- PAT/BAR acetyltransferase: 2 → 3 → 91 of 93, all in 7 lineage groups. A selection failure in both sequence arms.
- chloramphenicol acetyltransferase: 38 → 37 → 69 of 80

Every other family changed by at most five recoveries.

## Fusion subregions (exploratory)

The fusion set has 14 documented subregions in 13 constructs and 10 lineage groups. Nine are GFP-like. The predefined comparison is at four inputs per construct. Region-supported recovery requires a correct call from an input covering at least 80% of the documented region.

| Inputs at cap 4 | ESM-2 region / family | BLOSUM62 region / family | Inputs | Residues | Share of supporting input outside the region |
| --- | ---: | ---: | ---: | ---: | ---: |
| Intact Pyrodigal candidates | 10 / 10 | 14 / 14 | 48 | 17,318 | 39–45% |
| Pyrodigal + 300-residue windows | 10 / 11 | 10 / 11 | 52 | 15,479 | 22–25% |
| Whole documented CDSs (diagnostic) | 10 / 10 | 14 / 14 | 32 | 13,302 | 42–47% |
| Documented regions (diagnostic) | 14 / 14 | 14 / 14 | 14 | 4,052 | 0% |

ESM-2's four misses with intact inputs were scorer failures on an input that covered the region. Mean pooling over the fusion partner diluted the signal, and the exact region recovered all four. Windowing fixed some of those failures but lost others. Windows used the query budget before reaching two regions, cut across one, and left one scorer miss. ESM-2 therefore stayed at 10. Alignment needs no segmentation: BLOSUM62 found all 14 from intact inputs and fell to 10 with windows. A midpoint or fixed-window split is therefore not justified as the next method. The evidence instead favors local scoring of intact inputs, consistent with the [native-protein study](../benign_protein_pilot/NATIVE_FUSIONS.md). That study's residue-alignment scorer recovered 7 of 8 labels from intact proteins. Ten related groups, mostly GFP fusions, keep this exploratory.

## Low similarity

The construct cohort contains a predefined low-similarity tier: 250 corroborated components with less than 40% identity to every same-family reference at ≥ 70% coverage. They span 154 lineage groups but only about 100 distinct proteins, and about half are GFP-like. In this tier, Pyrodigal's preprocessing gain held with both scorers (+7.6 and +7.4 points). On identical inputs BLOSUM62 beat mean-pooled ESM-2 by 8.0 points (3.9–12.5) with Pyrodigal inputs and by 6.1 points with annotation-derived inputs. This agrees with the [known-protein comparison](../benign_protein_pilot/RELATED_FAMILIES.md), where normalized BLOSUM62 ranked 38 of 41 low-identity labels first, mean-pooled ESM-2 19 and ESM-2 residue alignment 30. On this evidence, mean-pooled ESM-2 did not outperform conventional alignment on distant components at the construct level. The tier's few distinct proteins, its GFP dominance and its 37-component insufficient-coverage subgroup are the main limits. In that subgroup ESM-2 was ahead by 15.0 points (3.3–28.3). Testing whether that subgroup reflects fusion context needs its own analysis.

## Anchor-only tier (exploratory)

The anchor-only tier has 33 components in 25 lineage groups. ESM-2 recovered 15 with minimal inputs, 18 with Pyrodigal and 27 with annotation-derived inputs. BLOSUM62 recovered 19, 23 and 32. The equal-group gains are +9.3 (0.0–21.3) and +13.3 (1.3–28.0). This tier is too small to change the conclusion.

## What remains open

- **Generalization beyond these vectors.** One NCBI vector query underlies the corpus. Common resistance markers dominate the families. A cohort with more eukaryotic, plant and reporter cassettes would test where a prokaryotic gene caller stops helping. For example, PAT/BAR genes were recovered almost only from annotation.
- **A confirmatory test** needs fresh held-out groups. All further analyses of these data, including any new selection rule motivated by the 412 below-cap misses, are exploratory.
- **Construct-level low similarity** has a defined tier here but few distinct proteins. A stronger test needs independently labelled distant components from more families.

## Reproduction and provenance

- **Code:** locked protocol and runner at `deaf88b`; the holdout ran with that revision.
- **Inputs:** holdout FASTA `238684dc…`, references `893a7db7…`, corroborated key `398d6816…`, freeze `1c10becd…`, public freeze `2eb4dada…`.
- **Selected-input fingerprints**, identical across scorers: minimal `28eea256…`, start ranking `c0682dfc…`, Pyrodigal `b00f633b…`, annotation-derived `b5df2c77…`.
- **Private reports** are in ignored storage under `results/paired-cohort-v2/`:

| File | SHA-256 |
| --- | --- |
| `holdout/analysis.corroborated.json` | `74afa039…` |
| `holdout/analysis.anchor_only.json` | `509dfd15…` |
| `holdout/exploratory.by_family.json` | `aeba66ef…` |
| `holdout-fusion/analysis.json` | `00563798…` |

- **Replay:** re-running the analysis from the saved traces reproduced `analysis.corroborated.json` byte for byte.
- **Versions:** Python 3.13.18, Biopython 1.88, PyTorch 2.6.0 CPU, Transformers 4.48.1, Hugging Face Hub 0.28.1, Pyrodigal 3.6.3, NumPy 2.5.5. ESM-2 `facebook/esm2_t6_8M_UR50D` at `c731040f…`, weights `24c5fa47…`.
- **Tests:** 306 passed with one historical optional skip.

The commands are in the [protocol's reproduce section](COHORT_STUDY.md#reproduce). The holdout arms are `minimal`, `start_ranked`, `gene_caller` and `annotation_reference`. The fusion arms are `gene_caller`, `gene_caller_windows`, and `annotation_reference` with `fusion/annotation-units.json` and `documented_regions` with `fusion/documented-regions.json`. Run each for both scorers, then run `cohort_analysis` (`--answer-key …/answer-key.anchor_only.json` for the exploratory tier) and `cohort_fusion analyze`.
