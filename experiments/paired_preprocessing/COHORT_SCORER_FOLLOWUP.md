# Exploratory follow-up: stronger ESM-2 configurations and the identity crossover

Exploratory. These analyses reuse the **public** cohort v2 holdout, whose answer keys are
committed, so they cannot serve as a confirmatory test; a confirmatory result needs fresh
held-out groups. They were run to answer a reader question: *does the cohort-v2 finding
(mean-pooled ESM-2 8M no better than BLOSUM62) mean foundation models do not beat sequence
alignment once a preprocessing toolkit exists?* The short answer: **not in general — the two
methods trade places by identity regime, and residue-level ESM-2 decisively beats alignment
once no close reference is available.** Motivation and the regime boundary come from
Wittmann et al., _Science_ 2025 and _Front. Bioeng. Biotechnol._ 2026; see the
[research design](../../docs/project/ALIGNMENT_VS_MODELS_PLAN.md).

## Scorers compared

All share the locked cohort-v2 selected inputs (fingerprints identical across scorers; the
8M mean-pooled traces reproduce the locked ESM-2 holdout traces byte for byte — predictions
and full rankings). Only the characterization step differs.

- **Normalized BLOSUM62** — the locked conventional comparator.
- **Mean-pooled ESM-2** at 8M and 150M — cosine of mean residue embeddings.
- **Residue-level ESM-2** at 8M and 150M — local alignment over the embedding cosine
  matrix, the related-family study's scorer and gap settings.

## Finding 1: on the whole cohort, no ESM configuration beats alignment

Documented-component recall at the primary cap (2 inputs), gene-caller inputs, 1,221
components in 422 lineage groups:

| Scorer | Recovered | vs BLOSUM62, equal-group gain (95% CI) |
| --- | ---: | ---: |
| Normalized BLOSUM62 | 746 | — |
| ESM-2 8M mean | 725 | −2.5 (−3.9, −1.2) |
| ESM-2 8M residue | 734 | −1.4 (−2.7, −0.2) |
| ESM-2 150M mean | 727 | −2.1 (−3.6, −0.7) |
| ESM-2 150M residue | 735 | −1.3 (−2.5, −0.1) |

Residue-level scoring adds about +1 point over mean pooling; a 4× larger model adds about
+0.2. Neither closes the gap to alignment. This extends the cohort-v2 result across model
size and representation.

## Finding 2: the methods barely overlap in their errors, except where alignment has no foothold

The oracle union of any ESM configuration with BLOSUM62 exceeds BLOSUM62 alone by at most
6 of 1,221 components. ESM-2 is nearly a subset of alignment here: at the primary cap it
uniquely recovers 5–6 components, while alignment uniquely recovers 16–26. Every component
ESM uniquely recovers lies in the **insufficient-coverage tier** — components with no
reference aligning over ≥ 70% of both proteins, the closest in-cohort analog of the
alignment-failure regime. There ESM strictly dominates: 34/37 (all four configurations) vs
29/37 for BLOSUM62, +5 and −0; all five are phage integrase, a divergent recombinase family.
So on this cohort the gap is not closeable by ensembling — the regime where the model wins is
barely present.

Reproduce: `scorer_comparison.py` writes the recall table, tier breakdown, oracle union and
complementarity (`scorer-comparison.json`).

## Finding 3 (headline): the gap reverses below ~35% identity, and residue-level scoring is what realizes it

The cohort hides the regime because most components have a close reference. The
leave-close-references-out ablation (`cohort_ablation.py`) manufactures the low-identity
regime from the same benign data: for each component's validated protein, strip every
reference at or above a falling identity threshold, then recompute each method's top-family
recall over the components that still retain a same-family reference. Equal-group recall with
95% whole-group bootstrap intervals:

| Highest retained-reference identity | n (groups) | BLOSUM62 | ESM-2 8M mean | ESM-2 8M residue | ESM-2 150M mean | ESM-2 150M residue |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| any (< 101%) | 1221 (422) | 96% | 96% | 97% | 95% | 97% |
| < 50% | 1211 (421) | 96% | 93% | 97% | 92% | 97% |
| < 40% | 1188 (419) | 89% | 87% | 97% | 90% | 97% |
| < 35% | 1102 (403) | 88% | 91% | **97%** | 86% | **97%** |
| **< 30%** | 639 (264) | **68%** | 83% | 91% | 74% | **96%** |
| < 25% | 71 (46) | 92% | 100% | 100% | 76% | 100% |

At **< 30% identity**, residue-level ESM-2 recovers **90.9% [87.6, 93.9]** (8M) and
**95.6% [93.3, 97.5]** (150M) of families, against BLOSUM62's **68.2% [63.1, 73.2]** —
intervals do not overlap, a reversal of 23–27 points. At < 35% the residue scorers are
already clear of alignment (97% vs 88%).

Two sub-findings matter for *how* to close the gap:

1. **Representation, not scale, is the lever.** Residue-level scoring is decisive;
   mean pooling is not. Scaling mean-pooled ESM-2 from 8M to 150M makes it *worse* in the
   low-identity regime (74% vs 83% at < 30%), because averaging discards the local signal the
   larger model carries. Residue-level local alignment over the embeddings recovers it.
2. **This explains why residue-level looked useless on the full cohort.** Its advantage only
   appears once close references are removed. Above ~40% identity every method is near
   ceiling, so the residue scorer's benefit is invisible in the pooled cohort number.

![crossover](crossover.png)

## What this does and does not establish

- It **does** show, on independently reused benign data, that the alignment-vs-model verdict
  is regime-dependent, with a crossover near 30–35% nearest-reference identity, and that
  residue-level ESM-2 (not mean pooling, not model scale) is what reverses it. This reconciles
  cohort v2 with the biosecurity literature's thesis.
- It **does not** establish a confirmatory result: the holdout keys are public, the query is
  the documented unit rather than the full preprocessed pipeline, labels still derive from
  anchor alignment, and ESM-2 saw these families in pretraining, so a share of the low-identity
  recall may be memorization. The [research design](../../docs/project/ALIGNMENT_VS_MODELS_PLAN.md)
  addresses each: independent HMM/structure labels, a profile-HMM comparator (stronger than
  pairwise BLOSUM62, which may itself lift the alignment curve in this regime), fresh held-out
  groups, and the full preprocessing pipeline rather than the documented unit.

The immediate consequence for the main study: the cohort-v2 conclusion ("mean-pooled ESM-2
8M does not beat alignment") stands **for that cohort and that configuration**, but it should
not be read as "foundation models do not help." The help is concentrated below ~30% identity
and requires residue-level scoring — exactly the regime the next phase is built to reach.

## Reproduce

Model environment as in the [cohort study](COHORT_STUDY.md), with ESM-2 150M at revision
`a695f604…` (weights `c3f1da8a…`) also cached.

```bash
# Re-score the locked selected inputs with each configuration (one model per stage)
python -m experiments.paired_preprocessing.cohort_scorers collect --out results/paired-cohort-v2-scorers/proteins.json
python -m experiments.paired_preprocessing.cohort_scorers score --model esm2_8m   --proteins results/paired-cohort-v2-scorers/proteins.json --references results/paired-cohort-v2/prepared/references.private.json --out results/paired-cohort-v2-scorers/rankings.esm2_8m.json
python -m experiments.paired_preprocessing.cohort_scorers score --model esm2_150m  --proteins results/paired-cohort-v2-scorers/proteins.json --references results/paired-cohort-v2/prepared/references.private.json --out results/paired-cohort-v2-scorers/rankings.esm2_150m.json
# then `trace` each scorer, and `scorer_comparison` for the recall/complementarity table.
```

```bash
# Identity crossover (independent of the traces)
python -m experiments.paired_preprocessing.cohort_ablation score --model esm2_8m  --out results/paired-cohort-v2-scorers/ablation.scores.esm2_8m.json
python -m experiments.paired_preprocessing.cohort_ablation score --model esm2_150m --out results/paired-cohort-v2-scorers/ablation.scores.esm2_150m.json
python -m experiments.paired_preprocessing.cohort_ablation analyze --scores results/paired-cohort-v2-scorers/ablation.scores.esm2_8m.json results/paired-cohort-v2-scorers/ablation.scores.esm2_150m.json --out results/paired-cohort-v2-scorers/ablation.crossover.json
```

Committed sequence-free artifacts (hashes only, no sequence): `scorer-comparison.json`,
`ablation.crossover.json`, `ablation.scores.esm2_8m.json`, `ablation.scores.esm2_150m.json`,
`crossover.png`.
