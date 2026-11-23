# Grouped construct cohort v2: preprocessing gain at a matched query budget

**Question.** When component annotations are held out, does Pyrodigal preprocessing help an existing protein scorer identify more documented components than the minimal six-frame conversion at the same allowed query budget? How large is the gain, where does it occur, and does it persist when BLOSUM62 alignment replaces ESM-2?

This document was committed as the protocol lock (`deaf88b`) before any holdout inference. The results are in a separate [results report](COHORT_RESULTS.md).

## Checkpoint A: readiness

| Item | State |
| --- | --- |
| Code | `paired-research-pilot` (crossed comparison, corpus audit; tip `c18b110`) merged with `main` (grouped reporting, related-family and native-fusion studies; tip `7c5b87b`) on the session branch at `da3e55b`. Neither branch contained the other. Only `STATUS.md` conflicted. |
| Environment | Python 3.13.18, Biopython 1.88 (earlier runs used 1.86), PyTorch 2.6.0 CPU, Transformers 4.48.1, Hugging Face Hub 0.28.1, Pyrodigal 3.6.3, NumPy 2.5.5. ESM-2 8M at revision `c731040f…` with weight SHA-256 `24c5fa47…`, matching the pinned hash. |
| Private artifacts at handoff | None were present in this workspace. The reporter sources, model cache, corpus snapshots and PR #19 audit were all reacquired. |
| Reporter panel | Reacquired byte-identically: input `b73d9309…` and references `4a54303a…` match the frozen protocol. |
| Frozen reporter result | Reproduced exactly at the primary cap: minimal 3/7 (10 queries, 3,344 residues) and Pyrodigal 6/7 (8 queries, 2,486 residues) with both scorers. The selected-input fingerprints match across scorers for all five arms. |
| New runner | `cohort_study.run` reproduces the old runner's selected-input fingerprints, queries, candidate counts and predictions for all four reused arms and both scorers on the reporter panel. Caching identical inputs changes cost, not results. |
| Corpus | All 50 batch files (1,000 records) are byte-identical to the frozen manifests. The audit reproduces PR #19: 2,338/2,338 units and labels, 1,311/1,311 anchor regions re-translate, 792 eligible records in 63 shared-protein groups (largest 706), and 14 documented fusion subregions in 13 records across 7 groups. `check_corpus_audit` runs clean. |
| Missing prerequisites | The cohort `select`/`prepare` stages named in `corpus_batch.py` were not implemented. They are implemented here as `cohort_study freeze/prepare`. No other prerequisite was missing. |
| Tests | 306 passed with one historical optional skip (`tests/test_cost_quality.py`, fictional plot fixture). |

## Checkpoint B: cohort

**Eligibility** follows the PR #19 shortlist rule without changes. A component is an exported CDS whose deposited translation the toolkit validates against the DNA. Its family must be supported by a reviewed UniProt anchor (identity ≥ 0.90, anchor coverage ≥ 0.80) in a family with a reference panel. It must not be a near-copy of a panel reference, and its deposited name must not name a different family. The five reporter-pilot records are excluded. Undocumented CDSs are outside every denominator.

**Evidence tiers** keep translation evidence and family evidence separate. Every eligible component has a toolkit-validated translation. The **corroborated** tier has an anchor match and a deposited CDS name or note naming the same family. It is the primary denominator. The **anchor-only** tier has an anchor match while the deposit is silent or names only a partner. It is reported separately as exploratory.

**Groups** were fixed from evaluator metadata before scoring:

- *Lineage (primary analysis groups).* Records are linked when they share a direct-submission event (authors and submission text), a published reference title, or an unpublished or catalog reference title. They are also linked when they are near-copy constructs: at least 90% of each record's sampled canonical 25-mers occur in the other. Groups are connected components.
- *Conservative (sensitivity).* Lineage links plus any identical eligible protein shared by two records.
- *Partition.* Development and holdout are assigned by conservative component. No eligible protein or lineage crosses partitions. The giant conservative component is in holdout. Of the other components, the lowest 20% by `sha256(20261123:first accession)` form development, together with the five reporter constructs.

Transitive linkage through common markers and backbones (TEM-1 β-lactamase, kanamycin resistance, the pUC backbone) chains almost every vector into one component. Before deciding the rule, a metadata-only comparison was run. It found a largest component with 52% of components when DNA containment ≥ 0.8 was added, and 97% when containment ≥ 0.5 was added. Submission-level lineage gives several hundred groups. The conservative view is therefore reported, but it cannot support inference.

| | Development | Holdout |
| --- | ---: | ---: |
| Records | 15 (+5 reporter) | 777 |
| Eligible components (corroborated / anchor-only) | 14 / 3 | 1,221 / 33 |
| Lineage groups with corroborated components | 11 | 422 |
| Largest lineage group, share of corroborated components | 14% | 20% |
| Effective lineage groups (inverse Simpson on component share) | 9.8 | 23.9 |
| Conservative groups (largest share) | 9 (21%) | 39 (96%) |
| Families with suitable references (of 28 scorable) | 7 | 27 |

Holdout records by archetype: 355 multi-CDS, 144 backbone-heavy, 128 single-label-other (one labelled CDS plus unlabelled CDSs, labelled share ≥ 15%), 99 single-CDS and 51 tagged or fused. Corroborated components by archetype: 750, 142, 125, 98 and 106. Class A β-lactamase (361) and aminoglycoside phosphotransferase (309) account for over half the components.

**Low-similarity tier (defined before scoring).** For each component, compute a local BLOSUM62 alignment against each panel reference, using the anchor region for fusion-class components. The tier is *low* when the best same-family identity is below 40% with ≥ 70% coverage of both proteins. These are the related-family study's thresholds. The *intermediate* tier is below 60%, and *high* is 60% or more. In holdout: 260 low (250 corroborated), 432 intermediate, 516 high, and 46 with insufficient coverage. The low tier spans 159 lineage groups but only 106 distinct proteins, and GFP-like components make up 124 of the 260. No component has a competing-family reference at least as close as its best same-family reference. None reaches ≥ 0.90 identity at adequate coverage to any reference.

**Repeated relationships.** The 3,693 near-copy construct pairs are inside lineage groups by construction. One reporter-development protein recurs identically in a holdout construct. No eligible protein crosses development and holdout.

**Fusion subregions.** All 14 documented subregions (13 records, 10 lineage groups; 9 GFP-like) are in holdout. Development has no documented subregion, so the fusion code path was integration-tested on six development anchor-region boundaries only.

**Can the cohort resolve five points?** The equal-group-weighted gain averages 422 lineage groups. Assuming a per-group standard deviation of paired differences near 0.4 (an assumption, not an estimate from holdout outcomes), its bootstrap half-width would be roughly ±4 points. A lower bound of five points is therefore reachable only if the gain is roughly nine points or more. A smaller gain would be inconclusive. Pooled recall has about 24 effective groups and is less precise. Under the conservative rule the holdout is effectively one group. That analysis is descriptive.

## Checkpoint C: development comparison

Fourteen corroborated components were scored at caps 1/2/4. Each cell gives recoveries, followed by actual inputs and input residues.

| Arm | ESM-2 cap 1 / 2 / 4 | BLOSUM62 cap 1 / 2 / 4 | Inputs at cap 2 (allowed 30) | Residues at cap 2 |
| --- | --- | --- | ---: | ---: |
| Minimal | 4 / 6 / 5 | 5 / 8 / 9 | 30 | 13,218 |
| Start ranking | 4 / 7 / 5 | 5 / 9 / 9 | 30 | 13,183 |
| Pyrodigal | 4 / 7 / 9 | 5 / 8 / 11 | 30 | 11,698 |
| Annotation-derived | 11 / 13 / 13 | 12 / 14 / 14 | 17 | 5,970 |

Failure breakdown at the primary cap:

| Finding | Minimal | Start ranking | Pyrodigal | Decision |
| --- | ---: | ---: | ---: | --- |
| No compatible candidate in the pool | 0 | 0 | 1 | Pyrodigal called no gene over one eukaryotic-context GFP fusion locus. This is a candidate-generation limit of a prokaryotic gene caller, not an implementation error. It is stated as a method limit. |
| Compatible candidate below the cap | 6 | 5 | 5 | This was the dominant failure for every sequence arm. Start ranking moved one into the budget. Pyrodigal calls about seven genes per development construct, so length order still pushes documented genes past two queries. This is a selection effect that the matched-budget comparison is meant to measure, so the caps stay. |
| Covering input reached the scorer, label missed | ESM-2 2, BLOSUM62 0 | 2 / 0 | 1 / 0 | Best ranks of the correct family were 4–9 (ESM-2 only). The annotation-derived input recovered one of them, so clean boundaries help ESM-2 in a minority of cases. |
| Unexpected input, budget or denominator differences | none | none | none | Fingerprints and candidate pools matched across scorers, and every arm used the same 14 components. |

**One verified implementation fix was made before the lock.** The annotation-derived arm dropped a documented 1,060-residue GFP fusion that exceeded the model limit, while the sequence arms split long candidates into 1,022-residue chunks. The diagnostic now applies the same chunking rule (`chunk_unit`), so that component is recovered. The sequence arms are unchanged.

**A noted interaction, unchanged.** With a fixed limit of five returned family calls and turn-taking across queries, more queries can displace a correct second-ranked call. ESM-2 minimal recall falls from 6 to 5 between caps 2 and 4. The output limit stays fixed across arms as specified, and budget curves are reported as observed.

## Locked protocol

- **Inputs.** Development FASTA `aa1dab79…`, holdout FASTA `238684dc…`, holdout fusion FASTA `3d6a39d2…`, references `893a7db7…` (134 reviewed UniProt sequences, 28 families; the component-study panel). Every private file hash is in [`cohort_v2.freeze.json`](cohort_v2.freeze.json) and [`cohort_protocol.json`](cohort_protocol.json).
- **Arms.** Minimal six-frame conversion (baseline); start ranking (control: same pool, reordered); Pyrodigal 3.6.3 in metagenomic mode (preprocessing test); annotation-derived documented units (diagnostic, using the same model length rule).
- **Scorers.** Mean-pooled ESM-2 8M cosine and reference-self-normalized local BLOSUM62, both taking the best reference per family. The same reference panel, limit of five calls, turn-taking rule, 80% same-strand locus overlap rule and one-to-one matching apply in every arm. The two scorers must receive identical selected inputs and candidate pools, checked by fingerprint.
- **Budgets.** Caps of 1, 2 and 4 protein inputs per construct, with 2 as the primary cap. Report allowed and actual inputs, input residues, candidate pools, ESM-2 encoder calls, alignment comparisons, alignment cells and time.
- **Primary endpoint.** Recall of corroborated eligible components at cap 2.
- **Primary quantities at cap 2.** (1) Pyrodigal − minimal with ESM-2; (2) the same with BLOSUM62; (3) ESM-2 − BLOSUM62 on identical inputs for each arm. Also report the difference in preprocessing gain between scorers.
- **Estimators.** The equal-lineage-group-weighted paired gain is primary. Pooled gain, leave-one-group-out ranges and the gain without the largest group are secondary.
- **Uncertainty.** Resample whole lineage groups with paired outcomes: 4,000 replicates, seed 20261123, 95% percentile intervals. Sensitivity checks are the same analysis on conservative groups (descriptive) and a two-way bootstrap that resamples lineage groups and families independently.
- **Decision rule (target five points), applied to the equal-group gain.**
  - A **practically meaningful benefit** requires an interval lower bound ≥ 5 points and a leave-one-group-out minimum ≥ 5 points.
  - A **positive signal, not established** has an interval lower bound > 0 while the interval still spans 5 points.
  - A **target-size gain is unsupported** when the interval upper bound is < 5 points.
  - Otherwise the result is **inconclusive**.
  - **Persistence with alignment** applies the same rule to BLOSUM62. A **scorer effect** is reported when its interval excludes zero.
- **Strata.** Archetype and similarity tier, including the predefined low-similarity tier, each with its lineage-group denominator. The anchor-only tier is analyzed separately as exploratory.
- **Failure categories.** No compatible candidate in the pool; compatible candidate below the cap; covered input with the label missed (with the correct family's rank); recovered. Report transitions from baseline misses to each arm.
- **Fusion (exploratory, 14 subregions across 10 lineage groups).** Compare intact Pyrodigal candidates with fixed windows: candidates over 300 residues become 300-residue windows at stride 150, with the last ending at the candidate end, and windows count against the cap. Documented regions and whole documented CDSs are diagnostics. Caps are 1/2/4/8, with 4 as the predefined comparison. Report region-supported recovery, family recovery, the miss category (region not reached, split cut across the region, scorer missed with a covering input), queries, residues and the share of each supporting input outside the region.
- **Holdout discipline.** The holdout is run once with this commit's code. Any later method or threshold is labelled exploratory. A confirmatory retest would need fresh held-out groups.

## Reproduce

Use the model environment above with `HF_HOME` pointing at the pinned ESM-2 cache. Stages write only to new paths under ignored `results/`.

```bash
python -m experiments.paired_preprocessing.corpus_batch acquire
```

```bash
python -m experiments.paired_preprocessing.corpus_batch audit
```

```bash
python -m experiments.paired_preprocessing.cohort_study freeze
```

```bash
python -m experiments.paired_preprocessing.cohort_study prepare
```

```bash
python -m experiments.paired_preprocessing.cohort_study fusion-units --partition-dir results/paired-cohort-v2/prepared/holdout
```

Each arm and scorer is one stage. First write the protocol for a partition and scorer (`holdout`, `holdout_fusion` or `development`):

```bash
python -m experiments.paired_preprocessing.cohort_study protocol --inputs-key holdout --scorer esm2 \
  --out results/paired-cohort-v2/holdout/esm2/protocol.json
```

```bash
python -m experiments.paired_preprocessing.cohort_study run --method gene_caller \
  --protocol results/paired-cohort-v2/holdout/esm2/protocol.json \
  --inputs results/paired-cohort-v2/prepared/holdout/method_inputs/constructs.fasta \
  --references results/paired-cohort-v2/prepared/references.private.json \
  --out results/paired-cohort-v2/holdout/esm2/gene_caller.json
```

```bash
python -m experiments.paired_preprocessing.cohort_analysis --run results/paired-cohort-v2/holdout \
  --answer-key results/paired-cohort-v2/prepared/holdout/evaluator/answer-key.corroborated.json \
  --construct-map results/paired-cohort-v2/prepared/holdout/evaluator/construct-map.private.json \
  --out results/paired-cohort-v2/holdout/analysis.corroborated.json
```

```bash
python -m experiments.paired_preprocessing.cohort_fusion analyze --run results/paired-cohort-v2/holdout-fusion \
  --key results/paired-cohort-v2/prepared/holdout/fusion/subregion-key.json \
  --out results/paired-cohort-v2/holdout-fusion/analysis.json
```
