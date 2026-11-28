# Q2 reporter replication — frozen protocol

**Frozen before scoring.** This document, the panel manifest, the reference manifest, the
grouping, and the analysis plan are committed before any arm is scored. The research
question, eligibility, endpoint, and sensitivity checks are fixed here and are not
changed after seeing results.

## Research question

At the same query budget, does gene-calling preprocessing improve identification of
documented components relative to minimal six-frame preprocessing, on **fresh public
benign reporter constructs** that are independent of the published 777-construct cohort?
Specifically, how much of the published cohort gain (+9.2 equal-group / +12.4 pooled
points with ESM-2 at a cap of two) carries over to a more diverse set of independently
documented reporter lineages, where two bacterial marker families did not dominate by
construction?

This is an independent **replication at pilot scale** using the project's existing
verifier, candidate generators, and alignment scorer. It is self-review on public,
benign material; it establishes family identification, not biological function or risk.

## Data and provenance

- **Source:** NCBI GenBank nuccore records fetched by EFetch (`curate.py`), from eight
  committed reporter-focused esearch queries. Source sequences live only in the ignored
  results tree; audits are sequence-free.
- **Independence:** every candidate accession is checked against the 792 cohort
  accessions (`exclude_cohort.txt`, from `cohort_v2.freeze.json`) and the 5 reporter-pilot
  accessions; matches (including version-stripped) are dropped before fetching.
- **Translation validation:** each documented CDS is validated against its deposited DNA
  by `iri_constructs.process_inputs`; only records whose every CDS is `exported`
  (translation agrees with coordinates) are admitted.
- **Labels:** the family of each component is taken from the deposited `/product` (or
  `/gene`) qualifier via a fixed lexicon (`curate.LEXICON`) — documentation, independent
  of the scorer under test. **Limitation:** labels rest on depositor annotation, a weaker
  corroboration than the cohort's reviewed >= 90% anchor; see HANDOFF for the stronger
  corroboration a confirmatory study would add.

## Eligibility (admission) rules

A construct is admitted only if, after toolkit validation:
1. it is not in the exclusion set;
2. every documented CDS is `exported` and maps to a family in the lexicon;
3. at least one exported CDS is a reporter that is model-eligible (30–1022 aa); and
4. every component family has an independent reviewed reference in the panel.

Rule 4 excludes the engineered **NanoLuc** family (no suitable independent reviewed
reference) and therefore its 13 constructs; this shortfall is reported, not worked around.

## Reference set

Reviewed UniProtKB/Swiss-Prot entries, one or more per family, chosen as canonical family
representatives independent of the deposited constructs (`references.manifest.json`;
`freeze.REFERENCE_ACCESSIONS`). `gfp_like` carries green (Aequorea GFP) and red (DsRed,
eqFP611) representatives because the deposited "fluorescent protein" family spans both
GFP-fold subtypes at low mutual identity.

## Arms and scorer

Reusing `experiments.paired_preprocessing.reporter_pilot` unchanged:

- `minimal` — length-ranked six-frame translation blocks (baseline).
- `start_ranked` — the same pool reordered by start evidence (control for reordering).
- `gene_caller` — Pyrodigal metagenomic-mode predicted proteins (candidate).
- `annotation_reference` — validated documented proteins (diagnostic ceiling).

Scorer: **normalized local BLOSUM62 alignment** (`characterizer: blosum62`), the project's
conventional comparator. It runs fully offline. The pinned ESM-2 arm is optional and is
not required for this replication (no cached weights in this environment); it would check
whether the preprocessing gain transfers across scorers.

## Endpoint and target

- **Primary endpoint:** documented-component recall at a query cap of **two** protein
  inputs per construct, weighted by **equal lineage group**. A recovery is a correct
  family call from an input covering >= 80% of the documented component (verifier rule).
- **Practical target:** +5 percentage points (the cohort's approved effect size), judged
  by the lower bound of the 95% whole-group bootstrap interval.
- Pooled recall (every component equal) is reported alongside, with explicit denominators.
- Caps **1** and **4** are reported as budget sensitivity.

## Grouping rules

`analyze.py group`, computed and frozen before scoring:

- **Lineage group (primary weight):** single-linkage clusters of components that share
  >= 90% protein identity *within a family*. 169 components → 21 lineage groups
  (14 reporter, 7 marker).
- **Family group (stricter sensitivity):** one group per family (9 groups).
- A large near-identical series (e.g. 39 EGFP-type `gfp_like#1`) therefore contributes
  one lineage group, limiting how much a single ancestry can drive the headline.

## Analysis plan (fixed before scoring)

1. Report pooled and equal-lineage-group recall per arm at each cap, with denominators.
2. Primary effect: `gene_caller − minimal` at cap 2, equal-group and pooled, with 95%
   whole-group bootstrap intervals (4000 replicates, seed 20261124), for **all**,
   **reporter-only**, and **marker-only** components.
3. Sensitivity: omit the dominant lineage group; stricter family-group weighting;
   caps 1 and 4; reporter-only grouping.
4. Coverage funnel per arm: components with no covering candidate, covering-but-missed,
   recovered. Resources: queries, input residues, candidate pool, scorer runtime.
5. Interpretation: whether the gain carries to reporter lineages, whether it survives
   removing the dominant lineage and the budget sensitivity, and which failure mode
   (selection vs characterization) drives it. Any post-result experiment is versioned
   separately from this frozen replication.

## Denominators (frozen)

86 constructs; 169 components (93 reporter, 76 marker); 21 lineage groups (14 reporter,
7 marker); 9 families; archetypes — multi_cds 33, backbone_heavy 27, single_cds 23,
tagged_fusion 3.
