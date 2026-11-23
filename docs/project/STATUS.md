# Current project status

Updated November 23, 2026. The [research scope](SCOPE.md) and [plan](PLAN.md) define the active direction. This is the current status ledger; historical reports preserve their original findings and conditions.

## Established work

- The annotation-aware `iri-constructs` toolkit is implemented. It validates and exports documented coding translations, and does not infer components from unannotated DNA. See the [release assessment](RELEASE_ASSESSMENT.md) and [contract](TOOLKIT_CONTRACT.md).
- The [component study v1](../results/component-study-v1/README.md) and [extension](../results/component-study-v1/extension-1/README.md) found greater construct coverage when querying more annotated units, with more model calls and no change to per-unit retrieval. An exploratory oracle segment recovered the evaluated fusion family slots; it did not demonstrate a deployable segmentation method.
- The prior high-identity study does not establish model performance on distant homologs or the effect of preprocessing unannotated constructs. Earlier plans and records are indexed in [history](HISTORY.md).

## Research pilot

The paired verifier and two exploratory pilots are ready on `paired-research-pilot`. At two queries per construct, start-aware conversion recovered 4 of 7 documented components versus 3 of 7 for minimal conversion, with ten model calls in each arm. At four queries per construct both recovered 6 of 7. On the separate known-protein pilot, local alignment ranked the correct label first on 7 of 8 low-identity queries versus 5 of 8 for ESM-2.

The reporter comparison v1 tested two more preprocessing arms. At the primary cap, Pyrodigal recovered 6 of 7 components with eight calls. Reordering the full candidate pool by start evidence recovered 4 of 7 with ten calls. Both results are exploratory. [Comparison and reproduction commands](../../experiments/paired_preprocessing/REPORTER_COMPARISON.md)

The crossed reporter comparison scored the same candidates with ESM-2 and normalized BLOSUM62. Both recovered the same components at every cap. At the primary cap each recovered 6 of 7 with Pyrodigal and 3 of 7 with the minimal arm. [Results and reproduction commands](../../experiments/paired_preprocessing/CROSSED_CHARACTERIZATION.md)

The related enzyme family comparison now includes 41 low-identity queries. ESM-2 residue alignment ranked 30 labels first versus 19 for mean-pooled ESM-2 and 38 for normalized BLOSUM62. Both ESM-2 arms used the same encoder outputs. [Results and reproduction commands](../../experiments/benign_protein_pilot/RELATED_FAMILIES.md)

## Grouped construct cohort v2

The fresh cohort comparison is complete. Its protocol was locked at `deaf88b` before holdout inference. The cohort has 777 held-out public benign vectors with 1,221 corroborated documented components, spanning 422 lineage groups and 27 families.

At two protein inputs per construct, Pyrodigal preprocessing beat minimal six-frame conversion with both scorers. With ESM-2 the equal-lineage-group gain was **9.2 percentage points (95% whole-group bootstrap 6.6–11.8)**. With normalized BLOSUM62 it was **9.0 (6.5–11.6)**. Pooled recall rose from 46.9% to 59.4% with ESM-2. Pyrodigal used slightly fewer inputs and 7.6% fewer input residues. Both gains meet the locked five-point criterion.

The gain comes from candidate selection: scorers gained almost identically, with a 0.2-point difference. On identical inputs, BLOSUM62 recovered about 2.3 points more than mean-pooled ESM-2. That gap is largest for fusions and for the predefined low-similarity tier, where it is 8.0 points.

The claim is qualified in four ways. About 82% of the extra recoveries are β-lactamase or kanamycin-resistance genes, and the two-way group-by-family interval is 2.8–23.6. The conservative shared-protein view is effectively one group. Tagged or fused constructs are inconclusive. At one input per construct the gain is 4.0 points.

In the fusion subset (14 documented regions in 10 groups, exploratory), fixed windows did not help. ESM-2's misses were scorer failures on intact covering inputs, which points to local scoring rather than splitting. [Results](../../experiments/paired_preprocessing/COHORT_RESULTS.md) · [Protocol](../../experiments/paired_preprocessing/COHORT_STUDY.md)

The grouped results and sequence-free traces are on `main`. The holdout answer keys are public, so a confirmatory test needs fresh held-out groups. [Grouped cohort PR](https://github.com/khgng026/iri-construct-screening-2026/pull/26)

## Exploratory follow-up: scorer configurations and the identity crossover

A reader asked whether cohort v2 means foundation models do not beat alignment once
preprocessing exists. Re-scoring the locked holdout with residue-level ESM-2 and ESM-2 150M,
and a leave-close-references-out ablation, shows the verdict is regime-dependent. On the whole
cohort no ESM configuration beats BLOSUM62 (residue +1 pt over mean pooling, 150M +0.2 pt over
8M; all still below alignment). But as close references are stripped, the methods cross: at
**< 30% nearest-reference identity**, residue-level ESM-2 recovers **95.6% [93.3, 97.5]** (150M)
of families against BLOSUM62's **68.2% [63.1, 73.2]**. The lever is residue-level scoring, not
model scale — mean pooling gets worse when scaled. This reconciles cohort v2 with the
biosecurity literature. It is exploratory (public holdout keys, alignment-derived labels,
pretraining overlap). The [research design](ALIGNMENT_VS_MODELS_PLAN.md) and
[follow-up results](../../experiments/paired_preprocessing/COHORT_SCORER_FOLLOWUP.md) set out the
confirmatory version: independent labels, a profile-HMM comparator, and fresh held-out groups.

**Next work, all exploratory or new.** Test a selection rule for the 412 documented genes Pyrodigal leaves below the cap. Evaluate local or residue-level scoring of intact fusion candidates. Build a cohort with more eukaryotic and reporter cassettes and more distinct distant components. The [pilot handoff](HANDOFF.md) records the earlier runs.

## End-to-end test on distant reporters

Tested 13 deposited reporter DNA records using the same query budgets for minimal preparation and gene calling. On seven below-25% cases from two reporter ancestries, gene calling delivered a covering input for six at a cap of two. Raw alignment, HMMER, and local ESM-2 recovered the same six. All three identified every case from its documented protein. The experiments extend the low-identity comparison to whole constructs and show how candidate selection and baseline choice affect the result. Independent low-identity construct families with documented failures of the strongest conventional scorers are still needed. [Study and paired results](../../experiments/paired_preprocessing/DISTANT_REPORTERS.md)
