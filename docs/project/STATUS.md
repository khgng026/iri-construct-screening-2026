# Current project status

Updated November 23, 2026. The [research scope](SCOPE.md) and [plan](PLAN.md) define the active direction. This is the current status ledger; historical reports preserve their original findings and conditions.

## Established work

- The annotation-aware `iri-constructs` toolkit is implemented. It validates and exports documented coding translations, and does not infer components from unannotated DNA. See the [release assessment](RELEASE_ASSESSMENT.md) and [contract](TOOLKIT_CONTRACT.md).
- The [component study v1](../results/component-study-v1/README.md) and [extension](../results/component-study-v1/extension-1/README.md) found greater construct coverage when querying more annotated units, with more model calls and no change to per-unit retrieval. An exploratory oracle segment recovered the evaluated fusion family slots; it did not demonstrate a deployable segmentation method.
- The prior high-identity study does not establish model performance on distant homologs or the effect of preprocessing unannotated constructs. Earlier plans and records are indexed in [history](HISTORY.md).

## Research pilot

The paired verifier and two exploratory pilots are ready on `paired-research-pilot`. At two queries per construct, start-aware conversion recovered 4 of 7 documented components versus 3 of 7 for minimal conversion, with ten model calls in each arm. At four queries per construct both recovered 6 of 7. On the separate known-protein pilot, local alignment ranked the correct label first on 7 of 8 low-identity queries versus 5 of 8 for ESM-2.

The related enzyme family comparison now includes 41 low-identity queries. ESM-2 residue alignment ranked 30 labels first versus 19 for mean-pooled ESM-2 and 38 for normalized BLOSUM62. Both ESM-2 arms used the same encoder outputs. [Results and reproduction commands](../../experiments/benign_protein_pilot/RELATED_FAMILIES.md)

The larger construct-level comparison and exploratory fusion results are in [draft PR #26](https://github.com/khgng026/iri-construct-screening-2026/pull/26). The [handoff](HANDOFF.md) records the earlier runs and available artifacts.

## Low-identity follow-up

- **Broad structure labels.** All 72 domains were below 25% global identity to their reference panels. Local ESM-2 8M outperformed sequence methods on two folds, and raw pairwise alignment led on the third. The comparison includes profiles built from the same references. [Study and results](../../experiments/benign_protein_pilot/STRUCTURAL_FAMILIES.md)
- **Exact enzyme reactions.** On 22 natural proteins below 30% identity to their same-reaction references, pairwise alignment and local ESM-2 150M both identified the same 19 labels. HMMER identified 18 and local ESM-2 8M identified 9. The benefits depend on the prediction task and reference set. The low-identity test from construct DNA remains open. [Study and results](../../experiments/benign_protein_pilot/EXPERIMENTAL_REACTIONS.md)
