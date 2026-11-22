# Current project status

Updated November 22, 2026. The [research scope](SCOPE.md) and [plan](PLAN.md) define the active direction. This is the current status ledger; historical reports preserve their original findings and conditions.

## Established work

- The annotation-aware `iri-constructs` toolkit is implemented. It validates and exports documented coding translations, and does not infer components from unannotated DNA. See the [release assessment](RELEASE_ASSESSMENT.md) and [contract](TOOLKIT_CONTRACT.md).
- The [component study v1](../results/component-study-v1/README.md) and [extension](../results/component-study-v1/extension-1/README.md) found greater construct coverage when querying more annotated units, with more model calls and no change to per-unit retrieval. An exploratory oracle segment recovered the evaluated fusion family slots; it did not demonstrate a deployable segmentation method.
- The prior high-identity study does not establish model performance on distant homologs or the effect of preprocessing unannotated constructs. Earlier plans and records are indexed in [history](HISTORY.md).

## Research pilot

The paired verifier and two exploratory pilots are ready on `paired-research-pilot`. At two queries per construct, start-aware conversion recovered 4 of 7 documented components versus 3 of 7 for minimal conversion, with ten model calls in each arm. At four queries per construct both recovered 6 of 7. On the separate known-protein pilot, local alignment ranked the correct label first on 7 of 8 low-identity queries versus 5 of 8 for ESM-2.

The reporter comparison v1 tested two more preprocessing arms. At the primary cap, Pyrodigal recovered 6 of 7 components with eight calls. Reordering the full candidate pool by start evidence recovered 4 of 7 with ten calls. Both results are exploratory. [Comparison and reproduction commands](../../experiments/paired_preprocessing/REPORTER_COMPARISON.md)

The next work is to run the arms on a verified batch from the existing corpus and add a fusion segmentation arm when independent subregion labels are available. A stronger low-similarity panel can proceed in parallel. The [pilot handoff](HANDOFF.md) records the earlier runs and available artifacts.
