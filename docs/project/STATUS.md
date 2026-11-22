# Current project status

Updated November 22, 2026. The [research scope](SCOPE.md) and [plan](PLAN.md) define the active direction. This is the current status ledger; historical reports preserve their original findings and conditions.

## Established work

- The annotation-aware `iri-constructs` toolkit is implemented. It validates and exports documented coding translations, and does not infer components from unannotated DNA. See the [release assessment](RELEASE_ASSESSMENT.md) and [contract](TOOLKIT_CONTRACT.md).
- The [component study v1](../results/component-study-v1/README.md) and [extension](../results/component-study-v1/extension-1/README.md) found greater construct coverage when querying more annotated units, with more model calls and no change to per-unit retrieval. An exploratory oracle segment recovered the evaluated fusion family slots; it did not demonstrate a deployable segmentation method.
- The prior high-identity study does not establish model performance on distant homologs or the effect of preprocessing unannotated constructs. Earlier plans and records are indexed in [history](HISTORY.md).

## Active next work

The scope defines four parallel streams: evaluation and verification, sequence-only preprocessing, baseline and model characterization, and low-similarity evaluation. The first integration goal is a budget-aware paired pilot using withheld annotations. Low-similarity evaluation is a required research outcome. These streams are planned, not yet completed under the new scope.

The existing annotation-aware toolkit can support the answer key and annotation-derived comparison; it is not the sequence-only preprocessing method. Use new versions for new experiments, and preserve frozen historical artifacts.
