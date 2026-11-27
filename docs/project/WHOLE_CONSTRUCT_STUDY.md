# Whole-construct characterization study plan

Test whether gene calling followed by ESM-2 recovers more documented protein-coding components than either minimal translation followed by the same ESM-2 scorer or BLASTx searching raw DNA. The primary endpoint is whole-construct characterization with the full candidate pools. Prioritize mid- and low-identity components and keep a higher-identity comparison for context. The benchmark uses existing, public benign research records and experimentally characterized redesigns when available.

## Starting evidence

At two inputs per construct, the historical ESM-2 run recovered 573 of 1,221 documented components with minimal translation and 725 with gene calling. At four it recovered 1,014 and 1,099. Those are capped results. The archived GFP whole-DNA search was a localized diagnostic on a narrow panel. Neither result settles the full-pool, multi-family comparison. [Historical preprocessing validation](PREPROCESSING_VALIDATION.md) · [Archived translated search](PREPROCESSING_EXPLORATION.md)

## Build the cohort before scoring

1. Inventory existing archived sequences, corroborated annotations, reference files, frozen identities, and pinned model weights. Recover reproducible inputs from their original public sources if the local copies are missing. Time one representative uncapped job before scaling. Preserve source release or commit, hashes, assay or annotation evidence, and inclusion dispositions.
2. Prioritize documented benign protein families with multiple independent lineages and constructs containing multiple verified coding features. For the redesign question, look for already published parent–redesign pairs with experimental evidence and complete deposited construct records. A protein sequence or a design repository alone can support an isolated protein diagnostic. A whole-construct case requires the corresponding deposited DNA and verified locus. Record a frequently repeated protein, vector backbone, or design series as related evidence.
3. Freeze a development subset and a held-out subset split by related proteins and construct lineage, with family-level sensitivity where the sample permits. Keep the test answer key evaluator-side. Classify identity from the full documented protein against the fixed same-family references using the existing coverage-qualified rules. Report insufficient coverage separately. Model-training exposure is a separate variable, with exact and near-relative overlap reported only where the training release can be audited.
4. Count independent families, protein lineages, constructs, and documented occurrences in each tier. Aim for at least eight families and 30 independent protein lineages in the combined mid/low primary cohort. Treat this as a planning floor and report intervals and family sensitivity even if it is met. A single source, duplicated reporter series, or a few accepted case studies can inform feasibility. A population claim needs enough diversity and precision after accounting for related records.

## Run the comparison

| Arm | Input and search space | Purpose |
| --- | --- | --- |
| Whole-DNA BLASTx | Raw construct and the shared protein library. Retain eligible local hits across the construct. | Practical sequence-based baseline. |
| Minimal translation + ESM-2 | All qualifying candidates from the existing minimal six-frame translator. | ESM-2 baseline. |
| Gene calling + ESM-2 | All qualifying candidates from the existing gene caller. | Test the effect of preprocessing while holding the ESM-2 scorer fixed. |
| Documented protein + ESM-2 | Verified protein input with known locations. | Diagnostic for candidate generation, boundaries, and family ranking. |

Use the same reference sequences and family labels. Retain the pinned ESM-2 checkpoint and cosine retrieval for the primary analysis. Search complete candidate pools and keep per-window model-length handling and recorded cost. Candidate count and returned prediction count are outcomes in the primary experiment. Preserve the earlier capped settings as secondary comparisons. The BLASTx and ESM-2 scores have different scales. Choose each workflow's call rule only on development data and compare a held-out recall versus call-burden curve. Record the aligned spans as well as inferred locus spans so the evaluator can distinguish a convincing local family hit from recovery of a documented full locus.

## Measurements and interpretation

- Primary contrasts are gene-calling ESM-2 versus minimal ESM-2 and versus raw-DNA BLASTx. Measure one-to-one documented family/locus recovery and report the two paired differences at an agreed level of unsupported calls. Report pooled and equal-lineage estimates, lineage-cluster intervals, per-family sensitivity, achieved precision on exhaustively annotated records, and wall time and input residues. On incompletely annotated archives, report unevaluable calls separately from verified errors. If the selected operating point cannot be calibrated, report the full curve descriptively.
- Keep the existing 80% documented-locus criterion for continuity. Also report a predeclared parent-level endpoint for calls supported by linked windows and a local, same-family hit diagnostic. Preserve the distinction between the diagnostic and the primary endpoint. Review discordant cases for failure to generate a usable candidate, boundary or window coverage, wrong family ranking, reporting threshold, and ambiguous ground truth.
- If gene calling beats minimal ESM-2 but trails raw BLASTx, report the benefit within the ESM-2 pipeline and where the end-to-end comparison favors BLASTx. If advantages concentrate in a family or only appear at high unsupported-call burden, report that scope. If recall improves over both baselines with stable lineage and family effects at comparable burden, prioritize an independent archived replication. If the full-pool preprocessing gap shrinks, distinguish ranking benefits at small budgets from full-pool benefits in candidate quality.
- If clean documented proteins also fail, inspect reference representation and scoring assumptions before attributing errors to pretraining. For the redesign subset, compare the preprocessing gain in verified parent–redesign pairs and record how many genuinely independent pairs support the estimate. Low identity to the retrieval library supports a statement about that library.

## Execution order and stopping conditions

Start model and reference access checks, the independent-data audit, and throughput measurement immediately after the plan is posted. Define and verify the common evaluator before running the final comparison. Defer large downloads, exhaustive inference, and secondary scorer variants until the eligibility audit and throughput check show that those runs can answer the intended question. Record a resumable pilot before scaling. Preserve sequence-free artifacts and summary counts, with source sequences and model files in local research storage.

Before closing, check whether the analysis supports a clear answer for both comparisons on a sufficiently diverse mid/low cohort. If the dataset, weights, reference inputs, or compute make that infeasible, record the completed evidence, the exact dependency, and the commands or checks needed to resume. At less than 10% remaining task compute, checkpoint the analysis and write a handoff.
