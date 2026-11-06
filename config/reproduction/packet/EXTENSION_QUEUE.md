# IRI additional research queue

Use remaining session capacity for the following substantive experiments after shared core comparisons are verified. Do not stop merely because the initial query grid finishes. Before starting each extension, estimate cost from observed throughput, freeze its design, and leave time for independent verification and the final report. All extensions inherit the resource, access, data, and spending boundaries in OPERATIONS.md.

Extensions are planned after observing the earlier 8M pilot. Any design selected in response to overnight outcomes must be labelled exploratory. Record both completed and attempted work. Exhausting this queue, hitting the session boundary, or having no feasible useful job remaining are valid reasons to finish; elapsed time alone is not a scientific objective.

## Automatic segmentation

This is the first extension priority after core model/method and averaging/context comparisons. Move beyond supplied boundaries by evaluating segmentation using only information available at inference.

Implement a clearly documented alignment-guided segmenter inspired by Nishit's described iterative best-hit approach, without claiming his unavailable implementation or statistical calibration has been reproduced. It may consult the fixed reference catalogue but must never consult constructed component intervals or target labels during prediction. If a defensible stopping rule cannot be established, report the candidate generation separately and do not describe it as a calibrated detector.

Include fixed overlapping windows as a simple segmentation baseline. Choose widths, strides, minimum spans, overlap handling, segment caps, and candidate selection rules on a separate development panel. Partition by unordered protein group so reversals/permutations of one group cannot appear on both sides. This tests new combinations of the same proteins, not unseen-protein generalization. Freeze parameters before final evaluation; report how often caps bind. Keep known boundaries as the oracle control.

For each segmenter, record predicted spans, segment count, total queried residues, retrieved component coverage, runtime, and extraneous results. Compare end-to-end component recovery as well as interval overlap with known constructed components. A biological domain can differ from an engineered component, so exact boundary agreement is not the sole success criterion.

Separate label-free candidate generation from the evaluation's use of known target labels. Report two output budgets: the same total returned-hit budget per original query, and an explicitly larger per-segment budget with total review volume. Freeze label-free score aggregation or rank fusion on development data; do not select the best segment using the expected answer and call that deployable performance. Best rank over segments may be a labelled optimistic diagnostic only. With exact references present, alignment-guided segmentation is an easy control, not a demonstration of remote-boundary discovery.

## Catalogue expansion and reference removal

Prepare additional reviewed reference catalogues up to the protocol's 250/500/1,000-entry targets, then run the same frozen query comparisons. Data acquisition/review can proceed while another admitted model computes if resource limits allow. Keep source licensing/provenance checks and catalogue-family composition explicit.

Run exact-reference removal diagnostics using a common mask across all components/methods, as defined in PROTOCOL.md. If independent family/domain labels and surviving matches are defensible, run the labelled harder retrieval panel. Otherwise record why it is unscorable and continue another extension rather than manufacture labels.

## Representation sensitivity

Compare prespecified intermediate-layer and final-layer representations for a small number of admitted models, using the common suite. Set layer choices before outcomes; relative layer depths can differ across architectures. Compare mean pooling with the already defined length-weighted mixture and contextual component arms. Do not tune on the final benchmark or train a new model.

Run one layer at a time or extract only required states where supported. Avoid retaining all-layer tensors simply to make a broad sweep. A layer chosen retrospectively is exploratory and must be evaluated on an untouched partition before a claim of improvement.

## Failure analysis and stronger validation

Construct a reproducible case table for successes, persistent failures, and disagreements between alignment and embeddings. Include component fraction, length, order, related catalogue matches, clean-query behaviour, segmentation boundaries, and competing scores. Select examples by explicit rules covering best, worst, and discordant outcomes rather than attractive screenshots.

Review the closest published domain-segmentation and retrieval approaches, including DCTdomain, and identify which comparisons are actually missing. Run an existing implementation only if official installation, licence, input requirements, and resources permit. A documented incompatibility is useful evidence; it does not authorize spending hours forcing an unsuitable package to run.

If the core answers are stable and approved real annotated benign construct records are available, prepare a small real-input feasibility panel. Keep artificial and real results separate. If source/label review is incomplete, prepare the record inventory and parser checks without scoring unreviewed data. Do not claim source annotations are independently verified biological truth.

## Allocation and stopping

Before taking another expensive job, compare estimated runtime plus audit cost with remaining session time. Prioritize automatic segmentation, then credible catalogue/reference tests, then representation sensitivity. Keep an independent reviewer available for completed outputs. A blocked extension should lead to the next feasible queue item.

The lead should continue useful feasible work until the agreed reporting reserve, then finish validation and a morning report. Do not wait idly because the core is done, but do not add redundant runs, widen scope without justification, or conceal resource failures to fill the night.
