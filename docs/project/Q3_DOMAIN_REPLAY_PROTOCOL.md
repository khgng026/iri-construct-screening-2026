# Matched-domain robustness replay

This retrospective check uses the recorded decisions from PR #28 at
`4171d33bc0f7404b652f17c45c599d67e5c6a23f`. The original labels, model sizes,
reference libraries, and headline outcomes are already public. This is a
sensitivity analysis of the same experiment.

Use the original below-25%-global-identity domain cohort and pair each domain
with the wider-reference decision by its source ID, fold, fine family and
true structural superfamily. Re-check that all retained domains qualify under
both reference policies. Repeat the summaries on the subset without a
coverage-qualified same-label local hit under either policy. Report each
method's hit counts, the associated counts by fold, and exact counts of
discordant outcomes. Compare local ESM-2 150M against raw BLOSUM62,
reference-normalized BLOSUM62, and each HMM individually. Report the model's
gain over the best conventional scorer fixed across all domains, and an
optimistic conventional union of correct calls. The latter uses identified
labels to select an answer and does not estimate an available policy.

Record the paired domain-level agreement of the advantage across reference
policies. Compare the number of unique protein domains the model gets right
where every conventional scorer misses with the reverse count. Summarize
these counts within folds and count repeatable rescues. Report the
single-scorer and optimistic comparisons for both reference policies, and
include leave-one-fold-out ranges to show the effect of the small number
of folds. Treat all numbers as descriptive for these domain data.

Because the input is a known protein domain and one of a few closed-set
labels is always eligible, this replay does not estimate end-to-end gain,
open-set specificity, or biological activity. Agreement across the two
libraries uses the same queries and is a robustness check, not replication.
