# Historical scorer combination replay

Use the stored cohort v2 Pyrodigal traces at the fixed two-candidate cap.
They contain paired per-family ranks for normalized BLOSUM62 and mean-pooled
ESM-2 8M on identical candidates. The historical answer key is already
visible. This is a retrospective test of model information available in
those ranks, with a much weaker model and baseline than the planned Q3
comparison. The traces do not have comparable confidence scores or verified
out-of-family targets.

- Fix an equal-weight rank sum over the shared family ranking for each
  query. For an exact tie, keep the order of the conventional rank. Use one
  returned label per query, preserving the candidate and prediction budgets.
  Require identical selected inputs and a complete, shared family ranking
  from both scorers. Use the existing locus-aware one-to-one verifier.
- Report paired gains and group-bootstrap intervals for the full holdout,
  its original low-similarity tier, and the separately reported
  insufficient-coverage tier. Report the number of known-locus assignments
  that are incompatible with a single overlapping corroborated label.
  Classify all queries without one unambiguous documented component as
  unknown for this measure. The label error rate is a diagnostic on this
  subset, not an open-set false-assignment rate.

After the Open Reporters audit reported zero exact full-input matches,
check a prespecified explanatory distinction on those same selected
candidates. Count cases where a covering candidate contains the complete
documented reporter peptide as a contiguous substring, versus cases with
position overlap alone. Mark this as a post-hoc diagnostic on the published
reporter prototypes. It changes neither eligibility nor the frozen
reachability comparison.
