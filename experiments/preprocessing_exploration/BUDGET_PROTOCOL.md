# Exploratory allocation of a shared candidate budget

This retrospective experiment uses the committed sequence-free traces of
public, benign cohort v2. It asks whether some constructs can usefully receive
more candidate evaluations while others receive fewer at the same total cost.
The within-construct candidate order, selected input representation, reference
panel, and scorer are held fixed. The answer keys were already known to the
project and there is no new holdout.

## Comparison fixed before running

- Study minimal six-frame and Pyrodigal candidate orders separately. Read
  candidate loci, strand, and residue counts from the committed metadata. The
  allocator receives these features only. Compare caps of 1, 2, 3, and 4.
  Cap 3 is reconstructed from saved cap-4 traces.
- For each uniform cap, use the actual total selected candidates and residues
  as simultaneous ceilings for all allocations on the same candidate order.
  Reserve a first candidate for every nonempty construct. A policy can spend
  less than either ceiling. Report both realized costs. All alternative
  policies can assign at most four candidates to one construct because family
  rankings were saved only through cap 4. Cap 4 is consequently a saturation
  control. Caps 2 and 3 are the informative comparisons within this constraint.
- Compare five deterministic priorities for the next candidate in each
  construct: cheapest residues first, longest first, largest full-pool demand
  per candidate already allocated, demand divided by the next candidate cost,
  and least redundant locus in the within-construct prefix. Compare two
  conservative priorities that partially spread effort: rank-balanced
  cheapest-first and rank-balanced demand. Use a fixed hashed pseudonym as
  the final tie-break. The allocator never changes within-construct order.
- First score documented-region coverage by the original 80% same-strand
  overlap. Then replay the saved ESM-2 8M mean-pooled and normalized BLOSUM62
  family rankings using the same five returned calls per construct and the
  original one-to-one assignment verifier. Check that replays reproduce all
  saved uniform cap 1, 2, and 4 calls. Report pooled and lineage-weighted
  outcomes and a paired lineage-bootstrap interval for differences. Coverage
  describes locus reach, and does not confirm reading frame or a correct label.
- Use an evaluator-only oracle to choose one of the first four prefixes per
  construct while respecting both ceilings and the first-candidate floor.
  Separately maximize coverage and each saved scorer's family recovery under
  the original output rule. Solve this diagnostic as a binary linear program.
  If optimality is unresolved, report the demonstrated feasible objective and
  solver bound. The oracle measures room to rearrange the existing prefixes,
  without proposing a policy that can know a component label at inference.
- Report uniform full-pool coverage at caps 8 and 16 and without a cap.
  These tiers use recorded coordinates only. They cannot support a claim
  about family calls beyond the stored four-candidate traces.

All shared outputs contain aggregate metrics and input hashes. No DNA or
protein sequences, per-record predictions, or per-record outcomes are released.
