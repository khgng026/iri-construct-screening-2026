# Ranking candidates by sequence-free geometry

## Question and scope

Can geometry alone prioritize useful candidates in an existing topology-blind
translation pool for annotated benign archived constructs? This study changes the
ordering of each existing pool. Each policy uses only the pool's original order,
translated length, strand, and coordinates. The original three generators supply
the inputs separately: minimal six-frame, start-ranked six-frame, and Pyrodigal.
The evaluator checks each ordering against the previously documented component
regions. Component coordinates, families, and groups do not enter selection.

This is a retrospective exploration on the publicized cohort v2 traces.
Development is reported separately as a pipeline check. These settings were
written before inspecting the exploratory results. The historical holdout is
already known to this project, so no claim of a new independent test follows.

## Policies

- `recorded` keeps the order in the original trace.
- `defer_ceiling` postpones candidates at the 1,022-residue model length ceiling
  within each generator's recorded order.
- `length_center_250`, `length_center_450`, and `length_center_650` prefer
  lengths closest in relative terms to 250, 450, and 650 residues, respectively.
  Quality is the smaller of length divided by center and center divided by length.
- `novel_span_any` and `novel_span_same` maximize the number of newly covered
  residues, counting existing coverage from either strand or only the same strand.
  Span coverage is calculated from the union of selected intervals.
- `soft_diversity_any` and `soft_diversity_same` favor recorded rank with a
  penalty for overlap with prior selections. Quality is reciprocal square root of
  original one-based rank, multiplied by one minus half the overlap fraction.
- `balanced_diversity` applies the same modest all-strand overlap penalty to
  the 450-residue length quality.

All ties are broken by the recorded order. The policies keep the same candidate
pool and make no new predictions or sequence changes. Candidate caps are 1, 2, 3,
4, 6, 8, and 12. Separate residue-budget curves use 512, 1,024, 2,048, and 4,096
residues per construct, taking an eligible candidate in preference order and
skipping ones exceeding the remaining allowance. Candidate count and residue
count are cost proxies. The study logs ranking time separately.

## Endpoints

The primary endpoint is documented-region reach, with the historical verifier's
strand and overlap rule. For single-span regions, report coordinate-frame
agreement as a sensitivity check. It compares the first coding coordinate
modulo three on the same strand. It neither verifies the translated amino-acid
sequence nor reconstructs a phase for split features. Report pooled component
reach, equal-lineage-weighted reach, offered candidates and residues, and a
stratified count by annotated layout and similarity tier. These are descriptive
analyses of correlated cases.

For the separate Q2 replay, restrict each generator to the candidates that
already have BLOSUM62 and ESM-2 scores in its historical four-candidate trace.
Reorder within that pool and replay the original five-assignment policy at caps
1, 2, 3, and 4, including the recorded ordering as a control. Use the existing
one-to-one verifier to count correct family assignments at the annotated
position. Unknown biological components remain outside this recall endpoint.
This replay can measure whether a ranking policy helps on the restricted
scored pool. Full-pool reach gains require fresh fixed-scorer inference before
they can be interpreted as identification gains.

## Follow-up after inspecting the length ranking

The 250-residue policy gained most on high- and intermediate-similarity
components, while the insufficient-coverage tier lost recoveries. Compare
one mixed policy that keeps the first recorded candidate and moves
candidates closer to 250 residues up in the remaining slots. Replay both
policies and the original order on the same historical inputs. This is a
post-hoc trade-off check.
