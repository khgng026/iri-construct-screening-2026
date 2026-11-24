# Whole-DNA translated search pilot

## Question

How does a direct translated protein search compare with gene calling and
minimal translation on the same archived DNA, when scoring uses the same
protein reference library? This is an end-to-end comparator for the original
preprocessing result. It also measures candidate reach separately from
assignment on already selected candidates.

## Inputs and boundaries

Use the archived [Open Reporters Collection](https://github.com/Reclone-org/Open-DNA-Collections/tree/55d4b762194055c8027fd6c780068d8a78a03282/Open%20Reporters%20Collection)
and the [independent natural GFP panel](https://doi.org/10.7554/eLife.75842)
as frozen by the previous [Open Reporters audit](https://github.com/khgng026/iri-construct-screening-2026/pull/32).
Reuse its source hashes, verified donor matches and seven non-anthozoan GFP
references. None of the reference proteins exactly equals one of the assessed
reporters. Scorers receive plain DNA or protein FASTA and the reference library.
Use the independently corroborated 19 GFP-like protein regions as positives.
The five other reporter labels in the catalogue can support a descriptive
out-of-family check. Those five labels lack independent sequence
corroboration. Four uncorroborated GFP-like constructs remain excluded from
the positive and out-of-family summaries.

Record strand and coding-frame agreement. A correct region hit requires one
aligned region with at least 80% overlap with the documented coding region
after excluding its terminal stop. Also show results at 50% overlap to expose
partial matches. These are region-localization results for one family in a
shared cloning background. Unmatched hits away from the reporter locus have
unverified labels under this benchmark.

## Methods

- Whole archived DNA directly through BLASTx and DIAMOND blastx.
- Minimal six-frame candidate translation and Pyrodigal, followed by BLASTp
  or DIAMOND blastp with the same references and family label.
- For each candidate method, evaluate actual selected prefixes of 1, 2, 4
  and 8 candidates and its full pool. Separately record compatible selected
  candidate coverage and the fraction of those regions receiving a correct
  localized hit. The historical candidate generators and order are unchanged.
- Use an E-value ceiling of 1e-5 in the primary runs, with one CPU thread.
  Rank returned hits by bit score. Collapse redundant aligned regions in the
  same frame if they overlap by at least 80% of the shorter hit. Report
  localization at hit-report caps of 1, 2 and all passing distinct locations.
  These result caps apply across raw and preprocessed search. Candidate caps
  are a separate constraint in the explicit preprocessing arms.
- Use one reference database per engine and report end-to-end wall time
  separately from deterministic accuracy and work counters. Timing is a
  small-batch observation on this machine.

The panel already served as an independent source audit. This evaluation
is exploratory and uses existing archived benign annotations. It does not
estimate generalized construct-annotation performance or screening efficacy.

## Follow-up comparing search sensitivity

The default DIAMOND configuration recovered few targets in the first run.
Re-evaluate using its very-sensitive mode, keeping the references, candidate
generation, coordinate verifier and score cutoff fixed. Run both DIAMOND
configurations with separate output files. Report the difference as a
configuration sensitivity analysis on these same reporter records.

## Follow-up applying the historical length ranking

The retrospective ranking replay favored candidates near 250 residues at
the tight cap. Apply that same locked policy to each complete candidate pool
of the existing reporter archive, with BLASTp and very-sensitive DIAMOND
blastp. Compare at all specified caps. This checks transfer of the ranking to
different records and one narrow class of reporter lengths. The archive was
already inspected by this project.
