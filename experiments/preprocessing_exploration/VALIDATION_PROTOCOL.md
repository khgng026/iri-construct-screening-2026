# Preprocessing-only validation

Compare minimal translation, the current Pyrodigal baseline, and changes
to candidate handling. Use only benign archived records. The primary
sequence-only endpoint is an exact documented peptide at its original,
in-frame locus. For proteins longer than one model input, also count
whether the selected windows jointly contain all documented residues.
Run family identification only where the historical ESM-2 scores are
available. Keep its model, reference panel and five-call output rule fixed.

## Panels

- The new panel is a fixed set of 27 records from the Reclone Open Enzyme
  catalogue at the revision in `validation_panel.json`. Selection uses
  the catalogue, the target CDS name, and the availability of archived
  bytes. Target translations must pass the existing annotation validator.
  Report exclusions and distinguish displayed catalogue groups from
  independent lineages. The vectors share a cloning background. The
  collection's readme describes OpenMTA distribution for the physical
  parts. This study reads the public digital records locally and shares
  summary statistics and source metadata.
- Reuse the Open Yeast and LAMPRA panels as sensitivity checks. They were
  already inspected and do not provide fresh confirmation.
- Replay the existing historical cohort with the ESM-2 arm alone. It
  is retrospective and can use only candidates scored in a saved run.

## Candidate policies and budgets

Use the original candidate generators and keep all historical artifacts
unchanged. Compare the recorded gene-call order, the prior mixed ranking
that keeps the first gene call and favors candidates close to 250 residues
in later slots, the resolved-translation arm, and two fallback policies.
For the fallback policies, retain all ranked gene calls first, then fill
spare slots with distinct, fully resolved translation candidates. On the
archive records with unresolved bases, also test offering one distinct
fallback after the first two ranked gene calls. Distinctness uses the
existing same-frame locus overlap rule. The uncertainty flag comes from
the input DNA alone.

Include two ways to retain a long gene's peptide windows together in the
shortlist. One preserves the original window boundaries and processes
siblings consecutively. The other divides the same peptide into the
minimum number of similarly sized windows. Record the number of model
inputs and total residues for every cap. Align lineage-level claims to the
data source and avoid extrapolating from catalogue groups.

Compare caps 1, 2, 4, 8 and 16 for candidate reach, including a full-pool
diagnostic. Report realized input calls and residues. Repeat at the
per-record residue ceilings used by original gene calling at each cap.
The ESM-2 replay compares caps 1, 2 and 4 because extra rankings are
unavailable. Show the historical low-similarity and insufficient-coverage
tiers as descriptive strata. Keep them distinct from the archive's
peptide-reach endpoint.

## Exploratory check added after seeing the archive results

Keep the initial mixed-rule candidate and all its remaining recorded
windows at the head of the queue. Compare exact peptide availability and
cost. Parent-gene membership is available in the archive run, so this
policy stays in that run. The same yeast case motivated and checked the
rule.
