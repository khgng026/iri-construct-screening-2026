# Component-coverage study v1 — extension 1 (replication sample and fusion units)

Written 2026-11-17 after the v1 results and before any extension record was
fetched or scored. Agent self-review only. v1 files and results are not changed.

## Why

v1 left two things thin: the parent-level pattern was seen on one sample, and
the observation that a mean-pooled embedding of a whole fusion protein loses its
documented family rested on ten sequences. This extension adds a second,
disjoint sample to check both. It does not add methods, models or labels.

## Fixed from v1 (unchanged)

Model, pooling, cosine, alignment scorer, `families.v1.json`, anchors, label
thresholds, the **same 134-reference panel**, K = 5, the four conditions, the
component-attributed primary endpoint, eligibility (exported, <= 1022 residues),
deduplication and the eight-cluster rule for population statements.

## New input

A second seeded draw of 800 records from the same NCBI query, seed 20261118,
excluding every record of the v1 sample. Processed by the same toolkit command.
A unit labeled in family F that is a near copy (>= 0.90 identity, >= 0.80
coverage) of an F reference of the frozen panel is excluded from scoring and
counted, because the panel's lineage exclusion was computed for v1 units only.

## Registered analyses

1. **Toolkit behavior** on the 800 records: processed/failed records, exported
   and withheld CDSs by reason, warnings. Reported whatever it shows.
2. **Parent-level replication**: the v1 table (four conditions, both methods,
   attributed recall, distinct unit sets, per-unit ceiling) on the new parents.
   Parents whose eligible-unit sequence set already occurred in v1 are counted
   and removed, so the replication is on unit sets not seen before.
3. **Per-unit retrieval by label class** for distinct supported sequences that
   did **not** occur in v1 (`full`, `embedded`, `multi_family`), both methods:
   hit@1, hit@5, and for `multi_family` the family slots recovered at 5. The
   fusion classes are also totalled over v1 plus extension, with each part shown.

Statement fixed in advance: if the extension's fusion-class hit@5 for ESM-2 is
not clearly below that of `full` units, the v1 fusion observation is withdrawn
as a small-sample artifact.

## One exploratory diagnostic (labelled as such in the results)

For each fusion-class unit, embed only the segment that aligns to its anchor and
score it against the same panel. The cut uses the label's own alignment, i.e.
the answer key: this is an **oracle upper bound on what within-protein splitting
could recover for mean pooling**, not a decomposition method, not a toolkit
feature and not evidence that such boundaries can be found without labels. It
is run only to tell whether the fusion loss is caused by pooling over the extra
sequence. Units of class `full` trimmed to their anchor-aligned segment are the
control.

## Bounds

CPU; each scoring stage under 30 minutes and 6 GiB; private inputs stay in
ignored storage.
