# T07 scoring audit and T17 preparation

The E2 parent scorer, E3 reference-tier scorer and E4 cost/quality aggregator
already implemented most of the predeclared design. This audit reused those
modules rather than creating a second benchmark. All numerical checks described
here use fictional software fixtures, not biological model measurements or
scientific-clearance decisions.

## Corrections

### Correlated CDS queries were counted as independent observations

The shared pilot/E3 paired-comparison code previously used only the protein's
reviewed `independence_group`. Eight CDS queries from one construct could have
eight different protein groups and incorrectly reach the eight-cluster
confidence-interval threshold. E3's E4 measurement records similarly counted the
protein groups rather than independent construct observations.

The scorer now freezes conservative connected query clusters using shared source
parent, reviewed independence group and exact query-protein fingerprints.
Indirect links are joined transitively. Original reviewed groups remain visible
for source/reference exclusion and audit; they are not replaced by a guessed
biological family. Paired differences and E3/E4 independent-group counts use the
same resulting clusters. All original reference guards and the saturated-budget
exclusions remain in effect.

The new regression demonstrates that eight distinct protein groups from one
parent produce **one** independent cluster and a null interval. A separate
regression checks indirect parent/group/exact-copy linkage and deterministic
cluster identities. Tier reports distinguish raw reviewed-group counts from
connected-cluster counts. An archived panel without the additional cluster
metadata is interpreted from its original frozen fingerprints; historical
artifact files have not been edited.

### Numerical pilot results could be mislabeled as scientific results

The legacy B/C/D executor previously set `scientific_scores_produced` whenever
any score row existed, including unclassified inputs and fictional numerical
adapters. E2/E3 already had explicit evidence classifications.

The legacy executor now records `scores_produced` separately and sets the
scientific flag only for explicit `evidence_kind=admitted_research`. Missing
classification remains `unclassified`; fictional results remain fictional.
Unsupported classifications block before an adapter is loaded. Existing
execution statuses still describe whether computation completed and do not
constitute scientific admission. Structural flags cannot establish that an
external review actually occurred; the T02/T07 source and scientific-clearance
requirements still apply.

### E4 accepted internally contradictory denominators

E4 now rejects coverage values that contradict recorded scored/requested
counts, and independent-group or saturated-observation counts exceeding scored
observations. Missing coverage remains null. The legacy-report adapter no
longer treats a raw per-protein group count as proof of independent constructs;
without recorded conservative clusters it preserves the count as unknown.

## Verified invariants

- E2 chooses the longest eligible documented CDS with a stable ID tie break.
- Every E2 condition uses the same frozen eligible references. Scores merge by
  maximum per reference, then return five total unique hits per parent.
- E2 excludes the source parent, all retained documented component fingerprints
  and reviewed component groups, including components outside the baseline.
- Repeated components in one family do not enlarge the distinct-family
  denominator. Full documented and independently supported denominators remain
  separate; unknown labels remain explicit.
- Missing/nonfinite scores produce abstention rather than a smaller reference
  catalogue. Failed additional CDSs do not erase a successful baseline result.
- Reference tiers preserve mandatory parent/copy/group guards. Unsupported
  strict policies remain blocked; valid earlier tiers are retained.
- Methods share identical frozen panel identities and five-hit budgets.
  Returning the whole eligible catalogue remains descriptive only.
- E4 keeps source, protocol, hardware/software, cache policy and scored-support
  identities separate. Incompatible rows cannot silently enter a model ranking.
  Shared physical work is labeled rather than summed across conditions/tiers.

## Validation

Focused regression commands:

```bash
.venv/bin/python -m unittest discover -s tests -p 'test_pilot_scoring.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_pilot_experiments.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_extension_execution.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_reference_tiers.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_multi_cds.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_cost_quality.py' -v
.venv/bin/python -m unittest discover -s tests -p 'test_pilot_alignment.py' -v
.venv/bin/python scripts/verify_pilot_release.py
```

The original pBR322 source, translations, historical E1 measurements and saved
blocked E2–E4 artifact digests still verify. No frozen historical protocol or
configuration was altered. Existing CI and the integrated Week 2 validation
record provide the final repository-wide result after all workstreams merge.

## Status and limits

T17 has tested scoring infrastructure and corrected uncertainty/evidence
accounting. It is not a completed real-data benchmark. T18–T21 still depend on a
cleared independent panel and actual paired experiments. This engineering
regression is not a post-experiment T20 preprocessing improvement and is not
claimed as T21 scientific before/after evidence.

No new biological inference, scientific alignment scores, E1–E4 experiments,
thresholds, family assignments or delegated clearances were produced by this
audit. Connected provenance clusters are conservative dependence controls, not
proof that the remaining clusters are biologically or pretraining-independent.
