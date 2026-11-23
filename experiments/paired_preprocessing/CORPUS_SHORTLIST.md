# Corpus audit and shortlist for the next paired batch

Stream 1 (evaluation and verification). This pass audits the existing component-study
corpus and produces a metadata-only shortlist of independently reviewed benign records.
It contains no preprocessing method, comparator choice, or model scoring; those are
developed separately. All biological payloads and the shortlist itself stay under the
ignored run directory `results/paired-corpus-batch-v1/`.

## Audit integrity

- **Source availability.** All 50 frozen batch files (1,000 records) reacquired from NCBI
  are byte-identical to the recorded hashes, with no missing records. All 40 UniProt
  family anchors and 134 panel references match their frozen sequence hashes.
- **Toolkit agreement.** 1,000/1,000 records processed; 2,338/2,338 exported translations
  match the frozen unit hashes and 2,338/2,338 recomputed family labels match the frozen
  labels. 1,311/1,311 anchor regions map back to DNA and re-translate to the same protein
  segment.
- **Per-occurrence location fix.** The audit previously cached protein-level family
  evidence and then wrote construct-specific DNA spans into that shared dictionary, so a
  protein repeated in another construct could overwrite the location saved for an earlier
  one while the tallied counts still passed. Location evidence is now assembled per CDS
  occurrence in `occurrence_anchor_evidence()` and never written back into the cache.
  `tests/test_corpus_batch.py` is a regression: two records that share one protein at
  different coordinates must keep distinct DNA spans, and the shared cache must not be
  mutated. It fails on the previous behaviour and passes now.
- **Metadata checker.** `check_corpus_audit.py` runs clean on the corrected audit (no
  findings, exit 0). It was extended to summarise programmed-frameshift CDSs
  (`join(a..b,b..c)`, two records) instead of aborting the run on their overlapping source
  segments. It independently reproduces the shared-protein diagnostic: **193 groups across
  1,000 records, largest group 773**.

## Shortlist (metadata only)

Eligible component: an exported CDS whose translation is validated against the DNA, whose
family is supported by a reviewed UniProt anchor (identity ≥ 0.90, anchor coverage ≥ 0.80)
that also has a reference panel, that is not a near-copy of a panel reference, and whose
deposited name does not name a different family. Deposited names are recorded separately
from the anchor evidence. The five reporter-pilot records are excluded. Missing or
unsupported annotations are outside the denominator.

### Records reviewed and eligible components by construct type

| Construct type | Records with an eligible component | Related groups (shared protein) |
| --- | ---: | ---: |
| Tagged or fused | 56 | — |
| Multi-CDS | 355 | — |
| Backbone-heavy | 149 | — |
| Single-CDS | 102 | — |
| Single label, other layout | 130 | — |
| **All independently reviewed records** | **792** | **63 (largest 706)** |

Archetype is assigned from the evaluator-side evidence: tagged/fused if any eligible
component is a fusion-class CDS (anchor covers < 80% of the CDS); else multi-CDS if two or
more CDSs carry an eligible label; else backbone-heavy if the single labelled CDS is < 15%
of the construct and another CDS carries no eligible label; else single-CDS if the record
has one CDS. "Single label, other layout" is one labelled CDS with extra unlabelled CDSs
at ≥ 15% labelled fraction; these are eligible but fit no clean archetype quota.

### Fusion subregions: documented boundaries vs anchor-only

A fusion component always counts in the ordinary whole-construct comparison. It is marked
segmentation-eligible only when the record itself documents the subregion: a deposited
sub-feature inside the fusion CDS that names the family, is codon-aligned, spans ≥ 30
residues, and overlaps the anchor region. That boundary comes from the deposit,
independent of the anchor alignment used to label the whole protein.

| Fusion subregion category | Components | Records | Related groups (largest) |
| --- | ---: | ---: | ---: |
| Documented boundary (segmentation-eligible) | 14 | 13 | 7 (4) |
| Anchor-only match (whole-fusion comparison) | 53 | 45 | 23 (20) |

The 14 documented boundaries are all deposited `misc_feature` subregions whose residue
span matches the anchor region to within the fusion's start methionine (e.g. `AB684655.1`
deposited gfp region 148–388 vs anchor 149–387). `MW879734.1` deposits both of its
components (`gst_mu` and `solute_binding_protein_1`). Both tomato histone–fluorescent
fusions flagged in the handoff (`AB684655.1`, `AB684657.1`) are in this set. A stricter
exact-boundary test (both-way containment, in the checker) counts 6 such regions; the
looser overlap test here counts 14. Verification of these boundaries does not require a
segmentation method to exist.

### Dependence, not record count

Shared common markers connect most vectors: the whole corpus forms 193 groups with a
largest group of 773 of 1,000 records. The 792 independently reviewed records collapse to
**63 related groups** (largest 706), and the 13 documented-boundary fusion records to
**7 groups**. Record counts therefore overstate the number of independent cases, and
setting aside a fixed number of fusions does not by itself establish an independent
holdout. Final evaluation groups should be chosen from the group structure, not the record
count.

## Reproduce

Stages are separate; only `acquire` uses the network. Outputs go to the ignored run
directory. A model environment is not needed for any stage here.

```sh
python -m experiments.paired_preprocessing.corpus_batch acquire     # network: NCBI, UniProt
python -m experiments.paired_preprocessing.corpus_batch audit        # toolkit + local alignment
python -m experiments.paired_preprocessing.corpus_batch shortlist    # metadata only
python -m experiments.paired_preprocessing.check_corpus_audit \
  --acquisition results/paired-corpus-batch-v1/acquisition.private.json \
  --audit results/paired-corpus-batch-v1/audit.private.json
python -m pytest tests/test_corpus_batch.py tests/test_corpus_audit_check.py
```

Artifacts (ignored): `acquisition.private.json`, `audit.private.json`,
`shortlist.private.json`, and the source snapshots under `sources/`. The opaque-ID
method inputs and the evaluator-only answer key are a later preparation step and are not
produced in this pass.

## Scope boundary

This pass covers audit and shortlist verification only. Preprocessing methods, candidate
selection and ordering, comparator choices, and model scoring are handled separately. The
anchor alignment is used here to label and locate components for the answer key; it is not
offered to any method under test.
