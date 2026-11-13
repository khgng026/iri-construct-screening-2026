# Week 2 parser and representation QA

T10–T13 now have a tested, reusable preparation pipeline. It reads actual
permitted preparatory sources, verifies their exact identities and hashes,
normalizes source annotations into unchanged T03 v0.1 objects, and exports only
validated source-provided CDS translations. T14 has a recorded stratified
technical QA and schema binding. Scientific source admission, reviewed family
labels and experimental clearance remain separate prerequisites for a scored
benchmark; this is not a new E1 result.

## Implementation and corrections

- [Shared ingestion](../../src/iri_constructs/bench_inputs.py) detects FASTA and
  GenBank from content, rejects a declared-format mismatch, validates the exact
  record ID and DNA/protein alphabet, and checks original bytes before and after
  parsing. Malformed text, invalid UTF-8, duplicate IDs and unsupported annotations
  produce controlled errors without raw source content in logs.
- [Preparatory QA](../../src/iri_constructs/week2_pipeline.py) reuses the existing
  normalizer instead of introducing a second schema. It writes normalized
  feature tables, exact original qualifier sidecars, stable parent/feature/unit
  provenance, and a local annotated-protein export. No raw DNA enters a protein
  model, no unrelated proteins are concatenated, and no ORF is inferred.
- Multiple related-lineage proposals now survive curation and QA. Previously
  they were present only in the review ledger, leaving the QA lineage map empty.
  The regression verifies shared and legacy links without treating them as
  independent samples.
- Registry feature completeness now requires exact integer page 1 and total
  count. Boolean `true` cannot masquerade as page/count 1. Unverified coordinate
  conventions retain original annotations but withhold normalized locations.

## Actual-source QA results

The [main QA report](qa/main-report.json) records **18 prepared records, 271 source
features, 18 validated provided translations, and zero processing failures**.
The main pool comprises eight GenBank and ten iGEM sources. Every exported protein
passes the existing 1022-residue sequence preflight; this does not imply model
inference, accepted family labels, or scientific admission.

The [separate full-fusion QA](qa/fusion-report.json) records **one prepared record,
27 original features and one 825-aa provided translation**. It remains outside
the main pool. The existing pBR322 demo also processes all 50 original features
and three provided translations. Preparatory feature counts include source
bookkeeping entries; they are not the earlier E1 observable-feature denominator
of 48 and must not be compared as if the endpoint changed.

Two fresh main runs and two fresh supplemental runs have byte-identical
normalized objects, qualifier audit, feature/unit provenance, protein exports,
and source-free reports. Exact paths, input hashes, content hashes and artifact
digests are in [schema-freeze.json](schema-freeze.json). The unchanged
[T03 schema](../week1/T03_archetypes_and_schema.md) remains version 0.1; no schema
amendment was needed. This freeze is technical evidence for a subsequent review,
not a benchmark/source clearance.

## Stratified inspection checklist

| Check | Inspected evidence | Finding |
| --- | --- | --- |
| Source identities, versions and original bytes | All 18 main inputs and MK448012.1 supplement | IDs and checksums agree; each raw annotation response is separately bound |
| Backbone-heavy reverse CDS | pBR322 beta-lactamase | Original `[3292:4153](-)` and its strand/segments survive; provided translation validates |
| Multiple documented CDSs | AF298787.1 pUG35 | URA3 `[416:1220](+)`, EGFP `[2270:2987](-)`, AmpR `[4656:5517](-)` remain separate units |
| Single provided CDS | JX560321.2 pSEVA111 | One source CDS becomes one protein; noncoding annotations remain metadata |
| Tagged/fusion architecture | MK448012.1 pCHX50 | Source product `3xUbVR-tTA-Dendra2` is one complete 825-aa CDS; annotated subfeatures remain metadata |
| Misleading labels | MK448012.1 marker notes | `NeoR_CDS`/`AmpR_CDS` occur as `misc_feature` notes; no extra CDS proteins are invented |
| Missing translations | All ten iGEM entries | Source CDS roles retained, zero exported proteins; unknown boundaries remain ambiguous |
| Fusion-ready architecture | AF298787.1 | GFP fusion vector is not relabeled a completed target fusion |
| Parent/feature/unit linkage | Every normalized row | Source feature counts, audit linkage, ordered segments and source-provided translation strings checked |
| Shared sources and lineages | Main source ledger and QA report | Exact parent/protein duplicates plus five proposed lineage groups are visible |
| Circular joins, malformed inputs and length cap | Explicit fictional regression fixtures | Ordered origin-crossing segments preserved; invalid/unsupported inputs abstain; over-cap provided proteins remain in coverage but are not model eligible |

No selected real source supplied a joined CDS needed to establish real-data
origin-crossing performance. That behavior is regression-tested on a clearly
fictional circular fixture and is not presented as observed biological coverage.
Missing/ambiguous annotations remain limitations, not fabricated negative labels
or a perfect annotation-preservation score.

## Runnable checks and outputs

```bash
.venv/bin/python -m unittest discover -s tests -p test_week2_pipeline.py -v
.venv/bin/python -m unittest discover -s tests -p test_dataset_curation.py -v

# One-command clean-checkout source demonstration; no network or scientific scoring:
.venv/bin/python scripts/validate_week2.py \
  --inputs data/manifests/week2/demo.preparatory.json \
  --output results/week2/new-demo
```

Fourteen focused parser tests and seventeen curation tests pass. The parser tests
cover valid and malformed formats, sequence-only FASTA, exact identity/version
and hashes, duplicate aliases and parents, stable deterministic outputs, shared
lineage links, missing translations, circular joined CDSs, model caps, and
conservative Registry handling. Structured per-record logs contain status and
failure class rather than sequence data. Output directories are never reused.

Sequence-bearing outputs and full qualifier sidecars remain ignored. Committed
[feature provenance](qa/feature-provenance.json) and
[unit provenance](qa/unit-provenance.json) contain IDs, coordinates, counts,
fingerprints and eligibility states without source sequences or embeddings.
