# Example inputs

## `synthetic/` — fictional software fixtures

Every sequence here was produced by a seeded random generator
(`python examples/make_synthetic_examples.py`, seed 20261117, Biopython 1.86).
They have no biological source, so there is no third-party right to review and
they are safe to commit. They test *implementation behavior*; they are not
biological observations and say nothing about real annotation quality.
`EXPECTED.json` states, by construction, the reason code each feature must get;
`tests/test_preprocess.py` compares the toolkit against it.

| File | What it exercises |
| --- | --- |
| `single_cds.gb` | Promoter, one forward CDS, terminator: the plain path |
| `edge_cases.gb` | Circular record with a reverse CDS, an origin-spanning join, a reverse two-segment join, translation table 4, missing / mismatching / non-canonical translations, a CDS without stop codon, `/pseudo`, a 5'-partial CDS, identical duplicate features, a between-bases site, `order(...)`, a remote location and a feature beyond the sequence end |
| `long_protein.gb` | A valid 1,100-residue protein: exported, flagged only against a model limit |
| `multi_record_one_bad.gb` | Three records; the middle one has an unparseable location and fails alone |
| `no_accession_export.gb` | A tool-style export with `LOCUS Exported` and no accession |
| `sequence_only.fasta` | Nucleotide FASTA: a construct with zero features and zero invented proteins |
| `not_a_sequence.gb` | Not a sequence file: an input-level failure |

Processing the whole directory gives 7 processed records, 2 failures, 11
exported proteins and 6 withheld CDSs, exit status 2.

## Real public records

- `data/approved_public/J01749.1.gb` (pBR322, NCBI GenBank) is the one real
  record in Git; its source, hash and sharing basis are in
  [`data/approved_public/README.md`](../data/approved_public/README.md). The
  test-suite uses it as an integration check (50 features, 3 of 3 CDSs exported).
- A seeded sample of 200 NCBI vector records was processed for the
  [component-coverage study](../docs/results/component-study-v1/README.md):
  200 of 200 records processed, 488 of 528 CDSs exported, 36 withheld for
  partial boundaries and 4 for exceptional annotations. Those records are held
  in an ignored private cache with a recorded acquisition manifest (accessions,
  hashes, retrieval date, NCBI terms note); they are not redistributed here.
  Source-level review: NCBI GenBank cloning/expression vector records, used
  privately for research with attribution; no record-specific conflicting
  notice was looked for beyond excluding the patent division. This is the
  agent's own review, not an independent one.
