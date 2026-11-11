# Week 2 annotation semantics verification

FASTA contains an identifier and sequence, not a GenBank feature table. The
fictional FASTA regression deliberately uses a misleading `.gb` filename:
content detection still recognizes FASTA, reports zero observable annotations,
and exports no invented CDS. This is a software-format check, not a biological
example or an annotation-preservation score.

The inspected pBR322 GenBank record has 50 original features and three explicit
CDSs. Source GenBank intervals are displayed as 1-based inclusive; Biopython and
T03 preserve them as 0-based half-open intervals. Its beta-lactamase CDS is
represented as `[3292:4153]` on strand `-1`, corresponding to the source's
complemented interval. Reverse strand is a valid annotation, not a parser error.
The original provided translation agrees with extraction using the documented
frame and translation table.

AF298787.1 contains separate URA3, GFP and AmpR CDSs, including reverse-strand
features. Their original products and translations remain traceable through the
qualifier audit and checksums. The construct is fusion-ready; that name does not
establish a completed fusion with an unspecified target.

MK448012.1 instead provides a complete `3xUbVR-tTA-Dendra2` CDS and source-provided
825-aa translation. The annotated subcomponent locations remain source evidence;
no new isolated protein was reconstructed. Notes named `NeoR_CDS` and `AmpR_CDS`
are still `misc_feature` annotations and cannot become coding units merely
because the note contains `CDS`.

The iGEM responses supply original feature UUIDs, roles and raw locations, but
neither a verified coordinate-indexing convention nor source translation
qualifiers. Their raw annotations are retained, normalized locations remain
unknown/ambiguous, and no protein is synthesized to remove that uncertainty.

The [QA checklist](PARSER_QA.md) records the inspected sources and the distinction
between real-source findings and fictional malformed/joined regression cases.
Original full qualifier sidecars and sequences are private local artifacts;
shareable provenance records their fingerprints and locations.
