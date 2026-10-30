# Data: original fictional examples and local T05 pilot

T05 adds `manifests/t05_benign_catalogue.json`: five reviewed public protein records with versions, licence evidence, inclusion notes and checksums. This metadata manifest contains no sequences. The user explicitly authorized this narrow protein-only feasibility test; it does not approve a T06/T08 DNA benchmark. See `docs/week1/T05_design.md` for the scope and review authority.

Original API responses live in ignored `data/raw/t05/`; canonical FASTA files, artificial query FASTAs and schema v0.2 mappings live in ignored `data/processed/t05/`. Both directories remain outside Git. Source data: UniProt Consortium, UniProtKB/Swiss-Prot, CC BY 4.0; four yeast records are served by EMBL-EBI Proteins API. Canonical sequences are unchanged; artificial concatenations are identified in the query design and mappings. This branch does not redistribute sequences or the confidential presentation.

## Original T04 handoff (historical)

The `manifests/` and `examples/` folders contain only invented metadata for schema and testing demonstrations. No real accessions, sequences, evidence of biological function, or model outputs are stored here.

T02 data-governance policy v1.0 was approved on 2026-10-28 and documents source selection, benign eligibility, licensing, and provenance checks. This policy approval does not approve any individual record. Add real public benign records only after the record-specific T02 checklist, permission, and safety reviews have been completed in T06/T08. Downloaded `.gb`, `.fasta` or related sequence data belong outside Git in ignored local directories. No user/customer/proprietary data are allowed.
