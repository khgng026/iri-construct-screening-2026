# T03 → T05: minimal v0.2 extension

The original construct → feature → unit IDs and foreign keys are retained. The
existing checker accepts both v0.1 and v0.2; all three original fixtures are unchanged.

| Addition | Meaning |
| --- | --- |
| `sequence_ref` on a parent or ready unit | `{path, format, record_id, sha256, length, alphabet}`. Repository-relative path inside ignored `data/raw/` or `data/processed/`; one-record FASTA; SHA-256 of the exact uppercase sequence bytes without header/newline. No inline sequence. |
| Parent `molecule_type` | Protein-only parents in this pilot; original v0.1 parents still use DNA coordinates. |
| `sequence_length_aa` | Protein parent length; `sequence_length_nt` is null when no DNA record exists. |
| `coordinate_system` on parent/features/units | Explicit `protein_0_based_half_open` for this pilot. Protein strands are null. Never reinterpret these numbers as nucleotide positions. |
| Parent `primary_archetype: unknown` | Avoid assigning a DNA-construct archetype to an artificial protein string. |
| Parent `artificial`; feature `source_record_id` | Distinguish constructed inputs and retain linkage to the frozen source manifest. |

Mappings are generated in `data/processed/t05/normalized/`. Each component has a
stable feature ID and an interval in the artificial parent. The whole unit references
all features; the oracle unit references only GFP and its exact interval. Each unit
points to its own local FASTA and checksum. The source manifest separately stores
the downloaded JSON file checksum, sequence version and annotation version.

For `gfp_pgk`, protein parent length is 654, with GFP at `[0,238)` and PGK1 at
`[238,654)`. For `pgk_gfp`, GFP is `[416,654)`. These boundaries are arithmetic
from the known recipe, not a discovery by the model. The same reviewed GFP is
reused in all cases, so these are not independent biological observations.

`ready` here means a reviewed canonical protein string that passes this model's
input checks. It does not mean verified biological function, wet-lab feasibility,
or safety classification. Future DNA features still need an actual annotation/
translation pipeline before receiving this status. Unsupported noncoding features
remain metadata-only; none are sent to ESM-2 in T05.
