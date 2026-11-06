# T03 | Construct archetypes and normalized output schema

**Week:** 1  |  **Budget:** 2 hours  |  **Status:** Approved schema v0.1 baseline (researcher sign-off 2026-10-29)  |  **Prerequisite:** T01  |  **Calendar source:** `Weekly To-Dos!A7:K7`.

## Archetype definitions (selection labels, not inferred hazards)

| Archetype code | Inclusion description | Key annotation limitation |
| --- | --- | --- |
| `single_cds` | A source record with one principal annotated coding sequence. | Backbone and regulatory annotations can still be present. |
| `multi_cds` | One record with at least two distinct annotated CDS features, e.g. a benign multi-gene cassette. | Preserve each feature and its orientation, not one merged protein. |
| `tagged_fusion` | A construct with an annotated fusion, tag or linker associated with a coding product. | Subunit boundaries may not be known; do not infer them. |
| `backbone_heavy` | A plasmid/construct whose annotated vector/backbone and noncoding elements are substantial for its interpretation. | No automatic protein interpretation of noncoding regions. |

**Multiple labels can apply.** Store a `primary_archetype` for stratified counts plus `secondary_archetypes` for overlays. This avoids misclassifying a tagged fusion occurring in a multi-CDS construct. Archetypes may be assigned only from supported record metadata and reviewer notes; `unknown` and `ambiguous` are first-class states.

## Normalized representation: v0.1 (three linked objects)

The approved v0.1 JSON examples in `data/examples/` are **hypothetical schematic benign metadata only**: no real accessions, sequences, or verified biological findings.

### `construct` (one per source record)

| Field | Type and meaning |
| --- | --- |
| `construct_id` | Stable unique local ID used by features and units. |
| `source` / `source_url` / `accession` / `accession_version` | Original public source and immutable version; required for real records. |
| `sequence_length_nt` / `topology` | Positive length, `linear` or `circular` (or `unknown` if missing). |
| `primary_archetype` / `secondary_archetypes` | Controlled category values above; include rationale in review notes. |
| `annotation_status` | `known`, `unknown`, `ambiguous`, or `pending_review`. |
| `license_note` / `review_note` | Required audit fields before production use; fictional examples marked as such. |

### `features[]` (zero or more source-annotated regions)

| Field | Type and meaning |
| --- | --- |
| `feature_id` / `construct_id` | Stable feature ID and parent foreign key. |
| `feature_label` | Controlled value: `CDS`, `promoter`, `terminator`, `tag`, `linker`, `selection_marker`, `backbone`, `other`, `unknown`. |
| `display_name` | Original provided label or `null`; not a model-generated functional guess. |
| `location` | `{start0, end0, strand}` in **0-based half-open** parent coordinates when confidently available. `strand` is `1`, `-1` or `null` if unknown; null bounds when absent. |
| `segments` | Optional list of parent-coordinate spans for origin-spanning, joined or discontinuous annotations. Do not flatten without preserving segments. |
| `annotation_source` | Source record/feature qualifier or reviewer evidence trail; never silently drop original annotation. |
| `annotation_status` | `known`, `unknown`, `ambiguous`, `pending_review` with an optional `ambiguity_note`. |

### `units[]` (candidate model input representations)

| Field | Type and meaning |
| --- | --- |
| `unit_id` / `construct_id` | Stable child ID linked to its source parent. |
| `feature_ids` / `source_segments` | List of originating feature IDs and original spans so output can be audited. |
| `molecule_type` | `protein`, `dna`, or `metadata_only` (**intended** input domain, not an inferred assay result). |
| `model_input_status` | `pending_model_selection`, `ready`, `unsupported`, or `needs_review`. Before T05 all examples remain `pending_model_selection` or `needs_review`. |
| `representation_note` | Why a unit was created or withheld. Never use speculative reconstruction or completion. |

**Invariants:** every feature/unit references an existing construct; every unit feature reference exists; coordinates satisfy `0 <= start0 < end0 <= sequence_length_nt` where known; unknown feature location remains explicit; `feature_ids` are ordered and stable; IDs do not depend on processing order; raw public feature annotations remain recoverable. For circular-origin crossing, keep multiple original spans instead of fabricating one continuous interval. The model may only receive unit types its **selected** interface supports; metadata-only output is retained for unsupported features.

## Three illustrative benign examples

1. **`example_single_cds.json`** — *fictional reporter plasmid*: one annotated reporter CDS, parent coordinates, one tentative protein unit. Demonstrates one-to-one provenance and the baseline control.
2. **`example_multi_tagged.json`** — *fictional two-CDS research cassette*: two coding annotations plus one labeled affinity-tag segment. Its primary archetype is `multi_cds` and secondary archetype is `tagged_fusion`. Demonstrates one-to-many output and a shared/overlapping feature without inventing a new protein product.
3. **`example_backbone_heavy.json`** — *fictional routine cloning vector*: one benign selectable-feature annotation, a large vector/backbone element and a promoter whose location is not provided. Its promoter is `unknown` and its input is `needs_review`, illustrating abstention instead of speculative coordinates.

All three examples can be checked with `python -m unittest discover -s tests -v`, which validates IDs, parent/child links, basic field types and coordinate sanity; **these are structural tests, not biological validation**.

## Handoff and decisions

**T03 planning sign-off — completed 2026-10-29:** The researcher approved the four archetypes, normalization schema v0.1, and three fictional examples for the Week 1 deliverable. This does **not** claim biological accuracy or validation on real records. T03 can serve as the baseline for T06/T10/T11/T12; the schema must still be checked against a benign annotated pilot and any revisions must be versioned. T05 must decide the model interface before any example is labeled truly `ready`. T10/11 must preserve raw GenBank feature semantics; FASTA-only inputs typically lack feature labels and should receive `unknown` unless a separate approved annotation source exists. Additional annotation keys may be added after a Week 2 pilot, with a documented version bump and migration.

**Source basis:** IRI proposal, pp. 1–2 (mixed constructs and four archetypes); four-week workbook `Weekly To-Dos` T03, T10–T12 and `MVP Criteria` MVP-02/03/04.

## Local T05 extension: v0.2 sequence references

The original v0.1 examples remain unchanged and valid. See [T05 sequence-reference extension](T05_sequence_refs.md) for the small v0.2 addition: local storage references, hashes and explicit protein-coordinate units. The protein-only artificial examples do not invent nucleotide coordinates or imply that a DNA parser has been implemented.
