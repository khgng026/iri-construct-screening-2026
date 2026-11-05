# T05 annotation parser: separate audit contract

T03 v0.1 and its three accepted metadata examples remain the normalized baseline.
The parser emits a list of those objects in local `normalized.json`. This document
defines **separate** `t05-annotation-audit-1` sidecars and `t05-inputs-1` experiment
configuration; neither changes T03 fields/enumerations or implies schema approval
for a future production pipeline.

`python scripts/parse_annotations.py --inputs data/processed/reviewed-inputs.json --output results/annotation-demo/run01`

The command requires Biopython and reviewed local DNA inputs. See
[input contract and dependency setup](../../experiments/README.md). It produces:

| Local artifact | Contents |
| --- | --- |
| `normalized.json` | T03 v0.1 construct/features/units, no sequence fields |
| `annotation-audit.json` | Objects keyed by construct ID; source SHA256, original record ID/annotations and feature maps |
| `units.fasta` | Only source CDS `/translation` qualifiers verified against their exact locations/frame/table |
| `report.json` | Provenance digests, versions, counts, timing/resources and limitations |

Each audit feature is keyed by its normalized feature ID and preserves original
feature type, ID, textual location, compound operator and a deep copy of **all**
qualifiers, including translation, product, gene, notes and evidence. Biopython
reference objects in record-level annotations use their text representation.
T03 `annotation_source` links to the corresponding sidecar feature. Original
source bytes remain in the reviewed local file identified by SHA256. Sidecars
and exports are local, ignored artifacts; permissions still apply to them.

Feature/unit IDs use parent ID and a content digest of original feature
semantics/qualifiers, preserving identity when feature order changes. Identical
duplicate feature signatures fail for manual disambiguation. Changes to source
annotations yield new feature IDs. Units retain parent/feature foreign keys and
ordered 0-based half-open nucleotide segments. Protein experiment IDs additionally
prefix the local source `input_id`, as documented in `t05-inputs-1`.

Reverse-strand and exact joined/origin-spanning CDS use Biopython's strand-aware
feature extraction **only to validate an existing annotated translation**. Joined
features preserve segment order and use null bounding `location`; no misleading
continuous interval is fabricated. `codon_start` and `transl_table` are honored.
The complete CDS must pass Biopython's CDS/start/stop checks and exactly match the
provided translation. A single annotated terminal `*` is removed for comparison;
the original qualifier remains intact in the audit.

Missing translations remain `unknown`. Conflicting/multiple qualifiers, sequence
disagreement, fuzzy/remote/missing boundaries, unknown or mixed strands,
unsupported compound operators, pseudo/exceptional CDS, incomplete CDS and
noncanonical proteins remain `ambiguous`/`needs_review` with a reason; they are
not sent to models. Invalid/out-of-bounds coordinates fail. Fuzzy coordinates are
not rounded; any exact segments are retained alongside the original location in
the audit. No ORF inference, completion or remote sequence reconstruction occurs.

Other features are retained as unsupported metadata. FASTA DNA produces an
`unknown` construct annotation state and empty features/units. Archetype labels
come from the reviewed manifest, never sequence guessing. Valid source
translations remain `pending_model_selection` in T03 output; structural/translation
checks do not approve source annotations or complete T05. The experiment runner
separately checks the selected adapter's alphabet, cap and token alignment.

Generated in-memory test records are explicitly fictional software fixtures.
They are not public source annotations, approved reference labels, or experimental
evidence. Real-record manual spot checks remain T06/T08/T14 work.
