# E1: Annotation quality and appropriate abstention — v1.1

This prospective extension implements the researcher-approved E1 direction.
It preserves T02, T03 v0.1, the original T06 inventory and frozen T07 v1.0.
It does **not** admit source records, approve source annotations or establish
protein function. No real-data quality measurement is claimed without admission.

The [implementation](../../src/iri_constructs/annotation_quality.py) reuses
`pilot.prepare_pilot` and `audit_annotation_preservation`. All T02 source,
annotation, safety and permissions gates, exact source version/digest and existing
T03 checks apply before the executable reads a biological source. Candidate
selection, discovery metadata and the pending review queue confer no admission.

## Source strata fixed before model outcomes

| Stratum | Rule from reviewed original annotations | Interpretation |
| --- | --- | --- |
| `well_annotated` | At least one documented CDS; no unresolved substantive source feature or CDS-translation issue | Source annotation is sufficient under the declared extraction policy, not independently proven biological truth |
| `partially_annotated` | At least one documented CDS plus missing, conflicting, fuzzy, unsupported or uncertain substantive annotation | Faithful withholding is appropriate; incomplete source evidence is not automatically parser error |
| `annotation_limited` | No documented CDS, including feature-free FASTA, source-only bookkeeping or noncoding-only records | No invented protein unit; report any genuine noncoding annotations as metadata |

Exact reverse-strand and joined/multipart annotations can be well annotated.
Multipart and reverse counts are reported separately and are not defects by
themselves. Missing/fuzzy/remote locations, unknown/mixed strand, exceptional CDS,
missing translation, contradictory qualifiers and failed source-coordinate/table
translation validation are explicit concerns. Strata are determined from original
source evidence before assessing the parser's output or any model outcome. A
parser losing a valid source translation does not make the source low quality.

`source` bookkeeping features are excluded from the **observable** feature
preservation denominator. A source-only record reports
`no_observable_annotations`, null observable preservation/input rates and a
justified record-level abstention. A record with documented noncoding features
but no CDS reports `no_observable_cds`; its noncoding metadata can be preserved
without claiming protein-model coverage. The reused structural audit is also
retained with an explicit bookkeeping-inclusive scope.

## Audits and denominators

Every rate is an object with `numerator`, `denominator` and `value`; a zero
denominator produces null, never an invented zero or perfect score.

| Metric | Numerator / denominator |
| --- | --- |
| Observable feature preservation | Retained substantive original feature IDs / substantive original features |
| Observable qualifier preservation | Unchanged original qualifiers on substantive features / substantive original features |
| Confirmed source-aligned translation rate | Exported proteins matching the original validated annotated translation / documented source CDSs |
| Useful model-input rate | Correctly exported, known canonical source translations within the declared residue cap / documented source CDSs |
| Justified model abstention | Correctly withheld or explicitly incompatible CDS units / CDSs requiring annotation or cap abstention |
| Source uncertainty preservation | Original uncertain substantive features retained as unknown/ambiguous / original uncertain substantive features |
| Missing/conflicting CDS annotation | CDSs with the respective source issue / documented source CDSs |
| Record input/abstention rates | Records with usable input / audited records; justified abstaining records / records requiring abstention |

Per-record reports retain the existing independent preservation checks for stable
feature/unit IDs, exact original location text and ordered segments, strand,
qualifiers, source-aligned translations, source checksum and parent-child links.
Invention counters distinguish new/duplicate feature identities, unsupported
feature-label changes, and unsupported/untraceable protein exports; their sum is
`erroneous_invention_count`. These are parser-audit errors, not biological labels.
Metadata-only noncoding units never become proteins.

The residue cap defaults to 1022 and may be reduced by the frozen E1 config. It
does not affect source-quality classification: a long, well-annotated protein can
be incompatible with the selected input policy. Its annotated local export is
not a model-ready assertion. E1 runs no pretrained inference and assesses the
appropriate abstention policy explicitly.

Withholding reasons and unit states are recorded by stable IDs, without sequence
or original qualifier payloads in the quality report. Aggregates are available by
archetype and source-quality stratum. Failed source parsing or preservation audits
remain visible, including invention counts, but are excluded from successful-rate
aggregates. If all observable audits fail, `scientific_scores_produced` is false.
An empty/unadmitted run produces no quality measurements.

## Reproduce

Use the existing parsing dependencies. The shared v1.1 CLI freezes its E1 design
before invoking the module:

```powershell
python scripts/run_evaluation_extension.py --config config/extensions/E1.json --inputs config/pilot/inputs.pending.json --output results/extensions/E1/pending-run01
python -m unittest discover -s tests -p test_annotation_quality.py -v
```

`inputs.pending.json` contains no admitted source records; its run must remain
`no_inputs`/blocked with no scientific metric. For approved local sources, replace
it with the reviewed `t05-inputs-1` mapping. Local outputs include the unchanged
prepared T03/audit artifacts under `prepared/`, `quality-audit.json` when records
can be audited, and `report.json`. The shared runner owns protocol freezing and
top-level configuration/provenance; raw qualifiers and exports remain ignored.

API: `classify_quality(source_counts)`,
`audit_quality(record, document, audit, extracted, max_residues=1022, source_sha256=None)`
for an already admitted record, and
`run_annotation_quality(inputs_path, output_path, max_residues=1022)` for the fully
gated pipeline. Low-level audit helpers never grant source approval.

Tests use explicitly fictional in-memory/temporary records and numerical
corruption controls. They cover source-only/no-CDS, missing/conflicting/fuzzy,
reverse joined features, lost valid translations, invented features/exports,
checksum mismatch, empty admission and independent cap handling. Test gates are
software assertions, not biological approvals or experimental results.

## Versioned candidate selection

[Shortlist v1.1](../../data/manifests/pilot_shortlist_v1_1.json) records approved
candidate identities/order only. The current
[pending queue](../../data/manifests/t08_pilot_review_queue.json) includes
`addgene-108228` and excludes `ncbi-L09137.2` from the ten priorities. pUC19c is
retained as an optional, unadmitted annotation-limited control hypothesis, not a
validated low-quality source record. All source/annotation/licensing/safety gates
remain pending/uncertain. The [v1.0 queue](../../data/manifests/t08_pilot_review_queue_v1_0.json)
is archived byte-for-byte; the original twenty-record discovery inventory is intact.

```powershell
python scripts/inventory_candidates.py --inventory data/manifests/t06_candidate_inventory.json --shortlist data/manifests/pilot_shortlist_v1_1.json --check-pilot-queue data/manifests/t08_pilot_review_queue.json --output results/extensions/shortlist-review/run01
```

The optional `--shortlist` changes candidate review order only. Without it, the
existing API still reproduces v1.0. Source approval fields are forbidden inside
selection decisions, including optional controls. The revised queue binds both
the original inventory and selection manifest by canonical JSON SHA256.
