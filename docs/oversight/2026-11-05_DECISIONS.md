# Scientific oversight decision record — 2026-11-05

**Authority:** Researcher instruction in the IRI Fall 2026 project conversation on 2026-11-05. This document records decisions; it does **not** alter the already approved T02 governance policy, change T03 schema v0.1, or approve an unreviewed biological source.

**Status:** APPROVED for candidate *shortlist selection* and T07 *evaluation design in principle*. **Not approved** for biological source admission, source-file redistribution, reference-family labels, or inference on records lacking the recorded T02 review gates. No scientific improvement has yet been demonstrated on real constructs.

## Decision D-001 — ten-construct pilot shortlist approved for review

The researcher approves the following **ten candidate identities for the prioritized pilot shortlist** (not as ten admitted, scored dataset records). This supersedes the prior *selection order* in T08 where pUC19c was prioritized; preserve the original inventory and record this selection change transparently.

| Priority | Candidate record ID | Name | Primary intended archetype | Known review concern |
| --- | --- | --- | --- | --- |
| 1 | `ncbi-J01749.1` | pBR322 | backbone-heavy | Check original versioned source CDS/qualifiers, terms and checksum |
| 2 | `addgene-212936` | pGL3 Basic Vector | single principal reporter | Confirm exact whole-vector file and luc+ CDS |
| 3 | `addgene-18964` | Luciferase-pcDNA3 | single principal reporter | Confirm whole-vector annotation/version and source rights |
| 4 | `addgene-87121` | pcDNA-RLuc8 | single principal reporter | Renilla family/reference independence |
| 5 | `addgene-18673` | pYIC | multi-CDS | Conflicting tag descriptions; confirm separate products |
| 6 | `addgene-64034` | pGL4Luc-RLuc | multi-CDS | Multiple sequence versions and exact CDS annotations |
| 7 | `addgene-108228` | pGreen2x35SLUC | multi-reporter, provisional | Do not mistake an insertion site for an inserted ORF |
| 8 | `addgene-21073` | pEGFP-LC3 | fusion | Missing insert ATG note; confirm full-fusion boundaries |
| 9 | `addgene-51970` | Nano-lantern/pcDNA3 | fusion | Full-vector version vs insert; Rluc8 relatedness |
| 10 | `addgene-255006` | pcDNA3-Renilla-FLAG-firefly | fusion | Exact tag/fusion bounds; luciferase relatedness |

**pUC19c (`ncbi-L09137.2`) remains in the 20-record inventory as an annotation-limited backbone control/reserve, NOT in the approved ten-identity priority shortlist.** Its original structured record appears insufficient for scored protein-family retrieval; assess independently rather than invent missing CDS annotations. All other candidates remain available as documented reserves, with original statuses unchanged until T02 review.

The list is deliberately **not ten independent protein-family examples**. Group linked vectors, shared reporters, duplicate sequences and homologous references conservatively; a scored family retrieval benchmark will need sufficient separately reviewed, independent reference proteins.

**Data-release gates stay intact:** exact source/version and original bytes, license/reuse scope, safety rationale, source annotation review, source-file SHA-256, reviewer/date, independently supported labels and relatedness group. Until recorded, `review_status`, `annotation_status`, `licensing_status` and `safety_review` remain pending/uncertain as appropriate. Never infer approval from this D-001 shortlist decision.

## Decision D-002 — T07 evaluation design approved in principle

The researcher approves the scientific comparison design documented in PR #2:
[`docs/week1/T07_benchmark_protocol.md`](https://github.com/khgng026/iri-construct-screening-2026/blob/codex/real-construct-pilot/docs/week1/T07_benchmark_protocol.md).

Retain:
- Separate **annotation preservation**, **documented-family retrieval**, and **actual function** claims.
- Protein-only model compatibility; noncoding DNA stays metadata.
- Matched source records and common eligible reference panel; exact duplicate, source-parent and relatedness-group exclusions.
- One fixed **total** hit budget (default K=5), never K per CDS or window; separate runtime/model-call budgets.
- Source/protocol/panel hash freezes before scoring, status/coverage/abstentions, and no unsupported labels or fabricated scores.
- Cluster-aware paired uncertainty; no population-wide inference from a small pilot.
- Four planned experiments A–D, with multi-CDS parent aggregation explicitly still prospective.

**Approval of design is not a claim of validated metrics, model performance, approved biological records, or T05/T07 scientific completion.** Codex should record the approval while keeping implementation-scientific-validation statuses separate.

## Decision D-003 — four additional experiments/extensions authorized

Authorize implementation of the four *prospective* extensions in `EVALUATION_EXTENSIONS_V1_1.md`. This must be a **versioned addendum** (v1.1), not a silent edit to the frozen T07 v1.0 protocol. Preserve v1.0 for provenance. No current results are available from these extensions.

## Decision D-004 — future larger dataset

Record the intention to expand beyond the ten-identity pilot **after** the pilot confirms annotation quality, suitable independent references and method feasibility. A practical later target is approximately **150–250 reviewed benign constructs if quality permits**, consistent with the project's original Week 2 aspirational target; it is **not an immediate acceptance criterion** and must not incentivize weak labeling or review shortcuts. Predefine inclusion, archetype/independence targets, versioning and staged QA before expansion.

## Decisions the assistant may make on the researcher's behalf

The researcher delegates ordinary implementation and process-direction decisions *unless existing project rules expressly restrict them*. The assistant may approve technical housekeeping, non-sensitive task sequencing, conservative evaluation-quality guardrails, and document/version coordination. **Do not use this delegation to bypass T02 record-level human review or assert permissions, benign-status verification, independent family labels, or license rights without actual evidence.** A verified reviewer may record a specific T02 approval with supporting evidence; missing/conflicting evidence requires HOLD and escalation. No operational screening, hazardous datasets, model training, or protected data.

## Coordination references

- [PR #1 — draft parallel research (read only)](https://github.com/khgng026/iri-construct-screening-2026/pull/1)
- [PR #2 — pilot implementation (do not overwrite while being developed)](https://github.com/khgng026/iri-construct-screening-2026/pull/2)
- [Existing T02 governance](https://github.com/khgng026/iri-construct-screening-2026/blob/main/docs/week1/T02_data_governance.md)
- [Oversight message bridge](COMMUNICATION_BRIDGE.md)
