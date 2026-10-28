# T02 | Public-data and safety rules

**Week:** 1  |  **Budget:** 1 hour  |  **Status:** Approved — finalized policy v1.0 (2026-10-28)  |  **Prerequisite:** T01  |  **Calendar source:** `Weekly To-Dos!A6:K6`.

**Approval record:** User approved the written governance policy on 2026-10-28. This version is the current baseline for source selection, provenance, license and biosafety review. Future amendments should be versioned and documented; this approval does **not** approve any specific sequence, accession, dataset, source license, or release. The checklist below remains unchecked because each entry must be applied independently to each future candidate record in T06/T08.

## Purpose and permissible data

Use **only public, benign synthetic construct records** for an evaluation of annotation-aware preprocessing. Candidate sources named in the research proposal are **Addgene** and **NCBI**; neither public visibility nor a downloadable sequence automatically authorizes redistribution. Capture the terms specific to each source and record before using it. No proprietary synthesis orders, customer identifiers, unpublished laboratory submissions, protected accessions or purchased screening datasets.

### Inclusion gate: a record is eligible only if all apply

- [ ] An accessible public source page and durable accession/record ID (including version where applicable) exist.
- [ ] The stated purpose and annotations are consistent with **ordinary benign research** (e.g., non-hazardous reporter, cloning, or routine expression research); reviewer records the basis for that decision.
- [ ] The record has enough documented features to support one of the four archetypes without guessing from raw sequence alone.
- [ ] Source terms, licensing/reuse information, and any redistribution restrictions have been checked and recorded; unclear terms mean **hold**, not implicit permission.
- [ ] A reviewer confirms that the construct is suitable for the benign-only project. Any uncertain case remains **pending** until reviewed.
- [ ] Identifiers, retrieved date, input checksum, version, source annotations and review notes are retained for reproducibility.

### Exclusion / escalation gate

- [ ] Exclude records associated with hazardous agents or toxins, suspect pathogenicity, explicit evasion/obfuscation, or unreviewed dual-use concerns; do not assemble a library of hazardous sequences.
- [ ] Exclude private or customer-facing data even if obtainable through credentials, a company account, or collaborators.
- [ ] Exclude records where provenance or source permissions cannot be established; mark as `rejected` or `uncertain`, with reason.
- [ ] Do not treat absent, conflicting, inherited or incomplete annotations as verified ground truth. Record `unknown`/`ambiguous` and review.
- [ ] If any record raises biological safety or release concerns, stop processing it and seek project oversight before inclusion or disclosure.

## Required provenance and review fields

**One row per source record** in a versioned manifest (`data/manifests/`):

| Field group | Required fields / rule |
| --- | --- |
| Identity | `record_id`, `source`, `source_url`, `accession`, `accession_version` |
| Reproducibility | `retrieved_on`, `checksum_sha256` (once downloaded), record/annotation version notes |
| Construct summary | `archetype`, `sequence_length_nt`, `feature_count` (derived only after review) |
| Review gates | `annotation_status`, `review_status`, `licensing_status`, `safety_review` (values: `pending`, `approved`, `rejected`, `uncertain`) |
| Audit explanation | `notes`: explicit scope of permission, missing annotations, reviewer, date, inclusion/exclusion rationale; never sensitive free text |

**Feature-level provenance** is specified in T03: a source feature ID, type/label, coordinates and strand when present, original qualifier/annotation source and evidence status. `checksum_sha256` identifies the downloaded local file without storing the contents in Git.

## Storage, sharing, and version-control conventions

1. The repository contains **only fictional metadata examples, schema specifications, code, and tests** until source and safety checks are approved. Raw DNA/protein files, restricted source copies, personal/customer data and secrets are excluded by `.gitignore`.
2. Use `data/raw/` and `data/processed/` locally for approved inputs; record hashes and source metadata in a shareable manifest only when source terms allow. No raw sequence in application logs, exception messages, issues or screenshots.
3. Document whether a derived artifact may be shared; default to **not redistributing** downloaded sequences or source annotations until licensing and safety reviews are complete.
4. Keep `pending`/`uncertain` entries out of the scored reference set. Distinguish unclear data from negative examples. Report exclusions and sample size changes.
5. For duplicates and shared backbone records, track linkage and avoid presenting highly overlapping entries as independent observations.
6. Freeze reviewed manifest version and label policy before downstream evaluation; save change history and reasons for rejected cases.

**T02 policy milestone — done (2026-10-28):** The written rules and manifest fields have been approved as v1.0. **Operational application remains pending:** T06/T08 must apply every applicable check to each proposed record before any candidate enters a scored dataset. No real accession, source permission, biological classification, or record-level inclusion decision has been verified here.

**Source basis:** IRI proposal, p. 2 (public Addgene/NCBI examples; no proprietary/customer data); four-week workflow `Overview`, `Weekly To-Dos` T02/T06/T08 and `MVP Criteria` SCOPE-01. The papers provided with the project discuss screening research; they are not a license to obtain or publish their restricted datasets.