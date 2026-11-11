# T02 | Public-data and safety rules

**Week:** 1  |  **Budget:** 1 hour  |  **Status:** Approved v1.0 (2026-10-28); v1.1 storage amendment merged 2026-11-06; v1.2 delegated-review amendment proposed 2026-11-10  |  **Prerequisite:** T01  |  **Calendar source:** `Weekly To-Dos!A6:K6`.

**Approval record:** User approved the written governance policy on 2026-10-28 and reaffirmed T02's approval with T01–T04 on 2026-10-29. Version 1.0 remains the baseline for source selection, provenance, licensing and safety review. The subsequent v1.1 storage amendment below supersedes only the repository-storage restrictions after individually verified approvals. Future amendments should be versioned and documented; this approval does **not** approve any specific sequence, accession, dataset, source license, or release. The checklist below remains unchecked because each entry must be applied independently to each future candidate record in T06/T08.

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

## T02 v1.1 — user-authorized public-pilot acquisition and cloud/GitHub storage amendment (2026-11-05)

**Authority and scope.** The researcher explicitly authorizes acquiring the ten selected publicly described, ordinary-benign pilot candidates and permits their source files to be stored in project cloud storage or a private GitHub repository **once individually reviewed and cleared for the relevant access and storage rights**. This section supersedes the earlier blanket local-only/raw-file Git exclusion above for the approved eligible records only; all other v1.0 inclusion, exclusion, safety, provenance, annotation and independent-label checks continue to apply.

**Key distinction:** Being publicly viewable is not itself a grant of permission to upload or redistribute a provider's files. An NCBI accession, an Addgene plasmid catalogue ID, and a downloadable full-plasmid GenBank file are not interchangeable proof of licensing or complete annotation. Private GitHub storage still involves a cloud copy and access by repository collaborators. Confirm the exact record's use and hosted-copy rights before pushing it, and record that finding. Do not upload records obtained under access terms that forbid such storage.

**Acquisition authorization:** The researcher approves obtaining the ten prioritized pilot candidates from their authorized public sources, including any necessary ordinary permitted access. Acquisition does not admit them into the scored dataset. For each downloaded record, preserve exact source URL, complete record/sequence version, date, SHA-256 of the original downloaded bytes, observed feature table/annotation status, safety review basis, documented rights and reviewer/date.

**Two explicit gates per source:**
1. `source_admission` / `scoring_admission`: Existing T02 review statuses, correct source annotation and independent endpoint labeling must justify inclusion in real scored experiments. An incomplete/ambiguous CDS may be held even if downloading is allowed.
2. `cloud_storage_approval`: An explicitly recorded rights determination must permit uploading the exact source content to the chosen private cloud/GitHub location, with attribution and any required restrictions. `cloud_storage_approval` can be true independently of scientific score eligibility, but not without source permission and benign-scope review.

**Storage after these gates:** Individually cleared files may be added under a dedicated, clearly documented `data/approved_public/` path (or a configured, access-restricted cloud store) with a source-specific metadata/provenance manifest. This is an exception to the earlier .gitignore/local-only defaults only for exact reviewed, redistributable files. Avoid changing .gitignore to allow all `.gb`, `.gbk`, `.fasta` or `data/raw/` recursively. Do not commit unapproved files, credentials, nonpublic records, restricted annotations, or model embeddings. Keep data storage private by default, and treat a public repo release as a separate publication decision.

**Current pilot state as of authorization:** Ten candidate identities are researcher-selected, and acquisition/cloud-storage eligibility is authorized in principle. **Zero exact downloaded biological records have completed the independent T02 source/version/checksum, annotation, permission, and safety checks, so zero files are currently admitted for scientific scoring or GitHub upload.** Do not change their manifest `pending`/`uncertain` statuses to `approved` until the corresponding evidence exists. This is not a denial of the researcher's permission; it is the implementation of the evidence conditions that remain in force.

**Review evidence requirements and decision record:** See `docs/oversight/T02_V1_1_PUBLIC_PILOT_STORAGE.md`. Any future modification of T02 scoring or content eligibility must be versioned separately.

## T02 v1.2 — researcher-initiated, delegated experimental-data admission (2026-11-10)

The researcher authorized a **new approval workflow** in which ChatGPT may conduct and sign a documented, evidence-based review of the selected benign public constructs, protein-family labels, independent-reference panels, and readiness **after the researcher explicitly approves the scoped experimental run**. This user request authorizes the governance *mechanism*, **not the actual experiment execution**. Its complete signed-off decision criteria and audit record are in [T02 v1.2 delegated-review amendment](../oversight/T02_V1_2_DELEGATED_REVIEW.md). Earlier T02 source and reuse rights, safety exclusion, provenance, evidence standards and narrow cloud-storage controls remain mandatory. No absence of evidence can be changed into an approval.
