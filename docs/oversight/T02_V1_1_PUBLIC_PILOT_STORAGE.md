# T02 v1.1 storage authorization and admission status

Date: 2026-11-05 (researcher authorization). Original T02 v1.0 approved 2026-10-28.

## Decision
The researcher authorizes downloading the ten previously selected public benign construct candidates and **changes the project's categorical local-only storage restriction**. After exact-file review, cloud/GitHub storage is permitted for a record if the source's terms allow that hosting and the record passes benign-scope checks. The permission applies to the reviewed record only, not all files from that database or source. Repository access stays private by default. This is a new option, not a requirement to commit biological sequences.

## Distinct statuses
- `candidate_selected`: the ten named identities are approved as a prioritized discovery list.
- `acquisition_authorized`: the researcher authorizes legitimate download of those candidate records under their providers' terms.
- `downloaded_locally`: evidence of an actual exact source file, date and original-byte checksum.
- `cloud_storage_approved`: reviewed terms specifically permit project cloud/GitHub hosting of this exact source file, with attribution and any limitations recorded.
- `source_admitted_for_scoring`: all T02 source, safety, annotation, version/checksum and relevant independent-label requirements are met. Protein-family scoring additionally requires independent references and labels.

**Do not set the last three statuses based only on this authorization.** No exact ten-record pilot release has yet been created or scientifically admitted.

## Prioritized candidates
1. `ncbi-J01749.1` — pBR322
2. `addgene-212936` — pGL3 Basic Vector
3. `addgene-18964` — Luciferase-pcDNA3
4. `addgene-87121` — pcDNA-RLuc8
5. `addgene-18673` — pYIC
6. `addgene-64034` — pGL4Luc-RLuc
7. `addgene-108228` — pGreen2x35SLUC
8. `addgene-21073` — pEGFP-LC3
9. `addgene-51970` — Nano-lantern/pcDNA3
10. `addgene-255006` — pcDNA3-Renilla-FLAG-firefly

pUC19c is a prospective annotation-limited control, not one of the ten.

## Required follow-through
1. Acquire only permitted full-file versions and capture source URL, access date, version and SHA-256.
2. Inspect the exact source's original feature table, uncertainties and valid protein translations. Do not infer missing annotations.
3. Document case-specific data-use and cloud/GitHub hosting permission; log reviewer, date, attribution, and chosen private storage/access policy.
4. Review benign research suitability. Exclude or hold any concerning/uncertain record.
5. For scientific scoring, separately confirm documented families, independent references and pre-frozen panels. Keep unreviewed fields pending.
6. Stage cleared files individually under `data/approved_public/` with a provenance manifest, using narrow allowlisting, then commit only after all required checks pass. Do not broadly unignore `data/raw/`.
7. Retain linked/relatedness grouping, conservative exclusions and a complete audit trail.

The researcher may later expand toward 150-250 reviewed public benign constructs if pilot evidence and data rights support it; this authorization does not approve bulk acquisition or automatic admission.
