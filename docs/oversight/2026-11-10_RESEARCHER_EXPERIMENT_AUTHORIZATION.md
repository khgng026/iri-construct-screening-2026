# Researcher authorization of scoped E1-E4 experimental program

**Decision date:** 2026-11-10 (America/New_York).  
**Authority:** Researcher's explicit approval in project conversation.  
**Decision status:** `RESEARCHER_SCOPE_APPROVED_PENDING_POLICY_MERGE_AND_EVIDENCE_CLEARANCE`.  
**Scope reference:** proposed [PR #9](https://github.com/khgng026/iri-construct-screening-2026/pull/9), T02 v1.2 / T07 v1.2 delegation workflow.

## Researcher's authorization (verbatim)

> I approve the proposed experimental authorization for E1–E4 on the eligible subset of the 11 GenBank and 10 iGEM constructs, subject to all stated evidence, safety, rights, resource and protocol conditions. I authorize ChatGPT to perform delegated evidence review and clear eligible stages once PR #9 is merged.

## Permitted scope once conditions are fulfilled

- E1 annotation preservation, E2 multi-CDS parent comparisons, E3 pretrained-model versus alignment retrieval, and E4 computation-versus-quality evaluation, **only** for a separately admitted and prospectively frozen eligible subset of the eleven NCBI/GenBank plus ten iGEM researcher-selected candidates.
- ChatGPT may review the documented source-by-source data provenance, annotations, permitted analysis/storage rights, benign scope, source identities and sequence hashes, supported family labels and independent reference/exclusion groups, and issue a dated, explicit clearance or HOLD for each stage. The researcher need not sign off on each routine record when ChatGPT has defensible evidence.
- Codex may prepare data and evidence, run non-scored software/readiness tests, freeze proposed evaluation panels, and run the **cleared** E1-E4 stages only after the authorization/clearance gates below. Use the established model protocols, existing computational resources, and private research outputs; **no additional datasets, external spending, or public release** is authorized.
- The nine previously selected Addgene records remain retired from the active plan. Preserve historical approvals and previous pBR322 E1 data. The prior frozen pilot and model/hit-budget hypotheses cannot silently expand to all 21 candidates; any new scored panel must be prospectively versioned.

## Separate gates (not waived by this authorization)

1. **Policy:** PR #9 must be merged so the delegated-review governance amendment is in force. A draft/open PR does not make its rules active.
2. **T02 admission:** Exact source bytes/version/sha256, source-specific rights/attribution, benign suitability, annotation/translation validity, and distinct local analysis versus cloud-hosting permissions must be checked for each eligible record. Exclude/HOLD any unsupported entry; do not manufacture data.
3. **T07 reference and scoring:** Per-protein family label evidence, independent positive and competing references, conservative exact-copy/shared-parent/lineage exclusions, protocol/manifest/model freeze before examining scores, and acceptable resource budget must be documented. Only E1 is independent of protein reference labels; E2/E3 are not.
4. **Delegated approval:** ChatGPT must inspect the actual evidence packet and issue `CLEAR_FOR_E1`, `CLEAR_FOR_E2`, `CLEAR_FOR_E3`, and/or `CLEAR_FOR_E4` with explicit panel fingerprint and review date. No clearance may rely solely on this authorization or the existence of successful acquisition/download checks.
5. **Execution:** Codex must execute only cleared stages and accurately report run status, blocked experiments and costs. E4 requires actual comparable E2/E3 measurements. No post hoc change to frozen labels or thresholds.

## Current stage

Researcher authorization is **recorded**. As of this decision, PR #9 remained draft/open and the newly selected NCBI/iGEM records had not been fully admitted or assigned reviewed independent protein-family reference panels. **No new scored experiment has yet been cleared or executed under this authorization.**

**Next action:** Merge PR #9, prepare a detailed source/label/reference/compute approval packet, and request ChatGPT's documented delegated evidence review. Codex's preparation instructions are in `docs/coordination/CODEX_PREPARE_EXPERIMENT_AUTHORIZATION.md`. The execution clearance is a separate future record, not this authorization statement.
