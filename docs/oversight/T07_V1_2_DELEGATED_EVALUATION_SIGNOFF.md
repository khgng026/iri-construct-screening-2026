# T07 v1.2 | Delegated experimental-scoring authorization

**Authorized policy-change request:** 2026-11-10. **Actual run approval remains pending.** This addendum changes the evidence-review workflow; it does not alter the predeclared E1–E4 scientific endpoints, baseline models, shared five-hit budgets, comparison controls, or the original protocol freeze. Retain T07 v1.0 and its v1.1 extension as historical scientific versions. See [T02 v1.2](T02_V1_2_DELEGATED_REVIEW.md).

## Gate sequence

1. Researcher explicitly approves a clear scoped experiment request (named eligible GenBank/iGEM subset or versioned panel, E1–E4 stage(s), model and compute budget, source/output storage scope). This governance request and previous candidate selection approval alone are insufficient.
2. Codex prepares record-level T02 evidence, original source hashes, original annotated features, accepted CDS units, provisional family/reference proposal and a **prospective** scored-panel/protocol freeze before any new inference. Exact-copy, shared-parent, lineage and near-duplicate exclusions remain conservative and documented.
3. ChatGPT reviews that packet as a delegated evidence reviewer, with the researcher's run authorization already recorded. For each case, writes `APPROVE`, `HOLD`, or `REJECT` with factual evidence, assumptions, reviewer label and dated scope. An explicit `CLEAR_FOR_E1`, `CLEAR_FOR_E2`, `CLEAR_FOR_E3`, `CLEAR_FOR_E4` determination may be issued **only** for eligible experiment stages. ChatGPT cannot authorize missing independent families/refs by opinion alone. A contested family assignment or insufficient evidence remains unscored.
4. Codex executes only cleared stages and manifests, and reports inference was performed only if it actually ran. Stage failures remain blockers, not manufactured scientific outcomes.

## E1–E4 readiness requirements

- **E1:** real source bytes, per-file T02 admission, traceable original annotated features and explicit ambiguity handling. No protein family reference requirement for source-to-parser annotation agreement.
- **E2:** validated independent protein families and competing references for each scored parent, conservative shared-parent/copy/lineage exclusions, same reference panel, longest-CDS vs all-documented-CDS comparison, common five-hit budget. No saturation-based advantage claims.
- **E3:** all E2 source/independence evidence, frozen tier memberships, exact pinned model snapshots and environment validation, same queries and reference sets across methods. Unsupported strict tier 3 can remain blocked if a non-post-hoc exclusion rule cannot be justified; this does not imply approval of it.
- **E4:** only real measured E2/E3 outcomes with matching hardware/software/panel/cache accounting; incompatible records remain separated.

## Amendments, limitations, stop conditions

- Candidate pool has expanded to 11 NCBI + 10 iGEM. The earlier frozen pilot was approximately 8–12 records. The researcher may authorize a new **prospective** panel from the eligible pool; do not silently substitute into, or change the denominators of, a completed frozen protocol. This policy is not automatic approval to score all 21.
- Genuine protein fusions must be source-documented and have reviewed amino-acid boundaries. iGEM cassettes are not interchangeable with complete plasmids; annotations are not biological-function proof.
- Retain alignment baseline and no training/fine-tuning/raw plasmid input to a protein model. No public data redistribution without separate rights and human release approval.
- Missing rights, hazardous annotations, unsupported reference independence, model-access limitations or resource uncertainty must be explicitly logged and held. Reviewer cannot waive these scientific gates.

**Experiment permission record:** `docs/oversight/EXPERIMENT_RUN_AUTHORIZATION.md` is to be created only **after** a researcher explicitly responds affirmatively to a scoped approval question and the assistant validates evidence. Do not fabricate this file to make CI pass. The current status is `NOT_AUTHORIZED_TO_RUN_NEW_SCORED_EXPERIMENTS`.
