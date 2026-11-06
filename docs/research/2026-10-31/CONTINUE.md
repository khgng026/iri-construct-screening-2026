# IRI session continuation and recovery

Status: research closed under the amended eight-hour window; no workers or keep-awake processes remain. Final code revision is `1406697b55680c262764f98bc1216138f0a7401c` on `rudra/iri-session-2026-10-31`. No push or publication occurred. Read [REPORT.md](./REPORT.md) before choosing follow-ups.

## Safe inspection and saved-output reruns

Use the existing interpreters without modifying either environment:

- T05/ESM-2: `/Users/Rudra2.0/Desktop/IRI/.venv-t5/bin/python`
- ESMC/SciPy/plots: `/Users/Rudra2.0/Desktop/IRI/.venv-esm-readiness/bin/python`

The following commands only validate saved evidence. From the repository directory:

```sh
cd /Users/Rudra2.0/Desktop/IRI/iri-construct-screening-2026
/Users/Rudra2.0/Desktop/IRI/.venv-t5/bin/python scripts/session_independent_audit.py \
  results/iri-session-20261031-032506/jobs/esmc-600m-cpu-full-v1-8ae67f23a1d8 \
  --output /tmp/iri-esmc600-saved-audit.json
/Users/Rudra2.0/Desktop/IRI/.venv-esm-readiness/bin/python -m unittest discover \
  -s tests -p test_session_dct_experiment.py
```

The first command independently rebuilds rankings/denominators and samples saved vectors; it does not run the model. The second verifies adapter code, not DCT scientific admission. The full supported test-module list and excluded drafts are in `audit/closeout_test_results_v1.json`; do not treat unreviewed drafts as part of the accepted regression suite.

Each result store has `identity.json`, immutable per-unit attempts, checksummed `complete.json` files and a job `completion-*.json` indexing the exact accepted units. Use that index, not a bare directory count or exit code. `src/iri_constructs/session_store.py` implements publication/validation. `audit/closeout_verify_v1.py` documents the final report/export checks but writes a dated audit path and must be versioned before a later run, preserving the original closeout record.

Original queue commands and frozen configurations are in `inputs/queue-*.json`. The completed closing reports used `queue-closeout-reports-v1.json` (catalogue succeeded; segment roster failed) and corrected `queue-closeout-segment-report-v2.json` (four-model segment report succeeded). Their exact commands regenerate reports from saved sources, with no new model scoring. For a new report, choose a new report name/output identity and verify exports afterward.

## Continuing experiments requires a fresh bounded session

**Do not launch a historical queue unchanged.** `launch.json` and `inputs/deadline_anchor.json` preserve the original 12-hour deadline; the later eight-hour amendment is in `audit/session_limit_amendment_8h_v1.json`. They are historical records, not a current execution allowance. Preserve them. A future authorized session needs a new session directory, immutable launch/deadline, input/source bindings, resources and compute ownership.

Recheck AC, normal memory pressure, swap, free disk and competing jobs using the unchanged prepared guard. Keep one heavy worker, batch one, finite timeouts and the handoff thresholds. Guard SHA-256: `a5435e14ce302c8d7a78844baed1e015ce526931965a33f960552716d439673d`. Last closeout had approximately 26 GiB free before the final recovery copy, so do not assume the launch's 55 GiB remains available. New official weight payloads total 6,837,035,672 bytes of the original 12,000,000,000-byte authorization; reuse verified local weights rather than redownloading. No ESM3 weights were acquired.

Reference completed core jobs and the 1,000-catalogue manifest through immutable hashes. Do not change old labels, pooling, model revisions, precision or thresholds to make checks pass. New designs get new identities and affected reruns. `session_dispatch.py` catches optional failures and validates completion records; `session_guard_entry.py` supplies tested signal cleanup without changing the prepared guard. A future process resuming an incomplete accepted identity must validate completed units, retain corrupt/staging evidence, and recompute only incomplete units.

## Deferred work: exact status

| Work | State at closeout | Requirement before execution |
|---|---|---|
| ESMC 600M automatic segmentation | Designed, not run | New guarded session and current resources; preserve frozen settings |
| Five-model segmentation report in queue 12 | Prepared, incompatible with current comparison roster | Current adapter explicitly expects four original core models plus two SW arms; develop/test/review a new five-model report version first |
| ESMC 600M total-budget oracle and family diagnostics | Prepared, not run | New session; unchanged masks, budgets and denominators |
| Seven-fragment fresh numerical replay | Prepared, not run | Use frozen seven-case design over 1,931 combined attempts / 1,117 distinct cases; do not confuse with already-passed stored-data audit |
| Expanded family scoring | Labels reviewed/frozen; no scoring | Use revised source-admission-v2 design/queue; no expanded near-90% scoring without the frozen global-identity calculation |
| Expanded near-90% identity table | Feasibility only | Implementation/review and resource forecast; no pair alignments completed |
| Midpoint-layer sensitivity | Design and audit implementation prepared | New layer probes and official-API verification before admission |
| pLM-BLAST | Original numerical probe passed; original runtime gate failed | Rerun both resource-v2 probes with the exact reviewed overlay/import path; validate forecast and a new time budget before admitting full comparison |
| DCTdomain | Build succeeded after repair; discontinuous fixture failed | Diagnose against original source/fixture without altering expected labels; new version and independent validation if repaired |
| ESM3 small | Metadata checked; unacquired/unprobed | Access/licence/interface/resources and download headroom checks; no new agreement or bypass |
| Real benign constructs | 15 candidates pending; zero molecular records admitted or scored | Record-level provenance/reuse rights, annotations, boundaries and independent references |
| Unit-normalized weighting and enriched casebook | Unreviewed, untested draft code | Source review, fictional tests, frozen identity and resources before numerical use |
| ESMC ranking compression | Approval record only; no plan/trial/apply | A separately verified plan/trial if later needed; not required to read current results |

The six untracked draft files are `scripts/session_casebook.py`, `scripts/session_unitnorm_controls.py`, `scripts/session_unitnorm_audit.py` and their corresponding tests. Preserve them as drafts; they are excluded from accepted tests and findings. Earlier `SESSION_STATE.md` notes describing these as pending next jobs are superseded by this closure.

## Evidence map

- Five-model core: `./jobs/compare-core-five-models-v1-fb9fecb900a2`.
- Four-model held-out segmentation and family tables: `./analysis/heldout-segment-four-model-closeout-v4-33b8c9f01052`.
- Five-model expanded catalogue plus alignment: `./analysis/expanded-catalogue-closeout-v1-cd8aa1563580`.
- Weighting: `./analysis/equal-component-weight-v1-e9b0c06e0305`; independent audit: `./audit/independent_weighting_audit_v1.json`.
- Context: `./analysis/context-common-repeat-v2-9dbff8491506`; fresh replay job `context-profile-8m-independent-replay-v1-d151fab81d1e`.
- Failure profiles: `./analysis/failure-profile-v1-3ff70216e42c`.
- Original T05 integrity: `./audit/closeout_prior_outputs_v1.json`.
- Attempts, failures and commands: `./job_ledger.jsonl`; compact counts in `audit/closeout_verification_v1.json`.
- Optional method/source records: `audit/alignment_methods_review.md`, `audit/independent_dct_reccut_index_review.md`, `audit/plmblast_resource_v2_overlay_binding_execution_recipe.md`, `audit/real_benign_genbank_followup_v1/README.md`. Historical reviews' “pending” language is superseded by the report's final statuses.

## Recovery copies and limitations

The first verified snapshot, `snapshots/core-verified-20261031T0913Z`, covers six earlier core stores. The final selective snapshot, `snapshots/closeout-20261031T1704Z`, adds ESMC 600M core, equal-weight controls, order/context/failure/competition reports, the five-model comparison, inputs, selected audits, code and a fixed ledger prefix. It copies only selected small segmentation exports and manifests; the large segmentation report and expanded ranking stores remain in their original locations. Its manifest explicitly lists coverage and hashes. This is not a complete duplicate of all results and does not protect against drive failure.

Inspect `snapshot_report.json` and `manifest.json` before restoring. Restore into a **new** directory, compare hashes, and preserve the existing source as evidence. Do not overwrite scientific stores or replace newer ledger records with the saved prefix. A snapshot of selected exports is not itself a resumable full UnitStore.

The native compression audit verified all 3,441 selected records and 6,745,153,107 logical bytes. Only redundant owned backup links were removed after verification; no scientific results were deleted. There are no pending compression transactions.
