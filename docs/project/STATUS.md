# Current project status

Updated November 16, 2026 (America/New_York). This is the only active status ledger. [Scope](SCOPE.md), [autonomy](AUTONOMY.md), [contract](TOOLKIT_CONTRACT.md) and [plan](PLAN.md) are authoritative. Historical reports retain the dates and limits of their original observations.

## Reset action items

| Action | Current disposition |
| --- | --- |
| Supersede old acceptance criteria | Done in SCOPE.md and root AGENTS.md: toolkit delivery is independent of E2–E4, a large dataset or proof of a performance benefit. Historical criteria are preserved. |
| Establish standing autonomy/resource policy | Done in AUTONOMY.md: routine in-scope review, dependencies, pinned ESM-2 download, suitable private caches, tests and reviewable PRs can proceed without repeated panel approval. |
| One current plan/status and clear history | Done: PLAN.md and this page. HISTORY.md identifies superseded planning, governance and experiment entry points. |
| Choose v0.1 input/output contract | Done in TOOLKIT_CONTRACT.md. It is a target interface; public CLI/API implementation is still pending. |
| Preserve local panel-v2 evidence | Prior local E1 and alignment artifacts/reviews/runners included in this update with a provenance index; no new experiment claimed. Historical paths/bytes retained where bound by hashes. |
| Make clean package installation the release gate | Done in contract and U5; installed CLI packaging work is still pending. |
| Defer extensions | Done in PLAN.md; existing research remains accessible, outside core release dependencies. |

## Existing implementation and evidence

The project already has source-guided normalization, translation validation, provenance, parser QA, test coverage and optional model/scoring modules. The current public-facing scripts still expect repository-specific input/review manifests and output roots. `pyproject.toml` does not yet expose the new installed CLI or parsing dependencies. Those are U1–U5 work, not completed by this policy reset.

The merged original pilot reports source agreement on seven new records and a tiny matched alignment/ESM-2 comparison; see [historical research report](../EXPERIMENT_STATUS.md). The separate [panel-v2 evidence index](../results/panel-v2-2026-11-13/README.md) records 43/43 observable features and 5/5 translations in its two-source E1 helper audit and 4/4 expected families with 2/20 wrong-family hits in alignment only. These prior runs do not prove a general preprocessing benefit. The panel-v2 model comparison was not run; its old no-download restriction is preserved in its frozen proposal and is not active authority for a separately versioned future run.

Validation of this reset: smoke passed; parsing-environment suite ran 251 tests, with 246 passed and five optional skips; historical pilot-release verifier passed. Active documentation relative links resolve. Preserved E1 review, alignment review and alignment freeze hashes match their saved result bindings; result JSON parses. No new packaging success or model inference is asserted here.

## Next work

Execute U1–U5 in PLAN.md as one bounded implementation stretch; U1 starts by separating the public parser API from benchmark approval metadata. A further permission question is unnecessary for in-scope implementation. U6 is optional once core release acceptance passes.
