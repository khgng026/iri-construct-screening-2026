# Current execution readiness and decision request (2026-11-10)

## Confirmed state of merged main

- PR #6 merged: pBR322 J01749.1 source reviewed and E1 annotation-preservation audit completed. This does not demonstrate model retrieval performance.
- PR #7 merged: 11 NCBI GenBank candidate files successfully retrieved/parsed in ephemeral Actions storage; the committed source-free audit identifies all records as `ACQUIRED_NOT_ADMITTED` and lacks per-record license/annotation clearance.
- PR #8 merged: 10 selected iGEM candidate parts returned published sequence bytes and nonempty feature responses; source-free audit is committed, but the raw responses were not staged; no part-specific T02 scientific admission, reference family approval or license/private-cloud review is complete. Some are cassettes or tagged/fusion constructs, not full plasmids.
- Original nine Addgene catalog records retired by researcher direction. Do not request signed-in Addgene download. Historical evidence remains unchanged.
- E2 `blocked_insufficient_parent_panel`, E3 `blocked_insufficient_panel`, E4 `blocked_no_measurements`. No real retrieval comparisons or confirmed ESM-C inference yet. Hosted CI checks passing do not establish scientific readiness.

## Current decisions

**Authorized by researcher already:** 11 GenBank candidate identities, ten iGEM candidate identities for further review, retire nine Addgene holds, and request a delegated-approval **policy mechanism**.

**Not yet authorized:** a particular real experimental run on a dated frozen panel. Researcher must explicitly approve named experiments and scope after the proposed rule-change question. The earlier selection permission is not execution permission.

**Evidence and engineering required from Codex before ChatGPT can clear a real run:**

1. Capture exact downloaded original files, timestamps and digests in safe, permitted project storage; inspect source terms and author attribution. Limit raw-file Git commits to individually cleared cloud-copy rights.
2. Record per-construct benign scope, annotated feature/CDS/fusion status, integrity/ambiguity/exclusions, and proposed archetypes. Retain a distinct original full-plasmid vs partial cassette classification.
3. Prepare stable, evidence-backed per-protein family labels and *independent* positive/competing references; exclude shared source, identical sequence, shared lineage, and tightly related vectors. Provide a reviewable exclusion map, no post hoc decisions.
4. Choose a defensible prospective pilot panel (do not silently score all 21), lock the experiment protocol, model revisions, query/reference set and scientific denominator before inference.
5. Run offline dependency checks and tests; assess ESM-2, ESM-C, and alignment availability, memory/time budgets. Produce separate E1 readiness, E2/E3 readiness and E4 measurement readiness statuses.
6. Bring the evidence packet to ChatGPT for its delegated review **after** an explicit scoped researcher run authorization. No new scored inference until both authorization and approval exist. If only E1 clears, run E1 and report other stages blocked.

The researcher does not need to review every small technical choice; only a clear scoped run-authorization request, and exceptional scope/safety/rights questions, should be brought back. The assistant can approve based on evidence, never on missing evidence.
