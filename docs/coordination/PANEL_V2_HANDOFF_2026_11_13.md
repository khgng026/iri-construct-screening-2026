# Panel v2 authorization and review handoff — November 13, 2026

## Researcher authorization

The researcher replied in the project conversation on November 13, 2026:

> I authorize panel v2. Run the cleared experiments and continue progressing on the project. After finishing, analyze your work and in a comprehensive scope and present your results. Include what the most important questions you think should be discussed next and if any other possible preprocessing approaches may be more effective.

This responds to the request for private, bounded research on the exact [panel v2](../../data/manifests/prospective-2026-11-12/panel.v2.json): annotation preservation for AF234299.1 and MK448012.1, and the proposed O58462 reference and round-robin comparison on the existing pilot queries. Its canonical proposal SHA-256, recomputed with `proposal_sha256` omitted, is `f9111dd17ec028fdc170552aacc8af7275fccaf46d1e9712e3f7a681bdc0d70b`. The existing panel bounds prohibit paid services and model downloads.

The researcher subsequently explicitly authorized Codex to verify the packet. The dated [AI-assisted review](PANEL_V2_AI_REVIEW_2026_11_13.md) records this instruction, its evidence and limits. It cleared an in-memory E1 annotation audit of the two exact sources. A subsequent [alignment-only E3 decision](PANEL_V2_ALIGNMENT_REVIEW_2026_11_13.md) cleared the changed reference panel after all source/reference hashes were verified. ESM-2, E2 and new E4 remain on hold in this environment. The user instruction supersedes the policy's prior role separation for this review; the reviewer does not claim independent human review.

## Review and execution handoff

The [current evidence and limitations](../EXPERIMENT_STATUS.md#independent-evidence-and-exact-prospective-review) remain relevant. E1 covers the two supplemental sources for source-annotation preservation only; E2 requests the existing supported parent and a separately labeled round-robin/component-attributed diagnostic with the original primary unchanged; E3 requests reference-composition sensitivity using the existing four queries and O58462. The proposed extra reference does not add an independent query. Fusion retrieval, ESM-C, Renilla and stricter tiers remain outside this proposal.

After any clearance, verify and restore the exact permitted private source/reference inputs and, for inference, the pinned ESM-2 environment and weights. This checkout has neither the original private caches nor the model stack. Missing or mismatched inputs remain blocked; obtaining fresh material is not equivalent to restoring the frozen bytes without verification. Prepare versioned execution manifests and report each stage separately. The existing launcher is scoped to the original clearance and must not be presented as authorizing panel v2.

Follow-up after the review: provider responses for all 17 previously reviewed UniProt reference records matched both their frozen original-byte and protein-sequence hashes, and the three original query-source records matched their frozen byte hashes. The added O58462 record was also verified in the review. These checks were performed in memory and did not restore the private caches. The [alignment-only execution freeze](../results/panel-v2-2026-11-13/E3.alignment.freeze.json) and [source-free result](../results/panel-v2-2026-11-13/E3.alignment.json) record the revised-panel run; the pinned model environment remains unavailable. The historical review text records what had been checked at the time of that decision.

## Current checks

The [in-memory E1 helper run](../results/panel-v2-2026-11-13/E1.in-memory.json) checked the exact provider hashes and retained 43/43 observable features and qualifiers and 5/5 source-aligned translations across the two sources, with zero detected inventions. It used existing normalization and quality-audit helpers, not the file-based production CLI; raw inputs and sidecars were not written. The [repeat script](../../scripts/panel_v2_e1_in_memory.py) verified the same inputs and produced identical report fields apart from timestamp. These are two source-agreement observations, not independent biological validation or evidence of retrieval performance.

E3 panel-v2 alignment found the expected family in 4/4 queries at five hits, with 2/20 wrong-family hits. The earlier 17-reference alignment had 4/4 and 3/20: adding O58462 displaced one wrong-family hit for the OMP query. Reference-self normalization gave 4/4 and 2/20 on the new panel, with different ranks. This is a descriptive catalogue-composition sensitivity on the same four queries, not an independent replication or evidence of general method improvement. The [runner](../../scripts/panel_v2_alignment_in_memory.py) rechecks source-byte and protein hashes and keeps source/protein contents in memory; it does not run the unavailable matching ESM-2 method.

A separate rerun on the same inputs produced the same canonical raw and normalized score hashes and identical scored results; timing varies. The parsing-environment suite passed again (251 tests run, 246 passed and five optional skips). This repeat demonstrates reproducibility, not new independent observations.

On the merged main snapshot, the existing source-free E4 analysis rebuilt byte-identical tracked outputs using the approved E3 evidence and E4 protocol hashes. This is a reproduction of a report, not new inference or a new biological observation. The parsing-environment suite ran 251 tests: 246 passed and five optional tests skipped. The system-Python suite ran 251 tests: 192 passed and 59 skipped. The fictional metadata smoke test passed. These software checks are not a fresh retrieval rerun.

## Research interpretation and next iteration

The original E1 experiment showed source-annotation agreement in seven new records; it did not establish independent biological truth. In the original E3 pilot, alignment and ESM-2 each recovered the expected family for 4/4 queries at five hits, with 3/20 wrong-family hits each. The queries represent only two conservative dependence clusters. The original E2 primary result on one parent was 2/2 supported parent families for the longest eligible component and 1/2 for all eligible components under maximum-score merging. Component attribution shows that the baseline's second family arose from off-family matches and all five all-component hits came from GFP; this points to aggregation and endpoint limitations, not proof that decomposition harms or helps generally. The full-family endpoint is unavailable because AmpR is unsupported. The observed ESM-2 work-time ratio is specific to the measured small run and is not a general cost ranking.

Fast candidate directions to freeze before evaluation:

1. Compare maximum-score merging against the already proposed round-robin ranking under the same global hit budget, reporting both parent-family and component-attributed recovery and wrong-family hits. Round-robin may force low-quality component hits or be sensitive to component order and duplicates; evaluate these tradeoffs.
2. Consider per-component rank normalization or calibrated scores to reduce dominance by one component. Calibration requires a separate training/development set and a held-out evaluation to avoid tuning on the pilot. Rank normalization sacrifices score magnitude and may elevate weak matches.
3. Report a small per-component result set alongside a parent summary; this improves interpretability but changes output volume or budget unless explicitly constrained.
4. Explore domain-aware segmentation only where independent source-supported boundaries exist. Keep full-protein and source provenance controls; invented fusion boundaries would undermine the evaluation.
5. A DNA-native model or multi-modal approach could eventually preserve noncoding context unavailable to the protein-only pilot, but requires its own validated input limits, references and matched baseline, and would add complexity and computational cost.

The highest-priority questions are whether the intended product must surface every supported component, what shared output budget and false-hit burden are acceptable, and which independently annotated multi-family constructs can provide genuinely independent queries. Broader sampling and frozen outcome definitions are needed before ranking approaches or claiming general benefit. Panel v2 remains descriptive and cannot resolve the population-level hypothesis by itself.
