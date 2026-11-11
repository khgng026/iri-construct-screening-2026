# Additional iGEM construct shortlist: audit and selection

**Decision:** Ten *additional* iGEM constructs/parts are conditionally approved for source-level T02 review under the researcher's 2026-11-10 delegation. This is **candidate approval, not scientific scoring admission**. Existing eleven GenBank approvals remain untouched.

## Verified download evidence

GitHub Actions [rate-aware audit](https://github.com/khgng026/iri-construct-screening-2026/actions/runs/42942853053) downloaded all 16 selected public part sequences and checked part-specific sequence-feature responses through the anonymous iGEM Registry REST API. **Ten had non-empty feature sets** and are selected; the other six had zero sequence-feature annotations in the API audit and are reserve/annotation-limited only. All 16 returned published statuses; API sequence and feature-response SHA-256 digests and metadata are in [report.json](report.json). This is an availability check, not an independently reviewed GenBank file, translation validation, family annotation, or biological truth.

**What was downloaded:** The original REST JSON for complete published part sequences plus source feature JSON. Raw sequence-containing responses were retained only in ignored temporary GitHub Actions storage and were not committed. The Registry advertises annotated GenBank, FASTA and SBOL exports through the part page; native GenBank download endpoints were **not** independently exercised by the audit. Codex can re-download original public JSON and features and convert conservatively, preserving original files/annotation source for traceability.

## Ten provisionally selected constructs

| Part | Documented context | Feature annotations | Important caveat |
|---|---|---:|---|
| [BBa_25FDWALN](https://registry.igem.org/parts/bba-25fdwaln) | multi_cds | 24 | Multi-component red/green fluorescence device; examine reported coding features, not just design text. |
| [BBa_252QAUN3](https://registry.igem.org/parts/bba-252qaun3) | multi_cds | 21 | Multi-regulator green reporter; related to 25B5R4ZL. |
| [BBa_26NKZH8J](https://registry.igem.org/parts/bba-26nkzh8j) | tagged_fusion | 5 | mCerulean3 with C-terminal His6 tag; confirm actual annotated fusion bounds. |
| [BBa_25NPQIRU](https://registry.igem.org/parts/bba-25npqiru) | backbone_heavy | 16 | GfasPurple full plasmid candidate, family-linked to 25UIP3DW. |
| [BBa_25ULGTE2](https://registry.igem.org/parts/bba-25ulgte2) | backbone_heavy | 20 | TyrBm pigmentation plasmid candidate, family-linked to 25ECYPMT. |
| [BBa_25B5R4ZL](https://registry.igem.org/parts/bba-25b5r4zl) | multi_cds | 25 | Optogenetic green reporter plasmid; shared YF1/FixJ lineage with 252QAUN3. |
| [BBa_25K8XDBD](https://registry.igem.org/parts/bba-25k8xdbd) | tagged_fusion | 8 | BslA-dCBM engineered coating fusion; manual benign-scope confirmation before admission. |
| [BBa_26BQEN31](https://registry.igem.org/parts/bba-26bqen31) | multi_cds | 4 | Nitrate-responsive sfGFP/lacZ insert; not a whole plasmid. Validate each CDS feature. |
| [BBa_25UIP3DW](https://registry.igem.org/parts/bba-25uip3dw) | single_reporter | 12 | GfasPurple expression cassette, related to 25NPQIRU; not independent protein reference. |
| [BBa_25ECYPMT](https://registry.igem.org/parts/bba-25ecypmt) | single_reporter | 12 | TyrBm tyrosinase expression cassette, related to 25ULGTE2; not independent protein reference. |

All ten were retrieved anonymously and have nonempty deposited feature responses. They are distinct Registry entries but some share the same target proteins and are **not independent reference samples**. Separate full plasmids from inserts/cassettes. Do not assume annotations are complete, that tags are correctly bounded, or that downloaded data are experimentally validated.

## Six accessible but not selected for this E1-oriented supplemental set

- [BBa_K4046704](https://registry.igem.org/parts/bba-k4046704): no deposited sequence-feature entries returned by the public API audit. A downloadable sequence exists, but annotation completeness does not meet this selection gate.
- [BBa_K3084001](https://registry.igem.org/parts/bba-k3084001): no deposited sequence-feature entries returned by the public API audit. A downloadable sequence exists, but annotation completeness does not meet this selection gate.
- [BBa_K1758354](https://registry.igem.org/parts/bba-k1758354): no deposited sequence-feature entries returned by the public API audit. A downloadable sequence exists, but annotation completeness does not meet this selection gate.
- [BBa_K2507008](https://registry.igem.org/parts/bba-k2507008): no deposited sequence-feature entries returned by the public API audit. A downloadable sequence exists, but annotation completeness does not meet this selection gate.
- [BBa_K823025](https://registry.igem.org/parts/bba-k823025): no deposited sequence-feature entries returned by the public API audit. A downloadable sequence exists, but annotation completeness does not meet this selection gate.
- [BBa_K5218004](https://registry.igem.org/parts/bba-k5218004): no deposited sequence-feature entries returned by the public API audit. A downloadable sequence exists, but annotation completeness does not meet this selection gate.

## Source-policy and research gates

1. **T02 source rights/safety:** Verify the exact per-part license label, attribution and local/private-cloud copy scope, original authorship, intended benign role and any ambiguous/source-derived material. A Registry-published/cleared status is supporting evidence, not a substitute for the project's record-level review. Do not commit raw files until that independent gate is documented.
2. **T02 provenance/annotations:** Re-download the part JSON and sequence-features using `scripts/audit_igem_public_parts.py` or the Registry GenBank export, save immutable source hashes, validate DNA and feature coordinate consistency and freeze a source manifest. Store originals only in permitted data paths.
3. **T07 family/reference controls:** Require reviewed per-CDS protein-family identities, independent positive and competing references and linkage exclusions. Source-sharing pairs here (e.g. GfasPurple, TyrBm, YF1/FixJ/GFP) must not be counted as independent retrieval positives.
4. **Pilot-size scope:** The frozen T07 pilot expects approximately 8-12 records. These are **ten supplemental candidates beyond the existing eleven GenBank candidates**, not an automatic extension of the scored pilot to twenty-one. Adding them to scored panels requires a versioned, pre-score protocol/panel revision, rather than modifying a prior freeze post hoc.
5. **Source-type scope:** Some iGEM published parts are expression cassettes or composite inserts rather than whole circular plasmids. They may be suitable for feature-level E1 if documented, but cannot stand in for whole-plasmid backbone-heavy comparisons.
6. **Previous Addgene records:** The nine abandoned Addgene HOLDs are historical, not part of this new selection. The supplied PX551/181490 sequence page is Cas9/AAV-based and remains excluded both for scope mismatch and Addgene login restrictions. Addgene's developer API requires approved token access.

## References

- [iGEM Registry](https://registry.igem.org/) and [2026 GenBank export announcement](https://registry.igem.org/about/whats-new)
- [Researcher-approved GenBank public archive PR #7](https://github.com/khgng026/iri-construct-screening-2026/pull/7)
- [T02 governance](../../docs/week1/T02_data_governance.md) and [T07 frozen protocol](../../docs/week1/T07_benchmark_protocol.md)

**Readiness:** Candidate selection complete; **no new scored E1/E2/E3/E4 run approved yet**. Next task for Codex is individual source admission with a reviewer and a new frozen experimental input panel.
