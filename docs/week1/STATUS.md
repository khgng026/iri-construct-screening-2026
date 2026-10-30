# Week 1 | Researcher approvals and task status

**As of 2026-10-29:** The researcher explicitly approved **T01–T04**. Status **Done** refers to acceptance of each Week 1 deliverable, not completion of downstream benchmark criteria or biological validation.

| ID | Task | Status | Decision and remaining boundaries |
| --- | --- | --- | --- |
| T01 | MVP question, scope, acceptance criteria | **Done — approved 2026-10-29** | Scope/criteria accepted as planning baseline v1.0. Model choice, paired benchmark, and results remain later tasks. |
| T02 | Public data and safety rules | **Done — approved 2026-10-28; reaffirmed 2026-10-29** | Written governance policy v1.0 accepted. No specific record, source license, or dataset is approved by this sign-off. |
| T03 | Archetypes and normalized schema | **Done — approved 2026-10-29** | Four archetypes, schema v0.1, and three fictional examples accepted. Real annotated pilot validation and future version bumps remain pending. |
| T04 | Repository and reproducible environment | **Done — approved 2026-10-29** | Private GitHub repo, configuration/logging, smoke check and 11 unit tests accepted; GitHub Actions CI previously passed. |

**T05 — In progress (not approved as done).** The researcher reported running ESM-2 `facebook/esm2_t6_8M_UR50D` locally on public GFP UniProt `P42212` (238 amino acids), obtaining an output shape `[1, 240, 320]`. This is evidence of a successful local representation-generation probe, **not yet a functional-label evaluation or full feasibility sign-off**. Runtime, reproducibility, interface limits and valid benchmark outputs remain to be documented; the experiment script has not been independently verified in this repository.

**Evidence and source of record:** [four-week workflow](../../planning/4-week-workflow.xlsx), [T01](T01_scope_and_acceptance.md), [T02](T02_data_governance.md), [T03](T03_archetypes_and_schema.md), [T04](T04_repository_setup.md), [private repository](https://github.com/khgng026/iri-construct-screening-2026).

**Scope guardrail:** User acceptance of the written work does not validate model predictions, biological records, source permissions, or subsequent MVP criteria. Continue to apply the approved T02 checklist record by record.
