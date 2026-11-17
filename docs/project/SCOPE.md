# Toolkit-first scope and acceptance

Effective November 16, 2026 (America/New_York). Authority: the researcher asked to retain the current preprocessing approach, prioritize a tractable toolkit, replace unnecessary review bottlenecks, and then instructed: “perform the key action items, then follow up with a plan for a large stretch of agent work for progressing towards working pre-processing toolkit that can help answer the key research question”.

## Deliverable

A researcher can install `iri-constructs`, process an annotated GenBank file or batch, and obtain traceable normalized annotations, validated protein FASTA, feature/unit tables and actionable diagnostics without model weights, benchmark labels, approval CSVs or historical project caches.

The research question remains whether preprocessing annotated benign constructs into traceable model-compatible units improves characterization against a valid minimally processed baseline. Toolkit readiness and scientific evidence of improvement are separate outcomes. A usable toolkit may be complete while the general benefit remains unresolved.

## Current approach

Reuse source annotations. Normalize construct/feature/unit relationships, exact strand and ordered segments. Export existing CDS translations only after validating them against the source DNA coordinates, frame and translation table. Missing, conflicting, fuzzy or unsupported evidence stays explicit; other valid units continue. Keep a full annotated fusion as one protein. Preserve noncoding features as metadata. Do not require independent family truth to extract an annotated protein.

## Release acceptance

- R1: wheel installs in a fresh supported Python environment; installed API and CLI work outside the checkout without `.git`, legacy project files or model packages.
- R2: `iri-constructs preprocess input.gb --out output/` and a documented batch form produce the contract outputs. These commands were the target interface at the reset and are implemented as of v0.1.0 (see [release assessment](RELEASE_ASSESSMENT.md)).
- R3: original feature semantics, qualifiers, source hash, parent links and ordered locations are preserved; feature-type normalization retains the original source type.
- R4: known translations are validated, unsupported units carry stable reason codes, and supported units are exported independently of benchmark label availability and model length limits.
- R5: single-record and multi-record files, reverse strands, joined locations, missing translations, ambiguous annotations and malformed inputs have meaningful automated coverage. Partial batches retain successful records.
- R6: repeated processing of the same bytes/configuration gives deterministic biological outputs and IDs; variable timestamps/resources are separate run metadata. Duplicate identities never silently overwrite outputs.
- R7: a compact example collection, source/provenance notes, offline regression tests, clean-install demo and concise limitations report are delivered.

Acceptance is based on covered behavior, not a target of 150–250 constructs or a minimum number of independent families. Synthetic fixtures cover edge cases; appropriate public records provide integration examples. No statistical model-performance claim follows from these checks.

## Existing E1–E4 work

E1 preservation/translation checks become core QA. E2 longest/all-CDS aggregation is an optional research demonstration. E3 method comparisons/strict reference tiers and E4 matched cost-quality analysis are deferred extensions. Simple parser runtime/failure reporting remains useful. None of E2–E4, a model download, independent family annotation, a fusion retrieval benchmark or a positive scientific outcome is a release gate.

## Deferred

Gene calling, inferred missing translations, protein/domain splitting, arbitrary windows, training/fine-tuning, ESM-C/ESM-3 expansion, large dataset acquisition, production deployment and public biological-data publication. A later scope decision may add them after the toolkit is delivered.

Historical scope and workbook criteria are preserved in place; this document supersedes their active acceptance/dependency rules. [Plan](PLAN.md) | [Contract](TOOLKIT_CONTRACT.md) | [Autonomy](AUTONOMY.md) | [Status](STATUS.md)
