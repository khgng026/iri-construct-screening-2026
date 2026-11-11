# Evidence packet for delegated review — no execution clearance issued

This packet completes evidence preparation for the current autonomous mission.
It accompanies actual source acquisition, parser QA, model-interface execution,
scoring corrections and reproducible tests. **It is not an approval record.**
The executor has not signed its own ChatGPT scientific clearance.

## Authority, exact scope and evidence binding

PR #9 merged as `7bd2bdb`; the T02/T07 v1.2 delegated workflow is active. The
[researcher's recorded E1–E4 authorization](../oversight/2026-11-10_RESEARCHER_EXPERIMENT_AUTHORIZATION.md)
and [current autonomous mission](../oversight/2026-11-11_AUTONOMOUS_MISSION_SCOPE.json)
authorize the preparation performed here without another routine approval.
The original source, rights, reference and scientific-review gates remain in force.

The [prospective program](../../config/week2/program.proposed.json) defines the
candidate/source subsets, exact configurations, snapshots, resource limits and
storage policy. The [review freeze](../week2/review-freeze.json) binds its actual
file bytes together with code, evidence, source originals and references. The
freeze's canonical `freeze_sha256` is the identifier for a separate review; its
own `scientific_clearance` field explicitly grants nothing. Verification refuses
missing or changed bytes. A fresh acquisition with changed annotations requires
a new packet/review rather than a silently refreshed scientific input.

Read the source artifacts in this order:

| Evidence | Exact versioned artifact |
| --- | --- |
| Source identity/rights/scope/annotations | [source_review.json](../../data/manifests/week2/source_review.json) |
| Preparatory input files and checksums | [inputs.preparatory.json](../../data/manifests/week2/inputs.preparatory.json) |
| Actual feature/unit QA and technical schema binding | [QA report](../week2/qa/main-report.json), [unit provenance](../week2/qa/unit-provenance.json), [schema freeze](../week2/schema-freeze.json) |
| Family definitions, citations, query identity and reference proposals | [reference proposal](../../data/manifests/week2/references/reference_proposal.json), [query specification](../../data/manifests/week2/references/query_specification.json), [query proposals](../../data/manifests/week2/references/query_proposals.json) |
| Copy/parent/ancestor/lineage exclusions | [exclusion evidence](../../data/manifests/week2/references/exclusion_evidence.json) |
| Exact offline model identities and interface measurements | [model freeze](../week2/model-freeze.json), [model runtime](../week2/MODEL_RUNTIME.md) |
| Implementation/results/handoff | [completion report](../week2/COMPLETION_REPORT.md), [research findings](../week2/RESEARCH_FINDINGS.md), [run evidence](../runs/week2-2026-11-11/README.md) |

## Source decisions proposed for review

All 22 main source originals were acquired anonymously: eleven selected GenBank
records, ten selected iGEM parts, and historical pBR322 as a separate control.
Eighteen meet the preparer's technical QA conditions, with 271 retained features
and 18 validated provided protein translations. This count is **not** scientific
admission or independent sample size. Exact source terms, attribution, retrieval
dates, annotation/sequence fingerprints and case-specific caveats are recorded
per source. Local analysis and private-hosting rights proposals are independent
of source/family/scoring approval.

NCBI's captured molecular-data policy is retained with its submitter-IP caveat
and exact-file notice inspection; it is not called CC0 or patent clearance.
Each iGEM part's license UUID resolves to CC-BY-4.0, with original source credits
and attribution retained. Raw source files, original qualifier dumps and exported
proteins have not been committed by this task.

Four NCBI records are HOLD and outside preparatory/experimental eligibility:
KJ541668.1 and KJ541670.1 have toxin-selection annotations and reporter translation
discrepancies; KJ170897.1 and PP098727.1 have publication-context questions under
the project's conservative source policy. The latter does not establish that
the exact reporter vector is inherently hazardous. No missing scientific or
specialist judgment is converted into an approval.

iGEM data remain annotation-only where original translations are absent and
coordinate indexing is unverified. Their feature UUIDs, roles and original raw
intervals are auditable, but model proteins are not inferred. These entries are
not in the initial scored E1/E2/E3 proposal.

A separate [native-export check](../week2/REGISTRY_EXPORT_CHECK.md) inspected two
official Registry GenBank exports. Their DNA matches the original selected parts;
all four CDSs still lack provided translations. The 33 export features reflect
two alternate format snapshots, not extra independent records or a replacement
of the frozen REST denominators. The observed coordinate mapping is evidence for
these snapshots only, not an unreviewed universal API convention.

The separate [MK448012.1 source](../week2/SUPPLEMENTAL_FUSION.md) demonstrates a
genuine 825-aa annotated fusion for technical QA. Its 27 features and full provided
translation validate. It remains outside the original 21-candidate scored scope;
no isolated-target case, scored fusion panel or scope expansion is assumed.

## Prospective stage proposals and unresolved conditions

| Stage | Proposed eligible study | Current decision / remaining evidence |
| --- | --- | --- |
| E1 | Eight exact GenBank sources named in the program, including historical pBR322 as a reproducibility control | **BLOCKED pending actual source-specific and E1 delegated clearance.** Seven are new proposed audit sources; repeating pBR322 adds no independent source. Family references are not required for E1. |
| E2 | AF298787's supported URA3/GFP subset, longest eligible CDS versus all eligible CDSs | **BLOCKED pending actual source/family/panel clearance.** Literal source-longest is unsupported AmpR, not eligible URA3; the full documented-family endpoint stays unresolved. A supported-subset result must not be presented as complete-family recovery. |
| E3 | Four proposed supported query components across U47295, U55762 and AF298787 against the same 17-entry reference catalogue | **BLOCKED pending actual source/family/reference/tier clearance.** Strict tier 3 remains blocked without its independently reviewed fingerprint-bound relatedness policy. At most three parents means null comparative bootstrap intervals. |
| E4 | Actual compatible E2/E3 quality and computation records only | **BLOCKED: no biological measurements.** Fixture timing cannot become scientific cost/quality evidence. |

The fifth query proposal, Renilla, has no independent positive and remains HOLD.
The four proposed supported queries have 15–16 eligible references after
exclusions, exceeding the common total K=5 budget. GFP-like colors are not
different molecular families. Shared reporter/marker fingerprints and vector
lineages remain explicit. Natural source ancestors P42212, P08659, P03962 and
P27652 are excluded from their respective engineered descendants; all components
of AF298787 share its parent-wide exclusions, including unsupported AmpR copies.
Species differences do not establish evolutionary, statistical or pretraining
independence.

The full documented component set and the independently supported subset remain
separate. Unknown marker labels are not silently dropped to improve the reported
denominator. No reference boundary, identity cutoff, family label or endpoint
was tuned using model outcomes; none has been produced.

## Fixed model/configuration and resource proposal

Separate new configuration files preserve the unchanged v1.0/v1.1 historical
files: [E1](../../config/week2/E1.proposed.json),
[E2](../../config/week2/E2.proposed.json),
[E3](../../config/week2/E3.proposed.json),
[E4](../../config/week2/E4.proposed.json). They retain the total five-hit budget,
1022-residue cap, deterministic seed, float32 policy, fixed alignment parameters
and uncertainty rules. ESM-C's mutable historical request is replaced only in
these prospective copies by the exact locally validated snapshot.

| Model | Exact snapshot | Weight SHA-256 |
| --- | --- | --- |
| ESM-2 8M | `c731040fcd8d73dceaa04b0a8e6329b345b0f5df` | `24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189` |
| ESM-C 300M legacy SDK | `7f10b20ae75017b2dbc884070e03434515709a8d` | `323dff9fbf3fef297a74f4f18b6528e6f2e599b0bcf72b6927516804015becea` |

Both execute offline on the fictional fixture; no biological inference occurred.
The ESM-C provider alias points to the retained legacy checkpoint, not a silently
substituted model. Provider terms must remain applicable to local noncommercial
research; the captured third-party notice alone is not complete model-weight
rights evidence. No weights are redistributed.

The prospective executor is bounded to existing CPU resources, sequential models,
two intra-op/three inter-op threads, no scoring-time weight downloads, at most
64 unique proteins and 65,536 processed residues per model/stage, a 1,800-second
wall bound per model/stage, and a 6-GiB process RSS bound. These are reviewable
limits, not measurements or a claim that legacy CLIs enforce OS resource limits.
The eventual executor must enforce/monitor them and block on failure. No external
spending or public release is authorized. Scientific reports remain source-free;
source/protein/qualifier/embedding/weight artifacts remain ignored/private unless
a separately documented source-copy determination permits otherwise.

## Required actual review record

[T02 v1.2](../oversight/T02_V1_2_DELEGATED_REVIEW.md) states:
“An editable audit cannot replace a real approval message.”
[T07 v1.2](../oversight/T07_V1_2_DELEGATED_EVALUATION_SIGNOFF.md) requires a dated
ChatGPT review of actual evidence before cleared stages execute. Codex agents
provided technical cross-review; none is represented as that separate reviewer.

A valid external decision must identify its original review message/reference,
date, **AI-assisted delegated review** role, this exact freeze/program hash,
source decisions and rights scope, accepted family/reference/lineage decisions,
stage-specific `CLEAR_FOR_E1/E2/E3/E4` or HOLD decisions, and resource/storage
bounds. It may clear E1 while withholding later stages. An approval string in
an editable JSON file or successful unit test is insufficient to establish
review authenticity or scientific adequacy.

The available connectors did not provide a scientific-review submission action.
The shared-browser attempt failed at its browser control bridge and native
accessibility provider; no reviewer message was submitted or received. The
executor therefore retains blocked scientific statuses and delivers the completed
technical study. No new execution-authorization file was manufactured. A separate
actual reviewer can now inspect this concrete packet rather than another list
of unperformed preparation tasks.
