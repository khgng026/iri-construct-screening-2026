# Week 2 implementation, source QA and research handoff

The scientifically achievable preprocessing work is implemented and validated on
actual public records. **Week 2's technical pipeline is complete for the prepared
subset; the overall benchmark exit gate remains incomplete because scientific
source/panel admission and the separate ChatGPT execution decision are pending.**
No new E1–E4 experiment is reported as completed.

The work starts from main `12d0ad5`, including merged PRs #6–#10. The active
researcher mission authorizes source preparation, implementation, multiple
agents and bounded local validation while retaining the delegated-review gate;
its [authority record](../oversight/2026-11-11_AUTONOMOUS_MISSION_SCOPE.json)
does not create scientific clearance. Historical sources, results and protocol
versions remain unchanged.

## Original acceptance criteria and task status

The [workflow acceptance extract](WORKFLOW_ACCEPTANCE.json) preserves all relevant
workbook cells and the original XLSX checksum: Weekly To-Dos, Dependency Map,
Technical Skills and MVP Criteria. The workbook itself was not rewritten.
[TASK_STATUS.json](TASK_STATUS.json) distinguishes completed, partially completed
and blocked tasks T01–T28. Approved T01/T02/T03 records provide the proposal/scope
baseline; the original proposal PDF was not available in this checkout, and this
report does not claim to have independently reread it.

| Task / original definition of done | Work and completion evidence | Remaining limit |
| --- | --- | --- |
| T09: versioned manifest, counts by archetype and limitations | [Source ledger](../../data/manifests/week2/source_review.json), exact original retrievals, permissions/attribution proposals, deduplication and four explicit HOLDs; [dataset report](DATASET.md) | Partially complete: 18 preparatory records, zero new scientific admissions; 150–250 expansion deferred under the workbook's quality condition |
| T10: representative inputs parse; invalid inputs fail cleanly | Content-based FASTA/GenBank detection, exact identity/version/checksum checks, controlled malformed/unsupported errors; [parser implementation](../../src/iri_constructs/week2_pipeline.py) and regressions | Completed for supported formats; Registry coordinate conventions remain explicitly unverified |
| T11: standardized pilot feature table | 271 main source features and 27 separate fusion-example features normalized with source UUIDs/ordinals, qualifiers, strand and ordered locations; [feature provenance](qa/feature-provenance.json) | Completed technical preparation; unknown Registry locations stay unknown, source agreement is not biological truth |
| T12: deterministic units and machine-readable provenance | 18 main source-validated proteins plus one separate 825-aa fusion; metadata units retained; [unit provenance](qa/unit-provenance.json), stable parent/feature/unit IDs, 1022-aa interface preflight | Completed technical preparation; no inferred ORFs, missing translations, raw-DNA protein inputs or artificial concatenations |
| T13: passing tests and diagnosable logs | Final clean core environment: 246 tests, 241 passes, five optional skips; structured source/failure events, exact-cache checks, deterministic repeats | Completed; separate model environment validates the optional interfaces; missing dependency skips are disclosed |
| T14: stratified QA checklist and recorded schema changes | [QA checklist](PARSER_QA.md), repeated byte-identical artifacts, [technical schema freeze](schema-freeze.json); original T03 v0.1 unchanged | Technical QA complete; scientific benchmark readiness still awaits admission. The genuine fusion is a separate supplementary source, not a silent panel expansion |

## Dataset and inspected outputs

The core acquisition inspected **22 exact sources**: eleven selected NCBI records,
ten selected iGEM entries and separately reacquired historical pBR322. All original
responses were captured anonymously with retrieval date, source/version, original
byte checksum and a checksum-verified cache. Eighteen were retained for
preparatory QA; four NCBI records remain source-scope HOLD. Their reasons and
the distinction between toxin-selection annotations and publication-context
questions are documented in the ledger. No restricted Addgene access was used.

The preparatory set has 11 multi-CDS, four single-CDS and three backbone-heavy
architecture proposals. Ten iGEM responses retain deposited features but provide
no source protein translations; their location indexing convention remains
unverified. They contribute metadata/abstention coverage rather than fabricated
protein units. Five declared vector/cassette lineage groups now survive from
curation into QA reports.

The separate [pCHX50 fusion example](SUPPLEMENTAL_FUSION.md), MK448012.1, provides
the fourth architecture for technical QA. Its full 825-aa annotated fusion
translation validates. `NeoR_CDS`/`AmpR_CDS` text on `misc_feature` annotations
correctly remains metadata instead of becoming invented CDS proteins. No fusion
subcomponent was exported or admitted as an isolated retrieval target.

A bounded [native Registry export check](REGISTRY_EXPORT_CHECK.md) then inspected
two documented public GenBank exports. Their DNA matches the existing selected
parts, and preparatory QA retained 33 export features. All four CDSs still omit
provided translations, so zero additional proteins were gained. This is stronger
format/coordinate evidence for those two snapshots, not a universal indexing
contract or two new biological records. The main input snapshot is unchanged.

All 18 main records and the separate fusion repeated deterministically. Original
normalized objects, qualifier sidecars and protein FASTA remain ignored/private;
only source-free evidence and fingerprints are committed. The source rights
proposals distinguish local analysis, private hosting and scientific admission.

## Models, scoring and downstream progress

Both fixed models load fully offline in isolated Python 3.11.17 with the original
supported direct pins and a documented missing `httpx` dependency correction.
ESM-2 and ESM-C encoded a fictional 20-residue alphabet fixture twice, confirming
float32 parameters, residue/token alignment, finite representations, mean
pooling, deterministic repeats and rejection of unsupported inputs. Their exact
snapshots and weight hashes are in [model-freeze.json](model-freeze.json).
These checks complete interface engineering, not a newly cleared biological T05
feasibility result or T15/T16 experimental outputs.

The [scoring audit](SCORING_AUDIT.md) fixed three substantive defects:

- Multiple CDS queries from one parent could reach the eight-group bootstrap
  threshold. Conservative connected parent/group/exact-copy clusters now prevent
  that; eight CDSs from one parent remain one cluster with a null interval.
- Legacy numerical fixture/unclassified results could set a scientific flag.
  Computation and scientific evidence are now separately classified.
- E4 could accept contradictory coverage/group denominators. Those records now
  fail validation rather than entering comparisons.

The existing longest-eligible versus all-eligible runner, maximum per-reference
merge, five-total-hit budget, alignment baseline and cost accounting are reused.
T15–T17 have validated software infrastructure; T18's paired scientific run is
blocked. T19 empirical model errors, T20's evidence-driven preprocessing change
and T21's scientific regression cannot be claimed without that run. Current
engineering fixes are not relabeled as a post-experiment T20/T21 study.

## Reference dataset and research interpretation

[Seventeen reviewed UniProt entries](REFERENCE_EVIDENCE.md) support prospective
GFP-like, coleopteran luciferase, OMP-decarboxylase and Renilla-lineage evidence.
Four source-bound query components across three parents have proposed independent
positives and competing references; Renilla has no independent positive and is
held. Exact parent/component copies, known natural ancestors and reviewed lineage
proposals are excluded before any score. Unsupported marker components remain
visible in full-denominator accounting.

Several vector records contain identical reporter proteins, so accession count
does not establish independent sample size. Fluorescent colors do not define
separate molecular families. AF298787 could support a two-family, reviewed-subset
decomposition comparison, but its literal longest source protein is an unsupported
marker, and the full documented-family endpoint remains unresolved. At most three
parent observations are proposed; the eight-cluster comparative interval cannot
be justified. See [research findings](RESEARCH_FINDINGS.md) for the analytical
limits and next experiments.

## Reproducibility, costs and final gates

[Run evidence](../runs/week2-2026-11-11/README.md) contains raw software logs,
environment inventories, original model fixture reports, verification costs and
artifact digests. The clean core environment passed all nine technical checks,
including two repeat runs for each actual-source QA set. A fresh committed-source
archive also passed the full suite and produced exactly matching QA reports,
bound to its archived code commit and archive checksum. Its recorded costs are
technical validation costs, not E4 quality-versus-computation observations.

The offline fixture run measured ESM-2 initialization at 3.12688 s and its first
fixture forward at 0.01831 s; ESM-C initialization at 8.16608 s and its first
fixture forward at 0.33019 s. Highest sampled process RSS was approximately
0.48 GB and 3.12 GB respectively. Models ran sequentially; filesystem cache was
uncontrolled, RSS sampling can miss peaks, and the fixture cannot establish
latency for biological proteins. GPU time and monetary cost remain unmeasured.

The [complete review packet](../coordination/EXPERIMENT_APPROVAL_PACKET.md) and
[prospective program](../../config/week2/program.proposed.json) identify exact
sources, proposed labels/exclusions, model snapshots, configurations and bounded
resources. The [review freeze](review-freeze.json) binds the actual code, evidence
and original input bytes. It cannot grant scientific approval. No new
`EXPERIMENT_RUN_AUTHORIZATION.md` was fabricated.

A separate actual ChatGPT delegated decision remains required by T02/T07 v1.2.
The available connector set did not offer a scientific-review submission tool;
the shared-browser attempt could not access a reviewer because browser controls
failed. No review message was sent or returned. The executor therefore reports
E1–E4 blocked while completing all permitted engineering and preparation.
Additional blockers remain explicit: four source HOLDs, unresolved full marker
families, Renilla independence and strict-tier relatedness policy. External
clearance may authorize eligible stages/subsets; it cannot manufacture those
missing scientific conditions.
