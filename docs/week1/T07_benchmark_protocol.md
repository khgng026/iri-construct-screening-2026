# T07: Real-construct pilot benchmark protocol v1.0

Design frozen on 2026-11-05, before any scored results on this panel. This is an
implemented research design awaiting researcher review and record admission; it
does not approve a source, reference label, pilot release or scientific finding.
The approved T01 scope, T02 policy v1.0 and T03 representation v0.1 remain intact.
The target is 8-12 suitable reviewed records from the T06 candidate inventory;
annotation quality and a valid independent reference panel take priority over n.

The primary hypothesis is that annotation-guided decomposition into supported
protein units improves recovery of documented component-family information over
an appropriate minimally processed representation. There is no training,
fine-tuning, sequence generation, inferred ORF assembly or safety classifier.

## Claims and reference annotations

Three claims must remain separate:

1. **Annotated component presence:** a source feature/qualifier documents a
   component in a parent record. A parser can preserve or lose that annotation;
   it does not independently establish the component's biological identity.
2. **Documented-family retrieval:** a method returns an independent reference
   with the same reviewed annotation family. This is the model endpoint here.
3. **Biological function:** an actual activity/function assay or validated
   functional predictor. Neither parser preservation nor embedding similarity
   establishes this claim, and no such endpoint is evaluated.

Reference truth is the frozen, independently reviewed source annotation, not a
model-generated label or a keyword heuristic. Record accession/version, source
URL, retrieved date, original qualifiers, exact location/strand/segments, source
file SHA256, reviewer/date, reuse scope and supporting family evidence. Keep the
original reviewed bytes and qualifier sidecar locally. Unknown, ambiguous,
conflicting and pending labels remain explicit and outside scored denominators;
they remain in coverage/exclusion accounting. A validated translation only checks
agreement with documented CDS coordinates/frame/table, not biological truth.

All T02 source, licensing, safety and annotation gates must be approved before
pilot source parsing/scoring. Candidate metadata and software test fixtures do
not confer approval. No real record has been admitted by this design document.
Family labels must be `known` and separately `approved`, with evidence and an
independence group linking related vectors, homologs or shared source products.

## Valid input conditions and archetypes

ESM-2, ESM-C and ESM-3 sequence interfaces accept proteins. Plasmid DNA is never
passed to them. A noncoding/backbone DNA feature remains metadata; lack of a
protein interface is a coverage limitation, not a negative functional prediction.

| Archetype | Minimally processed protein baseline | Annotation-guided condition | Interpretation |
| --- | --- | --- | --- |
| Single CDS | The same source-annotated complete CDS translation | The identical validated CDS unit | Identity/control comparison; no improvement expected from unchanged input |
| Multiple CDSs | Longest eligible annotated CDS, stable feature-ID tie break | All eligible source CDS units; merge scores per reference, then return only K total hits per parent | Tests information retained by a single-CDS reduction, not a protein model on intact DNA |
| Tagged/fusion | The complete documented fusion protein | Isolated annotated target; contextual residues at reviewed target bounds from the same full-fusion pass | Known-boundary characterization control; boundaries are annotation/oracle information, not automatic discovery |
| Backbone-heavy | Longest eligible annotated CDS, when present | Reviewed CDS units plus retained unsupported noncoding metadata | Report substantial unsupported parent content; no extrapolation to backbone function |

Do not concatenate separate CDS proteins to manufacture a fusion condition.
Mixed-strand, joined and circular CDS require preserved ordered parent segments;
missing/fuzzy/remote boundaries and exceptional/pseudo CDS are withheld for review.
FASTA DNA with no annotations has no observable feature-preservation denominator
and no invented CDS. A malformed or unsupported source remains a reported failure.

The executable prioritized B experiment implements **matched annotated fusion
target queries**, not the prospective multi-CDS parent-level comparison above.
The multi-CDS baseline/score-merge specification is frozen for a later bounded
implementation; its availability as a design must not be reported as an executed
construct benchmark. For B, virtual `case/<case_id>` query identities allow the
same target in different documented contexts without overwriting outcomes. The
same reviewed family expectation is used for all three representations.
Evaluation parent/archetype refer to the documented fusion construct; the
isolated source remains explicit lineage. Cases sharing any target/fusion parent
or reviewed target/fusion independence group are conservatively joined into one
bootstrap cluster rather than counted as independent targets.

## A priori experiments

**A - Annotation preservation.** Compare original reviewed feature semantics to
the normalized T03 objects and original-qualifier audit. Primary checks are
feature retention, source-coordinate/strand agreement and parent-child links.
Record qualifier preservation, source CDS usable-translation coverage, exported
translation agreement, missing translations, withheld CDS, ambiguous/unknown
features and failures by archetype. Retain the denominator for each rate. A
zero-feature FASTA record is `no_observable_annotations`, never perfect
preservation. Disagreement is a parser/source-agreement issue to inspect, not
proof the source is wrong. Run before any pretrained inference when admitted
records permit it.

**B - Whole versus decomposed fusion representations.** Use the same immutable
ESM-2 checkpoint, float32 policy, final adapter output and residue mean pooling
for isolated target, complete fusion, and the target region sliced after the
complete-fusion forward pass. Require exact reviewed amino-acid boundaries and
exact isolated/span agreement. Report within-model pooled cosine and contextual
residue differences descriptively; the primary retrieval endpoint is matched
documented-family recall@5 on the same eligible reference identities. Similarity
alone is not a function prediction or evidence of characterization improvement.
All targets/fusions' source parents are excluded from the common reference set;
this prevents target-only provenance from allowing a reference from the fusion's
own construct. Optional independent references are admitted separately, never
derived from the query's own annotations/sequences as exact lookup targets.
Exact sequence copies of any isolated target **or complete fusion**, including
copies stored under different parents, are removed from every B condition's
common reference set. Reviewer-linked target/fusion groups are also excluded;
all removed reference IDs and reasons are saved before panel freezing.

**C - Embeddings versus alignment.** Compare the same admitted query/reference
IDs and five total returns. Raw local alignment is the primary conventional
baseline: Biopython PairwiseAligner, BLOSUM62, gap-open -10.0, gap-extend -0.5.
Reference-self-normalized local score is a declared sensitivity analysis; it is
not a calibrated detector or an endpoint selected after seeing this pilot.
Compare within-method ranked family hits, not raw alignment scores to embedding
cosines. Duplicate sequences, parent sharing and linked independence groups are
excluded. This panel evaluates independent-family retrieval; any separate
exact-reference diagnostic must be labeled separately and cannot substantiate
unseen-protein generalization. No embedding advantage is presumed.

**D - Model cost sensitivity.** Default to ESM-2 8M and ESM-C 300M, with the same
reviewed panel, fixed cap and endpoint. ESM-3 is optional only after explicit
protocol revision, model access and resource validation. Freeze model revisions,
weight digests, library versions, precision, pooling and token alignment. Measure
loading and first-pass inference/pooling/transfer separately; run models
sequentially. Compare coverage, matched annotation-family recovery, latency and
memory costs. Never compare cosine values between different embedding spaces.
No model winner follows from an absent or insufficient independent panel.

## Frozen panels, budgets and leakage controls

Machine-readable designs are [A](../../config/pilot/experiment_A.json),
[B](../../config/pilot/experiment_B.json), [C](../../config/pilot/experiment_C.json)
and [D](../../config/pilot/experiment_D.json). The design status is
`design_frozen_pending_admission`; it means choices precede pilot scoring, not
that a real dataset is frozen or approved. The shared standard-library
[scorer](../../src/iri_constructs/pilot_scoring.py) uses existing reviewed-label
and leave-construct-out eligibility utilities rather than a second label system.

Before inspecting scores, save the exact input configuration/review-manifest
digests, admitted source checksums and model-input IDs. `freeze_panel` binds the
protocol's canonical SHA256, query/reference identities, parent/feature identity,
sequence/source digests, reviewed family expectations, independence groups,
archetypes and per-query eligible references. Its own SHA256 is persisted before
method scoring. Empty or insufficient panels have blocked statuses, not success.
Do not alter family labels, reference membership, expected components, budgets
or thresholds after inspecting scores. Amend version/config hashes and restart
the scored analysis if a substantive review correction is necessary.

Query and reference role lists may overlap for leave-construct-out evaluation,
but every query excludes itself, exact sequence copies, its source parent and
reviewer-linked independence group. B additionally excludes every target/fusion
parent represented by its cases. A positive independent reference for every
documented expected family and at least one competing family are required.
Source-related or highly homologous records may still remain confounded; review
those groups conservatively and disclose possible pretrained-data overlap.

K=5 is a **total** returned-hit budget for a matched target query or a parent
construct ranking, never five per CDS/window. If fewer than five eligible
references exist, every method returns the same available number. Independent
panels with at most K eligible references return every reference and cannot
distinguish method quality. These query rows are flagged `budget_saturated`,
reported descriptively, and excluded from method-advantage interval summaries.
A run with only such scored rows returns a descriptive-only status, not a
successful comparative result. Independent
per-CDS scores in the prospective multi-CDS condition must be merged by maximum
score per reference before truncating one parent list. Record every method's
number of forward passes/alignment calls, total processed residues, maximum
input length, encoded/cached unique sequences and total returned hits. Equal
return budgets do not imply equal compute budgets; report both explicitly.
Automatic windows are outside the prioritized implementation. Any future window
analysis needs settings frozen on independent development groups, a matched
whole baseline and the same total merged hit budget. Known boundaries remain a
separate oracle/annotation control rather than an automatic segmentation gain.

Rank ties break by stable reference ID. A complete fixed eligible reference
score row is required for each scored query. Missing/nonfinite reference scores
cause abstention; a method may not silently shrink its reference panel. Reject
duplicate, linked, out-of-panel or over-budget ranked hits. Compare models both
on original requested-panel coverage and exactly common successful query AND
reference IDs. Store panel identities, rather than merely matching sample counts.

## Metrics, uncertainty and failures

Report requested/admitted/eligible/scored counts, review exclusions, format/model
coverage, abstentions and failure reasons before agreement metrics. At protein
target level, report same-family hit@5. For a frozen set of documented expected
families, report the fraction recovered and all-families recovered; repeated
components sharing one family are not identifiable as separate biological
components using this endpoint. These are family-observation endpoints, never
direct component identity or actual-function accuracy. The scorer names them
`documented_family_recall_at_k` and `all_documented_families_recovered`.

Report micro family-observation recall, macro target-query recall, macro
construct recall and per-archetype counts/recall. Summaries across several
protein targets in a construct state the total returned hits; they are target
query summaries, not a claim of one-parent-budget performance. Missing supported
annotations are excluded from agreement denominators and visible in annotation
coverage. An empty denominator produces null, never an invented zero or 100%.

Paired method contrasts use identical scored query IDs and save omitted/unpaired
IDs. Primary uncertainty is a seeded percentile cluster bootstrap (2,000
replicates, seed 20261105) of right-minus-left query proportions, averaging first
within reviewed independence groups and then weighting groups equally. Repeated
CDSs, fusion targets and backbone-linked records are not independent samples.
With fewer than eight independent groups, report paired observations and point
differences but a null confidence interval and `insufficient_independent_clusters`.
Even with eight or more groups, the pilot CI is exploratory and does not justify
population-level claims. No p-value or multiple-comparison winner is specified.

Record hardware, process/thread settings, elapsed time, process RSS sampling
policy, allocator GPU peaks where available, model-load failures, inference
failures, cap/alphabet/tokenization failures and source-validation errors. First
passes include cold effects; cached sequence hits are separate. No unavailable
dependency, skipped experiment, rejected source or zero-score panel is silently
converted into successful experimental evidence.

## Commands and evidence categories

The annotation audit can be exercised independently after source admission:

```powershell
python scripts/prepare_pilot.py --inputs data/processed/reviewed-inputs.json --output results/pilot/A/run01
python scripts/run_pilot_experiment.py --protocol config/pilot/experiment_A.json --inputs data/processed/reviewed-inputs.json --output results/pilot/A/scored01
python scripts/run_pilot_experiment.py --protocol config/pilot/experiment_B.json --inputs data/processed/reviewed-inputs.json --output results/pilot/B/run01
python scripts/run_pilot_experiment.py --protocol config/pilot/experiment_C.json --inputs data/processed/reviewed-inputs.json --output results/pilot/C/run01
python scripts/run_pilot_experiment.py --protocol config/pilot/experiment_D.json --inputs data/processed/reviewed-inputs.json --output results/pilot/D/run01
```

Use the documented optional parsing/model environments from
[experiments/README.md](../../experiments/README.md). Input and review paths are
local, raw sequences and outputs remain ignored, and no command grants approval.
Shared scorer tests are numerical/fictional software checks. A source-agreement
audit on admitted real records is a parsing experiment. A completed pretrained
forward pass is model execution. A sufficiently independent frozen, reviewed
reference comparison is interpretable pilot evidence. None implies the next.
Currently the real candidate/review queue is pending, so this protocol claims
no real pilot metric, pretrained comparison, result figure or improvement.

## Reuse and PR #1 boundary

Existing `main` supplies T02/T03 validation, local record normalization, fixed
model adapters, fusion context extraction and within-model retrieval utilities.
This contribution adds the prospective protocol, frozen shared panel/scoring
contract, leakage checks and paired uncertainty. No PR #1 source code or result
tables were copied into this scorer, and that branch was only read.

PR #1 was inspected at
[`7b15f02`](https://github.com/khgng026/iri-construct-screening-2026/tree/7b15f02ac929cd06a8dee76ccec0541641844ccf):
[START_HERE](https://github.com/khgng026/iri-construct-screening-2026/blob/7b15f02ac929cd06a8dee76ccec0541641844ccf/START_HERE.md),
[REPRODUCE](https://github.com/khgng026/iri-construct-screening-2026/blob/7b15f02ac929cd06a8dee76ccec0541641844ccf/REPRODUCE.md),
and the [research report](https://github.com/khgng026/iri-construct-screening-2026/blob/7b15f02ac929cd06a8dee76ccec0541641844ccf/docs/research/2026-10-31/REPORT.md).
Its artificial known-protein combinations show that global pooling can miss
components, alignment is a strong exact-reference comparator, fixed windows can
improve embedding retrieval, and known boundaries/extra returned hits create
oracle or budget advantages. Those lessons motivate this design's independent
references, matched hit budgets and distinct controls. They are not real DNA
construct/function validation, and held-out combinations sharing proteins do
not establish unseen-protein generalization. No historical measurement is
transferred into the present pilot results.
