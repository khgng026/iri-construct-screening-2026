# T07 v1.1: approved evaluation extensions

This versioned addendum extends the unchanged
[T07 v1.0 design](T07_benchmark_protocol.md). The researcher approved v1.0 in
principle and authorized E1-E4 on 2026-11-05; that decision does not retrospectively
change its frozen text/configuration or claim any scientific outcome. T02 and
T03 v0.1 remain unchanged. Machine-readable designs are in `config/extensions/`;
each run freezes the chosen design before data admission and scoring.

Authority is the researcher instruction and the read-only oversight records at
[`4fd8256`](https://github.com/khgng026/iri-construct-screening-2026/tree/4fd825691f5441f88fd0419996689f8bb76a98d9/docs/oversight):
`2026-11-05_DECISIONS.md` and `EVALUATION_EXTENSIONS_V1_1.md`. Implementation
details below are deterministic engineering choices under that authority, not
new biological source approvals or post-result tuning.

## What is approved, pending, or historical

The approved ten-identity shortlist is pBR322; pGL3 Basic Vector;
Luciferase-pcDNA3; pcDNA-RLuc8; pYIC; pGL4Luc-RLuc; pGreen2x35SLUC;
pEGFP-LC3; Nano-lantern/pcDNA3; and pcDNA3-Renilla-FLAG-firefly.
pUC19c remains an optional annotation-limited control only after its own T02
review. Preserve the original 20-candidate inventory and prior queue as discovery
history; the new selection manifest records the explicit change in priority.

No sequence file, source permission, annotation, family label, relatedness claim
or scoring panel is approved by candidate selection. Source-byte checksums,
reviewer/date, safety/reuse scope, source annotations and independent references
remain admission requirements. Missing evidence produces HOLD/abstention rather
than inferred features or approval. Software fixtures, model execution and
interpretable scientific evidence remain distinct.

## E1: annotation quality and appropriate abstention

Classify admitted records using source annotations and caveats, never model
outcomes: `well_annotated`, `partially_annotated`, or `annotation_limited`.
Reuse source-versus-T03 preservation audits for exact coordinates, strands,
ordered multipart spans, original qualifiers/translations and parent-child IDs.
Joined/reverse annotations are not intrinsically errors; unsupported/fuzzy,
missing or contradictory evidence is retained and may withhold model input.

For every record and archetype report source CDS/features, observable evidence,
validated emitted units, withheld/unknown/ambiguous units by reason, unsupported
input counts, invention checks and justified abstentions. Retain every rate's
numerator and denominator. Null denominators remain null. A faithful source-only
record is annotation-limited and has no observable CDS evidence; it must not
become perfect feature agreement, a parser failure, or a fabricated CDS.
Source-agreement checks do not validate the source's biological correctness.

An admitted pUC19c control may test this limitation, but no source is currently
admitted. Automated tests use clearly fictional records; no real annotations
are deliberately corrupted to manufacture scientific sensitivity findings.

## E2: multi-CDS parent-level evaluation

Compare the longest eligible documented CDS translation against all eligible,
reviewed CDS translations from the same parent. Choose the longest by amino-acid
length, with ascending stable unit ID breaking ties. Preserve the number and
reasons for unsupported/missing translations or unreviewed family labels.
No model receives raw DNA, inferred ORFs or a concatenation of separate CDSs.

Freeze documented expected family sets per parent before scoring. Source
component multiplicity and distinct family sets are different denominators:
two CDSs with the same family are not two independent recovered families.
Independently supported subsets and full documented sets are reported separately;
unsupported families are not silently dropped from the full denominator.

Every condition uses one common eligible independent reference set, excluding
the source parent, all documented component sequence copies and reviewed linkage
groups. For the all-CDS condition, take the maximum score for each reference over
its eligible CDS units, deduplicate reference IDs, then truncate to **five total
hits per parent**. The single-CDS condition uses that same panel and budget.
At least an independent positive family and a competing family must remain;
missing evidence/scores yield explicit blocked subsets or abstention.

Report unit extraction/representation coverage, distinct-family recall and
complete-family recovery, returned hits, calls, residues, latency, memory and
failures. Keep matched per-unit controls separate: representing more CDS units
can recover more documented information without improving accuracy on an
individual protein. Reusing identical unit scores in both conditions is a
control, not an extra independent observation. Shared caches/load work must be
identified, not charged as two independently measured cold experiments.

## E3: independently reviewed reference difficulty

All scored tiers retain **every** v1.0 source/parent/exact-copy/independence-group
guard. The first tier is the original T07 guard, not a relaxed self-retrieval
panel. Tier 2 adds predeclared reviewed vector/linkage exclusions when available;
otherwise it is transparently redundant with tier 1. Tier 3 requires an
independently justified, reviewed, pre-score sequence-relatedness exclusion
policy, bound to exact admitted sequence fingerprints. No identity cutoff is
chosen automatically or after seeing results. Missing/invalid strict policy
blocks tier 3 without fabricating results or weakening supported tiers.

Freeze all tier identities, allowed references, removed IDs/reasons, rules and
policy/source/protocol hashes before producing any method result. Every method
uses identical query IDs, eligible reference IDs and five total hits within a
tier. Exact ESM model snapshots must resolve before the first model forward.
Cold loading, full catalogue encoding, scoring and reusable tier work are disclosed
separately; costs shared across tiers are not additive per-tier measurements.

Compare ESM-2, ESM-C where available, and raw local BLOSUM62 alignment. The
v1.0 reference-self-normalized alignment remains a labeled sensitivity, not a
new primary endpoint selected after outcomes. Report family retrieval and
wrong reviewed-family labels among returned hits with explicit denominators;
these are not biological-function or operational screening false-positive rates.

Report catalogue sizes, independent positives/competitors, removals, attrition,
coverage/abstentions and matched paired differences. No independent positive or
competing family means null/abstention. Cluster bootstrap intervals require the
predeclared number of independent groups; saturated budgets that return every
candidate remain descriptive only. Changed catalogue composition, common-family
coverage and unknown pretraining overlap limit generalization claims even with
strict exclusions. No unseen-protein or biological-function claim follows.

## E4: quality and measured computation

Aggregate only actual measurements with exact comparable panel identities,
protocol/task/budget, source configuration, precision, cache policy and matching
hardware/software conditions. Preserve incompatible records as separate rows
with reasons rather than pooling them into a model ranking. Unclassified or
fictional records cannot be promoted into scientific evidence by a cost plot.

The common `cost-quality-record-1` contract records evidence kind, admission
evidence, identities, coverage, documented-family/complete-family recovery,
requested/scored counts and independent groups. Costs include separate cold-load
and inference/pooling/transfer time, model/alignment calls, unique encoded inputs,
processed residues, process CPU and GPU time when actually measurable, sampled
process RSS, allocator GPU peaks and cache hits. Unknown measurements remain null;
wall time is not relabeled GPU execution time. Shared work is explicitly marked.

Reports contain tables and may produce a simple standalone cost-versus-quality
plot only for actual admitted, comparable observations. Optional plotting
dependencies are loaded lazily; their absence is an explicit plotting limitation,
not a fabricated figure. No module automatically selects a winning model.
Until reviewed measurements exist, the infrastructure is tested using clearly
labeled fictional numerical fixtures and scientific reporting remains pending.

## Reproduction, status, and changes from v1.0

```powershell
python scripts/run_evaluation_extension.py --protocol config/extensions/E1.json --inputs config/extensions/inputs.pending.json --output results/extensions/E1/pending01
python scripts/run_evaluation_extension.py --protocol config/extensions/E2.json --inputs config/extensions/inputs.pending.json --output results/extensions/E2/pending01
python scripts/run_evaluation_extension.py --protocol config/extensions/E3.json --inputs config/extensions/inputs.pending.json --output results/extensions/E3/pending01
python scripts/run_evaluation_extension.py --protocol config/extensions/E4.json --output results/extensions/E4/pending01
```

Use a new output directory each time. These pending commands must return
blocked/no-input statuses rather than successful science; no download or source
approval occurs. After independent T02 admission, use a separate local reviewed
configuration and the existing optional parser/model environments. Inference
is offline by default and sequential. E4 consumes explicit measurement-record
files; it does not infer comparable measurements from absent or ambiguous logs.

V1.1 adds quality strata/abstention, executable parent score merging, controlled
reference tiers and strict cost reporting. It does not alter T02, T03 v0.1,
v1.0 files/hashes or claims; it does not train models, generate biological
sequences, reconstruct fragments or implement operational screening decisions.
Larger-model experiments and expansion toward 150-250 reviewed benign constructs
are conditional on a useful pilot, annotation quality, independent references,
review effort and resources, not immediate acceptance targets.

The human-readable milestone handoff is
[CODEX_STATUS.md](../coordination/CODEX_STATUS.md). Specific remaining source and
policy decisions are in [DECISIONS_NEEDED.md](../coordination/DECISIONS_NEEDED.md).
That handoff requires requested human review; ChatGPT is not automatically
monitoring it. While PR #2 is unmerged, the stacked PR targets its branch; after
PR #2 is merged, the extension PR targets main. No branch is merged by this
implementation.
