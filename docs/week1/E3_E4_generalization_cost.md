# E3/E4: independent references and measured computational cost

These implementations extend [T07 v1.1](T07_v1_1_evaluation_addendum.md).
They have software validation using fictional numerical adapters and metadata,
with no admitted biological panel, pretrained inference or scientific result.

## E3: freeze before evaluation

`reference_tiers.freeze_reference_tiers` fixes query/reference roles, reviewed
family expectations, sequence/source fingerprints and exclusions before scoring.
Every tier retains the original T07 exact-copy, same-parent and reviewed-group
guards; tier 1 never reintroduces a source-parent match merely to make it easier.

| Tier | Additional rule | Missing evidence |
| --- | --- | --- |
| `tier1_t07_invariants` | Original T07 exact-copy, parent and reviewed-group exclusions | Query abstains without independent positives for every expected family and a competing family |
| `tier2_additional_linkage` | Reviewed vector/linkage pair exclusions beyond those existing guards | Explicitly redundant with tier 1 when no additional policy is supplied |
| `tier3_reviewed_relatedness` | Independently justified, predeclared sequence-relatedness pair exclusions | Blocked; no identity threshold is invented |

Local reviewed input JSON can name `linkage_policy` and `relatedness_policy`
files relative to itself. Each policy is a `reviewed-relatedness-policy-1` object
with `policy_kind: reviewed_pair_exclusions`, `policy_scope: vector_linkage` or
`sequence_relatedness`, an approved `review_status`, `declared_before_results: true`,
and nonempty `policy_id`, `reviewer`, ISO `reviewed_on` date and `evidence`.
`excluded_pairs` contains unique query/reference pairs from the frozen roles;
each has `query_id`, `reference_id`, `query_sequence_sha256`,
`reference_sequence_sha256`, `reason` and `evidence`. Fingerprints must match.
This validates the declaration and binding, not the scientific truth of a review.
An explicit empty list is meaningful only when the human review justifies it.

`run_evaluation_extension.py` freezes every tier before local BLOSUM62 alignment
or model scores. ESM-2 and ESM-C run sequentially; each actual revision/weight
identity is frozen before its first forward. Raw and reference-self-normalized
alignment reuse one physical alignment execution. Embeddings are computed once
on the frozen source catalogue, then reused across tiers. All methods receive
identical eligible references within a tier and return at most five total hits.

Missing scores cause abstention. An unavailable model leaves empty common scored
support rather than silently becoming a smaller successful comparison. The
underlying v1.0 alignment helper conservatively drops an entire query score row
when a base-panel reference fails; a stricter tier does not recover that row.
Report this attrition, rather than comparing it to a model on a different support.

Per-tier outputs retain removed-reference reasons/counts, panel sizes, positive
and competing counts, coverage, abstentions, documented-family retrieval and
common scored query IDs. Paired differences use matched independent clusters;
saturated panels are descriptive and too few clusters produce null CI bounds.
False-family hit fractions describe disagreement with reviewed family labels,
not biological-function false positives. No unseen-protein generalization claim
is justified by these software tests or by pretrained-model retrieval alone.

## E4: consume recorded observations

E2 and E3 write `cost-quality-records.json` with a `records` array of normalized
`cost-quality-record-1` observations. E4 never starts inference or approves data.

```powershell
python scripts/run_evaluation_extension.py --protocol config/extensions/E3.json --inputs config/extensions/inputs.pending.json --output results/extensions/E3/pending01
python scripts/run_evaluation_extension.py --protocol config/extensions/E4.json --measurement-records results/extensions/E3/pending01/cost-quality-records.json --output results/extensions/E4/pending01
```

Both commands intentionally return nonzero while the panel/measurements remain
pending. Use new output directories. The standalone `report_cost_quality.py`
also accepts explicit record batches via `--records`; its legacy A-D adapter
preserves missing context as null/unclassified and cannot infer admission.

A direct comparison requires actual admitted measurements and matching panel,
protocol, input/review manifest, total budget/scope, task/tier, exact scored IDs,
runtime/cache policies, cache state, precision, hardware and software context.
Model configuration identities are retained. Ambiguous reuse of one method label
for different model configurations is excluded before pairs, groups or plots.
At least eight independent groups and measured quality denominators are required
by the checked-in design; catalogue-saturated panels cannot establish a gain.
Rejected comparisons retain all descriptive rows and reasons.

Resource fields distinguish cold model initialization, inference/pooling wall
time, total wall time, successful model calls, alignment calls, unique encodings,
amino acids processed, RSS/CUDA peaks and cache hits. Unobserved CPU/GPU time or
memory stays null; wall time is never relabeled GPU compute time. Successful calls
are distinct from attempted/failing calls retained in execution events. Weight
files may already be cached; a fresh model instance does not imply cold disk/OS
caches. E2 baseline costs are the measured reference-plus-baseline prefix of one
shared execution, not an independent cold rerun. Shared execution IDs prevent
costs from being summed across conditions, alignment normalizations or tiers.

Outputs are JSON, a CSV table with explicit `NULL` cells and, when permitted,
one standalone SVG of measured inference/pooling seconds versus family recall.
Production plotting requires one comparable admitted group with at least two
measured points. Multiple groups require selection through the plotting API.
Matplotlib is an optional dependency loaded only when a plot is eligible; it is
absent in this session, so its fixture rendering test is skipped. The separate
test-only plotting option visibly labels fictional data. No figure is produced
for the current pending configuration and no automatic winner is selected.

## Before real execution

Supply exact permitted local files with complete T02 admission evidence and
independently reviewed family/reference groups. Declare `evidence_kind` explicitly
as `admitted_research` only for those reviewed data, provide the tier-3 policy,
and freeze the intended model/runtime configuration. Validate model dependencies,
access and memory separately, then run fixed checkpoints sequentially and preserve
all failures/denominators. E4 can summarize the resulting records only when its
strict comparability checks pass. See [DECISIONS_NEEDED.md](../coordination/DECISIONS_NEEDED.md).
