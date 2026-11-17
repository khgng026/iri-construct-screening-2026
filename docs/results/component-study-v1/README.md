# Component-coverage study v1 — results

Run 2026-11-17 under [the specification](../../../experiments/component_study_v1/SPEC.md),
which was committed with the [panel freeze](panel.freeze.json) (commit `17ea54b`)
before any score was computed. Agent self-review only. Everything here is
descriptive: the 127 multi-unit parents fall into **6** dependence clusters
(parents sharing a protein sequence), below the eight the specification requires
before any population statement. Counts are over *distinct* unit sets unless
noted, so repeated vectors do not inflate them.

## What ran

| Item | Value |
| --- | --- |
| Inputs | Seeded sample of 200 NCBI vector records ([acquisition](acquisition.json)); processed by `iri-constructs preprocess` |
| Model | `facebook/esm2_t6_8M_UR50D` @ `c731040f…5df`, weight SHA-256 `24c5fa47…4189` (same as the historical runs); mean-pooled, cosine, CPU float32 |
| Comparator | BLOSUM62 local alignment, raw score (the project's existing scorer) |
| Panel | 28 scorable families, 134 reviewed UniProt references, none ≥ 0.90 identical to a queried component or its anchor |
| Queries | 476 eligible units = 248 distinct sequences; 243 units (78 distinct sequences) carry a supported family label |
| Budget | 5 hits per parent in every condition |
| Work | ESM-2: 382 encoder calls (121,754 residues), 35 s, 556 MB peak RSS. Alignment: 33,232 scores, 25 s, 50 MB. Repeat ESM-2 run gave the identical score-matrix hash |

Records, reference sequences, score matrices and weights are in ignored private
storage. [results.json](results.json) holds every count below;
[scoring.esm2.json](scoring.esm2.json) and [scoring.alignment.json](scoring.alignment.json)
record model identity and resources.

## 1. Coverage: what the protein interface can and cannot see

| Of 200 sampled constructs | Count |
| --- | --- |
| Source features | 2,459 (1,931 noncoding, kept as metadata) |
| Annotated CDSs | 528 |
| Exported (translation validated) | 488 (92.4%) |
| Withheld | 36 fuzzy/partial boundary, 4 exceptional annotation |
| Exported but longer than the model limit | 12 (kept, flagged, not scored) |
| Eligible query units | 476 |
| With a supported family label | 243 (51% of eligible) — 230 exported units are `unknown`, 3 `fragment`, 7 labeled in a family with too few references |
| Constructs with no exported CDS | 13 (11 have no CDS at all) |
| Constructs with no supported component | 45 |
| Share of construct length inside exported CDSs | mean 34%, median 31% |

Two thirds of an average construct is never shown to the protein model, and half
of the proteins that are shown have no family label to score against. Neither is
dropped from the tables; neither is evidence for or against the method.

## 2. Per-unit retrieval (the model, not the preprocessing)

Identical in every condition by construction.

| Distinct supported sequences | ESM-2 8M hit@1 / hit@5 | Alignment hit@1 / hit@5 |
| --- | --- | --- |
| All (78) | 70 / 72 | 78 / 78 |
| Single-family, whole protein (`full`, 68) | 67 / 67 | 68 / 68 |
| Family plus extra sequence (`embedded`, 6) | 1 / 3 | 6 / 6 |
| Two families in one CDS (`multi_family`, 4 units, 8 family slots) | 3 of 8 slots at 5 | 7 of 8 slots at 5 |

For ordinary single-family proteins both methods are at ceiling, and this is not
because the references are near copies: the nearest same-family reference has a
median 38% identity, and an [exploratory check](exploratory.reference-distance.json)
that lets only references below 30% identity count still gives alignment 52/52
at rank 1 and ESM-2 47/52 (50/52 within five). The panel is small (134
references); it cannot separate the two methods on ordinary proteins.

The ten fusion/tagged sequences are where they differ. A mean-pooled embedding
of a whole fusion protein loses the labeled family in half the cases; local
alignment does not. This is the observable cost of keeping a fusion as one
protein. It is ten sequences and is reported as an observation, not a rate.

## 3. Parent-level comparison (110 distinct multi-unit parents, 185 supported components)

Components recovered through the component that owns the family
(the registered primary endpoint), five hits per parent:

| Condition | Query units | ESM-2 8M | Alignment |
| --- | --- | --- | --- |
| `single_longest` (baseline, label-blind) | 1 | 55 / 185 | 59 / 185 |
| `single_longest_supported` (label-aware, historical style) | 1 | 106 / 185 | 110 / 185 |
| `all_max` (historical merge) | all (mean 3.2) | 151 / 185 | 135 / 185 |
| `all_round_robin` | all (mean 3.2) | 178 / 185 | 184 / 185 |
| Per-unit ceiling (five hits *per unit*, a larger budget) | — | 180 / 185 | 185 / 185 |

Mean per-parent recall, in the same order: ESM-2 0.32, 0.69, 0.85, 0.96;
alignment 0.34, 0.71, 0.80, 0.99. Over all 127 parents (not deduplicated) the
picture is the same (57, 123, 181, 208 of 215 for ESM-2).

Reading it:

- **More queries, not better representations.** A one-unit condition can recover
  at most one component per parent (110 of 185 here). Everything `all_*` gains
  over `single_longest_supported` comes from querying components the baseline
  never asks about. It costs 349 instead of 110 query-side model calls (3.2×)
  and 100,457 instead of 42,479 residues (2.4×); the 134 reference encodings are
  shared.
- **The label-blind baseline is weak for a mundane reason.** In 51 of 110 parents
  the literal longest protein has no supported label (a replication protein, a
  cargo gene), so `single_longest` scores zero there. Which single protein is
  "the" baseline changes the baseline from 55 to 106.
- **The merge rule decides how much of the gain survives.** Maximum-score merging
  loses 27 (ESM-2) and 49 (alignment) of the components that round-robin
  recovers. One unit takes all five hits in 19 of 110 parents with ESM-2 and 62
  of 110 with alignment, whose raw score grows with protein length. Round-robin
  reaches essentially the per-unit ceiling without comparing scores across units.
  Paired by parent, round-robin versus maximum: better in 26, equal in 82, worse
  in 2 (ESM-2); 41 / 68 / 1 (alignment).
- **Decomposition can also cost.** Among the 50 parents with exactly one supported
  component, the label-aware single unit recovers 50, while `all_max` recovers 45
  and round-robin 48 with ESM-2: unlabeled units take budget. Under a fixed
  budget, adding units is not free.
- **Budget dependence** ([exploratory](results.json), added after seeing K = 5):
  the maximum-merge loss is 44 components at K = 3, 27 at K = 5 and 14 at K = 10
  for ESM-2 (47, 49, 8 for alignment). The panel holds five references per
  family, so a dominant unit can fill a five-hit budget with one family; the
  size of the loss depends on that design choice, its direction does not.

Consistency control: in all 28 single-unit parents the four conditions return
identical results, as they must.

## 4. The historical off-family credit did not recur

The original E2 baseline was credited with a second family through off-family
hits. Here the origin-agnostic family endpoint and the attributed endpoint agree
almost exactly: unattributed recoveries are 4 for `single_longest` and 0–1 for
every other condition (ESM-2), 1 and 0 (alignment). With 28 competing families
an accidental match to another expected family is rare; with the four families
of the historical panel it was likely. The earlier complication was a property
of a very small catalogue, and attribution remains the safer endpoint.

## Extension 1 (added later the same day)

A [pre-registered extension](extension-1/README.md) repeated the comparison on
800 further records with this panel unchanged. On 363 unit sets not seen here
the proportions are nearly identical (ESM-2: 33%, 57%, 83%, 96% of 614
components for the four conditions), the round-robin gain over the best single
protein survives a conservative cluster bootstrap (+0.12 [0.04, 0.20], 20
clusters) while the maximum-merge gain does not, and the fusion observation
holds on new sequences (ESM-2 28 of 53 family slots over both samples,
alignment 50 of 53). An exploratory oracle cut recovers all 53, which locates
the loss in mean pooling over unrelated sequence.

## Answer supported by this run

Keeping every validated annotated protein lets a fixed model retrieve the
documented family of more components under the same five-hit budget: 178 of 185
with ESM-2, against 106 for the best label-aware single protein (1.7×) and 55
for the label-blind longest protein (3.2×) — provided per-unit rankings are
merged by turn-taking rather than by raw maximum score. That benefit is one of
**coverage**: it comes from asking about more components and costs proportionally
more model work. Nothing here shows that preprocessing improves how the model
represents any individual protein; per-unit retrieval is unchanged by
construction. On this panel the conventional alignment baseline equals or
exceeds ESM-2 8M per protein and under round-robin merging; it is worse only
under maximum-score merging (135 against 151), because raw alignment scores grow
with protein length. The study therefore gives no reason to prefer the small
language model for this task. Run times here (25 s against 35 s) describe this
134-reference panel only: alignment work grows with queries × references,
embedding work with the number of sequences.

## Limits that matter

- Six dependence clusters; common markers recur across most vectors. Descriptive
  only. The 78 distinct supported sequences are variants of just 30 curated
  anchor proteins in 24 families (the 68 single-family ones: 27 anchors, 21
  families), so per-unit counts overstate the number of independent observations.
- Labels require ≥ 0.90 identity to a curated anchor, so only well-known
  components are scored; homologs of declared families outside that threshold
  and all undeclared families are `unknown`. The scored set is easier and
  narrower than the constructs themselves.
- Labels and the comparator both rest on sequence alignment. Alignment's perfect
  per-unit score partly reflects that shared basis.
- The panel is small and curated. Ceiling results here say little about
  database-scale search, where remote hits compete with many more candidates.
- All records are well-formed GenBank submissions. Tool exports with missing or
  inconsistent translations are not represented.
- Fusion observations rest on ten sequences.
- The family table, thresholds and conditions were fixed before scoring; the
  budget and reference-distance sensitivities were added afterwards and are
  marked exploratory.

## Reproduce

```bash
python experiments/component_study_v1/fetch_ncbi_sample.py          # network: NCBI
iri-constructs preprocess data/raw/component_study_v1/ncbi \
  --out results/component_study_v1/bundle --max-model-residues 1022 \
  --provenance data/raw/component_study_v1/ncbi/provenance.json
python experiments/component_study_v1/build_panel.py                # network: UniProt
python experiments/component_study_v1/score.py --method alignment
HF_HOME=results/component_study_v1/model-cache \
  python experiments/component_study_v1/score.py --method esm2       # model environment
python experiments/component_study_v1/analyze.py
```

A later provider snapshot of a record or UniProt entry is a new input: the
freeze records every hash, and `score.py` refuses a panel that no longer matches.
The original private bytes of this run are not in Git.
