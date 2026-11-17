# Component-coverage study v1 — specification

Written 2026-11-17 before any retrieval score for this study was computed. Agent
self-review only; no independent human review is claimed. This is a new,
separately versioned study. It does not modify or re-score the historical E2/E3
runs, whose rules and results stay as recorded.

## Question

For a construct with several source-annotated coding sequences, does keeping
**all** validated proteins characterize more of its documented components than a
principled reduction to **one** protein — with the same fixed model, the same
reference panel and the same total number of returned hits — and how much of
any difference is due to (a) simply asking more queries, (b) the rule used to
merge per-protein rankings?

"Characterize" means only: a reference from the component's documented protein
family is returned. It is not a function or safety assessment.

## Inputs

- **Parents.** A seeded random sample (seed 20261117, n = 200) of the 11,211
  records matching one fixed NCBI nuccore query for complete cloning/expression/
  shuttle vector sequences (`fetch_ncbi_sample.py`; identities in
  `docs/results/component-study-v1/acquisition.json`). Drawn before any record
  content was inspected. No record is added or removed by hand.
- **Units.** The `iri-constructs preprocess` bundle of those records: every
  source CDS whose `/translation` validates. A unit is *eligible* as a query when
  it is exported and has at most 1022 residues (the model limit). Withheld and
  over-length units stay in the coverage tables and are never truncated.
- **Model.** `facebook/esm2_t6_8M_UR50D` at revision
  `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`; final-layer residue states, mean
  over residues (special tokens excluded), cosine similarity; float32, CPU. The
  conventional comparator is the project's existing BLOSUM62 local alignment raw
  score. Both use identical queries, references, conditions and budget.

## Family labels (what counts as documented)

Families are declared in `families.v1.json`: curated UniProtKB/Swiss-Prot
evidence (Pfam cross-references plus the reviewed "Belongs to …" statement where
one exists). The table was written after tabulating source product/gene names of
the sample and before any retrieval was computed.

A unit is labeled with family F when a local BLOSUM62 alignment to one of F's
anchors has identity >= 0.90 over the aligned region and covers >= 0.80 of the
anchor. Class `full` if the alignment also covers >= 0.80 of the unit, otherwise
`embedded` (the unit carries additional unlabeled sequence, e.g. a tag or fusion
partner); `multi_family` if anchors of two or more families are contained. A
unit matching an anchor at >= 0.90 identity but covering < 0.80 of it is recorded
as `fragment` and left unlabeled. Everything else is `unknown`. Labels therefore
rest on near-identity to a curated protein, which is independent of the model
under test; they are not experimentally verified function. Source product names
are tabulated beside the label as a consistency check and do not create labels.

## Reference panel (one panel for every parent and condition)

For each family, candidates are reviewed entries satisfying the family rule
(UniProt REST, first 500 by accession), 50–1022 residues, canonical residues, not
fragments. In the order of `sha256("20261117:" + accession)` take the first five
that are (i) not anchors, (ii) below 0.90 identity / 0.80 coverage to every anchor
and every labeled unit of that family, and (iii) below the same threshold to
references already chosen. This removes exact copies and direct lineage copies of
the queried components. A family is *scorable* with at least two references.

## Conditions (total budget K = 5 hits per parent in every one)

| Condition | Query units | Merge |
| --- | --- | --- |
| `single_longest` (baseline) | the longest eligible unit, chosen without looking at labels | its top 5 |
| `all_max` | every eligible unit, labeled or not | best score per reference over units, top 5 (the historical rule) |
| `all_round_robin` | every eligible unit, ordered by length (longest first) | units take turns adding their next unseen reference until 5 |
| `single_longest_supported` (sensitivity) | the longest *supported* unit (the historical style of baseline) | its top 5 |

Ties: higher score first, then reference ID. Unlabeled units take part in the
`all_*` conditions because a real user does not know labels in advance; their
hits use budget and are counted as unscorable, not as right or wrong. Only these
two merge rules are evaluated; none is tuned.

## Endpoints

- **Primary — component-attributed recall@5.** Per parent: the fraction of its
  supported components (eligible, labeled, family scorable) for which a hit
  *produced by that component* belongs to its family.
- **Historical-style — parent family recall@5**, where a family counts whichever
  unit retrieved it; reported with the number of recoveries that attribution does
  not support, to show how much credit that endpoint gives by accident.
- **Per-unit retrieval** (each distinct supported sequence queried alone):
  hit@1, hit@5, reciprocal rank of the first same-family reference. It is
  identical in every condition by construction; it measures the model, not the
  preprocessing.
- **Coverage accounting** per parent: source CDSs, exported, withheld by reason,
  over the model limit, labeled, supported; noncoding features; share of the
  construct length inside exported CDS segments.
- **Cost:** query-side model calls and residues per condition.

Parent analysis sets: *multi-unit* (>= 2 eligible units and >= 1 supported
component) for the comparison; *single-unit* parents as the consistency control
(all conditions must coincide). Results are reported over all scored parents and
over distinct parents after collapsing those with identical eligible-unit
sequence sets. Parents sharing any unit sequence are dependent; the number of
connected dependence clusters is reported and no population inference is drawn
when it is below eight. Per-unit results are over distinct sequences.

## Interpretation fixed in advance

- Under attribution a one-unit condition can recover at most one component, so
  some advantage for `all_*` is structural. It is reported as the effect of
  asking more queries, with its cost, and is not described as a better model
  representation.
- The informative contrasts are `all_max` versus `all_round_robin` (aggregation
  under a shared budget), both versus the per-unit ceiling (what the budget
  costs), and attributed versus historical-style recall (endpoint artifact).
- `embedded`/`multi_family` units are whole annotated proteins. Their results
  describe the limit of not splitting proteins; they are not a test of
  within-protein decomposition.
- Noncoding content is outside what a protein model can characterize and is
  reported as uncovered, not dropped.

## Resource bounds

CPU only; each scoring stage under 30 minutes and 6 GiB; one model; weights
under 1 GiB. Private inputs (records, reference sequences, embeddings, weights)
stay in ignored storage; committed results are source-free.
