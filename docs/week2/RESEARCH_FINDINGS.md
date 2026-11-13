# Findings from source preparation and technical validation

The central retrieval hypothesis remains **unanswered by new model experiments**.
No biological embedding, alignment ranking or scored E1–E4 comparison was run.
The completed source work nevertheless identifies concrete constraints and a
smaller defensible next study, rather than treating successful downloads or
software tests as scientific performance.

## Annotation quality limits the available experiment

The main preparatory set contains 271 source features and 18 source-validated
protein translations. All ten inspected iGEM entries have deposited feature
responses, but no provided protein translations in those responses. Coordinates
are retained as original metadata while their indexing convention remains
unverified. Downloadability and a `CDS` role therefore do not establish a
model-ready protein or a reference-family label.

The bounded native GenBank export follow-up succeeded for two selected Registry
parts. It confirms source DNA equality and usable native coordinate semantics
for those snapshots, but all four formal CDSs still omit provided translations.
Thus the missing-protein limitation persists in both inspected representations;
this is not a claim about every Registry part or export format.

Fusion-ready architecture similarly does not establish a genuine fusion. The
separate MK448012.1 example supplies a complete documented fusion CDS and validates
the parser's full-protein/provenance handling. It does not establish automatically
discoverable boundaries, isolated-target effects or scored fusion performance.
Descriptions of marker CDSs on `misc_feature` rows must remain metadata.

This favors a smaller source-reviewed dataset over expansion toward 150–250
records. More copies of uncertain annotations would increase denominator size
without resolving the core evidence problem.

## Distinct constructs can encode identical proteins

Exact protein fingerprint inspection found that U47295.2, U47296.2 and AY738222.1
encode the same 550-aa reporter translation. Their parent DNA and annotation
records are distinct, but they are not three independent protein queries. Common
marker proteins create further dependence across vectors. This finding uses
exact equality, not an alignment threshold selected after retrieval results.

The scoring correction also prevents multiple different CDSs from one parent
from being counted as independent bootstrap groups. Any apparent confidence
gain from increasing CDS count would otherwise be misleading. Conservative
connected provenance groups control known dependence; they do not prove that
the remaining groups are evolutionarily or pretraining-independent.

## More components need not mean more families or better predictions

GFP-like green, red, cyan and nonfluorescent proteins belong to a proposed
molecular homolog family; emission color is not a separate family definition.
Two represented fluorescent components may therefore add CDS coverage while
leaving distinct-family recall unchanged.

AF298787 provides a prospective example with documented URA3, GFP and AmpR
translations. The proposed independently supported subset covers URA3 and GFP;
AmpR remains unsupported in the current reference proposal. The longest eligible
supported CDS would be 267-aa URA3, whereas the literal source-longest CDS is
286-aa AmpR. This must be disclosed rather than silently treating the filtered
baseline as the complete construct. A future result recovering GFP after adding
its CDS could demonstrate retained component information; it would not by itself
show better per-protein prediction or biological-function validation.

## Reference composition is part of the claim

The prospective catalogue has 17 reviewed UniProt entries. Four query proposals
have positive and competing references after copy/ancestor/lineage exclusions;
Renilla has no independent positive and remains held. Renilla and firefly
luciferases are not merged because they share a light-producing function.

All proposed supported rows have more than five eligible references, avoiding
the trivial saturated-budget case. Nevertheless, the four query components come
from only three parents. They cannot support the predeclared eight-independent-
cluster confidence interval. A strict relatedness tier also remains unavailable
without an independently justified, fingerprint-bound policy. No homology cutoff
was invented to make that tier pass.

## Model feasibility is now demonstrated at the interface level

Both fixed ESM-2 and ESM-C checkpoints execute offline on a fictional software
fixture, with validated float32 residue outputs and deterministic pooling.
Initialization costs differ substantially, and their representation widths are
320 and 960. These facts establish runnable interfaces, not comparative retrieval
quality. Fixture timing cannot answer when added compute outweighs improved
family retrieval, and cosine coordinates cannot be compared between models.

The ESM SDK's omitted import dependency and the adapter's missing explicit-SHA
match guard were corrected. Passing dependency metadata checks alone was not
sufficient to establish an importable model runtime; an actual offline forward
and exact-snapshot verification were needed.

## Next valid experiments

1. Obtain actual delegated decisions on the exact source/label/reference packet.
   E1 can be cleared independently of protein-family reference completeness.
2. If cleared, run the source-bound E3 supported queries on identical reference
   sets with alignment, ESM-2 and ESM-C. Report descriptive paired outcomes and
   null uncertainty intervals at the current parent count.
3. Run AF298787's supported-subset E2 only if explicitly cleared, retaining its
   full unresolved component/family denominator and eligible-versus-literal
   longest-CDS distinction.
4. Feed only actual matched E2/E3 costs and quality observations to E4, separating
   loading, shared catalogue work and sequence-cache hits.
5. Expand independent benign parent/family coverage or repair missing source
   annotation evidence before pursuing a larger population-level comparison.

These recommendations preserve the predeclared endpoints. None asserts a
decomposition benefit, winning model, biological-function prediction, completed
paired benchmark or scientific clearance that has not occurred.
