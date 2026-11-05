# IRI component retrieval protocol

## Research questions

1. Does retrieval change across ESM-2 model sizes and other compatible protein representations?
2. How does recovery vary with component count, component fraction, order, and catalogue competition?
3. Does local alignment recover components missed by whole-query embedding retrieval?
4. How much of the change arises from averaging components together, and how much from encoding them in joined context?

The completed pilot used ESM-2 8M, 105 reference accessions (101 distinct sequences), 16 target proteins, and 230 ordered two-component queries. Its 52.4% component recall@5 is an exact-accession result for this selected dataset. It is neither general protein identification accuracy nor evidence that existing tools fail. The existing report and manifest are the source of those counts.

## Model admission and priorities

| Priority | Candidate | Purpose | Admission condition |
| --- | --- | --- | --- |
| Core | Existing ESM-2 8M | Preserve and extend the baseline | Reproduce the existing benchmark first |
| Core | ESM-2 650M, `facebook/esm2_t33_650M_UR50D` | Same-family checkpoint comparison and connection to Nishit's reported checkpoint | Exact revision pinned; official loading route; resource probe passes |
| First extension | ESMC 300M | Different protein representation | Hardware handoff identifies supported checkpoint and environment; residue mapping and licence verified |
| Further extension | ESMC 600M | Capacity comparison within ESMC | 300M interface verified and resources permit |
| Further extension | ESM3 small, `biohub/esm3-sm-open-v1` | Sequence-only representation comparison | Sequence-only hidden states verified; extra tracks consistently absent; access and resource gates pass |
| Optional bridge | ESM-2 150M | Intermediate same-family capacity point | Core comparisons complete and marginal download/runtime justified |

Candidate names are not claims of local compatibility. Resolve immutable revisions and file hashes at launch. For ESMC, reconcile checkpoint format with the hardware chat's chosen official SDK; do not interchange legacy and current checkpoint IDs silently. No 3B/6B/15B model, DNA model, structure generation, sequence generation, or training is included by default. Adding one requires a documented research question and feasible resource/access assessment.

Use each model's last-layer actual-residue hidden states with mean pooling as the primary comparison. Exclude padding and special tokens. Use a fixed, documented representation for ESM3 without supplying structures or function labels that other models lack. Embedding dimensions may differ; compare retrieval metrics, not raw cosine magnitudes across models. Do not select a best layer after viewing test outcomes. Any layer/pooling exploration is secondary and recorded separately.

Changing checkpoints also changes architecture and training, not just parameter count. A checkpoint comparison does not isolate a causal effect of model size.

## Fixed shared suite

Use the existing reviewed catalogue and 16 unique-sequence targets. The prepared `challenge_manifest.json` carries original query IDs, component order, exact boundaries, length/fraction, sequence hashes, exclusions, repeatability IDs, and deterministic selection rules. It contains no sequence strings. At run initialization, verify all source file hashes against it and recheck unchanged data-review status.

The shared operational cap is 1,022 residues, inherited from the pilot. This is an experimental policy, not a universal architectural limit. If a model supports less, declare its coverage and report both its full supported set and the cross-model intersection. Never truncate silently. Longer-input tests may be a separate extension after verifying each interface.

| Suite | Selection | Purpose |
| --- | --- | --- |
| Clean | All 16 existing targets | Self-retrieval, mapping, and clean-query checks |
| Two components | Preserve all 230 existing ordered queries | Direct comparison with prior results |
| Three components | Up to 90 unordered groups, balanced by minimum fraction; all six orders | Component count and internal-position tests |
| Four components | Up to 45 unordered groups, balanced by minimum fraction; cyclic orders and their reversals | Additional count and position tests at bounded cost |

Select groups by deterministic hash within fixed minimum-fraction bins before new inference. Use full distinct canonical proteins; do not mutate, repeat, trim, or insert designed linkers. Record all considered groups and length exclusions, as well as eligible groups omitted by the workload cap. The builder defines the exact seed and tie rules.

This suite is designed after observing the earlier pilot. It is a planned extension, not a blinded external test. More permutations are correlated observations, not additional independent proteins. Count effects and fraction effects covary: report stratified descriptions and matched subsets where available, without claiming isolated causality from their aggregate association.

## Catalogue competition extensions

The primary catalogue remains the fixed 105 accessions. Do not change it during a model run. The first competition sensitivity test uses three independently hashed nested subset selections from existing unique-sequence groups: always retain all 16 target sequences and sample other groups by fixed hash to reach 25, 50, and 101 unique sequences. The full 101-sequence endpoint is identical across the three selections and needs scoring only once. Compare against the deduplicated 101-sequence catalogue so accession duplication does not masquerade as catalogue size. Freeze subset IDs before scoring; retain the original accession analysis separately. These sets contain different competitors as well as different counts, so they are not a pure causal estimate of catalogue size.

If the shared core, controls, and at least two model comparisons are complete, prepare larger nested catalogues targeting 250, 500, and at most 1,000 entries. Apply T02's record-specific checks and source terms. Include related competitors, document the selection mechanism, and report achieved counts rather than filling with unreviewed records. Cache each catalogue once per model. This extension must not displace the shared comparisons or audit.

Do not treat relatedness composition as size. Report closest-neighbour similarity and exact duplicates, and keep separately sampled random versus related-competitor extensions distinguishable.

## Alignment comparison

Use exhaustive Smith-Waterman local alignment against the same references, with BLOSUM62 and a documented affine gap convention (open 10, extend 0.5; verify the library's first-gap convention). Pin the implementation and settings. Validate scores and aligned coordinates on short hand-checkable examples and a representative subset using an independent implementation.

Rank references by descending raw alignment score and accession for exact ties. Record local query/reference spans, aligned length, identities, coverage of both sides, and ties. A raw score favours longer high-scoring matches; report that limitation and any prespecified normalization as a secondary result. Do not tune the ranking to match the known components. Known boundaries may be used for evaluation only, never to cut alignment queries in the primary whole-query arm.

Identical canonical components are in the catalogue, making this an easy but essential alignment control. A win here does not demonstrate remote-homology identification. The alignment arm provides retrieval scores, not calibrated hazard or function decisions.

## Averaging and context controls

For a query with components of lengths L_i, keep one fixed reference catalogue and representation layer. Use raw, unnormalized residue/mean vectors until after aggregation; cosine normalization happens only for ranking.

| Arm | Construction | What the comparison can reveal |
| --- | --- | --- |
| Isolated component | Encode each component alone; mean its residues | Clean reference control, expected self-retrieval |
| Isolated mixture | Length-weighted mean of separately encoded component vectors, sum(L_i v_i)/sum(L_i) | Mixing summaries without joined context |
| Contextual component | Encode the joined query once; mean only the residues within each known component span | Change from joined context without mixing other spans into the final mean |
| Whole query | Mean all residues from that same joined forward pass | Combined context and averaging |

Verify numerically that the whole-query vector equals the length-weighted mean of contextual-component vectors, within the dtype tolerance. Verify isolated-mixture invariance to component order. The contextual-component arm still uses supplied boundaries and is not automatic segmentation.

Whole versus isolated mixture compares context under whole-query aggregation. Contextual versus isolated component compares context at a component. Isolated mixture versus isolated component illustrates mixing. These are controlled representation contrasts, not an additive causal decomposition of rank loss; ranks are nonlinear and depend on all competing scores. An equal-component-weight mixture may be a labelled secondary control, not a general solution requiring no boundaries.

Run all controls on 8M. Repeat them on 650M and admitted alternative models where feasible. Preserve per-component and per-query outputs, including cases whose ranks improve.

## Reference removal and generalization

Removing an exact reference makes exact-accession recall undefined for that target. Never score the absent accession as a retrieval miss. Do not label the most similar remaining protein a correct functional match merely because it is most similar.

Two distinct optional extensions are allowed after freezing their definitions:

1. Exact-reference removal diagnostics: remove every duplicate of the exact target, record surviving neighbour identities, rankings, and similarity changes. Do not report identification accuracy without acceptable-match labels.
2. Annotated-family retrieval: define acceptable reference sets from reviewed family/domain evidence independent of the methods under test, remove exact/near-duplicate matches under a documented identity/coverage rule, then measure family recovery. Mark targets with no valid surviving reference as unscorable and report the denominator. Family-level and domain-level labels are not equivalent to whole-protein function.

Freeze similarity convention, masks, labels, and exclusions before model scoring. A catalogue holdout does not remove proteins from pretrained model weights. If credible labels or enough reference relatives are unavailable, deliver the diagnostic extension and explain why generalization accuracy was not measured.

For each multi-component query, construct one common catalogue mask removing exact-sequence equivalents of every evaluated component, not just the component currently being scored. Apply that same mask to every model/method and arm. Freeze any near-duplicate removal and acceptable surviving labels together. Mark unscorable components individually, retain their parent-query membership, and report the number of queries with all components scorable. All-components accuracy is undefined for a query with an unscorable component; report these exclusions against the original query count instead of silently reducing the denominator.

## Metrics and evidence

Primary outcomes are component recall@5 and @10 and all-components-hit@5 and @10. Also report any-component hit, component ranks/MRR, and top-1 with its denominator. With m components and a top-k list, per-query component recall cannot exceed min(1, k/m), and all-components-hit is impossible when k < m. Do not compare raw metrics across component counts as if their ceilings were fixed.

Report exact-accession and duplicate-equivalent outcomes, exact/near ties, per-target macro means, query counts, distinct group counts, fraction bins, position, and common supported subsets. All-components recovery gets harder as the number of required components increases; show component-wise recovery alongside it. Family outcomes require independently defined labels. Top-k hits are not calibrated positives, so do not call all non-component hits false positives or claim precision/specificity without an appropriate label definition.

Report fractions below .25, .25–<.40, .40–<.60, .60–<.75, and >=.75, with per-cell counts. Use descriptive analysis, matched order comparisons, and per-target summaries. Avoid p-values treating permutations as independent samples. Do not infer mechanism solely from correlations of length fraction and rank.

Store complete rankings, derived vectors where storage permits, IDs/spans, model revision/settings, source hashes, code commit, device/dtype, timings, memory measurements, skipped cases, exceptions, and cache metadata. Raw sequences stay in existing ignored storage. Retain residue states only as needed for controls; stream summaries instead of saving every model's full hidden-state tensors by default.

## Verification and completion

Before full execution, run the shared end-to-end trial defined in OPERATIONS.md. At completion, run existing checks plus necessary adapter, alignment, masking, multi-component metric, cache, and control-identity checks. Repeat a fixed representative subset in a fresh process. Use the existing CPU float32 1e-6 absolute tolerance for matching configurations; accelerator/reduced-precision runs require a separately justified tolerance and explicit rank-change/tie report. Do not require exact GPU bitwise agreement or silently loosen tolerances.

Independently recompute a sample of scores and every reported summary denominator. Preserve failed jobs in the ledger. Evaluate all models on a common subset before comparing them. A partial run is not a full-model result.

The morning report answers each of the four questions with observed evidence or a specific unresolved reason, supplies comparative tables/plots, identifies unexpected counterexamples, and recommends one next research step. Explain what established methods solve before claiming an unsolved gap. Attribute the conceptual experimental approach to Nishit; do not redistribute his confidential slides or claim to reproduce his code/results.

## Sources

- Existing `docs/week1/T05_expanded_report.md`, `T05_expanded_design.md`, and manifests in the local experiment repository.
- [ESM-2 650M official model card](https://huggingface.co/facebook/esm2_t33_650M_UR50D).
- [ESMC 300M official model card](https://huggingface.co/biohub/ESMC-300M).
- [ESM3 small official model card](https://huggingface.co/biohub/esm3-sm-open-v1).
- [DCTdomain research paper](https://genome.cshlp.org/content/early/2024/09/30/gr279127124): relevant existing segmentation/retrieval work for the literature comparison.
