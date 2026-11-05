# T07 evaluation extensions v1.1 — prospective specification

**Status:** User-approved experimental directions on 2026-11-05; **implementation and scientific validation pending**. This addendum extends, but does not retrospectively replace, frozen [T07 v1.0 on PR #2](https://github.com/khgng026/iri-construct-screening-2026/blob/codex/real-construct-pilot/docs/week1/T07_benchmark_protocol.md). Maintain design/config hashes, versioning and a clear change log. The source/label panels remain *unadmitted*, so experiment runs must be gated.

## E1. Annotation-quality sensitivity / appropriate abstention

**Question:** Does source annotation completeness determine how many model-ready protein units can be safely extracted, and does the parser preserve uncertainty instead of inventing features?

- Pre-register strata: reviewed well-annotated, partial/ambiguous, and source-only/no observable CDS. Classify from independently reviewed source features and their documented caveats, **not** from model outcomes.
- Compare originals versus normalized T03 source features, audit sidecar and eligible protein units; verify stable IDs, locations, strands, qualifiers, translation checks, and provenance. Count withheld/unknown/ambiguous units by reason.
- Report: source annotation coverage (with denominators), feature/qualifier preservation, confirmed source-aligned translation rate, erroneous invention count, useful-model-input rate, and justified abstention rate. Show category-level counts and sample records with non-sequence metadata.
- Source-only controls (e.g., pUC19c if separately admitted) must report `no_observable_annotations`, **not** 100% agreement. Bad/missing source annotations must not be scored as parser errors when faithfully preserved; parser may never hallucinate absent CDS.
- No review/permission shortcuts, ORF inference, gene completion, or manufactured ambiguous source mutations.
- **Acceptance:** machine-readable per-record audit, aggregate tables with eligible denominators, tests for no-CDS, missing/conflicting/joined/reverse/fuzzy locations, and no real-data measurement claimed if review gate is missing.

## E2. Multi-CDS parent-level comparison

**Question:** Does representing every *documented* CDS improve the range of documented protein families retrieved at the **parent construct** level compared with a valid one-CDS reduction?

- Implement T07's previously prospective multi-CDS score merging. Baseline = single eligible source-annotated CDS (largest translation; stable tie-break); decomposed condition = all eligible annotated CDS units. Neither is raw DNA passed to a protein-only model. No sequence concatenation.
- Pre-freeze expected **documented family set** per parent and independent reference eligibility. Merge per-CDS reference scores by maximum within each model; deduplicate references, then return exactly `K=5` **total** hits per parent, not K per CDS.
- Report at minimum: fraction of distinct expected documented families retrieved; proportion of parents with all eligible documented families retrieved; parsed CDS coverage; abstentions; additional model calls, processed amino-acid count, wall time, peak memory; full versus matched supported subsets.
- Explicitly distinguish *more units represented* from *higher per-unit retrieval accuracy*; the one-CDS baseline intentionally loses coverage. Where feasible also report per-unit performance under matched query conditions as a separate control. Do not claim function prediction.
- Require at least one independent same-family reference and a competing-family reference for a query to count as scored, plus exclusions from T07. Families lacking independent examples remain unscored with transparent denominator changes.
- **Acceptance:** deterministic score merge and top-K tests, multi-CDS/duplicate-family/missing-CDS/empty-parent tests, matched-panel attribution, resource logs, and clear 'blocked' behavior when no admitted panel exists.

## E3. Controlled independent-reference generalization

**Question:** Do embedding representations maintain useful documented-family retrieval when easy, overly related references are excluded, relative to raw local alignment?

- Declare difficulty tiers **before scoring**, for instance: (1) remove exact sequence duplicates, (2) exclude same parent, shared vector and independently reviewed linkage groups, (3) stricter sequence-relatedness exclusion where a defensible independent homology/identity policy and enough same-family references remain. No arbitrary threshold chosen after seeing results.
- Evaluate **the same** query IDs and eligible reference IDs for ESM-2, ESM-C (when feasible) and conventional alignment, with the **same total K**. Keep raw BLOSUM62 local alignment as primary baseline; normalized alignment exploratory as already specified in T07.
- Freeze and report all exclusions, relatedness definitions, valid positive/negative reference counts, family/query/parent independence, reference catalog size, source-version hashes and per-tier sample attrition.
- If a difficulty tier has no independently supported positive or competing reference, report **abstention/null**, never failure or artificial success. Avoid claiming unseen-protein generalization when proteins or families appear in development sources or pretraining exposure is unknown.
- Report retrieval and **false-family-positive** rates where defensibly measurable, coverage, per-tier paired differences/uncertainty, and explicit confounding by selection and catalog size. No sequence obfuscation, fragment reconstruction, hazardous screening, or synthesis.
- **Acceptance:** deterministic tier partition tests, reference-leakage tests, panel-identity parity, abstention and denominator accounting, and accurate limitations.

## E4. Quality–compute trade-off

**Question:** Are improvements, if any, worth the extra inference work?

- Combine measured E2/E3 and existing A–D run records on **identical eligible panels**. Do not generate fabricated speed/cost numbers or pool incompatible model scales.
- Report separately: coverage, documented-family recovery, complete-family recovery, latency (cold loading vs inference and pooling), total CPU/GPU usage where known, process memory, number of model calls, unique encoded sequences and residues, cache policy, and hardware/platform/version.
- Show a transparent cost/quality trade-off (tabular and/or simple Pareto visualization) with clear NULL for unavailable measurements; mark exploratory, incomparable-hardware or insufficient-independent-panel cases. Do **not** automatically pick a model winner.
- Include a simple low-resource baseline where possible (source annotation-only or alignment) and the same total hit budget. Record throughput as observations per appropriate unit, and flag unique-data-versus-cache effects.
- **Acceptance:** report aggregation that rejects mismatched panels/configs/hardware as direct comparisons, tests for absent metrics and denominators, and reproducible output artifacts.

## Scientific and procedural rules common to E1–E4

1. **T02 gating is independent:** approved candidate shortlist ≠ approved raw source files/labels. Do not score or run pretrained sequences until the required exact-record review and permissions are recorded.
2. Keep the original T07 v1.0 version as the historical baseline and version proposed extensions as v1.1. Protocol and pilot source/label panel must freeze before first scored result. Breaking design changes require a new version/reason, not quiet edits.
3. Test with clearly fictional *software fixtures*; they are not biological source evidence. No fabricated functional annotations, automatic family approval, or generated sequences.
4. Same-model within-space metrics only; embedding cosine spaces differ across ESM model types.
5. Keep overlap between the ten pilot constructs and the independent reference set explicit. Ten constructs are not ten independent targets. If the dataset is underpowered, report descriptive feasibility rather than statistical significance.
6. Preserve separate Git history and **PR #1 read-only**. Reuse useful T05 parser and scorer modules from PR #2 instead of inventing a second framework.
7. Apply E1 and E2 first. E3 depends on credible independent references; E4 depends on actual measured E2/E3 or prior T05 results.

## Future scale

After a defensible pilot, plan a versioned, stratified expansion toward **150–250 reviewed benign public constructs if quality permits**, with sufficient independent protein-family references. Include review workload, diversity, exclusion counts, data rights and separate development/test cohorts. Do not treat a ten-record pilot or this target as proof of population-wide performance.
