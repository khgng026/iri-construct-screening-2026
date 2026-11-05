# Verified lessons and unresolved findings

- Existing ESM and T05 environments have incompatible pinned Torch/Transformers versions; preserve both and launch their absolute interpreters.
- AC, normal memory pressure, zero swap and sufficient free disk confirmed at launch. A short probe does not certify a sustained model run.
- All target and reference review and source hashes must be revalidated before scoring.

- External SIGTERM of unmodified hardware guard bypasses its finally cleanup; tested signal-safe entry wrapper fixes ownership cleanup without changing resource thresholds.
- Initial runner failed before query output due to Path/string precedence; repaired, failed evidence retained. Real model recovery subsequently passed.
- Baseline code inventory includes every src module; new unrelated modules make strict dictionary comparison fail. Original module subset and all numerical outputs were independently checked unchanged.
- Toy durability checks are separate from actual model end-to-end recovery; both now passed.

- 650M MPS shortprobe passed and was faster, but nativecatalogue run triggered warningpressure duringload despite60snormalpreflight. Guard terminatedownedchild; CPUfull/acquisition appropriatelyblocked. Pressurelaternormalized; validatedCPUfallback chosen without repeatingMPSallocation. MPS processRSS excludes muchGPU/driverallocation and must not be compared as totalmemory to CPU.

### 2026-10-31T08:36:27.853012+00:00 — larger model and dependency routing
650M improved whole retrieval materially, yet exactcomponent pooling improved much more; fullquery averaging remains limiting. Establishednormalizedalignment already solves this exactreference suite. Route optionalSciPy-dependent tests to existing ESMCenv instead of altering frozenT05. FullESMCprobe now verifies actual raw forwarddtype/device, not merely converted arrays.

## 09:17UTC — budgets, family transfer and local context

Fixed total5 contextual-oracle recall is materially lower than5percomponent; report both. RRF60 remains a negative result. Family removal reverses650M/ESMC300 ordering butnormalizedalignmentstillleads; preserve mask-specificdenominators.8Mcontextchangesconcentratenearjoins including517matchedcontrasts, notcausalranklossproof. DCTsourcefailure remains an admissionblock. IndependentpLMreviewcaughtdescendantcleanupgapbeforeexecution;repairandtestrequired. Sixcore snapshotsverifiedsourcehashes; same diskonly.

### 2026-10-31 10:10 UTC — comparisons and provenance
- At equal five total returns on the same 950 held-out combined queries, alignment-guided segmentation ties exploratory self-normalized SW on exact matches and is worse on source-labelled families after reference removal. This is endpoint-specific: at ten returns, exact-removal segmentation recovers 2,088/2,088 components and276/276 fully scorable queries, versus normalized SW2,074/2,088 and270/276. Do not generalize five-return ordering to every budget.
- Completed common32 context comparison now explicitly binds the frozen design, numerical/storage policy, quantiles, sequence hashes and profile implementation. v2 scientific report equals v1 exactly; only acceptance checks strengthened. Both models'19 matched group-target pairs show higher far-join cosine and lower relative L2. This is descriptive evidence, not a causal decomposition or calibrated cross-model scale.
- A second guarded process correctly failed the single compute lock. Saved-scalar reporting can use the explicitly coordinated lightweight slot with a finite timeout; this does not authorize parallel model/alignment workers or a weakened guard.
- Dependency tree integrity must bind the actual runtime import path, not only an approved folder elsewhere. Independent review caught a sibling-overlay gap before pLM resource-v2 probes. New runner928e6a0a…2fc5ec binds parent and worker paths;25 focused tests passed, scientific functions/design unchanged. Original v1 negative admission and both earlier v2 freezes remain preserved. Rerun both affected probes under the new source identity.
- Fragment audit coverage must distinguish all966 held-out queries from950 combined queries. The original1,933 fresh-fragment attempts included2 clean-sanity attempts; prospective replay samples1,931 combined attempts/1,117 distinct fragment-mask cases. Seven fixed cases cover six policy×suite strata plus the longest906aa fragment; no outcome-driven substitutions.

## 11:08 UTC — catalogue competition, order, and acceptance checks

- All five models completed the same frozen1,000-reference panel. Combined whole-query recall@5 at101→1,000 unique references changed from.292614→.263636(8M),.359659→.316761(150M),.395170→.375568(650M),.424432→.395739(ESMC300),.449432→.407670(ESMC600). Contextual known-component pooling with five returns per component remained.999432/1.0 forESMC300/600 at1,000; isolated oracles all1.0. Catalogue size and composition change together, and the oracle has a larger budget. Exact-accession105 and deduplicated101 baselines are separate configurations.
- ESMC600 joined whole-query pooling can beat its length-weighted isolated mixture on aggregate; do not present contextual encoding as uniformly harmful. Known-boundary weighting controls are prospective secondary diagnostics, not an automatic method or additive causal decomposition.
- Dedicated saved order analysis passed all5×1146 mixture invariance checks within tolerance with zero rank inversions. Whole-query top5 hit status varied within27/37/34/39/31 of680 target-group cells. Four-component groups cover8of24 permutations but balance ordinal positions. Correlated group/target observations do not support IID significance or a causal position claim.
- Independent review caught that a clean-vector tolerance check alone did not enforce the design’s separate cosine/rank sanity requirement for new equal-weight controls. Producer97c84f27… now records both matched clean full-catalogue cosine checks under unchanged1e-12/2e-12 policies; original code preserved and outcomes remain unscored pending tests/dependencyfreeze.
- A real ten-record native compression trial preserved exact logicalbytes and metadata. Fresh-process verification covered all3441 units in the frozen three-store plan before removing ten redundant owned backup links. Allocated reduction7,577,600bytes is distinct from actual free space. Full compression is separately root-admitted and running; no full-scope success claim yet. macOS rejectedRLIMIT_DATA; the implementation accurately uses observed aggregateRSS stop checks, bounded files/streaming/one child, and unchanged prepared resource thresholds.

## Closeout — accepted results and execution limits

- Equal-weight controls and their independent scalar checker completed. Equal weighting does not generally repair pooling: every model lost four-component recall under equal-context weighting, despite small gains in some aggregate/query endpoints. Producer and checker shared an author; separate-agent static review is complementary, not independent empirical replication.
- All five models and both alignment arms completed the 1,000-reference comparison. Raw/normalized alignment stayed at82.56%/100% component recall@5 while whole-query model recall fell with added competitors.
- Four-model held-out segmentation comparison passed: whole→windows recall@5 was28.68→49.06%,35.64→59.77%,39.39→64.73%,42.27→68.21% for8M/150M/650M/300M. Normalized alignment already reached100% on that panel.
- A closing report invocation failed because it supplied a fifth core model to an explicitly four-model adapter. The corrected four-model invocation passed; the failure remains recorded. Prepared queue12 is not execution-ready until its five-model reporting adapter is separately revised/tested/reviewed.
- Native compression finished all3,441 records with full logical-byte verification; allocated blocks fell by5,115,174,912bytes, distinct from filesystem free-space change. No pending transaction remains.
- Usage limits interrupted oversight; completed queued work does not imply uninterrupted supervision. The amended eight-hour window prohibits new experiments on resumption. Closure uses saved evidence only. New experiments need a fresh bounded session, not altered historical deadlines.
- Six casebook/unit-normalization draft files remain unreviewed and untested. Do not infer readiness from their existence, a handoff description or fictional tests that have not run.
