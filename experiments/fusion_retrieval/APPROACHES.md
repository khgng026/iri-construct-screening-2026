# Approach inventory and missing comparisons

Decision review: 22 November 2026. Research plan, not a new result or a frozen
execution specification. Both upstream heads were checked: working pilot
8251868aa8785de1f025ef74b72fa7b261dcaa40; shared main c05559a.

## Already tested — reuse rather than restart

| Approach | Evidence | Remaining qualification |
|---|---|---|
| Annotation-guided CDS extraction | Shared released toolkit and component studies | Not gene discovery or fusion splitting |
| Sequence-only ORF calling | Ten-record diagnostic: longest one 7/21 exact CDSs, longest three 18/21; expanded start set 14/21 | No matched end-to-end unknown-plasmid retrieval comparison |
| Whole-protein mean embeddings | Five-model October comparison and later 650M studies | Keep as baseline, not a new experiment |
| Fixed overlapping protein windows | October 650M component recall 39.39%→64.73% at five total hits on held-out combinations; later 300/150 windows also tested | Shared proteins across splits; not independent-family performance |
| Motif/tag/linker cuts | Original pilot auto/auto_per_segment arms | Does not split motifless fusions; no new novelty claim |
| Alignment-guided segmentation | Original BHS; local fixed-budget studies | Hybrid using alignment; not a pure AI boundary method |
| Supplied correct boundaries | Original and local oracle diagnostics | Answer-assisted upper-bound diagnostic |
| Extra context around segments | 300-residue floor: RR 30→25, max 16→21 of 51 | This tests one context policy, not every possible collar definition |
| Encode separately then mix; local pooling inside joined embeddings; equal component weighting | October pooling controls | Known-boundary controls, not automatically discovered segments |
| Max / round-robin / whole-half output merging | Existing executable code; primary RR30 versus whole-half26 /51 | No universal winner; whole-half failed primary objective |
| Distinct-sequence refill | Historical replay adds54 slots, zero changed hits | Does not solve allocation/representation failure |
| Residue-level cosine alignment | Original positional implementation and calibration work | Not the published EBA normalization; no demonstrated general AI win |
| pLM-BLAST adaptation | 12-pair ESM150M probe succeeded; full run not admitted on cost | Original paper uses ProtT5; probe is not full retrieval benchmark |
| DCTdomain | Native coordinate fixture failed | No valid biological comparison completed |
| Cheap three-residue composition | Both parents4/6; supplied parts6/6 on public panel | Learned-model comparison still pending |

The 51-query figures use identity cutoff0.25, K10 and one designated target per
query; six targets are unavailable. They cannot establish all-component recovery.
Full original result summaries remain available through the pinned upstream source.

## Worth testing, in priority order

1. **Complete the missing matched multi-component comparison.** Existing guided
   and window results are useful, but they do not supply a prospective comparison
   on the public panel with all known components scored. Freeze preprocessing and
   hold the scorer, reference set, tie rules and final budget fixed.
2. **Consensus / weighted merging.** No completed held-out comparison found.
   Test against existing max and RR using complete saved score matrices first.
   Separate score weighting from earlier equal-weighted *embedding* averages.
   Overlapping windows are correlated; repeated votes are not independent evidence.
3. **Pfam-profile domain cuts.** Candidate additional preprocessing control:
   use predicted domain coordinates, preserve uncovered regions explicitly, then
   apply the same scorer/merger. Domain predictions are not truth labels. Existing
   Pfam/InterPro metadata audits are not an executed profile-guided retrieval arm.
4. **Reference-side windows.** Latest upstream code implements this. No verified
   local performance result found. Compare whole versus windowed references using
   original protein identities as the final budget unit, with length/opportunity
   controls: long reference proteins otherwise get more chances to match.
5. **Ankh residue-level local alignment.** Conditional next representation test if
   appropriate local queries still rank the correct reference poorly. Reuse the
   existing positional interface, retain traceback/coordinates, and distinguish
   localization from retrieval. Not an exact reproduction without a confirmed
   article checkpoint/configuration.
6. **Multi-size or adaptive windows; center-only pooling with a context collar.**
   Plausible development variants, not validated here. Fixed-width/stride sweeps
   and known-boundary local pooling already exist. Define what changes: input
   context, averaged residues, window count or final budget. Do not change all at
   once or call a new arbitrary width a new method.

Defer structural methods (Merizo, structure-based PLASMA), DCTdomain repair,
3Di/Foldseek scaling, learned attention pooling and a new model leaderboard until
an observed failure justifies their extra inputs, training, or infrastructure.
PLMSearch and EBA are relevant references, not interchangeable implementations of
our pooled-cosine search. Motif decomposition and plain/normalized alignment are
existing controls rather than new candidates.

## Exact next work

1. Verify cache identities and recover full chunk/reference score matrices for
   existing **development** cases. Replay max and RR and classify misses as absent
   reference, no useful local signal, or useful signal lost during merging. The
   meeting example alone is not an evaluation set and its handout lacks raw scores.
2. Preserve the sealed prospective v1.2 experiment unchanged. Its three methods
   are whole, BHS-RR, BHS-whole-half; primary K10 accession budget, secondary K5 and
   sequence-budget sensitivity. It has **no window or oracle arms**. Execute only
   after the existing runtime/cache/resource checks; do not rewrite its primary
   endpoint because whole-half previously lost.
3. Specify a separate comparison extension: whole, fixed windows, guided segments;
   initially one fixed merger, RR, to isolate preprocessing. Reuse max and the
   already measured context condition as sensitivities. Add a newly defined
   collar only if its distinction from the completed 300-residue condition is
   explicit. Keep alignment and supplied-boundary diagnostics separate.
4. The existing 12-query panel has already-inspected composition results and no
   calibration pool. New variants on it are descriptive development/diagnostic
   evidence, not an untouched confirmation. Choose methods on development cases,
   then freeze a candidate before genuinely new, family/carrier-aware evaluation.
5. Report correctly attributed components, all-supported-components recovery,
   losses and rescues, source-supported wrong matches, unresolved cases and cost.
   Do not turn unlabeled proteins into negatives or compare raw cosine magnitude
   directly with alignment score magnitude.

Proposed coordination: methods work produces versioned predictions and parameters;
verification checks source labels, scoring rules and reproducibility independently.
The exact meaning of the collaborator's verification commitment remains unconfirmed.

## Sources and scope of review

Historical reports: Oct31 five-model/pooling study; Nov7 fixed-budget adaptation;
Nov11 campaign/continuation; Nov16 literature review; current shared component
study and extension. Local original paths are in MANIFEST.json and the IRI workspace.
This is source reconciliation, not new inference or an exhaustive literature review.

Fresh primary-source search checked these additional options:
- Pfam sequence/profile searches and domain coordinates: https://pfam-docs.readthedocs.io/en/latest/searching-pfam.html
- Ankh local alignment: https://doi.org/10.1016/j.csbj.2025.12.002
- DCTdomain: https://doi.org/10.1101/gr.279127.124
- Merizo (primary abstract; full-page fetch failed): https://doi.org/10.1038/s41467-023-43934-4
- PLASMA (primary abstract, not full method audit): https://arxiv.org/abs/2510.11752

These sources motivate candidates, not measured improvements for this project.
