# Q2 reporter replication — handoff

## Observed result

On 86 fresh, independent public benign reporter constructs (169 validated components; 93
reporter, 76 marker), gene-calling preprocessing beat minimal six-frame preprocessing at a
cap of two inputs by **+11.8 pooled / +13.0 equal-lineage-group points** (BLOSUM62). The
gain carries to reporter components (reporter equal-group **+15.9**, >= marker +7.1) and
survives removing the dominant lineage (+13.1). Mechanism: candidate **selection**
(zero characterization failures at cap 2). Full numbers in [RESULTS.md](RESULTS.md).

## Independent lineages (the binding constraint)

- **21 lineage groups** total (>= 90% within-family identity), of which **14 are reporter**:
  gfp_like 7, firefly_luciferase 4, beta_galactosidase 2, renilla_luciferase 1.
- Two clusters dominate (gfp_like#1 = 39 EGFP-type, firefly#1 = 28), so the effective
  diversity is well below the component count.
- Consequence: wide bootstrap intervals. The all-component gain is +13.0 **[+2.4, +26.3]**;
  reporter-only is +15.9 **[0.0, +34.7]**. Every **point** estimate clears +5; no reporter
  **lower bound** does. This is a verified pilot, not a certified replication.

## What reproduced

- Published cohort denominators and gains reproduced from the tracked summary
  (`results/paired-cohort-v2/holdout/summary.corroborated.txt`): ESM-2 equal-group +9.2
  [6.6, 11.8], BLOSUM62 +9.0 [6.5, 11.6] at cap 2.
- This study's headline recomputed from the saved verifier output independently of
  `analyze.py` — exact match (+11.8 pooled, +13.0 equal-group, +15.9 reporter equal-group).
- Toolkit validated all 169 components against deposited DNA. Existing paired/reporter/
  grouped tests pass (14), plus 5 new regression tests for this study.

## Remaining limitations

1. **Labels rest on depositor `/product` annotation**, not an independent reviewed anchor.
2. **Few distinct reporter proteins** (renilla 3, beta_gal 2; gfp/firefly dominated by one
   near-duplicate cluster each).
3. **BLOSUM62 only** — no ESM-2 cross-scorer here (no cached weights); the cohort's
   cross-scorer agreement is therefore not re-checked on this panel.
4. **NanoLuc excluded** (13 constructs) for want of an independent reviewed reference.
5. Only **3 tagged/fused** constructs — fusion behaviour untested.

## Next decision

The effect's direction and magnitude replicate on reporters; the open question is precision.
To convert this pilot into a replication **claim**, the next agent should, in priority order:

1. **Grow independent reporter lineages to >= ~30** (target ~40+ lineage groups overall).
   Sources not yet mined: Addgene-deposited reporter plasmids, eukaryotic expression vectors,
   plant reporter cassettes, and more red-FP / luminescent families. Reuse `curate.py`
   (broaden `SEARCH_QUERIES`) and keep the exclusion set. The limiting factor is distinct
   proteins, not records — dedup by lineage, not accession.
2. **Add independent label corroboration** (profile-HMM/InterPro or structure) so families do
   not rest on depositor text; this is the same upgrade Codex's low-identity workstream needs.
3. **Add the pinned ESM-2 arm** once weights are cached, to check scorer transfer (protocol
   already supports `characterizer`).
4. **Re-freeze and re-run** `analyze.py` unchanged; the +5 target is met if the reporter
   equal-group **lower bound** clears 5.

Exact resume point: everything through scoring is committed and reproducible from the frozen
panel (`panel.json`) + fetched sources. To extend, add accessions via `curate.py`
discover/fetch/audit, rebuild with `freeze.py panel`, re-group, re-freeze, re-score. Do not
modify the frozen `7f41a8d` panel; version an extension as `panel-v2`.
