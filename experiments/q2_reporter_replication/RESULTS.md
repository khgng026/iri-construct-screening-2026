# Q2 reporter replication — results

**Self-review; public, benign material; family identification only.** Protocol, panel,
references, and grouping were frozen at commit `7f41a8d` before any arm was scored. The
headline was independently recomputed from the saved verifier output. See
[PROTOCOL.md](PROTOCOL.md) for the frozen design and [HANDOFF.md](HANDOFF.md) for the
resume point.

## Answer

**The preprocessing gain replicates on fresh, independent reporter constructs, and — unlike
the published cohort — it is not an artifact of the two bacterial marker families.** At a
cap of two protein inputs, gene-calling preprocessing beat minimal six-frame preprocessing
by **+11.8 pooled / +13.0 equal-lineage-group points** (normalized BLOSUM62), matching the
cohort's +12.4 pooled / +9.0 BLOSUM62 equal-group. The gain **carries to reporter
components specifically** (reporter equal-group +15.9, point estimate larger than markers'
+7.1) and **survives removing the dominant lineage** (+13.1). The mechanism is the cohort's:
candidate **selection**, not characterization.

**But this is a pilot and the +5-point target is not certified at 95% for reporters alone.**
With 21 lineage groups (14 reporter), the whole-group bootstrap intervals are wide. The
all-component gain is +13.0 **[+2.4, +26.3]** — positive, but the lower bound is below the
5-point target; the reporter-only interval **[0.0, +34.7]** touches zero. The direction and
rough magnitude replicate; the precision does not yet confirm the effect size for reporters.

## Panel (frozen)

86 fresh public benign reporter constructs from NCBI GenBank, independent of the
777-construct cohort (792 accessions excluded) and the 5-construct pilot. 169 toolkit-
validated components: **93 reporter, 76 marker**. 9 families; 21 lineage groups (>= 90%
within-family identity): **14 reporter, 7 marker**. Archetypes: multi_cds 33,
backbone_heavy 27, single_cds 23, tagged_fusion 3. NanoLuc (13 constructs) excluded for
lack of an independent reviewed reference.

Reporter families: gfp_like (55 components, 7 lineage groups), firefly_luciferase (33, 4),
renilla_luciferase (3, 1), beta_galactosidase (2, 2). Marker families: beta_lactamase (43),
aminoglycoside_phosphotransferase (22), hygromycin_phosphotransferase (9), others (2).

## Primary table — documented-component recall at cap 2

Pooled recall / equal-lineage-group recall. Candidate pool and residues are totals over 86
constructs; the query budget is identical across arms.

| Arm | All (pooled/eq) | Reporter (pooled/eq) | Marker (pooled/eq) | Pool | Residues |
| --- | ---: | ---: | ---: | ---: | ---: |
| Minimal (baseline) | 72.8 / 63.4 | 80.6 / 65.8 | 63.2 / 58.4 | 9,012 | 66,124 |
| Start ranking (control) | 77.5 / 65.4 | 79.6 / 66.0 | 75.0 / 64.1 | 9,012 | 66,030 |
| **Gene calling** | **84.6 / 76.3** | **89.2 / 81.7** | **78.9 / 65.5** | **320** | **59,963** |
| Annotation-derived (ceiling) | 88.2 / 82.7 | 93.5 / 84.1 | 81.6 / 79.8 | 169 | 51,260 |

### Gene calling − minimal (cap 2), 95% whole-group bootstrap

| Components | Equal-group gain | Pooled gain | Groups | Meets +5 (lower bound) |
| --- | ---: | ---: | ---: | :---: |
| All (169) | **+13.0 [+2.4, +26.3]** | +11.8 [+3.1, +26.9] | 21 | no (LB +2.4) |
| Reporter (93) | **+15.9 [0.0, +34.7]** | +8.6 [0.0, +38.5] | 14 | no (LB 0.0) |
| Marker (76) | +7.1 [0.0, +18.5] | +15.8 [0.0, +34.8] | 7 | no (LB 0.0) |

## Mechanism — selection, not characterization

At cap 2, **every arm has zero "covered-but-missed" components**: whenever a covering
candidate reached the BLOSUM62 scorer, the family call was correct. All misses are
selection misses (no covering candidate within the two-input budget).

| Arm | Recovered | No covering candidate | Covered but mislabelled |
| --- | ---: | ---: | ---: |
| Minimal | 123 | 46 | 0 |
| Start ranking | 131 | 38 | 0 |
| Gene calling | 143 | 26 | 0 |
| Annotation-derived | 149 | 20 | 0 |

Gene calling moved 20 components' covering candidates into the two-input budget (46 → 26
misses) with a **28× smaller candidate pool** (320 vs 9,012) and fewer input residues at the
same query budget. Reordering the minimal pool alone (start ranking) recovered part of it
(+4.7 pooled); gene calling recovered more by generating better candidates. Characterization
never failed here because the deposited reporters and markers sit at high identity to their
family references — so the alignment scorer always wins once the locus is in budget. This is
the cohort's "selection drives the gain" finding, in starker form.

## Sensitivity

| Check | Equal-group gain (all) |
| --- | ---: |
| Primary (lineage groups, 21) | +13.0 [+2.4, +26.3] |
| Omit dominant lineage (drop beta_lactamase#1, n=42) | +13.1 [+2.0, +27.0] |
| Stricter family-group weighting (9 groups) | +6.7 [0.0, +15.4] |
| Reporter-only lineage groups (14) | +15.9 [0.0, +34.3] |
| Budget cap 1 | +10.7 [-0.5, +25.5] |
| Budget cap 4 | +10.1 [0.0, +23.0] |

Removing the dominant (marker) lineage leaves the gain essentially unchanged — the opposite
of the cohort, where two marker families carried 82% of the extra recoveries and their
removal cut the pooled gain. Here the reporter point estimate equals or exceeds the marker
one. The family-group view (each of 9 families weighted once) softens to +6.7 with a lower
bound at 0, the expected loss of power from only 9 units.

## Interpretation

1. **Direction and magnitude replicate** on an independent, reporter-diverse panel with a
   different scorer (BLOSUM62 vs the cohort's ESM-2): ~+12 pooled / +13 equal-group points.
2. **The gain carries to reporters.** Reporter equal-group gain (+15.9) is at least as large
   as the marker gain (+7.1) and survives removing the dominant lineage, so the cohort's
   marker-family concentration does not explain the effect on this panel.
3. **Same mechanism:** candidate selection. Gene calling gets the right locus into a tiny
   budget; the scorer then rarely errs at these identities.
4. **Pilot-limited precision.** 14 reporter lineage groups give wide intervals; the +5-point
   target is met by every point estimate but not by the 95% lower bound for reporters alone.
   The result confirms the effect's direction and rough size, not the certified effect size.

## Limitations

- **Depositor-annotation labels.** Families come from `/product` qualifiers (independent of
  the scorer, but weaker than the cohort's reviewed >= 90% anchor). A confirmatory study
  needs independent corroboration (profile-HMM or structure).
- **Few distinct proteins.** gfp_like#1 (39) and firefly#1 (28) dominate their families;
  renilla (3) and beta_gal (2) are thin. Equal-group weighting mitigates but cannot
  manufacture diversity.
- **Scorer and families.** BLOSUM62 only (no ESM-2 cross-check in this environment); NanoLuc
  excluded; only 3 tagged/fused constructs, so fusion behaviour is not tested here.
- Self-review; no independent specialist review.

## Reproduce

```sh
run=results/q2-reporter-replication
dir=experiments/q2_reporter_replication
# 1. (optional) re-derive the panel from NCBI — needs network; see curate.py / freeze.py.
# 2. Validate + build inputs from the frozen panel and fetched sources:
python -m experiments.paired_preprocessing.reporter_pilot prepare \
  --sources $run/sources --panel $dir/panel.json --out $run/prepared
# 3. Score the four arms with the offline BLOSUM62 characterizer (pyrodigal==3.6.3):
for m in minimal start_ranked gene_caller; do
  python -m experiments.paired_preprocessing.reporter_pilot run --method $m --protocol $dir/protocol.json \
    --inputs $run/prepared/method_inputs/constructs.fasta --references $run/references.private.json --out $run/scored/$m.json
done
python -m experiments.paired_preprocessing.reporter_pilot run --method annotation_reference --protocol $dir/protocol.json \
  --inputs $run/prepared/evaluator/annotation-units.json --references $run/references.private.json --out $run/scored/annotation_reference.json
# 4. Verify and analyze:
python -m experiments.paired_preprocessing.verify --protocol $dir/protocol.json \
  --answer-key $run/prepared/evaluator/answer-key.json --baseline $run/scored/minimal.json \
  --candidate $run/scored/start_ranked.json --candidate $run/scored/gene_caller.json \
  --annotation-reference $run/scored/annotation_reference.json --out $run/scored/report.json
python -m experiments.q2_reporter_replication.analyze analyze \
  --report $run/scored/report.json --groups $dir/groups.json --out $run/scored/analysis.json
```

Committed sequence-free artifacts: `panel.json`, `protocol.json`, `references.manifest.json`,
`groups.json`, and under `results/q2-reporter-replication/`: `candidates.json`, `audit.json`,
`fetch_manifest.json`, `prepared/evaluator/answer-key.json`, `scored/report.json`,
`scored/analysis.json`. Source and reference sequences stay in ignored storage.

## Provenance

- Code frozen at `7f41a8d`; scoring ran after. Headline recomputed from `scored/report.json`
  independently of `analyze.py` (exact match: +11.8 pooled, +13.0 equal-group).
- Versions: Python 3.13, Biopython 1.88, NumPy 2.5.5, Pyrodigal 3.6.3. Scorer: normalized
  local BLOSUM62 (open −10, extend −0.5, reference self-score normalization).
- References: 12 reviewed UniProtKB/Swiss-Prot entries across 9 families
  (`references.manifest.json`).
