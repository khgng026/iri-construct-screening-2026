# Project research report (general audience)

`construct_screening_report.pdf` is a 10-page, plain-language summary of the completed
work: the construct-preprocessing comparison (cohort v2) and the exploratory
alignment-vs-model scorer follow-up. It is written for readers with little background
and assumes no biology or machine-learning knowledge.

It is a **derived summary**, not a primary record. The authoritative results, protocols,
denominators, uncertainty, and provenance live in:

- [`experiments/paired_preprocessing/COHORT_RESULTS.md`](../../experiments/paired_preprocessing/COHORT_RESULTS.md) — the matched-budget preprocessing comparison (Q1/Q2).
- [`experiments/paired_preprocessing/COHORT_SCORER_FOLLOWUP.md`](../../experiments/paired_preprocessing/COHORT_SCORER_FOLLOWUP.md) — the identity crossover (Q3, exploratory).
- [`docs/project/STATUS.md`](../project/STATUS.md) and [`docs/project/SCOPE.md`](../project/SCOPE.md) — current status and the research questions.

If any figure in the PDF disagrees with those files, **those files are correct.**

## What the figures are drawn from

- `fig_crossover.png` is computed from the committed, sequence-free result file
  `results/paired-cohort-v2-scorers/ablation.crossover.json` (equal-group recall by
  nearest-reference identity threshold, 150M scorers).
- `fig_gain.png` reproduces the budget-curve values quoted in `COHORT_RESULTS.md`
  (pooled recall, ESM-2 arm, caps 1 / 2 / 4, plus the annotation-derived ceiling at cap 2).

## Regenerate

Requires `matplotlib` and `reportlab` (both pip-installable). From the repository root:

```bash
python docs/report/generate_figures.py   # writes fig_gain.png, fig_crossover.png
python docs/report/build_report.py        # writes construct_screening_report.pdf
```

The build embeds the bundled DejaVu Sans font shipped with matplotlib, so the PDF renders
the same across machines.

## Status and boundaries

Self-review only; not peer reviewed. Public, benign material only; defensive purpose. The
report states plainly that the headline preprocessing gain is concentrated in two protein
families and reflects candidate **selection** rather than characterization, and that the
low-similarity (Q3) crossover is **exploratory** — it reuses the now-public holdout keys
and alignment-derived labels — so the project's low-similarity goal remains formally open
pending the confirmatory study described in
[`docs/project/ALIGNMENT_VS_MODELS_PLAN.md`](../project/ALIGNMENT_VS_MODELS_PLAN.md).
Identifying a protein family is a measurable result; it does not by itself establish
biological activity, function, or risk.
