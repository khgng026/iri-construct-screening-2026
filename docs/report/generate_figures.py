#!/usr/bin/env python3
"""Generate the two general-audience figures for the project report.

Figures are drawn from committed, reproducible result files:
  - results/paired-cohort-v2-scorers/ablation.crossover.json  (identity crossover)
  - budget-curve values quoted in experiments/paired_preprocessing/COHORT_RESULTS.md

Run from the repository root:
    python docs/report/generate_figures.py

Outputs fig_gain.png and fig_crossover.png next to this script.
Requires matplotlib.
"""
import json
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
ABLATION = os.path.join(ROOT, "results", "paired-cohort-v2-scorers", "ablation.crossover.json")

plt.rcParams.update({
    "font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13,
    "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
})
INK = "#1b2430"; MUT = "#5b6b7b"
C_ALIGN = "#2f6f8f"   # alignment / conventional
C_MODEL = "#c05a2e"   # foundation model
C_ANN = "#7a8a99"     # annotation reference


def fig_gain():
    """Budget curve: minimal vs gene-calling preprocessing (pooled recall, ESM-2 arm)."""
    caps = [1, 2, 4]
    minimal = [22.9, 46.9, 83.0]
    pyrodigal = [26.5, 59.4, 90.0]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(caps, minimal, "-o", color=C_MODEL, lw=2.4, ms=7, label="Minimal preprocessing (baseline)")
    ax.plot(caps, pyrodigal, "-o", color=C_ALIGN, lw=2.4, ms=7, label="Gene-calling preprocessing")
    ax.scatter([2], [91.4], color=C_ANN, s=90, zorder=5, marker="D",
               label="Using the hidden answer key (ceiling)")
    ax.annotate("+12.4 points\nat the same budget", xy=(2, 59.4), xytext=(2.25, 41),
                color=INK, fontsize=11.5, arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.4))
    for x, y in zip(caps, minimal):
        ax.annotate(f"{y:.0f}%", (x, y), textcoords="offset points", xytext=(0, -18), ha="center", color=C_MODEL, fontsize=11)
    for x, y in zip(caps, pyrodigal):
        ax.annotate(f"{y:.0f}%", (x, y), textcoords="offset points", xytext=(0, 10), ha="center", color=C_ALIGN, fontsize=11)
    ax.annotate("91%", (2, 91.4), textcoords="offset points", xytext=(14, 2), ha="left", color=MUT, fontsize=11)
    ax.set_xticks(caps)
    ax.set_xlabel("Screening budget  (protein units examined per construct)")
    ax.set_ylabel("Documented components correctly identified")
    ax.set_ylim(0, 100); ax.set_xlim(0.7, 4.4)
    ax.set_title("Better preprocessing finds more of a construct's\nknown parts for the same amount of work", color=INK, loc="left")
    ax.legend(frameon=False, fontsize=10.5, loc="lower right")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_gain.png"), dpi=200); plt.close(fig)


def fig_crossover():
    """Identity crossover: equal-group recall vs nearest-reference identity threshold (150M scorers)."""
    d = json.load(open(ABLATION))
    order = ["0.25", "0.30", "0.35", "0.40", "0.50", "0.60", "0.70"]
    xs = [25, 30, 35, 40, 50, 60, 70]

    def series(key):
        out = []
        for t in order:
            r = d["by_threshold"][t]["recall"]
            out.append(100 * r[key]["equal_group_recall"] if key in r and r[key] else None)
        return out

    blos = series("blosum62"); res150 = series("esm2_150m_residue"); mean150 = series("esm2_150m_mean")
    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    ax.axvspan(20, 30, color="#f0d9cf", alpha=0.55, zorder=0)
    ax.plot(xs, blos, "-o", color=C_ALIGN, lw=2.6, ms=7, label="Conventional sequence matching")
    ax.plot(xs, res150, "-o", color=C_MODEL, lw=2.6, ms=7, label="AI model, read residue by residue")
    ax.plot(xs, mean150, "--s", color="#d9a46a", lw=2.0, ms=6, label="AI model, whole-protein summary")
    ax.annotate("68%", (30, 68.2), textcoords="offset points", xytext=(16, -6), ha="left", color=C_ALIGN, fontsize=12, fontweight="bold")
    ax.annotate("96%", (30, 95.6), textcoords="offset points", xytext=(-2, 9), ha="center", color=C_MODEL, fontsize=12, fontweight="bold")
    ax.text(26, 86, "Hardest cases:\nno close known relative", color="#8a5038", fontsize=10.5, va="center", ha="center")
    ax.set_xlabel("Similarity to the nearest known reference  (% — lower is harder)")
    ax.set_ylabel("Protein families correctly identified")
    ax.set_ylim(55, 103); ax.set_xlim(23, 72); ax.invert_xaxis()
    ax.set_title("With no close relative on file, the AI model pulls ahead of\nconventional matching — but only when read residue by residue", color=INK, loc="left", fontsize=13.5)
    ax.legend(frameon=False, fontsize=10.5, loc="lower left")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "fig_crossover.png"), dpi=200); plt.close(fig)


if __name__ == "__main__":
    fig_gain()
    fig_crossover()
    print("Wrote fig_gain.png and fig_crossover.png to", HERE)
