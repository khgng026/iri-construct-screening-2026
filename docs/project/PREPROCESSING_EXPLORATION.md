# Preprocessing exploration

The strongest change in the historical cohort was reordering candidates at a
tight budget. The clearest gain on a specific archived record came from
retaining resolved translations in an incomplete plasmid map. The whole-DNA
comparison also showed that direct translated search can clear the historical
six-frame baseline without a gene caller. Each result addresses a different
part of the proposal’s question about improving the inputs a scorer receives.

## What changed candidate selection

All historical tests reuse the 777-construct, 1,221-component cohort and its
public answer key. The selection rules receive candidate geometry and
metadata. Evaluator labels stay out of the policies. The identification
replay can select only from candidates that already have saved rankings.
The independent reach calculation can use the full candidate pools.

| Change | Result on the historical cohort | Trade-off |
| --- | --- | --- |
| Prefer candidates near 250 residues | At cap 2, mean-pooled ESM-2 recovered **926** components from gene-called candidates versus **725** in the recorded order. Lineage-weighted recall rose from **65.7% to 76.2%**. BLOSUM62 rose from **746 to 933** correct calls. | ESM-2 inputs fell from **551,358 to 415,989 residues**. At cap 4, ESM-2 recovered **1,095 versus 1,099**. The gain is concentrated at tight budgets on familiar historical data. |
| Combine modest length preference with locus diversity | At cap 2, ESM-2 recovered **753 versus 725** from gene-called candidates. | A smaller effect on the same retrospective test. |
| Mix gene-called and six-frame candidates | Keeping the top gene calls and filling spare slots preserved the cap-2 result. Reserving a slot for six-frame candidates reduced cap-2 recall. At cap 8, a distinct six-frame fallback raised coordinate coverage from **1,175 to 1,207**. | At cap 8 it used **6,216 versus 3,932 candidates**, and those additional candidates have no saved scorer outputs. The combined full pool covers **1,211/1,221**, versus **1,176/1,221** for the gene caller. |
| Reallocate a shared budget between constructs | At the cap-3 total budget, prioritizing the longest next available candidate raised ESM-2 recovery from **1,016 to 1,029**, with a **5.7-point** lineage-weighted gain (95% CI 3.9–7.6). | It used **2,072 versus 2,206 candidates**, with nearly identical total residues. None of the tested allocation rules improved ESM-2 at the cap-2 budget. |
| Estimate the allocation ceiling | An evaluator-only allocation of the existing gene-call prefixes could recover **1,000 versus 725** with ESM-2 under the cap-2 total limits. | The optimization uses the answer key. It measures headroom in the existing pool and supplies no inference policy. |

The frame diagnostic also found that same-strand coordinate overlap can
overstate usable reach. In the 1,213 single-interval components, at cap 2
minimal translation overlapped **609** loci and matched the expected frame
at **597**. Gene calling overlapped **751** and matched the expected frame
at **750**. Frame agreement is a coordinate-based check. The real-archive
evaluation below directly checks the documented peptide.

The 250-residue gain concentrated in the high and intermediate tiers. In
the original low-similarity tier it gained **5** recoveries (**137 to 142**)
and in the separate insufficient-coverage tier it lost **26** (**34 to 8**).
A post-hoc mixed rule retained the first recorded gene call and preferred
shorter candidates after it. It recovered **804** components overall,
including **152** in the low tier and **27** in the insufficient-coverage
tier. The gain depends on which components matter and how the remaining
candidate slots are allocated.

## What transferred to other records

**Open Yeast.** Of 17 archived vectors, 11 selected CDSs passed complete
translation validation. At cap 2, gene calling and start-ranked six-frame
translation each supplied the full peptide for **10/11**, versus **9/11**
with minimal translation. At cap 4 the original methods reached **10/11**.
Favoring 250-residue candidates added no recoveries in the gene-caller arm
and reduced the two-candidate minimal arm from **9/11 to 8/11**.
The archive shares a cloning background and offers a bounded transfer check.
The remaining target is **1,392 amino acids**, above the 1,022-residue
input ceiling. Neither method can fit it in one ordinary candidate. An
evaluator-side comparison found all **1,392/1,392** exact residues across
the first two windows of either method. That distinguishes a windowing
and scoring question from a failure to generate the protein region.

**LAMPRA reporter.** Minimal translation, start-aware translation and gene
calling all missed the documented reporter even in their full candidate
pools. Its CDS has three **adjacent** annotated intervals. The map includes
unresolved bases outside the target. Keeping resolved translated segments
recovered the **entire exact reported peptide** by cap 4, with no change
in the Open Yeast results. The reference-guided miniprot run localized
**99.6%** of the coding region but proposed a **238-residue** translation
for the documented **239-residue** protein. This is one example where
handling an incomplete archival input mattered more than ranking candidates.

## Compare to a direct whole-DNA search

The Open Reporters test used **19 independently corroborated GFP-like
reporters**, five catalogue-labeled other reporters, and the same seven
natural GFP references across all search arms. The loci and strand were
checked after scoring. Each number below is the count of the 19 positive
reporters recovered by the top reported regional match.

| Search | 80% coding-region overlap | 50% coding-region overlap |
| --- | ---: | ---: |
| Whole-DNA BLASTx | 13 | 19 |
| Gene calling, then BLASTp at cap 2 | 12 | 19 |
| Minimal translation, then BLASTp at cap 2 | 8 | 12 |
| Minimal translation, then BLASTp at cap 4 | 12 | 19 |
| Whole-DNA DIAMOND blastx, default | 2 | 4 |
| Whole-DNA DIAMOND blastx, very-sensitive | 12 | 19 |

Whole-DNA search uses its internal translation and search across the entire
record. The explicit preprocessing arms additionally cap the number of
proteins searched. A shared cap on reported regions helps compare output
quality, while runtime and input work give separate cost measures. This
small family-specific panel gives no general false-assignment estimate.
Existing workflows such as pLannotate use direct translated search and
assess a broader reference library. [pLannotate](https://doi.org/10.1093/nar/gkab374)

On the same archive, the historical length preference improved the very
tight cap. Gene calling's top candidate covered **0/19** reporters in the
recorded length order and **19/19** after ranking near 250 residues.
BLASTp recovered **12/19** at the stricter overlap from that one candidate,
the same count as the original two-candidate gene-caller run. This archive
contains one class of relatively short reporter peptides in a shared
backbone. Its results complement the mixed Open Yeast findings.

## Interpretation

Retain a simple gene-caller baseline and evaluate ranking at several
operating points. Length ranking is worth a fresh, family-balanced test.
The cohort result alone cannot justify a universal length preference.
Keep full-pool reach and scored identification separate. The fallback and
allocation runs show substantial candidate-selection headroom, but they
also show how easily a plausible rule can trade away known components.
Archive features such as ambiguous sequence require explicit validation
of the peptide rather than coordinate coverage alone.

These experiments use archived benign sequences. The historical ESM-2
comparisons reuse stored predictions. The new sequence runs measure
candidate reach or conventional alignment, and the work leaves the
stronger prospective PLM comparison open.

## Reproduction

The separate modules print or write sequence-free summaries. Run from this
repository root with the project's normal Python environment. The historical
allocation oracle uses SciPy. The fresh archive and search runs used
Biopython 1.86 and Pyrodigal 3.6.3. The source checkouts must contain the
commits in [the archive manifest](../../experiments/preprocessing_exploration/archive_coverage_panel.json)
and [the translated-search protocol](../../experiments/preprocessing_exploration/TRANSLATED_SEARCH_PROTOCOL.md).

```bash
run=results/preprocessing-replay
mkdir -p "$run"
python -m experiments.preprocessing_exploration.ranking --output "$run/ranking.json"
python -m experiments.preprocessing_exploration.ranking --split holdout \
  --policy recorded --policy length_center_250 --policy keep_first_then_250 \
  --output "$run/rank-mix.json"
python -m experiments.preprocessing_exploration.hybrid --out "$run/hybrid.json"
python -m experiments.preprocessing_exploration.budget_allocate --out "$run/budget.json"
python -m experiments.preprocessing_exploration.frame_audit --out "$run/frame.json"
python -m experiments.preprocessing_exploration.archive_coverage \
  --reclone /path/to/reclone --lampra /path/to/LAMPRA \
  --ranking-followup --output "$run/archive.json"

git fetch origin refs/heads/study/q3-evaluation
git show 74161519082af1d7e124eaeb70bb3efb27fab67c:results/q3-followup/open-reporters.json > "$run/prior-audit.json"
python -m experiments.preprocessing_exploration.translated_search \
  --audit "$run/prior-audit.json" --construct-repo /path/to/reclone \
  --protein-repo /path/to/Orthologous_GFP_Fitness_Peaks \
  --blast-bin /path/to/blast-bin --diamond /path/to/diamond \
  --diamond-sensitivity very-sensitive --ranking-followup \
  --work "$run/search-work" --out "$run/search.json"
python -m experiments.preprocessing_exploration.compound_reporter \
  --lampra /path/to/LAMPRA --protein-repo /path/to/Orthologous_GFP_Fitness_Peaks \
  --miniprot /path/to/miniprot --work "$run/compound-work" --out "$run/compound.json"
```

The native executables were BLAST+ 2.16.0, DIAMOND 2.1.11, and miniprot
0.18-r281. The LAMPRA GenBank map produces parser warnings on unrelated
features. The selected reporter annotation validates against its translation.
Keep FASTA, GenBank and work directories out of the committed result set.
