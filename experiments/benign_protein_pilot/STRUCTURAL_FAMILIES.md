# Recognizing distant relatives of known metabolic protein domains

Can a protein language model recognize a distant relative when sequence comparison has few obvious matches? Wittmann and colleagues describe why this question matters in defensive screening and why performance depends on the method and the input. Here we test a narrower question with known domains of ordinary metabolic proteins and a fixed set of plausible alternatives. The labels come from the curated SCOPe structure classification. [Science, 2025](https://doi.org/10.1126/science.adu8578) · [Frontiers, 2026](https://doi.org/10.3389/fbioe.2026.1858951) · [SCOPe 2.08](https://scop.berkeley.edu/astral/ver=2.08)

## Comparison fixed before scoring

- Use curated SCOPe 2.08 TIM-barrel domains from selected carbohydrate, nucleotide, amino acid and cofactor pathways in model organisms. Use contiguous domains from the source's experimental structure records. Select at most one record per taxon and fine family, with at most four distinct sequences per fine family. Exclude exact sequence matches to the prior housekeeping panels.
- Predict the broader structural superfamily. All alternative labels share the TIM-barrel fold, so recognizing the general shape alone is insufficient. Hold out complete fine families, never splitting the members of one family between query and reference within a comparison. Every scorer receives the same reference domains. The split uses classification identifiers and requires at least two distinct fine families for each evaluated superfamily.
- Evaluate correct first-ranked superfamily primarily on queries with **less than 30% global identity to every available reference**, counting gaps in the alignment length. Divide that subset at 25% for descriptive reporting. Also record how many have a same-superfamily local alignment covering at least 70% of both domains. All cutoffs define this benchmark's comparison.
- Compare raw and reference-normalized local BLOSUM62, HMMER 3 family profiles built from the same references with FAMSA, and mean-pooled and residue-local ESM-2 8M. Build one profile per reference fine family, using HMMER's single-sequence model for singleton families. Each method ranks superfamilies by the strongest member or fine-family profile score. A cross-label tie or a missing ranking counts as a missed identification.
- Report raw hit counts, accuracy giving each held-out fine family equal weight, paired wins and losses, and the breakdown by superfamily. Interpret the results descriptively because the number of broad groups is small and some folds reuse references.

The outcome is classification of known, experimentally resolved protein domains into a fixed set of structural groups. Members of a group can perform different reactions. The study compares characterization after boundaries are supplied. It supports no estimate of construct-level recovery or the performance of a screening system.

## First experiment: TIM-barrel domains

The [panel](structural_panel.freeze.json) and protocol were committed at `de79857` before scoring. It contains 42 domains from 20 fine families and six superfamilies. The nearest reference for each query had 15.7% to 21.2% global identity, so all 42 qualify for the primary comparison and its under-25% subset.

| Scoring method | Correct first label | Accuracy weighting fine families equally |
| --- | ---: | ---: |
| Raw BLOSUM62 | 18/42 | 31.3% |
| Reference-normalized BLOSUM62 | 13/42 | 28.3% |
| HMMER fine-family profiles | 21/42 | 40.0% |
| ESM-2 8M, one vector per domain | 29/42 | 61.7% |
| ESM-2 8M, local comparison | 31/42 | 70.0% |

Local ESM-2 identified 12 domains missed by the profile method, and the converse occurred twice. It led in three superfamilies, tied in two, and trailed in one. One superfamily accounted for five of the ten extra recoveries. Omitting any one superfamily leaves an equal-fine-family advantage of about 13 to 39 percentage points. The profile method ranked five of the six domains with broadly covering local alignments first, versus six for local ESM-2. Among the other 36 domains the counts were 16 and 25.

This is evidence that a model can help classify structurally distant domains on a panel with little reference data. The profile method used only the reference domains available to the other scorers, and many fine-family profiles used one or two examples. Domain structures may also have been represented in ESM-2 pretraining. The difference between the two ESM-2 methods was much smaller than their difference from the profile result.

## Follow-up fixed before scoring additional folds

- Strengthen the profile comparison by also fitting a single HMMER model per superfamily from all available reference domains in that label. Keep the fine-family profile as a separate arm. This is a sensitivity analysis for the TIM-barrel panel.
- Test the same scoring methods on independently selected Flavodoxin-like (`c.23`) and adenine nucleotide alpha hydrolase-like (`c.26`) domains. Use the same curation, holdout, identity, and scoring rules. Each panel has its own references and competitors. Report each fold separately because its number of possible labels differs from the first panel.
- Report paired errors and the local alignment coverage check in all three panels. Treat the comparison as a bounded test of structural classification. The papers concern different materials and tasks, including shorter inputs in the 2026 study. Performance on those tasks requires separate evidence.
