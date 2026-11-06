# T05 expansion: query composition with ESM-2 8M

**Observed result:** combining proteins usually reduced retrieval of each component,
especially components occupying a small fraction of the query. Identity and order
also mattered. This is a completed local digital pilot, not a generalization test.
The original five-reference experiment and its outputs are preserved.

## Fixed experiment

Selection was committed in `ad79bf4`; reviewed inputs, every pair, exclusions and
repeatability subset were frozen in `a4ccad0` **before expanded scoring**. There are
105 reviewed UniProtKB/Swiss-Prot accessions from baker's-yeast S288c, representing
101 unique canonical sequences. Sixteen targets span 67–642 aa and ordinary
ribosomal, chromatin, redox, metabolic and protein-maintenance functions.

All 120 unordered target pairs were considered. Five exceeded the existing
1,022-residue policy and were excluded without truncation: PGK1–SSA1 (1,058),
ENO1–SSA1 (1,079), HXK2–SSA1 (1,128), HXK2–HSP60 (1,058), SSA1–HSP60 (1,214 aa).
The **115 admitted pairs** produce **230 ordered concatenations**, 170–1,009 aa,
plus 16 clean queries. Each combined query contributes two component observations:
**460 whole-query observations**, with 460 separately encoded oracle controls.
These are 16 proteins reused across pairs, not 460 independent biological samples.

Full canonical sequences are joined directly, unchanged, without a linker. These
are artificial strings, not observed or validated fusion proteins. The original
encoder produces a 320-number vector for each residue, then averages actual residues
only, excluding padding and special tokens. Both true components are evaluated
against the same complete cosine ranking. The oracle slices each known interval
and independently re-encodes it; it is not automatic boundary discovery.

Model remains `facebook/esm2_t6_8M_UR50D`, revision
`c731040fcd8d73dceaa04b0a8e6329b345b0f5df`; CPU float32, four threads, seed 0,
eval and inference mode. No tuning, training or model change. Conceptual attribution:
Nishit Kalani's supplied presentation, methods pp. 11–14. His larger-model numerical
results have **not** been reproduced, and his code was not used.

## Results

All ranks below use the **105-accession catalogue**. Positive rank change means worse.

| Arm | Component observations | Recall@1 | Recall@5 | Recall@10 | Mean reciprocal rank | Mean rank |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Component alone | 16 | 100% | 100% | 100% | 1.000 | 1.00 |
| Whole concatenation | 460 | 33.5% | 52.4% | 59.8% | 0.428 | 22.88 |
| Extracted, re-encoded oracle | 460 | 100% | 100% | 100% | 1.000 | 1.00 |

Relative to clean queries, **306 ranks worsened, 154 stayed at 1, and none improved**;
improvement was impossible from the rank-1 baseline. Median whole rank was 4, range
1–105. All 460 whole-query component cosines fell; mean cosine was 0.944300,
mean change −0.055700. Oracle vectors/scores matched their clean counterparts.
There were **zero clean-query failures**, retained as an explicit output field.
Tiny clean/oracle scores above 1 are floating-point roundoff, not probabilities.

For a two-component query, recall@k is the number of present components retrieved
divided by two, so maximum recall@1 is 0.5. At least one component was top-ranked
in **154/230 queries (67.0%)**. Both were in the top 5 in **19/230 (8.3%)**, and
both in the top 10 in **48/230 (20.9%)**. A hit on the other present component is
valid. Of the 76 top hits outside the pair, 25 shared ≥50% global sequence identity
with a present component under the documented alignment convention.

| Component fraction | Instances | Median rank | Recall@5 | Recall@10 |
| --- | ---: | ---: | ---: | ---: |
| <25% | 68 | 61 | 0% | 0% |
| 25–<40% | 98 | 36.5 | 3.1% | 16.3% |
| 40–<60% | 128 | 4 | 59.4% | 74.2% |
| 60–<75% | 96 | 1 | 95.8% | 97.9% |
| ≥75% | 70 | 1 | 100% | 100% |

Fraction versus rank had descriptive Spearman correlation −0.856. This is an
association, not a causal estimate: protein identity, partner and length covary.
Equal weighting of the 16 targets gave recall@5 **54.4%**, versus 52.4% across all
instances; length exclusions leave different numbers of partners per target.

Reversal comparisons hold proteins and fraction fixed. Across 230 protein/pair
comparisons, C-terminal placement ranked better in 47, N-terminal in 67, and 116
were unchanged. Mean absolute rank movement was 2.29; mean signed C-minus-N change
was only +0.013. There is no consistently preferred position. Examples illustrate
effects hidden by averages:

- TPI1 with TDH3 improved from rank **54 to 18** when moved from first to second.
- TRX1 with RPS28A ranked **12 or 14**, despite TRX1 occupying **60.6%** of the query.
  ENO1 with HXK2 ranked **1 in both orders** at just **47.3%**.
- SOD1 with HHT1 moved from rank **4 to 20** even as its cosine rose from
  **0.960761 to 0.965068**. Rank depends on every competing score, not just its own.

## Catalogue, ties and provenance

This was not an almost entirely unrelated distractor set: 58/105 entries and 10/16
targets have another catalogue entry at ≥50% global identity; 28 entries have a
neighbour at ≥90%. All 5,460 pairwise identities are saved. These alignment
thresholds describe similarity; they do not certify homology or function. Six
targets lack a ≥30% neighbour, and one yeast species remains a narrow sample.

Four duplicate-sequence pairs (RPS30A/B, RPL23A/B, RPL12A/B, RPS18A/B) create exact
score ties in every ranking. Accessions break exact ties alphabetically; complete
tie counts and component rank intervals are saved. No target had an exact tie;
one PRE1 observation had another score within 1e-6 (ranks 25–26). Targets have no
duplicate sequences, so duplicate-equivalent hit metrics match accession metrics.
A labelled **post hoc sensitivity** collapsing duplicate sequences to 101 entries
improved the numerical rank of 162/460 observations by up to four places, but changed
none of the recall@1/@5/@10 hit counts. Primary results remain the frozen catalogue.

Every admitted record has a source URL, entry/sequence versions, original response
and extracted-record hashes, canonical-sequence hash, review rationale and licence
evidence in `data/manifests/t05_expanded_catalogue.json`. UniProt Consortium data
are used under CC BY 4.0, retaining its other-rights disclaimer. Review was by Codex
under the user's T05 instruction, not institutional sign-off. Uncertain functional
annotations remain qualified. One alternative isoform, P17695-2, was excluded;
HHT1/HHT2 and HHF1/HHF2 each resolve to one accession. No post-score substitutions.
Sequences and source annotations remain in ignored local storage.

## Runtime and verification

Apple M4, 16 GiB RAM. Python 3.12.14; torch 2.15.0, transformers 5.22.0, numpy 2.5.4.
Original numerical packages are unchanged; additional plotting/alignment packages
are separately locked. No new model weights were downloaded.

| Timing | First full run | Full cached rerun |
| --- | ---: | ---: |
| Model setup/loading (checkpoint portion) | 3.102 s (0.076 s) | 1.560 s (0.077 s) |
| Catalogue encoding/cache write | 1.464 s | 0 s; cache read 0.002 s |
| Query encoding, including clean and oracle | 11.342 s | 10.713 s |
| Cosine ranking | 0.294 s | 0.280 s |
| Total process work | 17.633 s | 14.024 s |
| Peak process resident memory | 481 MiB | 475 MiB |

Totals include numerical imports and output work but exclude source acquisition,
shell startup and plotting. The first run overlapped a 5.36-second sequence-alignment
audit; these are descriptive timings, not a hardware speed benchmark. The isolated
fresh-process subset, including catalogue re-encoding, took 4.021 s and peaked at
511 MiB. The initial three-input timing check peaked at 384 MiB.

All 28 tests at the frozen implementation commit passed; the final suite has
**30 passing tests**, including two additional descriptive-analysis checks. Real padding checks passed
(max difference 3.43e-7); independent NumPy float64 recalculation checked all
**74,130 scores** (max difference 5.26e-7, tolerance 1e-6). The prespecified fresh
process repeated 40 evaluations/4,200 scores and all 105 catalogue vectors exactly.
The full cache-hit rerun reproduced all 706 evaluations/74,130 scores and 707 saved
arrays exactly. No inference failures or truncations. Acquisition encountered a
UniProt timeout and an EBI request-length error, resolved before freezing without
changing candidates. Final output/storage checks verified 31 original files
unchanged, all ranking rows and recall denominators, exact oracle/clean vector
agreement, ignored storage, and zero full sequences in Git candidate files.

## What this establishes, and next step

The expanded pipeline works and exposes composition-dependent retrieval changes
that the original one-target pilot could not characterize. It does not establish
general performance: clean and oracle inputs are exact catalogue entries, training
overlap is unknown, all pairs reuse 16 selected proteins, and artificial adjacency
is not biological construct validation. No function prediction, automatic
segmentation or reproduction of Nishit's results is claimed. Re-encoding extracted
components changes their context as well as the vectors being averaged.

**Next justified experiment:** keep this model and these frozen pairs, and compare
whole-query embeddings against length-weighted averages of separately encoded
components, plus component-only pooling from the whole-query hidden states. That
would help separate averaging dilution from contextual changes before changing
models. Treat family-held-out generalization as a later, separate question.

## Files

[Rerun instructions](T05_expanded_rerun.md); [frozen design](T05_expanded_design.md).
Numerical evidence and PNG/PDF/SVG figures:
`results/t05-expanded/expanded-01/` and its `analysis/` folder. Figure names:
`fraction_retrieval`, `target_partner_heatmap`, `position_and_identity`.
All per-target/bin/position tables, reversal comparisons, duplicate sensitivity and
complete rankings are saved there. Original results remain in `results/t05/`.

Official sources: [ESM-2 model card](https://huggingface.co/facebook/esm2_t6_8M_UR50D),
[UniProt licence](https://rest.uniprot.org/help/license),
[EMBL-EBI Proteins API](https://www.ebi.ac.uk/proteins/api/doc/),
[Biopython alignment documentation](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html).
