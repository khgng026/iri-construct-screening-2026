# Research design: when do foundation models beat sequence alignment, given a preprocessing toolkit?

Draft 23 November 2026. A forward-looking design, not a completed study. It defines the
experiment that would give a more definitive answer than the [cohort v2 result](../../experiments/paired_preprocessing/COHORT_RESULTS.md),
which found mean-pooled ESM-2 8M no better than normalized BLOSUM62 at family identification.
All work here is benign and defensive; see the safety boundary below.

## Why the cohort v2 result is not the final word

Cohort v2 compared mean-pooled ESM-2 8M with normalized BLOSUM62 on 1,221 documented
components and found alignment slightly ahead, even in its "low-similarity" tier. Three
features of that study keep it from answering the general question:

1. **It did not reach the regime where alignment is known to fail.** The low-similarity
   tier was 25–40% identity to the reference panel, in 28 well-separated families, with no
   confusable competitor. The literature locates alignment's failure *below* this band.
2. **Labels were alignment-derived.** Eligibility required a reviewed anchor at ≥ 90%
   identity, and the tier was defined by identity to the scored references. Using alignment
   to define the ground truth tilts an alignment-vs-model comparison toward alignment.
3. **The model configuration was the weakest reasonable one.** Mean pooling discards the
   local signal, and 8M is the smallest ESM-2. Residue-level scoring and larger checkpoints
   are the configurations the field actually proposes.

## What the biosecurity literature says about where alignment breaks down

Two studies from the same group bound the failure regime. They concern screening proteins
of concern; we take from them only the *methodological* boundary, and work entirely on
benign proxies (safety boundary below).

**Wittmann et al., _Science_ 2025 (adu8578), "Strengthening nucleic acid biosecurity
screening against generative protein design tools."** AI protein-design tools (ProteinMPNN,
EvoDiff) generate "synthetic homologs" that keep predicted function (TM-score > 0.5,
ΔpLDDT near 0) while dropping to low sequence identity. Best-match screening, including
HMM-based and neural-enhanced alignment, misses functional variants as identity falls, and
DNA obfuscation (fragmenting and shuffling) worsens detection. The authors note current
design tools are "largely incapable of producing the extremely divergent (< 30% sequence
identity) sequences that would necessitate" non-sequence methods, and point to "homology
detection through high-dimensional learned embedding spaces" as the needed alternative.

**Wittmann et al., _Front. Bioeng. Biotechnol._ 2026 (1858951), "The limits of
sequence-based biosecurity screening tools in the age of AI-assisted protein design."**
Fragment-level analysis (fragments 25–200 nt) over ~700k fragments of 5,505 benign-proxy
and concern sequences. The decisive quantitative result: flagging rate falls to near zero
once fragment identity drops below **~30%** *and* the longest common subsequence (LCS) with
the wild type is below **~20 amino acids**; one tool held out until identity < 50% with
LCS < 20. Shorter fragments fail sooner. The authors state that below this boundary,
in-silico structure (a higher-resolution function proxy than sequence) or "direct
prediction of sequence function may be necessary," and that structure "will not work for
shorter oligonucleotide sequences in isolation." Ground truth was a function proxy
(predicted structure/TM-score), independent of sequence alignment.

**The boundary we must cross.** Alignment is reliable above ~30% identity and fails below
it, worst for short fragments and for divergent-but-functional sequences. A study that
stays above 30%, or labels by alignment, cannot distinguish the methods where it matters.

## The decisive question

> Given a preprocessing toolkit that turns an unannotated construct into candidate protein
> units, does a foundation model identify a unit's documented family **below ~30% identity
> to any available reference** better than the strongest conventional sequence search, at a
> matched query budget — and at what identity does the better method change?

The headline is a **crossover identity**: the nearest-reference identity below which the
foundation model's recall exceeds conventional search, with a confidence interval, or a
statement that no crossover is observed down to the lowest identity the benign data reach.

## Design principles that make the answer definitive

1. **Labels independent of the method under test.** Family membership from curated HMM/
   profile evidence (Pfam/InterPro) or structure, never from pairwise alignment to the
   scored reference panel. This removes the circularity in cohort v2. The toolkit still
   verifies that the translation matches the DNA; the family label comes from the HMM.
2. **Identity-stratified reporting, pre-registered bands.** Recall as a function of the
   query's identity to its nearest same-family reference: < 20%, 20–30%, 30–40%, 40–60%,
   ≥ 60%. The comparison lives in the < 30% bands. Group-weighted, with whole-group
   bootstrap intervals, as in cohort v2.
3. **A reference-availability control to manufacture the low-identity regime honestly.**
   For each benign query, remove references within an identity threshold (leave-close-
   homologs-out), forcing the < 30% condition from naturally occurring benign homologs.
   Report recall as a function of the retained best-reference identity. This reaches the
   failure regime without designing divergent sequences.
4. **The strongest realistic alignment comparator, not only BLOSUM62.** Conventional
   screening uses profile methods. Add a profile-HMM / profile-search baseline (e.g.
   `pyhmmer` phmmer/hmmscan, or a PSSM) built from the family, which is far more sensitive
   than pairwise BLOSUM62 in the twilight zone. Keep pairwise BLOSUM62 as the weak baseline
   and the cohort-v2 link. A win over BLOSUM62 but not over a profile HMM is not a win.
5. **The strongest reasonable model configuration.** Residue-level ESM-2 (local alignment
   over embeddings, already implemented) and a larger checkpoint (150M now; 650M fits the
   1 GiB weight budget and is a stretch goal). Mean-pooled 8M stays as the cohort-v2 link.
6. **A fragmentation axis.** Evaluate on full units and on fragments (down to the papers'
   25–50 nt regime). This is where sequence methods fail hardest and where a preprocessing
   toolkit that reconstructs or scores units could matter most. It connects the question to
   the construct-preprocessing value proposition.
7. **Matched budget and honest cost.** Same query cap, reference panel and output limit per
   method. Report encoder calls, profile evaluations, alignment comparisons and time, as in
   cohort v2. A learned method that wins only with far more compute is reported as such.

## Safety boundary

This program is benign and defensive. It uses only public benign proteins — housekeeping
enzymes, reporters, selection markers and their natural remote homologs — with family
labels from curated public databases. It does **not** generate synthetic homologs of
toxins or other proteins of concern, does not design sequences to evade screening, and does
not reconstruct obfuscated hazardous sequences. The "divergent" regime is reached by
withholding close references from naturally occurring benign homologs, not by designing new
divergent sequences. The question is methodological — which characterization method
recovers a known benign family at low identity — and the result informs defensive screening
without producing an operational evasion method. Any extension toward hazard proxies would
need explicit human review and is out of scope here.

## Milestones, ordered by feasibility

**P0 — Stronger model configurations on the existing holdout (running now).** Re-score the
locked cohort-v2 selected inputs with residue-level ESM-2 8M, mean-pooled ESM-2 150M and
residue-level ESM-2 150M. Does a stronger configuration change the alignment-vs-model
verdict, overall and in the existing low-similarity tier? Exploratory (the holdout keys are
public), but the first direct evidence. Code committed at `072aa30`.

**P1 — Independent labels and a profile-HMM comparator, identity-stratified.** Re-label the
existing cohort components by Pfam/InterPro HMM membership rather than anchor alignment, add
a profile-HMM search baseline, and report recall vs nearest-reference identity for every
method. This reuses the frozen constructs and reaches ~25–30% honestly. Deliverable: the
identity-stratified table and the first crossover estimate, or its absence, at matched
budget. New dependency: `pyhmmer` (pure-Python/C, installable).

**P2 — A benign remote-homolog panel reaching < 25% identity.** Assemble benign families
with naturally occurring members spanning full identity range, label by HMM/structure, and
apply the leave-close-homologs-out control to force < 25% conditions. This is the panel the
papers' regime requires and that cohort v2 lacked. Deliverable: the crossover identity with
an interval, on data built for the question.

**P3 — Fragmentation stress test.** For the P2 panel, recall vs fragment length (25–200 nt)
for each method, with and without the preprocessing toolkit's gene calling and segmentation.
Deliverable: whether preprocessing plus a foundation model recovers families from short
benign fragments where alignment has dropped to near zero.

## What would count as a conclusion

- **Foundation models help where alignment fails** if, below ~30% identity, a model
  configuration beats *both* BLOSUM62 and the profile-HMM baseline at matched budget, with a
  group-bootstrap lower bound above zero, on independently labelled benign data, and the
  effect strengthens as identity falls.
- **They do not** if the profile-HMM baseline matches or beats every model configuration
  across the < 30% bands. Cohort v2 is weak evidence in this direction; it would be
  strengthened by P1's fair baseline and independent labels.
- **Inconclusive** if the benign data cannot reach < 30% with enough distinct families and
  groups, or if pretraining contamination cannot be ruled out. State which, and what data
  would resolve it.

## Threats to validity

- **Pretraining contamination.** ESM-2 saw many of these families. A model "win" may reflect
  memorization, not generalization. Mitigate by reporting per-family results and by favoring
  families and divergence ranges least likely to be well represented; note it as a standing
  caveat.
- **Label leakage.** HMM labels and the profile-HMM comparator share lineage. Keep the
  label HMM and the comparator profile from different sources, or hold out the query's own
  subfamily when building the comparator.
- **Few distinct proteins per low-identity band.** Cohort v2's low tier had ~100 distinct
  proteins across 154 groups, half GFP. P2 must add genuinely distinct families.
- **Budget confound.** Residue-level scoring and larger models cost more. Report matched
  query budgets and the actual compute, and separate "better representation" from "more
  compute."
