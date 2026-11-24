# Distant protein recognition with same-fold alternatives

This study addresses the characterization part of [Q3](../../docs/project/SCOPE.md) on experimentally resolved natural domains involved in ordinary cellular functions. The prior structural study showed different results across three folds. Here the structural labels come from additional folds, and the competing labels share the query's broad fold. [Previous comparison](STRUCTURAL_FAMILIES.md)

## Protocol fixed before scoring

- Take contiguous domains from curated SCOPe 2.08 structure assignments in the model-organism taxa and fine families listed in `fold_study.py`. Include domains from 55 to 620 residues, using the archived boundaries. Give priority to distinct proteins within each fine family, take at most six unique domains per fine family, and exclude exact sequences in the earlier protein panels. Both the source and the label are fixed before evaluating the scorers.
- Predict the structural superfamily. Hold out entire fine families with the same deterministic three-way assignment as the previous study. Compare only within one fold, retaining folds with at least two candidate superfamilies. Every method receives the same reference domains, including the same-label relatives and all competing labels.
- The primary group is below **25% global identity to every available reference**, counting gaps in the alignment length. Also report the 25% to below 30% tier and a check for a same-label local alignment covering at least 70% of both sequences. The two measures distinguish low whole-domain identity from absence of broadly covering local alignment.
- Use raw and reference-normalized BLOSUM62, HMMER profiles per fine family and per superfamily, and both mean and residue-local ESM-2 at 8M and 150M. Score the first-ranked label, with ties treated as missed identifications. The main paired comparisons use residue-local 150M against raw alignment and each profile baseline. Count model rescues among cases missed by all four conventional methods.
- Summarize success per domain, by fold, and as an average giving folds equal weight. Within each fold give superfamilies equal weight, then fine families within their superfamily. Report paired differences, a descriptive whole-fold bootstrap interval, and sensitivity when dropping one fold. Treat the intervals as panel summaries. Several groups share reference data and the domain sequences may have appeared in ESM-2 pretraining.

The measured label is a structural superfamily, which can contain enzymes with different functions. Evidence that a scorer can recover this label after receiving the domain is one part of Q3. The separate construct experiments address whether a useful input reaches the scorer.

## Reference-depth check fixed before running it

The first comparison holds out several fine families at once. As a check on reference availability, give every method all the other fine families within the query's fold. Keep the entire query family out. Use the same 132 queries, candidate labels and scoring rules. Recalculate identity against the larger reference set and compare decisions on both the original below-25% cohort and the cases still below 25%. Some superfamilies have only two fine families, so this adds same-label references for only some queries. This check reuses the original sample.
