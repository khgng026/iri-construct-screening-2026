# Identifying the exact reaction of a distant metabolic enzyme

The structural comparisons found a benefit for protein language models on two folds and a conventional alignment advantage on a third. This follow-up asks whether a model can also identify a documented biochemical reaction when given an entire known enzyme and only distant examples of that reaction. [Structure study](STRUCTURAL_FAMILIES.md)

## Protocol fixed before scoring

- Use reviewed UniProtKB records for selected natural enzymes involved in carbohydrate and nucleotide metabolism. Require one numbered Enzyme Commission (EC) label in the catalytic activity comments, with a PubMed citation marked as experimental evidence for that reaction. Exclude exact sequence matches to the previous housekeeping and structural panels. The evidence is curated by UniProtKB.
- Use one protein per taxon and reaction and cap each reaction in each broad taxonomic group at six examples. Hold out one entire taxonomic group at a time. Keep reactions present in at least two groups, and compare them only within the same broad activity category. Reference sequences, answer labels and candidates for comparison are identical across scorers.
- The primary group contains queries whose closest reference for their documented reaction has **less than 30% global identity**, counting gaps in the alignment length. Report below 25% and the subset with a local alignment covering at least 70% of both proteins. The latter measures whether the relationship can be evaluated over most of each protein. Low identity to a same-reaction reference can also reflect distinct protein families that perform the same reaction.
- Predict the first-ranked exact reaction with raw and reference-normalized BLOSUM62, one HMMER profile per reaction built from the reference proteins, and mean-pooled and residue-local ESM-2 8M. The test ranks known reactions on a fixed list. Report hit counts, accuracy weighting reaction labels equally, and paired wins and losses. Treat the estimates as descriptive given the group sizes.

The [frozen panel](reaction_panel.freeze.json) contains 61 proteins with 14 experimentally supported reaction labels. Before model scoring, 22 proteins qualified for the primary group across 11 labels, and 14 were below 25%.
