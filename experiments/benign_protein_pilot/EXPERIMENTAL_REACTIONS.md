# Identifying the exact reaction of a distant metabolic enzyme

The structural comparisons found a benefit for protein language models on two folds and a conventional alignment advantage on a third. This follow-up asks whether a model can also identify a documented biochemical reaction when given an entire known enzyme and only distant examples of that reaction. [Structure study](STRUCTURAL_FAMILIES.md)

## Protocol fixed before scoring

- Use reviewed UniProtKB records for selected natural enzymes involved in carbohydrate and nucleotide metabolism. Require one numbered Enzyme Commission (EC) label in the catalytic activity comments, with a PubMed citation marked as experimental evidence for that reaction. Exclude exact sequence matches to the previous housekeeping and structural panels. The evidence is curated by UniProtKB.
- Use one protein per taxon and reaction and cap each reaction in each broad taxonomic group at six examples. Hold out one entire taxonomic group at a time. Keep reactions present in at least two groups, and compare them only within the same broad activity category. Reference sequences, answer labels and candidates for comparison are identical across scorers.
- The primary group contains queries whose closest reference for their documented reaction has **less than 30% global identity**, counting gaps in the alignment length. Report below 25% and the subset with a local alignment covering at least 70% of both proteins. The latter measures whether the relationship can be evaluated over most of each protein. Low identity to a same-reaction reference can also reflect distinct protein families that perform the same reaction.
- Predict the first-ranked exact reaction with raw and reference-normalized BLOSUM62, one HMMER profile per reaction built from the reference proteins, and mean-pooled and residue-local ESM-2 8M. The test ranks known reactions on a fixed list. Report hit counts, accuracy weighting reaction labels equally, and paired wins and losses. Treat the estimates as descriptive given the group sizes.

The [frozen panel](reaction_panel.freeze.json) contains 61 proteins with 14 experimentally supported reaction labels. Before model scoring, 22 proteins qualified for the primary group across 11 labels, and 14 were below 25%.

## First results and model-size follow-up

With the protocol and panel locked at `bbad5cb`, raw BLOSUM62 identified 19 of 22 low-identity reactions. The profile identified 18. ESM-2 8M identified 13 from its mean representation and 9 from local comparison. Every correct prediction by local ESM-2 in this group was also correct with raw BLOSUM62. Below 25%, raw alignment and profile scoring each identified 11 of 14, versus 8 with mean ESM-2 and 5 with local ESM-2. Only six of the 22 queries had a broadly covering positive local alignment. On those six, raw alignment and mean ESM-2 each identified five, and the profile and local ESM-2 each identified four.

The gap motivates one check before interpreting the result. Repeat the same experiment with the larger ESM-2 150M checkpoint used in the earlier exploratory cohort study. Hold the proteins, reference set, answer labels, strata and scoring rules fixed. Report both model representations, and check that all classical decisions and identity measurements reproduce. Treat the larger checkpoint as a follow-up chosen after observing the 8M results.
