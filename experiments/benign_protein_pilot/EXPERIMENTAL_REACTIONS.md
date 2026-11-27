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

## Results with both checkpoints

The model-size check was committed at `e7d8279` before running the larger checkpoint. The classical correct/incorrect decisions and the identity measurements agreed for all 61 queries across the two runs. The table gives the number of correct first-ranked reaction labels. [Results and paired outcomes](reaction_results.json)

| Method | Below 30% (22 proteins) | Below 25% (14 proteins) | All 61 proteins |
| --- | ---: | ---: | ---: |
| Raw BLOSUM62 | **19** | **11** | **58** |
| Reference-normalized BLOSUM62 | 17 | 9 | 56 |
| HMMER reaction profiles | 18 | **11** | 57 |
| ESM-2 8M, mean | 13 | 8 | 48 |
| ESM-2 8M, local | 9 | 5 | 48 |
| ESM-2 150M, mean | 13 | 9 | 50 |
| ESM-2 150M, local | **19** | **11** | **58** |

Local scoring with the larger checkpoint improved from 9 to 19 correct labels in the primary low-identity subset. It succeeded on exactly the same examples as raw BLOSUM62 across the full panel, including in the low-identity subset. Mean-pooled 150M identified two of the three low-identity cases that raw alignment missed, while missing eight that alignment identified. This shows some complementary errors. It leaves open how to combine the evidence on fresh cases without knowing the answer.

The model-size check compared the same 61 encoder inputs. Encoding took 14.5 seconds for 8M and 101.4 seconds for 150M on two CPU threads. The shared evaluation considered 860 query/reference pairs and built 42 profiles. Most profiles had multiple references, ranging from one to seven members. A profile that uses a larger collection of homologs could behave differently.

## Interpretation and handoff

Taken together, very low identity can favor a pretrained model for broad structural relationships, as in two of the three [structure panels](STRUCTURAL_FAMILIES.md). Exact biochemical reaction prediction gave a different result. There the larger model reached parity with pairwise alignment, and its local scorer added no new correct labels beyond alignment. The two settings use different questions and reference panels. Neither defines a universal identity cutoff for choosing a method.

The next construct-level test needs fresh public benign records with independently documented distant components and credible competing labels. Preserve the same limits across preprocessing approaches, include pairwise and profile baselines, and use the known protein boundaries in a diagnostic comparison. That comparison can separate failures to find a protein candidate from failures to characterize it. The present cohorts start from known proteins. Their source evidence is curated and was not independently re-audited in the original publications, and overlap with model pretraining is unknown.

## Reproduce

Use the [offline model environment](../../docs/project/HANDOFF.md) with `pyfamsa==0.7.0` and `pyhmmer==0.12.3`. The larger checkpoint is [`facebook/esm2_t30_150M_UR50D`](https://huggingface.co/facebook/esm2_t30_150M_UR50D) at revision `a695f6045e2e32885fa60af20c13cb35398ce30c`. The runner checks both checkpoint hashes. For exact replay compare the UniProtKB release, entries and sequence hashes against the tracked freeze. Keep inputs in ignored storage. [UniProt terms](https://www.uniprot.org/help/license)

```sh
run=results/experimentally-annotated-reactions-replay
python -m experiments.benign_protein_pilot.reaction_study acquire "$run"
for size in 8m 150m; do
  OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
    python -m experiments.benign_protein_pilot.reaction_study score "$run" \
    --esm2-size "$size" --out "$run/report.$size.json"
done
```
