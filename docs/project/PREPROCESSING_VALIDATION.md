# Preprocessing validation with ESM-2 held fixed

Gene calling remains the best-supported change from minimal six-frame translation. The historical ESM-2 replay supports testing a different gene-call order when only two inputs can be scored. That order traded off recoveries across families and displaced one of the two windows of a long yeast protein. A separate public enzyme archive was too easy to confirm an improvement beyond gene calling.

## Identification with saved ESM-2 scores

The replay uses the same 777 archived benign vectors and 1,221 documented components as cohort v2. It keeps the 8M checkpoint (revision `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`), mean-pooled embeddings, reference library and cap of five family calls per construct fixed. An input is a translated peptide submitted to ESM-2. A family call is a label emitted for an input. The replay verifies coordinate, length, source and reference fingerprints and checks that the original four-input prediction schedules reproduce exactly. [Original cohort results](../../experiments/paired_preprocessing/COHORT_RESULTS.md)

| Maximum inputs per construct | Minimal translation | Gene calling | Gene calling with mixed ranking |
| --- | ---: | ---: | ---: |
| 1 | 279 (22.9%) | 323 (26.5%) | 323 (26.5%) |
| 2 | 573 (46.9%) | 725 (59.4%) | 804 (65.8%) |
| 4 | 1,014 (83.0%) | 1,099 (90.0%) | 1,099 (90.0%) |

The mixed rule holds the first gene-called input in place and orders the rest by how close their lengths are to 250 residues. At two candidates it found 234 components gene calling missed and lost 155 that gene calling found. The equal-lineage gain was 6.8 percentage points with a descriptive 95% lineage bootstrap interval of 3.1–10.2 points. Both used 1,537 model inputs, while mixed ranking used 525,567 residues versus 551,358 for gene calling. The mixed rule gained the same 79 components when each construct was subject to the original gene arm's residue ceiling. That ceiling is a sensitivity check. It can exclude a longer baseline translation that would otherwise fit the candidate cap.

GFP-like proteins and aminoglycoside phosphotransferases accounted for 95 net gains combined, with a net loss of 16 across the remaining families. In the historical low-similarity tier, recoveries rose from 121 to 137 to 152 of 250 across the three arms. The tier uses less than 40% identity to same-family references at 70% coverage or more. In the separate 37-component insufficient-coverage tier, recoveries were 31, 34 and 27. The archived scores were available for at most the original first four candidates from each method. The full mixed rule would request at least one unscored candidate in 35 constructs at cap two and 82 at cap four. On the 742 constructs where the full mixed selection could be replayed exactly at cap two, recovery rose from 713 to 783 of 1,162 components. The four-input result covers reordering only within the available shortlist. At two inputs the cohort's annotation-derived arm reached 91.4%, leaving 25.6 percentage points above the mixed replay. [Replay counts and diagnostics](../../results/preprocessing-validation/esm2-replay.json)

Offering distinct translation candidates after gene calls added 17 model inputs at cap two and 336 at cap four. It produced zero net ESM-2 recoveries at either cap. At cap four, one component was gained and one lost in each of the two fallback combinations. The five-call output cap can redistribute predictions when another input is added. The replay can use the saved start-ranked translations for this test. Resolved stretches and alternative windows from the archive analysis need fresh model scores.

## Candidate availability in public archives

These runs verify that the selected candidates contain the full documented peptide at the correct locus and frame. For long proteins, multiple selected windows can jointly cover the peptide. This endpoint comes before family identification. Source annotations serve as the answer key and enter only after candidate generation and ranking. [Archive results and source dispositions](../../results/preprocessing-validation/archive.json)

| Panel and input cap | Minimal translation | Gene calling | Mixed ranking | Mixed with first gene's windows kept together |
| --- | ---: | ---: | ---: | ---: |
| Open Enzyme, two inputs | 16/17 | 17/17 | 17/17 | 17/17 |
| Open Yeast, two inputs | 10/11 | 11/11 | 10/11 | 11/11 |
| LAMPRA, two inputs | 0/1 | 0/1 | 0/1 | 0/1 |

- **Open Enzyme.** The panel was selected from 27 public records. Seventeen passed the translation checks, nine lacked a deposited translation and one had an incomplete CDS. Each eligible record contributes one specified cargo target. Ten of the 17 were polymerases, and the vectors share a cloning background. At cap one, minimal translation reached 15 and gene calling reached 16. At cap four, both reached all 17. At cap two each used 34 inputs, with 18,552 and 16,761 residues respectively. This panel gives a fresh measurement of candidate availability for clean maps. A single, often prominent cargo gives little room to measure which extra step helps across a whole construct. [Pinned panel](../../experiments/preprocessing_exploration/validation_panel.json)
- **Open Yeast.** Eleven of 17 previously inspected cargos passed validation. The 1,392-residue target required two inputs. Its first two gene-called windows together covered the peptide, and the mixed rule displaced one. The exploratory first-gene guard restored it at the same 7,497 total residues as original gene calling. Keeping whole parent groups together and using two balanced windows also retained the target. This is one informative long-protein case. The guard was formulated after seeing the miss and was checked against the same case.
- **LAMPRA.** On this previously inspected incomplete map, all available gene calls missed the exact reporter. The resolved-translation arm exposed it by the fourth candidate after selecting 1,039 residues. Offering one distinct translation after the first two gene calls still missed it at four candidates. Gene-first followed by distinct translation exposed the reporter by eight candidates and 1,465 residues. Every arm missed under the 716-residue ceiling set by gene calling's full pool. This supports testing an explicit extra allowance on ambiguous archive inputs. The evidence is one map and one known reporter.

## What to carry forward

Keep minimal translation and original gene calling as the two anchors. The next informative test is a blinded fixed-ESM comparison on naturally archived benign vectors with multiple independently verified cargos, valid long CDSs, or incomplete records with verified local translations. The enzyme panel has little headroom for that test. Prioritize the two-input mixed rule with and without the long-window guard. For incomplete records, test a separate small resolved-translation budget. Prespecify family- and lineage-level reporting and count documented gains and losses alongside model calls, processed residues and unsupported predictions. Fresh embeddings are required for new candidates and for the 35 constructs where the full mixed rule exceeds the historical two-input trace. The pinned weight endpoints returned HTTP 403 from this workspace, so this study used only existing ESM-2 scores.

## Reproduce

Use Python with `biopython==1.86`, `pyrodigal==3.6.3`, and the local package on `PYTHONPATH`. Use Git checkouts of the Reclone and LAMPRA source repositories containing the revisions in the manifests. The runners read pinned record bytes from those revisions and verify their hashes. The two output paths must be new. [Protocol and panel provenance](../../experiments/preprocessing_exploration/VALIDATION_PROTOCOL.md)

```bash
PYTHONPATH=src python -m experiments.preprocessing_exploration.validation \
  --reclone /path/to/Open-DNA-Collections --lampra /path/to/LAMPRA \
  --out /tmp/iri-preprocessing-archive.json
PYTHONPATH=src python -m experiments.preprocessing_exploration.validation_replay \
  --out /tmp/iri-preprocessing-esm2.json
```
