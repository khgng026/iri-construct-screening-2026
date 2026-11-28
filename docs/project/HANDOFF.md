# Research pilot handoff

For the current low-identity results, open review branches, and next decision, start with the [Q3 evidence and handoff](Q3_EVIDENCE.md). The notes below describe the November 22 pilot checkpoint.

**22 November 2026.** Work on `paired-research-pilot`, based on the merged research scope in PR #16. Use the active [scope](SCOPE.md), [plan](PLAN.md), and [status](STATUS.md) for the project direction. The researcher approved a **five percentage point gain in documented-component recall at a matched budget** as the practical target for this pilot.

## What is ready

The paired verifier accepts traces from independently developed methods. It checks that the methods used the same constructs, characterization configuration, and budget grid. Sequence-only arms can be tied to the same input hash. It uses held-out component labels and coordinates, counts each eligible component once, and distinguishes a missing covering query from a covered locus with no correct call. Its tests cover repeated labels, joined coordinates, failed runs, and budget violations. The [verification interface and reproduction commands](../../experiments/paired_preprocessing/README.md) describe how to add an experimental arm.

The reporter pilot exercises the whole path. The annotation-aware toolkit validates the source translations and prepares the answer key and the annotation-derived inputs. Both sequence-only arms receive an identical FASTA with opaque IDs. Source record identities, feature labels, and references are pinned in the [reporter panel](../../experiments/paired_preprocessing/reporter_panel.json). The [protocol](../../experiments/paired_preprocessing/reporter_protocol.json) was committed before scoring (`2f46585`).

## Paired reporter pilot

Five public reporter constructs represent the four archetypes, with seven eligible coding components. The table reports recovered components out of seven, followed by the total actual model calls across the five constructs. Two calls per construct was the primary cap. All arms returned at most five family calls per construct.

| Query cap per construct | Minimal conversion | Start-aware conversion | Annotation-derived inputs |
| ---: | ---: | ---: | ---: |
| 1 | 1/7, 5 calls | 2/7, 5 calls | 5/7, 5 calls |
| **2** | **3/7, 10 calls** | **4/7, 10 calls** | **6/7, 6 calls** |
| 4 | 6/7, 20 calls | 6/7, 19 calls | 7/7, 7 calls |

At the primary cap, start-aware conversion recovered one extra component in the multi-CDS cassette. It used 2,979 input residues across its ten calls, compared with 3,344 for the baseline. The gain of 14.3 percentage points exceeds the target as a point estimate. One recovered component accounts for the entire difference, so this does not establish a general improvement. At cap four, the two arms shared five recoveries and each found a different sixth component. Both methods recovered one of two tagged or fused components at the primary cap. The eligible labels describe whole coding products. Subregion segmentation remains a separate experiment.

Every miss at the primary cap lacked a query covering the documented locus under the chosen strand and overlap rule. In a post-run check using the full candidate lists, the minimal arm had a covering candidate for all seven components and ranked one beyond the largest cap. Start-aware conversion moved a covering candidate into the primary budget for one component but had no covering candidate for another. This suggests testing candidate ordering and candidate preservation separately in the next version. Coverage is a location diagnostic. The successful calls also come from a small panel with two eligible target families and up to five returned calls per construct. A larger study should test characterization with stronger competing labels.

## Low-similarity protein pilot

Compared mean-pooled ESM-2 8M with local BLOSUM62 alignment on 56 known proteins against 26 separate reviewed references, spanning 15 housekeeping labels. Labels use the source EC and InterPro documentation. The full reference panel stayed fixed for each query. The [study](../../experiments/benign_protein_pilot/README.md) defines the tiers and comparators. The code and tier definitions were frozen before scoring (`c905380`).

| Query group | Queries | Alignment hit at 1 / at 5 | ESM-2 hit at 1 / at 5 |
| --- | ---: | ---: | ---: |
| Low identity, with adequate coverage | 8 | 7 / 8 | 5 / 7 |
| Intermediate identity | 41 | 41 / 41 | 39 / 41 |
| High identity | 1 | 1 / 1 | 1 / 1 |
| Insufficient same-label alignment coverage | 6 | 6 / 6 | 0 / 5 |
| All queries | 56 | 55 / 56 | 45 / 54 |

Reference-normalized alignment gave the same aggregate counts. In the low-identity tier, each method retrieved a top-ranked correct label missed by the other, with three alignment-only successes and one ESM-only success. Only one query had an adequately covered competing-label alignment with identity at least as high as its best same-label alignment. The sample gives an initial result for known housekeeping proteins and leaves the comparative benefit on a larger, more discriminating panel unresolved. It does not measure end-to-end characterization of distant components from DNA constructs.

## Where to resume

1. **Build the next answer key first.** Add independently documented benign constructs with more diversity and usable subregion labels for fusions. Choose eligibility and group related constructs before scoring. The verification work can proceed independently of preparation and scoring methods.
2. **Use the current arms as controls.** Separate whether a component is represented in the candidate pool from whether its candidate is queried at the matched budget. Compare new candidate selection and, where subregion labels support it, candidate segmentation on a fresh pilot. Keep the user-approved effect target explicit and use enough independent examples to interpret uncertainty.
3. **Strengthen the remote comparison in parallel.** Expand the independently labelled low-identity tier with genuinely confusable benign families and preserve adequate coverage and reference separation. Keep alignment as a comparator. Link this evaluation to construct-level results only when documented distant components are available.

A metadata check found two distinct in-frame regions annotated within each of two public tomato histone–fluorescent protein fusions. Treat the records as one related group if selected. Review the source context and assemble independent groups before a segmentation study. [AB684655.1](https://www.ncbi.nlm.nih.gov/nuccore/AB684655.1) · [AB684657.1](https://www.ncbi.nlm.nih.gov/nuccore/AB684657.1)

A separate actin–fluorescent reporter record documented a short actin fragment. It would contribute the fluorescent region alone under a complete-component endpoint. [DQ344484.1](https://www.ncbi.nlm.nih.gov/nuccore/DQ344484.1)

## Reproduction and artifacts

The run used CPU inference with Python 3.12.14, Biopython 1.86, PyTorch 2.6.0 CPU, Transformers 4.48.1, and Hugging Face Hub 0.28.1. This session used `/tmp/iri-benign-model` and `/tmp/iri-benign-hf`. Set `HF_HOME` to the cache path. For in-place commands, add `PYTHONPATH=src` or install the project from the checkout. Both pilots load the model offline. To prepare a missing cache, call `snapshot_download` from `huggingface_hub` for `facebook/esm2_t6_8M_UR50D` at revision `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`. Limit acquisition to the config, tokenizer files, vocabulary, and `model.safetensors`. The runners check the weight hash.

Local reports live at `results/paired-reporter-pilot-v0/report-with-diagnostics.json` and `results/benign-protein-pilot-v0/report.json`. The ignored run folders also contain the frozen source snapshots. The tracked panel files contain source identities and hashes. Exact reruns require the recorded input bytes. A changed provider snapshot should become a new version of the run. Keep raw sequences, annotation bundles, and model weights in ignored storage. Preserve the original reports and code of the historical studies.

The reporter acquisition and preparation commands reproduced the FASTA, answer key, and reference panel byte for byte from the pinned sources. The repository suite ran 277 tests with five historical optional skips. The dependency-free smoke check and preserved pilot evidence check passed.
