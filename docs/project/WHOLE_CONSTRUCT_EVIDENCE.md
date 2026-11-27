# Whole-construct benchmark checkpoint

The full-pool three-way comparison is still open. The retrospective audit shows that gene calling substantially reduces candidate work, while the minimal full pool overlaps more documented mid- and low-identity components. The uncapped BLASTx pilot identifies all verified low-identity reporters at the 50% locus rule. The pilot covers one family. Fresh ESM-2 inference and a diverse independent panel are still needed. [Study plan](WHOLE_CONSTRUCT_STUDY.md)

## Historical uncapped availability

For each documented component, checked every recorded candidate for the historical 80% same-strand documented-locus rule. For a single-span coding locus, also checked whether a covering candidate shared its expected frame. These are geometric checks. The saved first-four ESM-2 results give evidence about individual candidate rankings. The start-ranked arm reorders the original minimal pool, so identical geometry and residue count allow its saved rankings to extend coverage. The full pool and any duplicated rankings were verified before reuse. [Aggregate audit](../../results/whole-construct/historical-availability.json)

| Historical tier | Documented | Minimal pool covers | Gene-call pool covers | Minimal only | Gene call only |
| --- | ---: | ---: | ---: | ---: | ---: |
| Low | 250 | 248 | 237 | 11 | 0 |
| Intermediate | 426 | 421 | 421 | 1 | 1 |
| Intermediate + low | 676 | 669 | 658 | 12 | 1 |
| High | 508 | 506 | 483 | 23 | 0 |
| Insufficient reference coverage | 37 | 35 | 35 | 0 | 0 |
| All | 1,221 | 1,210 | 1,176 | 35 | 1 |

The intermediate/low counts come from 536 constructs, 143 distinct deposited protein hashes and 18 families. The two largest families account for 407 of 676 occurrences. Under the conservative linkage rule that joins constructs with shared proteins, one group contains 644 of 676. This panel is useful for paired descriptive analysis. Its size overstates the number of independent observations.

The 777 archived vectors supplied **101,407 minimal candidates (6,544,089 residues)** and **4,194 gene-call candidates (940,383 residues)** before cross-construct caching. That is about 24-fold fewer candidates and 7-fold fewer residues with gene calling. Saved rankings cover 3,270 minimal candidates and 2,772 gene-call candidates. No minimal pool and 336 gene-call pools have all rankings. Whole-pool scoring is a materially different compute job, especially for the minimal arm.

At least 535 of 676 intermediate/low components have a covering minimal candidate with the correct family ranked first in the saved traces. At least 599 have one in the gene-call traces. If every relevant remaining unknown rank were favorable, those figures could reach 612 and 631. This bounds the **gene-call minus minimal difference for the unfiltered per-component top-family diagnostic** between **−13 and +96** in the combined tier. In the low tier alone the range is **−23 to +22**. Intermediate alone ranges from **+10 to +74**. This arithmetic bounds potential family evidence when retaining every top call, without thresholding, locus consolidation, or one-to-one assignment. It provides neither a calibrated characterization estimate nor a comparison with BLASTx on this cohort.

## New uncapped BLASTx pilot on existing reporter records

Reused the independently corroborated Open Reporters archive after checking the pinned source hashes and exact deposited reporter translations. It has 19 verified GFP-like reporters with distinct translations and five other catalogue reporters, a shared cloning background, and seven natural GFP-like references. All 19 verified reporters fall below 40% sequence identity at 70% or more coverage of both proteins against this panel. The verified target records have only one scored cargo family. [Pilot results](../../results/whole-construct/reporter-pilot.json) · [Original archive and corroboration](https://github.com/khgng026/iri-construct-screening-2026/blob/study/q3-evaluation/docs/project/Q3_OPEN_REPORTERS_PROTOCOL.md)

Without a maximum number of reported loci, BLASTx recovered **15/19** reporters at 80% locus overlap and **19/19** at 50%, at both E-value thresholds `1e-5` and `1e-10`. At `1e-20`, it recovered **5/19** and **8/19**, respectively. The more permissive settings emitted 38 calls after consolidation, of which 23 remained unmatched under the strict 80% rule. The records have incomplete independent family adjudication, so unmatched calls are unverified. Some hits become supported when evaluated as local evidence.

Both complete translation pools contained **all 19 exact documented peptides** before scoring. Minimal translation produced 1,299 candidates (83,325 residues) across the 24 records and gene calling produced 48 (13,713 residues). On this machine raw BLASTx took about 0.5 seconds. Minimal candidate generation was about 0.04 seconds and gene calling about 0.10 seconds. These measurements cover this small archive. They do not estimate full ESM-2 runtime or identification accuracy.

## Redesign data and model access

| Source | Evidence and disposition |
| --- | --- |
| JointDiff fluorescent reporter validation | The paper describes experimentally tested designs, parent controls, and whole-plasmid sequence checks. The linked GitHub snapshot `9c481f69` did not expose an individually deposited complete plasmid sequence for an assayed pair. Retain as a candidate pending linkage to the primary whole-plasmid files. [Paper](https://doi.org/10.1002/pro.70340) · [Author repository](https://github.com/shen-Lab/JointDiff) |
| SAPP fluorescent protein redesign | The paper includes measured reporter variants and WT controls, with a code/data archive. The linked `70b695f5` GitHub snapshot exposed entry vector files. An entry vector does not identify the completed assay construct for an individual redesign. No validated full pair was obtained. [Paper](https://www.nature.com/articles/s41467-026-76740-9) · [Data archive](https://doi.org/10.5281/zenodo.20761616) |
| mSandy2 red fluorescent reporter | Identified an experimentally characterized variant with an Addgene full sequence map. The provider requires a signed-in account to download the bytes. It remains a candidate pending the matched parent and measured identity with the project panel. [Paper](https://doi.org/10.1039/D1SC05088E) · [Provider record](https://www.addgene.org/177760/) |
| LuxSit-i luciferase | Public deposition identifies the gene and the paper documents activity. Full constructs are referred to Addgene. Treat as a possible separate family only if an appropriate reference panel represents it. A gene-only deposition supports an isolated protein diagnostic. [Paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9946828/) |

The selected runtime enforces restricted networking. Direct HTTPS requests to the official Hugging Face checkpoint, the Meta weight endpoint, UniProt, NCBI, and Addgene returned HTTP CONNECT 403. The local workspaces contain no checkpoint or sequence-bearing copy of the frozen cohort inputs and 134-reference panel. The CPU machine and BLASTx 2.16.0 work, and the current research venv has the existing parser and Pyrodigal. It lacks the model stack. Package setup and scaled scoring were deferred until the exact weights and reference bytes can be obtained. The current task proceeded without an environment change request.

## Handoff for the next work session

1. Obtain approved local copies of the pinned ESM-2 8M cache and public input snapshots, or use a runtime with the necessary provider access. Confirm the model revision and weights against the existing constants. Follow the historical cohort reconstruction if rerunning that cohort. Keep its previously released answer key retrospective. [Cohort reconstruction](../../experiments/paired_preprocessing/COHORT_STUDY.md#reproduce)
2. Start with independently validated and archived benign redesign records. Freeze pair inclusion, family evidence, source-record hashes, identity tiers, and lineage partitions before looking at characterization results. Confirm the correct family is represented in the shared library. The source audit above has leads and no accepted redesign pair. Add independent natural archival cases for power and report redesigned and natural subsets separately.
3. Time full-pool ESM-2 on a representative development construct. Scale using measured distinct inputs, residues, wall time, and the saved workload above. The scorer consumes a FASTA for all three sequence-only arms and the same historical reference JSON convention (`accession`, `family`, `sequence`). Annotation-derived scoring consumes a construct-to-documented-unit JSON and stays diagnostic. Keep scoring outputs local.
4. Use the evaluator with separate development and holdout calls. It calibrates method-specific score cutoffs at predeclared unsupported-call burdens on development, checks fingerprints across comparators, and reports achieved held-out burden and paired lineage intervals by identity tier. Treat unsupported calls on partial records as unverified. Review the sensitivity of a calibrated answer to family and conservative group dominance. The primary endpoint still applies the 80% one-to-one rule. Add a separately declared linked-window analysis before interpreting long-protein misses.
5. Claim an overall improvement only if both held-out comparisons support it at similar burden with enough independent families and lineages. If the low tier yields mixed results, report that uncertainty separately from the intermediate tier. Confidence in the cached rankings or the small reporter pilot cannot fill in missing full-pool scores.

## Reproduce completed analyses

From this checkout with `PYTHONPATH=src` and a Python environment with Biopython 1.86 and Pyrodigal 3.6.3, write to new paths. The reporter runner uses the frozen Open Reporters audit from `study/q3-evaluation` and the two cited public source checkouts. Its assessment only writes summary counts and durations.

```bash
PYTHONPATH=src python -m experiments.whole_construct.availability \
  --out /tmp/full-pool-availability.json

PYTHONPATH=src python -m experiments.whole_construct.reporter_pilot \
  --audit /path/to/open-reporters.json \
  --construct-repo /path/to/Open-DNA-Collections \
  --protein-repo /path/to/Orthologous_GFP_Fitness_Peaks \
  --blastx /path/to/blastx --out /tmp/full-pool-reporter.json
```

Use `python -m experiments.whole_construct.scoring --help` and `python -m experiments.whole_construct.compare --help` for the uncapped runner and the held-out comparison. The new BLASTx arm ran on real archive records. The ESM-2 inference interface and evaluation were checked with a test double and deterministic fixtures. No fresh ESM-2 inference was possible here.
