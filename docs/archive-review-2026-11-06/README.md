# Public-archive acquisition audit — 2026-11-06

## Executed, not merely proposed

GitHub Actions [audit run](https://github.com/khgng026/iri-construct-screening-2026/actions/runs/42233935925) completed on the candidate branch using Python 3.11/Biopython 1.86. NCBI GenBank exact-version retrieval succeeded **11/11**; no ENA fallback was necessary. The job parsed each returned full flat-file record, required an exact accession-version match, and recorded original-byte SHA-256, file size, sequence length, feature types/counts and CDS counts.

The **source-free** results are in [report.json](report.json). Original biological files existed only in the job's ephemeral ignored `results/` directory and were **not** committed, uploaded as workflow artifacts, approved as inputs, or used for model inference. The downloader permits rerunning with a fresh documented checksum. No Addgene password, protected file or customer record was used.

All 11 are `ACQUIRED_NOT_ADMITTED` and require an exact source-record T02 licensing/annotation review before local/cloud storage or scored evaluation. Names are not proof that an NCBI vector is the same file as a selected Addgene plasmid. Preserve the **10/10 researcher-approved candidate sign-off** and all **nine Addgene HOLD records** until a separate explicit substitution/admission decision.

## Audited accession candidates (NCBI)

| NCBI accession | Annotated name from returned record | nt | Features | CDS | Scientific role / caveat |
|---|---|---:|---:|---:|---|
| U47295.2 | pGL3-Basic | 4,818 | 13 | 2 | Name-level match to Addgene #212936; exact-file identity not verified |
| U47296.2 | pGL3-Control | 5,256 | 15 | 2 | Alternative reporter construct; not the same plasmid as #18964 |
| AY864931.1 | pGL4.76[hRluc/Hygro] | 5,121 | 12 | 3 | Renilla reporter, not RLuc8 variant identity |
| U55762.1 | pEGFP-N1 | 4,733 | 4 | 2 | Fluorescent reporter and selection marker; not pYIC or EGFP-LC3 |
| AY738222.1 | pGL4.10[luc2] | 4,242 | 10 | 2 | Single firefly reporter, not verified dual luciferase |
| KJ170897.1 | pACYC184-ara-xylE | 6,369 | 5 | 2 | Cloning vector plus metabolic reporter; actual feature roles need review |
| AF298787.1 | pUG35 GFP fusion vector | 6,231 | 11 | 3 | Fusion-*ready* vector, not a completed GFP-target fusion |
| KJ541668.1 | pCellFree_G04 C-term eGFP 8xHis vector | 5,275 | 18 | 4 | Tagged fusion-ready architecture, not Nano-lantern |
| KJ541670.1 | pCellFree_G06 C-term mCherry 8xHis vector | 5,266 | 18 | 4 | Tagged fusion-ready architecture, not Renilla–FLAG–firefly |
| PP098727.1 | pETDuet-1-GFP | 6,104 | 2 | 1 | Name suggests dual-expression vector; only one annotated CDS; review before multi-CDS use |
| JX560321.2 | pSEVA111 | 2,094 | 9 | 1 | Modular backbone candidate |

The table is extracted from the actual acquisition report, not asserted from uninspected public catalogue descriptions. In particular, an available cloning/fusion site does not establish a completed fusion, and vector family names do not guarantee adequate source annotations.

## Next decisions and gates

1. Review/approve revised candidates in the separately versioned shortlist. Do not silently overwrite the previous nine holds or their researcher sign-off.
2. Audit whole-record annotations, source caveats and private-copy/reuse rights; record reviewer and exact SHA-256.
3. Fill coverage gaps with **documented benign multi-CDS and actual fusion constructs** if available without a login gate. The current alternative set includes reporter vectors and fusion-ready vectors but does not establish those archetypes as complete proteins.
4. Define reviewed component-family labels and independent references, keeping exact source/lineage exclusion groups separate. The public archive discovery step alone does not resolve E2/E3.
5. Only after review, stage permitted exact records in a separate T02-admitted manifest and run E1, then matched retrieval experiments. No model results or scientific accuracy claims were produced.

## Reproduction

```bash
python -m pip install biopython==1.86
python scripts/review_open_archive_candidates_v2.py --output results/open_archive_review_2026_11_06
```

Note: The initial `scripts/review_open_archive_candidates.py` was superseded by the corrected `_v2.py` entry point. Use **v2**. The GitHub Actions [metadata publication run](https://github.com/khgng026/iri-construct-screening-2026/actions/runs/42233935925) succeeded; the earlier verification-only workflow also succeeded.

## Source URLs

Original accessions are linked individually in [the candidate registry](../../data/manifests/open_archive_candidates_2026_11_06.json). NCBI records and provider terms, not search-result snippets, determine each source's admissibility.
