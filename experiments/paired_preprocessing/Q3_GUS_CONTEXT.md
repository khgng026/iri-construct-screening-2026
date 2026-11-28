# A distant plant reporter across construct layouts

This is a fixed end-to-end comparison for [Q3](../../docs/project/SCOPE.md). It asks whether the same reporter domain reaches the scorer and gets identified correctly across different archived construct layouts. Holding the reporter nearly constant helps distinguish candidate selection from scoring.

## Protocol fixed before scoring

- Use the seven dated NCBI records in `gus_context_panel.json`. They cover a compact moss reporter, an archived GFP–GUS fusion, a series of four long dual-reporter vectors, and one long single-reporter vector. Treat the four vectors from the same submission as one series. All cases represent one enzyme ancestry. Six have the same catalytic domain sequence. The compact case has one deletion and seven substitutions relative to the known structure.
- The answer key is the beta-glucuronidase catalytic domain. The archived coding sequence must independently reproduce the deposited translation and the domain interval must retranslate to the frozen protein. Map the experimental structure 3K46 chain A residues 274–601 to reviewed UniProt P05804 and then to each translated reporter. The structural label is SCOPe c.1.8. The embedded reporter in the fusion is validated as an in-frame segment of its complete deposited translation.
- Compare against the same six curated SCOPe 2.08 domains for every method. Two represent c.1.8 and four represent two competing superfamilies in the same TIM barrel fold. The structure used to establish the answer is outside the reference set. The primary tier has less than 25% global identity to every same-label reference. Also report whether any same-label local alignment covers at least 70% of both domains.
- Pass only DNA sequence to the existing minimal and Pyrodigal candidate generators. Rank by their existing rules and compare at a primary cap of two protein queries, with one and four as sensitivity checks. Count one correct component if a selected candidate covers at least 80% of its domain on the correct strand and its top-ranked label is correct. Report whether a covering candidate existed before the query cap, actual queries, and input residues.
- Apply raw and reference-normalized BLOSUM62, an HMMER profile per superfamily, and mean and residue-local ESM-2 at 8M and 150M to the same candidates. Run two diagnostics separately using the documented full reporter or fusion product and the independently mapped catalytic domain. Those diagnostics use the answer key and are not end-to-end conditions.

Report the seven records and four submission groups descriptively. Layout and submission differ together, and the enzyme ancestry is shared. This comparison measures recognition of a structural superfamily within a small fixed reference set.

[NCBI records and source hashes](gus_context_panel.json) · [Experimental structure](https://www.rcsb.org/structure/3K46) · [Reviewed protein](https://www.uniprot.org/uniprotkb/P05804/entry) · [SCOPe superfamily](https://scop.berkeley.edu/sunid=51445)

## Results

The sources and protocol were committed at `8bb9bbc` before scoring. The closest same-group reference was 16.99% to 17.06% identical to each verified domain by the global measure. The count of same-group alignments covering 70% of both domains was zero.

At one query per record, both preparation arms reached the documented domain on all seven DNA records. Raw BLOSUM62 and both ESM-2 8M methods identified all seven. The profile, normalized BLOSUM62, and both ESM-2 150M methods identified zero. The same seven records were reached and identified at caps of two and four. The table groups related submissions. [Paired decisions, query counts, and input sizes](gus_context_results.json)

| Layout | Independent submissions | Records | Reached at one query in both arms | Raw alignment | ESM-2 8M local | Profile HMM | ESM-2 150M local |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Compact single reporter | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| Fused reporters | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| Long dual reporter | 1 | 4 | 4 | 4 | 4 | 0 | 0 |
| Long single reporter | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| Total | 4 | 7 | 7 | 7 | 7 | 0 | 0 |

At the primary cap of two, the minimal arm used 14 queries and 7,216 residues. Gene calling used 13 queries and 6,717 residues. At one query the arms each used seven inputs, with 4,535 and 4,475 residues respectively. This held the input budget comparable. The models and conventional methods have different processing costs.

The documented full reporter or fusion gave the same correct and incorrect calls as the covering candidates for every scorer. With the verified catalytic domain supplied, the 8M methods still identified all seven and raw BLOSUM62 identified one. The six raw-alignment misses represent **one exact domain sequence** appearing in three submissions. The full product restored raw alignment on that sequence. This is a case where the model adds a correct judgment on an independently verified region and gives the same end-to-end recall once raw alignment receives the longer protein. The scoring difference is specific to a small reference set, one enzyme ancestry, and these two checkpoints. Estimating how often it occurs calls for more reporter ancestries.

## Reproduce

Use the [model environment](../../docs/project/HANDOFF.md) with Pyrodigal 3.6.3, PyFAMSA 0.7.0 and PyHMMER 0.12.3. Cache both immutable ESM-2 revisions in `cohort_scorers.py`. Acquisition verifies the deposited records, six SCOPe reference snapshots, and the RCSB mmCIF. The structure audit checks the direct mapping to UniProt and the domain translation of each source. The archived classification assigns 3K46 chain A residues 274–601 to superfamily c.1.8 as domain `d3k46a3`. A fresh source download and inference run reproduced the result JSON byte for byte. [SCOPe domain](https://scop.berkeley.edu/sunid=305651)

```sh
study=experiments.paired_preprocessing.distant_reporter_study
run=results/q3-gus-context-replay
python -m "$study" --panel gus-context acquire "$run"
python -m "$study" --panel gus-context audit "$run"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python -m "$study" --panel gus-context score "$run" --out "$run/report.json"
cmp "$run/report.json" experiments/paired_preprocessing/gus_context_results.json
```
