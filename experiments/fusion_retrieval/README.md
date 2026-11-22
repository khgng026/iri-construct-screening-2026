# Fusion retrieval: code and measured results

This adds the completed fusion-retrieval work alongside `component_study_v1`.
The toolkit extracts annotated CDSs; these experiments examine what happens
inside a protein containing multiple components. No toolkit behavior changes.

## Run now

From the repository root:

```bash
python experiments/fusion_retrieval/replay.py
python -m pip install -r experiments/fusion_retrieval/requirements.txt
python experiments/fusion_retrieval/aggregate.py experiments/fusion_retrieval/fixtures/scores.json --budget 2
python -m unittest discover -s experiments/fusion_retrieval/tests -v
```

The first command verifies source hashes and recomputes the example's target
positions from its actual saved output lists. It needs only Python. The second
runner accepts a complete chunk-by-reference score matrix and applies the
existing max-score or round-robin merger without target labels. The fixture is
fictional: max returns a,b; round-robin returns a,d. It is a software check,
not a new biological result.

## Results already measured

| Comparison | Whole ESM | Local ESM | Alignment | Meaning |
|---|---:|---:|---:|---|
| Expanded strict diluted cell, total K10 | 7/51 | guided RR 30/51 | guided RR 45/51 | Six targets absent; one designated target per query |
| Four deposited proteins, total K10 | 0/4 | guided 2/4 | 4/4 | Even supplied true slices miss two ESM targets |
| Six complete fusions, composition control | 4/6 both parents | supplied parts 6/6 | — | Cheap baseline; learned-model comparison pending |

Whole-half merging reduced the primary ESM result 30→26. Expanding context
reduced RR 30→25, while max improved 16→21. Neither is a general improvement.
Identity-deduplicated refill changed no hits. These studies use different
catalogues and must not be pooled with the shared 8M component study.

`results/` contains original summary tables and the existing fusion example.
`MANIFEST.json` gives their original paths and byte hashes. Source paths identify
archival locations; they are not required to run the packaged replay. The
example improves target rank 2→1, not a missed-target rescue. Raw chunk cosine
scores are not present in that example.

## Code and upstream work

`scripts/budget_compare.py` and `whole_preserving.py` are the existing total-budget
mergers, copied unchanged. Other scripts preserve the historical preparation,
scoring, reporting and audit implementation. They require original private
configs/caches and are **not** a fresh-clone end-to-end inference command.

The original working pipeline is a pinned Git submodule:

```bash
git submodule update --init experiments/fusion_retrieval/sources/nishit-working
```

It retains its code, public results, presentation, MIT license and citation
metadata. Pin `5dfef4734f07b981fa6305383fee170d6af28a6d` is the version used by
the local adaptation. `SOURCES.json` also records the newer reviewed revision;
its reference-window builders and length diagnostics have not been validated
in these local experiments. Do not silently upgrade a historical run.

`patches/conditional-macro.patch` is a separately tested historical repair:
exclude undefined payload ranks from conditional mean rank, while retaining
missing cases in recall/MRR. Not applied silently to the submodule.

The upstream 144/144 best-segment headline is a diagnostic with potentially
multiple K10 lists, not one fixed-budget final answer. Its exact historical
40,750-sequence catalogue and eight canonical raw prediction files are missing
locally. No exact canonical reproduction is claimed. Original benchmark:
https://github.com/nk-008120/plasmid-pilot-benchmark

## Next experiment

[Approach inventory and exact missing comparisons](APPROACHES.md) distinguishes
completed work, code-only candidates and the pending frozen study.

1. Recover complete existing chunk scores and reproduce max/RR from them.
2. Define consensus and one weighted rule on separate development cases.
3. Freeze one rule; compare whole, max, RR and alignment with the same references
   and five total distinct returned identities on untouched cases.
4. Score each component through its supporting region and report whether all
   supported components were found. Keep supplied boundaries oracle-only.

Reuse the released parser and existing scorers. Do not rerun the component
coverage study. Keep the prepared 12-query 650M comparison unchanged and pending.
No new inference was performed when preparing this contribution.
