# Reporter preprocessing comparison v1

Test candidate ordering separately from candidate generation on the five public reporter constructs in the original pilot. Reuse their seven eligible components, the same frozen sequence-only inputs and reference panel, and the same ESM-2 characterization. The annotations are available to the verifier and the annotation-derived diagnostic arm. [Reporter pilot protocol](reporter_protocol.json)

## Comparison fixed before scoring

| Arm | What the model receives |
| --- | --- |
| `minimal` | The original length-ranked translation blocks. |
| `start_aware` | The original length-ranked blocks shortened to their first methionine, with blocks that lack one dropped. |
| `start_ranked` | The full `minimal` candidate pool, with the same lengths, sequences, and loci. Rank by length minus the distance to the first methionine. Cap the penalty at the existing minimum query length, including when there is no methionine. Break ties in the original order. |
| `gene_caller` | Length-ranked translations from Pyrodigal 3.6.3 in metagenomic mode with default open sequence ends and no masking. This comparator uses pretrained prokaryotic models. |
| `annotation_reference` | Validated annotated proteins for the separate diagnostic. |

The translated candidates use the original model length rules. Report documented-component recall at query caps of one, two, and four, along with actual calls and input residues. Use two queries per construct as the primary comparison and five percentage points as the existing practical effect target. Some arms may exhaust their pools below a cap. Report those counts and distinguish missed loci from loci present in a query with no matching family call.

This is an exploratory comparison on the already studied reporter panel. A single component accounts for about fourteen recall points. The fusion labels here describe whole coding products. A separate experiment with independently documented subregions is needed to evaluate segmentation.

## Results and next decision

The comparison rules and code were committed before inference at `65e1abc`. Each table cell reports components recovered out of seven, followed by actual model calls across the five constructs.

| Query cap per construct | Minimal | Start-aware filtering | Start ranking | Pyrodigal | Annotation-derived |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1/7, 5 calls | 2/7, 5 calls | 2/7, 5 calls | 4/7, 5 calls | 5/7, 5 calls |
| **2** | **3/7, 10 calls** | **4/7, 10 calls** | **4/7, 10 calls** | **6/7, 8 calls** | **6/7, 6 calls** |
| 4 | 6/7, 20 calls | 6/7, 19 calls | 6/7, 20 calls | 7/7, 11 calls | 7/7, 7 calls |

At the primary cap, Pyrodigal recovered three more components than the minimal arm. It used 2,486 input residues compared with 3,344 for the minimal arm. Pyrodigal generated 12 candidates across the panel, compared with 262 for the minimal arm. Its smaller pool was sufficient for all seven components at the largest cap. A recovery means a correct family call supported by the protocol's locus overlap rule. This endpoint does not require exact coding boundaries or reconstruction of a joined protein.

Start ranking matched start-aware filtering on recall at every cap. It retained the original 262 candidates. In a post-run check of coverage before applying a cap, the filtered pool contained a covering candidate for six documented loci, while the minimal, ranked, and Pyrodigal pools covered all seven. A replay that counted only the first family call from each queried input gave the same recovery counts for every arm and cap. In this panel the primary differences came from which loci reached the model.

Carry Pyrodigal forward as the leading comparator and retain start ranking as the control for reordering the same candidate pool. The larger evaluation should use independently verified records to establish where those gains persist. The pilot's small reporter panel and Pyrodigal's prokaryotic assumptions limit generalization across construct types.

The [crossed characterization follow-up](CROSSED_CHARACTERIZATION.md) checks the same selected inputs with ESM-2 and normalized BLOSUM62.

All three original arms reproduced their queries, predictions, and counts with the extended runner. The repository suite ran 281 tests with five optional skips. Both the dependency-free smoke check and historical evidence check passed. The private run report and supplementary replay totals are under `results/paired-reporter-comparison-v1/`.

## Reproduce

First use the [original preparation commands](README.md#reproduce) to create the frozen inputs, or reuse the exact snapshots from a completed run. Use the same model environment with the pinned ESM-2 checkpoint cached. Install `pyrodigal==3.6.3` there for the gene-caller arm. Pyrodigal is an optional research dependency distributed under GPLv3. [Pyrodigal release](https://pypi.org/project/pyrodigal/3.6.3/)

```sh
source_run=results/paired-reporter-pilot-v0
run=results/paired-reporter-comparison-v1
protocol=experiments/paired_preprocessing/reporter_protocol.json
pilot=experiments.paired_preprocessing.reporter_pilot
mkdir "$run"
for method in minimal start_aware start_ranked gene_caller; do
  python -m "$pilot" run --method "$method" --protocol "$protocol" \
    --inputs "$source_run/prepared/method_inputs/constructs.fasta" \
    --references "$source_run/references.private.json" --out "$run/$method.json"
done
python -m "$pilot" run --method annotation_reference --protocol "$protocol" \
  --inputs "$source_run/prepared/evaluator/annotation-units.json" \
  --references "$source_run/references.private.json" --out "$run/annotation_reference.json"
python -m experiments.paired_preprocessing.verify --protocol "$protocol" \
  --answer-key "$source_run/prepared/evaluator/answer-key.json" --baseline "$run/minimal.json" \
  --candidate "$run/start_aware.json" --candidate "$run/start_ranked.json" \
  --candidate "$run/gene_caller.json" --annotation-reference "$run/annotation_reference.json" \
  --out "$run/report.json"
```
