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
