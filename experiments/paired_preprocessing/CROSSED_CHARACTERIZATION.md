# Crossed reporter characterization

The reporter pilot found higher recall with gene calling at the same query cap. This follow-up asks whether the gain depends on ESM-2 by crossing the existing preprocessing arms with two scorers. The five reporter constructs and seven eligible components are the same as in v1, so this is an exploratory replay. [Earlier comparison](REPORTER_COMPARISON.md)

## Comparison fixed before scoring

- Reuse the original inputs, reference panel, query caps of one, two, and four, and limit of five family calls. The primary cap is two queries per construct. Include minimal translation, start-based filtering and ranking, Pyrodigal, and the annotation-derived comparison.
- Compare mean-pooled ESM-2 cosine with local BLOSUM62 using the previous alignment settings and reference self-score normalization. Both scorers use the best reference score per family and the same turn-taking rule to assemble parent-level calls.
- Check the input and reference identities and a fingerprint of the selected protein inputs, source spans, and order before comparing recoveries. Evaluate both sets of traces with the same held-out answer key. Report differences in preprocessing gain and counts of components recovered only by one scorer. Query caps constrain the number of protein inputs. Relative runtime is a separate measurement.

## Reproduce

Use the cached inputs described in the [pilot handoff](../../docs/project/HANDOFF.md) and the model environment with Pyrodigal from the [v1 comparison](REPORTER_COMPARISON.md). The alignment arm uses Biopython without loading ESM-2. The traces and reports stay in ignored results storage.

```sh
source_run=results/paired-reporter-pilot-v0
run=results/paired-reporter-characterizers-v2
pilot=experiments.paired_preprocessing.reporter_pilot
protocols=experiments/paired_preprocessing
mkdir -p "$run/esm2" "$run/blosum62"
cp "$protocols/reporter_protocol.json" "$run/esm2/protocol.json"
cp "$protocols/reporter_alignment_protocol.json" "$run/blosum62/protocol.json"
for scorer in esm2 blosum62; do
  for method in minimal start_aware start_ranked gene_caller annotation_reference; do
    input="$source_run/prepared/method_inputs/constructs.fasta"
    if [ "$method" = annotation_reference ]; then
      input="$source_run/prepared/evaluator/annotation-units.json"
    fi
    python -m "$pilot" run --method "$method" --protocol "$run/$scorer/protocol.json" \
      --inputs "$input" --references "$source_run/references.private.json" --out "$run/$scorer/$method.json"
  done
done
python -m experiments.paired_preprocessing.compare_characterizers \
  --esm2-run "$run/esm2" --blosum62-run "$run/blosum62" \
  --answer-key "$source_run/prepared/evaluator/answer-key.json" --out "$run/report.json"
```
