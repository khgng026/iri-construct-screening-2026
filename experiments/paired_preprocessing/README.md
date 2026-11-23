# Paired preprocessing verification

The verifier takes a comparison protocol, a held-out answer key, and traces from the methods under test. Each arm uses the same construct IDs, candidate budgets, output cap, and characterization workflow. Sequence arms use the same input bytes. The annotation-derived arm is an optional diagnostic in which the characterizer receives the validated source proteins.

## Endpoint

Score the fraction of eligible documented components recovered at each query cap. A call must predict the documented family from a logged model input on the correct strand that overlaps the minimum fraction of the component's documented bases. Exact boundaries are unnecessary. Each prediction and each eligible component can be credited once. An unlabelled component is outside the denominator. Distinguish misses with no covering query from misses whose locus was covered. Locus coverage alone says nothing about whether the query translated the right protein. Report failures, number of calls, number of residues, and timing alongside recall. The same caps apply to each arm, and actual usage may differ if an arm has fewer candidates. Treat the effect target as a practical threshold for the point estimate, with uncertainty considered separately.

`verify.py` runs with the Python standard library. The [pilot protocol](reporter_protocol.json) illustrates the comparison settings. A key has `schema` and `constructs`, with `id`, `archetype`, `length`, and eligible `components` in every construct. A component has `id`, `family`, `strand`, and `spans`. `spans` are zero-based and half-open on the forward sequence, including all segments for a joined feature. A trace has `schema`, `method`, `characterization_id`, and `tiers`. Each tier contains `query_cap` and a `constructs` array. Each row lists `id`, `status`, `candidate_count`, `preprocess_seconds`, `characterize_seconds`, `queries`, and `predictions`. A query has `id`, `spans`, `strand`, and `input_residues`. A prediction has `query_id` and `family`.

Use [group-aware reporting](GROUPED_REPORTING.md) with the verifier output to summarize dependence on related constructs.

## Reporter pilot v0

Use five deposited reporter constructs covering a single CDS, a cassette with multiple CDSs, two fusion or tagged examples, and a construct with mostly sequence outside its documented reporter. The [panel](reporter_panel.json) records the source snapshots, the seven eligible coding components, and the reviewed UniProt references. The answer key comes from the source product descriptions and validated translations, with six reporter labels and one regulatory protein label. Treat each fusion as one documented coding product. The panel does not supply independently labelled fusion subregions.

- **Minimal conversion** uses stop-delimited translations in all six reading frames, with a fixed model input length rule and candidates ordered by length.
- **Start-aware conversion** uses the same rules but starts each eligible block at its first methionine. Orders and ties are deterministically resolved before characterization. This prototype follows linear coordinates and standard translation for the case study.
- **Annotation-derived** uses the separate validated protein bundle, with the same query caps and output rule. This measures characterization when the documented units are available.

All arms use the same ESM-2 8M model, reference panel, and family ranking with turn taking across the queried units. Compare at caps of one, two, and four queries per construct, with two as the primary cap. Return at most five calls. A gain of five percentage points is the practical target. With seven eligible components, a single extra recovery changes recall by about fourteen points. This pilot is descriptive.

The [v1 comparison](REPORTER_COMPARISON.md) adds candidate ranking and a Pyrodigal arm on the same inputs and caps.

## Reproduce

Use a model environment with the project installed, Biopython, and the pinned ESM-2 checkpoint cached. Run acquisition separately from inference. Exact reruns require the source snapshots recorded in the panel. Each command writes to a new or unused output path under ignored `results/`.

```sh
run=results/paired-reporter-reproduction
protocol=experiments/paired_preprocessing/reporter_protocol.json
pilot=experiments.paired_preprocessing.reporter_pilot
python -m "$pilot" acquire --out "$run"
python -m "$pilot" prepare --sources "$run/sources" --out "$run/prepared"
for method in minimal start_aware; do
  python -m "$pilot" run --method "$method" --protocol "$protocol" \
    --inputs "$run/prepared/method_inputs/constructs.fasta" \
    --references "$run/references.private.json" --out "$run/$method.json"
done
python -m "$pilot" run --method annotation_reference --protocol "$protocol" \
  --inputs "$run/prepared/evaluator/annotation-units.json" \
  --references "$run/references.private.json" --out "$run/annotation_reference.json"
python -m experiments.paired_preprocessing.verify --protocol "$protocol" \
  --answer-key "$run/prepared/evaluator/answer-key.json" --baseline "$run/minimal.json" \
  --candidate "$run/start_aware.json" --annotation-reference "$run/annotation_reference.json" \
  --out "$run/report.json"
```
