# Group-aware reporting for paired recall

The paired verifier counts each eligible documented component. This is the original pooled recall endpoint. When components or constructs recur in related records, a second view helps show how much of the result comes from one group. The group report uses the same verifier output and an evaluator-provided group manifest. Both inputs are metadata.

The manifest schema is `paired-recall-groups-1`. Its `constructs` list contains `id`, `group`, and `partition` for each construct. Use the same opaque IDs as the verifier. A manifest can include records from other partitions to check that a group stays together. Every scored construct needs a group. Each comparison uses one partition.

```sh
python -m experiments.paired_preprocessing.grouped_report \
  --report results/example/report.json \
  --groups results/example/groups.json \
  --out results/example/grouped-report.json
```

For each candidate and query cap, the report gives the original pooled recall difference, the mean of the differences with each eligible group weighted equally, and the share of eligible components in the largest group. It also gives the range of pooled differences after omitting each group in turn. That range shows sensitivity to group membership and is descriptive. It has no value when only one group contributes eligible components.

Group membership comes from independent evaluator evidence about shared components and construct lineage. This tool checks the supplied partition assignments and paired counts. The group counts and sensitivity result depend on the evidence used to assign those groups. Keep the group manifest with the evaluator artifacts.

## Reporter pilot replay

For a replay of the five-construct comparison, group the two StayGold entries together by their reported reporter lineage and assign the other records separate groups. This illustrates the calculation on four groups. It does not establish that those groups form an independent holdout. [Reporter annotations](reporter_panel.json) · [Source comparison](https://github.com/khgng026/iri-construct-screening-2026/pull/18)

At two queries per construct, Pyrodigal's 42.9 percentage point pooled gain becomes 41.7 points with equal group weights. Its gain ranges from 25 to 60 points after omitting each group in turn. Either start-based arm gains 14.3 points pooled and 8.3 points with equal group weights. Omitting the single group where those arms improve brings their observed gain to zero.
