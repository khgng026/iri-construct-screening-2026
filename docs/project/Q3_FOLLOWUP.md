# Q3 follow-up results

## Answer supported by the experiments

Protein language model scores added information on some distant, independently
labeled protein domains. On whole benign constructs, their end-to-end value
depends on getting a useful candidate into the scoring budget, on the reference
library, and on whether the task is structural or functional annotation. We
found an exploratory benefit from combining scores in the older whole-construct
cohort. The strictly distant whole-construct pilots remain too small to
establish a reliable incremental gain over strong conventional methods. The
results measure family or structural assignment on benign examples.

The original proposal's input-format premise is supported. A fresh archive
confirmed improved access to documented components at a fixed budget. The
proposal's expectation of general model superiority at low similarity is more
conditional. The most informative next comparison is between a strong
conventional scorer and the same scorer with model evidence on identical
whole-construct candidates. An isolated-domain result can then explain a
whole-construct result.

## 1. Available holdout and paired whole-DNA audit

The sequence-free [feasibility audit](../../results/q3-followup/evidence-audit.json)
found **217 distinct translations in 1,221** cohort v2 corroborated component
occurrences. The original low-similarity tier has **101 distinct translations
in 250 occurrences**. With the same coverage-qualified local identity metric,
the below-25% subset contains **seven occurrences of one translation**, even
though six original lineage groups contain them. The separate
insufficient-coverage tier has 21 distinct translations across 37 occurrences
and four families. These tiers are distinct from the global identity
definitions used in subsequent pilots. Exact protein fingerprints in the
original development and holdout sets did not overlap.

The three later whole-DNA panels have **14** occurrences from **three** reporter
ancestries below 25% global identity and without a sufficiently covering
same-label local hit. At two candidate inputs, gene calling reached 13 loci.
Raw alignment and local ESM-2 8M each recovered the same 13. Local ESM-2 150M
recovered six, as did the profile comparison. Neither model added a correct
whole-DNA call when both conventional methods missed. With known protein
boundaries, six model-only GUS occurrences represent one distinct domain
sequence. They are related examples, and the modelling conditions and
reference libraries vary by panel.

Even with genuinely independent, exchangeable units and zero observed model
rescues, 14 units would leave a one-sided 95% upper bound of 19.3% on the
rescue probability. It would take 59 such zero-event units to get that bound
below 5%. The existing three selected ancestries are more dependent than
that illustration. Use paired disagreements and independent lineages to plan
the eventual effect estimate.

## 2. Matched-domain robustness replay

The [replay](../../results/q3-followup/domain-replay.json) paired all 128 original
below-25% domains with their results under the wider reference policy from
[PR #28](https://github.com/khgng026/iri-construct-screening-2026/pull/28).
Local ESM-2 150M recovered 99 and 98, respectively. The best single conventional
scores were 86 and 85. The paired difference depends on whether the
comparator is a single method or the union of their correct calls.

| Same 128 domains | Original references | Wider references |
| --- | ---: | ---: |
| Correct calls among the union of all four conventional scorers | 111 | 109 |
| Model right when all four conventional scorers missed | 8 | 10 |
| Any conventional scorer right when the model missed | 20 | 21 |
| Correct calls among the union including the model | 119 | 119 |

Six model-only recoveries repeat across the two reference policies. The same
eight and ten recoveries remain in the 113-domain subset without sufficiently
covering same-label matches under either policy. These are robustness checks
on one set of domains. The optimistic unions use the labels to select the
answer. Actual performance of a combination requires a preselected decision
rule and independent evaluation.

## 3. Fresh source and candidate reach

The [Open Reporters audit](../../results/q3-followup/open-reporters.json)
validated coding-region translations on 28 public circular cloning vectors.
Nineteen GFP-like proteins across 11 source genera matched an independently
published natural protein exactly. Their exact fingerprints were absent from
the earlier cohort. All records belong to one collection and share the
reported cloning architecture.

At two candidates per record, Pyrodigal reached **19/19** of the independently
corroborated loci versus **12/19** for minimal six-frame preprocessing.
None of the covering inputs was exactly the reporter peptide. The follow-up
diagnostic found that every covered case contained that entire exact peptide
inside the longer input. This supports candidate reach and shows that
coordinate coverage and native protein boundaries are different diagnostics.
The five other catalogue-labeled reporters need independent validation
before they can serve as verified out-of-family controls.

There were **zero** qualifying strict low-similarity components with the
combined reference library. Two qualify in the copepod-only sensitivity
analysis, both reached by Pyrodigal and one by minimal preprocessing.
This source provides a useful additional construct background for
preprocessing and calibration. It cannot supply the primary Q3 test panel.
Eighteen of the exact-match donor accessions appear in Alieva et al. (2008)
and the remaining accession is indexed with Matz et al. (1999). All 19 were
publicly documented before the ESM-2 pretraining snapshot. This excludes
them from a temporal holdout, regardless of the vectors' release date.
Actual model exposure remains unknown.
[Open Reporters Collection](https://github.com/Reclone-org/Open-DNA-Collections/tree/55d4b762194055c8027fd6c780068d8a78a03282/Open%20Reporters%20Collection) · [Natural protein source](https://doi.org/10.7554/eLife.75842) · [Alieva et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2481297/#pone-0002680-t001) · [Matz et al.](https://pubmed.ncbi.nlm.nih.gov/10504696/)

## 4. Historical rank combination replay

On identical Pyrodigal candidates from the original 777 constructs, a fixed
equal-weight combination of stored mean-pooled ESM-2 8M and normalized
BLOSUM62 ranks improved recall over normalized BLOSUM62 alone. Both arms
selected at most one label per candidate. This output cap differs from the
original cohort v2 analysis.

| Cohort v2 tier | Components | Lineage-weighted difference | Paired 95% group-bootstrap interval |
| --- | ---: | ---: | ---: |
| All | 1,221 | +1.9 points | +0.5 to +3.4 |
| Original low similarity | 250 | +5.9 points | +1.3 to +10.9 |
| Insufficient same-family coverage | 37 | +35.0 points | +20.0 to +51.7 |

Under the documented-locus diagnostic, 736 of 1,537 assignments had a
single unambiguous documented component covering the queried locus by
the benchmark rule. In that subset, incompatible assignments dropped
from 28 to 11. The other 801 assignments were unknown under that definition.
The 37-component effect is descriptive and concentrated under conservative
grouping. The replay uses public answer keys, an older model representation,
and a single historical conventional comparator. It motivates testing combined
scorers. It does not measure open-set false assignments or establish benefit
over raw alignment with well-supported profiles.
[Sequence-free output and source fingerprints](../../results/q3-followup/rank-replay.json)

## Status and reproduction

The [primary comparison protocol](Q3_EVALUATION_PROTOCOL.md) is ready to
freeze once independent positive and negative groups and a calibration split
are available. Direct access to NCBI, UniProt and the pinned model distribution
was blocked by the current environment policy, so these runs used committed
model decisions or public records already distributed through their authors'
repositories. No new ESM-2 inference was run. The wider independently held-out
test remains open.

The evidence audit and the matched-domain replay require only Python and
the tracked JSON inputs. The rank replay uses NumPy. The reporter audit
used Biopython 1.86 and Pyrodigal 3.6.3. Reproduce from this repository root.
The two source checkouts for the reporter audit must be at the commits in
[the frozen manifest](../../experiments/q3_followup/open_reporters.freeze.json).

```bash
run=results/q3-followup-replay
mkdir -p "$run"
python -m experiments.q3_followup.evidence_audit --out "$run/evidence-audit.json"
python -m experiments.q3_followup.rank_replay --out "$run/rank-replay.json"
git fetch origin refs/heads/study/independent-low-identity
for file in fold_results.json fold_wider_reference_results.json; do
  git show "4171d33bc0f7404b652f17c45c599d67e5c6a23f:experiments/benign_protein_pilot/$file" > "$run/$file"
done
python -m experiments.q3_followup.domain_replay \
  --original "$run/fold_results.json" --wider "$run/fold_wider_reference_results.json" \
  --out "$run/domain-replay.json"
python -m experiments.q3_followup.open_reporters \
  --construct-repo /path/to/Open-DNA-Collections \
  --protein-repo /path/to/Orthologous_GFP_Fitness_Peaks \
  --out "$run/open-reporters.json"
```

The locked run descriptions are [holdout and whole-DNA audit](Q3_EVALUATION_PROTOCOL.md),
[domain replay](Q3_DOMAIN_REPLAY_PROTOCOL.md),
[Open Reporters audit](Q3_OPEN_REPORTERS_PROTOCOL.md), and
[historical rank replay](Q3_RANK_REPLAY_PROTOCOL.md).

All four output files reproduced byte for byte. The repository test suite
ran 322 tests with four skips.
