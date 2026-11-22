# Research plan

The active [scope](SCOPE.md) defines the goal and evaluation requirements. This plan organizes the next work without prescribing implementations. The annotation-aware toolkit is complete under its own [contract](TOOLKIT_CONTRACT.md); its release is distinct from the research study.

## Parallel workstreams

1. **Evaluation and verification:** Prepare public, benign annotated constructs, eligible held-out component labels, and reproducible matching and scoring. Agree on shared interfaces without exposing labels to methods under test.
2. **Preprocessing:** Develop and compare sequence-only approaches. Investigate segmentation as an early hypothesis motivated by the prior fusion result, without using oracle boundaries. Other approaches can proceed independently.
3. **Baseline and characterization:** Define a model-compatible minimally processed baseline and a consistent characterization workflow. Establish model constraints, budget measures, and relevant conventional comparators.
4. **Low similarity:** Prepare and execute a properly labeled low-similarity comparison against conventional methods. Plan for a separate benign dataset while assessing whether the construct corpus can support this question.

## Integration milestones

- Establish shared interfaces, eligibility and scoring rules, matched-budget or budget-curve comparisons, and a minimum practically interesting effect before evaluation.
- Prepare a small pilot while preprocessing, verification, and baseline work proceed in parallel.
- Compare baseline, candidate approaches, and annotation-derived units on the same eligible components. Diagnose missed components, characterization errors, and cost.
- Expand across the four construct archetypes as resources and eligible data allow, and complete the low-similarity evaluation. Report uncertainty, limitations, and negative or mixed findings.

The paired comparison depends on the first three workstreams; the low-similarity work can progress independently and connect to the end-to-end study where data allow. No unperformed study or positive result is implied by this plan. See [status](STATUS.md) for current progress and [history](HISTORY.md) for previous plans and studies.
