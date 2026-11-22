# Construct preprocessing for function-based screening

**Active research scope — 22 November 2026**

This document defines the active research direction and a high-level baseline for future planning, not a detailed implementation specification. It supersedes the toolkit-first research scope of 16 November 2026. The existing toolkit contract and release assessment describe the completed annotation-aware toolkit; historical results retain their original meaning.

## Purpose

Develop and evaluate preprocessing approaches that turn complex synthetic DNA constructs into units that existing biological foundation models can characterize. The project asks whether providing models with better-prepared units improves their ability to identify the documented components of a construct.

The motivation is defensive: sequence-similarity screening can be challenged when sequences are distantly related to known references. Foundation models may provide useful complementary information, but their value for complex constructs must be tested rather than assumed.

## Research questions

1. Can preprocessing approaches recover useful model-ready units from synthetic constructs without relying on their component annotations?
2. Does using those units improve identification of documented components compared with a clearly defined, minimally processed input at comparable query or compute budgets?
3. In the primary low-similarity evaluation, can model-based characterization identify components distantly related to its references, and how does it compare with conventional methods? How do results vary across construct types?

The project should report negative or mixed results as well as improvements. Identifying a protein or family is a measurable form of characterization; it does not by itself establish biological activity, risk, or a complete functional description.

## Scope of the initial research

The initial work focuses on protein-coding components in public, benign synthetic constructs. Constructs may contain multiple components and noncoding sequence, but evaluating noncoding functions requires its own methods and is not an initial acceptance requirement.

The intended eventual input is a construct whose contents are not already known to the preprocessing method. For evaluation, use constructs with documented components: retain their annotations as an answer key and withhold them from the methods being tested. The initial outcome is the fraction of eligible documented components identified. Define component eligibility and what counts as a correct identification before evaluating results.

Evaluate a few hundred public constructs if the available data and project resources support it, spanning the proposal's four archetypes: single coding sequence, multi-gene cassette, tagged or fused protein, and backbone-heavy construct. Begin with a smaller pilot and report any limitations in coverage or sample size.

The project integrates existing tools and pretrained models. Developing or training a new foundation model is not required. Specific preprocessing methods, models, baselines, matching rules, and datasets are intentionally left open for independent planning and evaluation.

## What the existing work contributes

The existing annotation-aware toolkit can help prepare traceable evaluation material: it checks documented coding translations against the corresponding DNA and exports supported protein units, while recording cases it cannot validate. The documented annotations remain an imperfect and potentially incomplete answer key.

Annotation-derived units can also provide a comparison showing how the downstream characterization method performs when documented component boundaries are available. This is an upper-reference condition for those annotations, not a complete ground truth or a sequence-only preprocessing solution.

The existing public corpus, provenance work, model evaluation harness, and prior fusion findings may be reused where appropriate. Existing results retain their original scope and do not establish the outcome of the new research questions.

The earlier component study found that increased coverage came from querying more units, with approximately 3.2 times as many model calls; per-unit retrieval was unchanged across those conditions. Its fusion analysis also found that mean-pooled embeddings missed some documented families, while providing the correct segment as an oracle recovered them. This motivates testing segmentation without access to the answer; the oracle result does not establish that a practical method can find the segment.

## Comparison requirements

The primary preprocessing comparison must use matched, predefined query or compute budgets, or report component coverage across budgets so that approaches can be compared at the same budget. Define the primary budget measure and candidate-selection rules before evaluation. Also report model calls, input size, preprocessing cost, and other relevant resource use: a call count alone may not capture actual cost. A matched comparison assesses the complete input strategy, including candidate selection, rather than isolating representation alone.

Specify in advance a minimum improvement that would be practically interesting at a comparable budget, together with the uncertainty and cost measures used to interpret it. Prior study estimates provide context, not an automatic threshold for this different experiment. Count the same eligible documented components in each condition; report unsupported predictions separately when assessed, without treating every undocumented component as necessarily incorrect.

## Workstreams and parallelization

The following streams can begin independently. Their teams should agree on minimal interfaces and evaluation definitions early, while retaining freedom to investigate different implementations.

### 1. Evaluation and verification

Select and assess documented public constructs; establish eligible component labels, provenance, and data splits; define how predictions are matched to documented components; and build reproducible scoring and reporting. Maintain the separation between the withheld answer key and the sequence-only inputs used by other streams.

**Output:** a documented evaluation set and scoring process that can assess any compatible approach. Verification of candidate approaches is independent of their development, although the final comparison requires their outputs.

### 2. Preprocessing approaches

Develop and compare candidate ways to derive model-ready units from unannotated construct sequence. Multiple approaches can be explored in parallel. Preserve enough information to trace each candidate back to its source construct and understand failures.

An early hypothesis is that segmenting candidate proteins may improve characterization of fused or tagged components, motivated by the prior oracle-cut finding. Test any segmentation approach without providing the known boundaries or reference-aligned answer. Keep its implementation open and compare it with other candidate approaches. Fusion constructs are a useful test of this hypothesis; backbone-heavy constructs can test behavior when much of the input is unrelated to the documented components, including candidate and unsupported-call burden.

**Output:** one or more independently testable preprocessing approaches with documented inputs, outputs, and limitations.

### 3. Baseline and model characterization

Define a minimally processed, model-compatible comparison and a consistent way to characterize its inputs and those produced by the preprocessing approaches. Identify and document any model input or length constraints. Explore relevant existing models and comparators without requiring a particular architecture or implementation.

**Output:** a reproducible baseline and a characterization workflow applicable across conditions.

### 4. Low-similarity evaluation

Develop an independently supported evaluation of the primary low-similarity hypothesis: that model-based characterization can identify components distantly related to the available references, compared with conventional methods. Establish how family or component labels are supported independently of the tested retrieval result, how similarity and coverage are measured, how close relatives are separated, and what claims the resulting evidence can support. Higher-similarity cases provide a sanity check, not a substitute for this evaluation.

Plan for a separate public, benign low-similarity dataset in parallel with assessing the construct corpus; do not assume the corpus contains sufficient eligible cases. First evaluate characterization on known units if appropriate, then connect it to the construct-level pipeline where the data allow.

**Output:** an executed low-similarity evaluation with documented labels, reference separation, comparators, results, and limitations. If it cannot be completed, report this research goal as unanswered.

## Integration and next steps

1. Establish shared input, prediction, and scoring interfaces so each workstream can proceed independently without exposing held-out annotations.
2. Prepare a small pilot dataset and a documented verification process while preprocessing approaches and the baseline are developed in parallel.
3. Run a paired pilot on the same constructs: minimally processed inputs, candidate preprocessing outputs, and annotation-derived units. Compare identification of the same eligible documented components at matched budgets or across budget curves, and record resources and candidate counts.
4. Examine where failures arise, distinguish missed units from characterization errors, and use those findings to guide independent follow-up work.
5. Expand the evaluation across construct archetypes and complete the low-similarity study. Report results, uncertainty, exclusions, and limitations.

The paired comparison depends on outputs from the verification, preprocessing, and baseline/model streams. The low-similarity stream can prepare its data and methods in parallel, then connect to the common evaluation where appropriate.

## Minimum research outcome

A completed initial study should provide a reproducible, budget-aware comparison showing whether one or more preprocessing approaches improve identification of eligible documented protein components relative to a defined minimally processed baseline. It should include an annotation-derived comparison, results across the construct types the dataset supports, and a clear account of candidate and computational costs. It must also execute the low-similarity evaluation against conventional methods; if this is not achieved, that goal remains incomplete rather than being counted as satisfied by an explanation.

The outcome may be positive, negative, or mixed. The research result is distinct from the existing toolkit's release readiness.

## Later extensions

After the initial study, possible extensions include incomplete or fragmented inputs, sub-segment reference methods, a shareable benchmark, and defensive robustness testing. Their designs and priority should follow the evidence from the initial work.

## Boundaries

Use public, benign material and preserve source provenance. Do not use proprietary or customer data. Any later robustness testing should use benign proxies or appropriately usable published benchmarks. Model outputs are research aids, not autonomous safety decisions. Preserve historical reports, frozen data, and their provenance as records of the work actually performed.
