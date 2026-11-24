# Q3 follow-up protocol

## Primary research question

Given archived benign construct DNA and a fixed candidate-input budget, what is the paired change in correct low-similarity component recall when a fixed residue-level ESM-2 scorer is added to a conventional comparison using raw local alignment and profile HMMs? The comparison concerns assignment to independently documented protein families. Structural superfamily recognition and experimentally documented exact reaction are distinct secondary endpoints.

## Stage 1. Audit the available evidence

Use the public, sequence-free exports on the project branches. This is retrospective. The underlying studies already reported their headline results and the answer keys are visible.

- In cohort v2, count occurrences, distinct protein fingerprints, declared lineage and conservative groups, families, and representation by construct archetype. Use corroborated components, and show the original low-similarity tier, cases below 25% identity within that tier, and the insufficient-coverage tier separately. The cohort uses a coverage-qualified local identity, unlike the later global-identity panels. A case without a sufficiently covering same-family match is missing that local identity by design. Count overlap with the original development set by exact fingerprint and report concentration under both existing groupings.
- For the whole-DNA pilots, read stored per-target decisions on cases below 25% global identity with no same-family local alignment covering at least 70% of each protein. Use the alternate reporter library, Renilla, and GUS as the three primary panels. The original reporter library is an additional sensitivity check on some of the same records. At the recorded primary cap of two candidates, measure candidate coverage and pairwise scorer wins and losses and compare with the known-boundary run. Show each panel separately and count shared reporter ancestries across the primary panels. Correct attribution across layout requires independent ancestries within a layout.
- Compare each fixed ESM-2 local scorer against raw BLOSUM62 and the profile separately. Also compute the number of cases it got right where both conventional scorers missed. A union of correct outcomes uses the answer key and is an optimistic diagnostic for any procedure choosing among the recorded assignments. It is not an achievable policy or a comparison at a measured false-assignment rate.
- Report how many independent units would be required to bound a missing effect. Illustrate the precision of zero observed rescues with a one-sided exact binomial bound under independent, exchangeable units. This is a planning illustration, not an interval for the existing selected reporter ancestries.

## Stage 2. Freeze a prospective comparison once a sufficient panel exists

- Acquire only archived public benign reporter or ordinary model-organism expression records, with source hashes and verified translations. Prefer independent published structural, functional, or source evidence for the label. Prefer exact sequences first appearing after the ESM-2 training snapshot, document relatedness to older public proteins when sequence history is available, and report unresolved exposure. Keep evaluator files apart from scoring inputs. Construct submission groups and protein ancestries must remain on one side of a development/test split.
- Fix the reference collection and family labels independently of either scorer, with plausible competing labels and documented out-of-family examples. Define similarity using a global identity band and a coverage-qualified local sensitivity analysis. Combine construct layouts only if the number and diversity of independent groups support the prespecified contrast.
- Fix candidate generation and budget, and give the same selected proteins to raw local alignment, profiles, and the model-assisted comparison. Pick one conventional policy and one model-assisted policy on development data, both with an option to abstain. Freeze all decisions and the analysis before final inference. Score a documented component as recovered only when both its family and the matching locus are correct. Count assignments to deliberately verified out-of-family examples separately from assignments to unannotated loci. Measure scoring cost and report which misses occurred before scoring.
- Use the paired change in lineage-weighted end-to-end recall as the primary estimate. Compare it at an agreed false-assignment budget and with a prespecified minimum practically relevant gain, supported by a precision calculation using planning data. Report the full paired interval and case denominators, both candidate-reaching and known-boundary diagnostics, and sensitivity to the weighting of lineages. On an underpowered panel report feasibility and descriptive evidence.

### Information needed to freeze the primary test

Record positive target counts, distinct protein and background group counts,
the number of plausible competing families, and the number of independently
verified negative records in each split. Include a column for the provenance
and earliest verified availability of each distinct protein. Keep accession
identity, sequence identity, and submission relatedness separate. Report
both availability of references and number of unrelated positives receiving
a usable candidate before selecting a decision policy.

For planning, let the paired outcome for an independent unit be +1 if only
the model-assisted method is right, −1 if only the conventional method is
right, and 0 otherwise. If the true mean difference is d and the chance of
a disagreement is q, its variance is q − d². For example, a 10-point true
gain and discordance of 10%, 20%, or 30% give approximate independent-unit
requirements of 139, 292, or 446 for a two-sided 95% interval with half-width
of five points. A 5-point practical margin is a planning example. These
numbers do not incorporate dependence and are not enrollment targets.
Estimate the variance of the planned per-lineage differences using
development results, account for related layouts and families, and lock a
capacity-supported objective before analyzing the held-out results.

Freeze the per-record definition of a false assignment before calibration.
Choose each method's abstention rule on the development split using the same
budget. Report false assignments and uncertainty on verified holdout negatives
alongside recall. If test false-assignment burden differs materially, report
the trade-off and withhold a claim of improved recall at matched burden.
Benchmark decisions and thresholds apply to the declared benign annotation
reference collection.

This protocol studies annotation on benign public data. The endpoint does not assess activity, risk, or operational screening performance.
