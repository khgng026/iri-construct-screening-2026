# Experimental results and research handoff — November 12, 2026

Authority: [actual ChatGPT delegated decision](https://github.com/khgng026/iri-construct-screening-2026/pull/11#issuecomment-6800553713), reviewed commit `32e0928`, evidence freeze `c90d063f…47c2e`. All original frozen inputs and criteria remain unchanged. Raw sequences, qualifiers, embeddings and full diagnostics remain private/ignored. This summary supersedes earlier readiness statements; their frozen audit trails are preserved.

| Executed experiment | Actual observation |
| --- | --- |
| E1, seven new sources | 67/67 observable features and qualifiers; 15/15 source translations; 74/74 original feature coordinates, strands and parent links agree. Zero parser failures/inventions. |
| E1, separate historical pBR322 repeat | 48/48 observable features/qualifiers, 3/3 translations; 50/50 original feature coordinates/strands agree. Not a novel independent source. |
| E3 descriptive partial-method pilot, tier 1 | Alignment and pinned ESM-2 each: 4/4 queries scored, family hit/recall@5 = 4/4, 20 returned hits, three wrong-family hits (3/20). Identical eligible panels and five total hits/query; different rankings, equal aggregate endpoints. |
| E2 descriptive supported subset, AF298787 | Longest eligible URA3: supported-family recall 2/2; all eligible URA3+EGFP3: 1/2. Five total parent hits in each. Full documented-family endpoint remains null because AmpR is unsupported. |

[E1](results/cleared-2026-11-11/E1.json), [E3](results/cleared-2026-11-11/E3.json) and [E2](results/cleared-2026-11-11/E2.json) retain source-specific denominators, rankings, panel fingerprints, environments, measurements and full local report digests. Source admission and the exact original reviewer decision are recorded in [the authorization record](oversight/EXPERIMENT_RUN_AUTHORIZATION.json).

## Interpretation and limits

E1 is agreement with deposited annotations, not biological truth. All seven new records are well annotated; uncertainty/justified-abstention denominators are zero and those rates remain null. Bookkeeping features are excluded from observable E1 denominators. Vector/protein duplicates do not become independent family observations.

E3 covers four components from three source parents, connected into two conservative dependence clusters. Confidence intervals requiring eight clusters remain null. Both methods reached the family-hit ceiling in this small catalogue; this does not establish equivalence, a winner, function prediction or generalization. Wrong-family hits are molecular-family label disagreements, not operational false-positive rates.

E2 preserves the frozen maximum-per-reference merge rule. GFP hits dominate the merged top five, displacing OMP-decarboxylase hits. The URA3-only baseline's two-family recovery includes off-target GFP-family matches that happen to belong to the parent's expected family set; it does not show that URA3 represents GFP. The per-unit URA3 control is unchanged. Thus this observation exposes a limitation of score merging and the family-recovery endpoint, rather than general evidence against decomposition. No scoring rule was changed after viewing the outcome.

Literal source-longest is unsupported AmpR (286 aa), while eligible-longest is URA3 (267 aa); EGFP3 is 238 aa. Two of three source CDSs have accepted supported labels. AmpR remains in source/unknown-family coverage and leakage exclusions. Full construct-family recall and population-level claims are unavailable.

## Exploratory E2 diagnosis and bounded improvement

[Component attribution](results/analysis-2026-11-12/E2.attribution.json) explains the unexpected result without changing E2: the baseline returns three correct URA3-family hits and two off-family GFP hits, while all five merged hits originate from GFP and match GFP. Each condition therefore recovers **one** family through its corresponding component, despite primary parent-family recovery of 2/2 versus 1/2. The actual representations still contain one versus two eligible documented proteins.

The eligible catalogue has six GFP, three OMP-decarboxylase, five beetle-luciferase and one Renilla reference. All six GFP references outrank the first OMP reference in the maximum-score merge; OMP first appears at **K=7**. This diagnoses the fixed five-hit budget and asymmetric score scales; changing K is an exploratory sensitivity, not a revised result. Cosine scores are already vector-normalized but are not calibrated across query proteins. Componentwise rank normalization recovers both corresponding families at K=5; it discards magnitude and can also promote poor matches. The predeclared alignment self-normalization changes two of URA3's top-five identities and none for the other three queries, demonstrating that normalization also changes method rankings.

One bounded improvement, [stable round-robin component merging](results/analysis-2026-11-12/balanced_merge.prospective.json), allocates the **same five total unique references** across component rankings. The algorithm and prospective status were recorded before its cached exploratory calculation. On this one parent it recovers both corresponding families; one later off-family hit remains. This is a mechanism check after observing the original outcome, not evidence of general benefit or a new confirmatory benchmark. Original maximum-score merging and all frozen scores remain unchanged. The focused regression checks shared budget, deduplication, stable ordering, component origin and invalid scores. Prospective primary evaluation requires delegated clearance.

Reproduce diagnostics with `.venv/bin/python scripts/analyze_cleared_pilot.py --stage E2`; exact cached E3/E2 matrices must pass the immutable evidence check. No new inference, reference substitution or model calls occur.

## Actual computation and concise E4 evidence

| Actual work | Measured cost |
| --- | --- |
| E1 seven new sources / historical control | Stage wall 0.377 / 0.358 s; OS peak RSS 36.2 / 35.8 MB |
| E3 alignment | 0.0855 s measured alignment work; 79 calls including 17 declared self-normalization sensitivity calls; 62 primary query/reference pairs |
| E3 ESM-2 | 4.006 s cached initialization, 1.888 s encoding/pooling; 21 unique model calls, 6,964 residues, zero sequence-cache hits; sampled inference RSS 612 MB |
| E3 whole process | 9.758 s wall; OS peak RSS 615 MB |
| E2 ESM-2 shared execution | 2.703 s initialization, 1.980 s inference/pooling; 17 unique calls, 5,670 residues; 15 eligible references encoded once, baseline once, extra CDS once |
| E2 inference phase breakdown | References 1.8821 s, baseline 0.0484 s, additional CDS 0.0493 s; no independent cold run per condition |
| E2 whole process | 6.411 s wall; OS peak RSS 568 MB |
| Validated E3-score reuse for E2 | Same two condition results with zero additional model calls; 1.640 s bounded process, 36.2 MB OS peak RSS, including exact evidence verification |

CPU, two intra-op/three inter-op threads, float32, deterministic seed 20261105, cap 1022 residues, no model downloads. All stages stayed below the 6-GiB/1,800-second/64-protein/65,536-residue bounds. Initialization includes cached loading, not an independently controlled cold-system benchmark. OS/file caches were uncontrolled. Sampled RSS can miss peaks; bounded-process reports also check the OS high-water mark. GPU execution time and monetary cost are unknown, not inferred from wall time.

**E4 completed under [recorded scientific clearance](https://github.com/khgng026/iri-construct-screening-2026/pull/11#issuecomment-6815404726).** [Versioned execution manifest](results/analysis-2026-11-12/E4.execution.json) and [results](results/analysis-2026-11-12/E4.json) bind the exact approved E3 and E4 protocol hashes. Both methods recovered 4/4 families with 3/20 wrong-family hits. Measured ESM-2 encoding/pooling took **22.07 times** the alignment work (1.888 versus 0.0855 s), plus 4.006 s separately measured cached initialization. Alignment's 79 calls include 62 primary pairs and 17 self-normalization sensitivity calls. Sampled process RSS was 211 MB during alignment and 612 MB during ESM-2 inference; these are absolute process levels with different prior imports/work, not incremental method memory. Monetary cost and GPU time remain unknown. No CI, controlled cold-cache, general efficiency or winning-method claim follows. E2 remains separate.

[Derived archetypes v1](results/analysis-2026-11-12/archetypes.v1.json) restores the source-supported multi-CDS classification for all four E3 queries without changing historical scores. Reproduce this source-free E4 analysis with `.venv/bin/python scripts/analyze_cleared_pilot.py --stage E4`.

ESM-C scored retrieval, tiers 2/3, Renilla, other non-cleared sources and the supplemental fusion remain HOLD.

## Independent evidence and exact prospective review

The [bounded audit](results/analysis-2026-11-12/independence-audit.json) rechecks all 22 acquired source hashes, groups identical provided translations and retains the ten iGEM translation gaps and four GenBank holds. Reviewed-UniProt Renilla search returned only P27652, the known excluded ancestor: **no independent positive gained**. Ordinary luciferase activity from another molecular lineage is insufficient to relabel it.

One additional natural, reviewed monofunctional OMP-decarboxylase reference, **O58462** (Pyrococcus horikoshii, 208 aa), has a distinct sequence and natural taxonomic origin. Its family evidence includes curated homology assertions, not an independently demonstrated new biological result. The reference is **proposed, unscored**, with original JSON and sequence hashes, attribution and private-cache provenance.

The targeted pCAMBIA-1303 source **AF234299.1** adds a genuine GUS/GFP fusion architecture and four source-validated translations. Together with previously prepared **MK448012.1**, it supports two new E1-only proposals. Neither record supplies complete independently reviewed protein-family/boundary truth for a new decomposition retrieval comparison. Shared reporter ancestry/backbones/markers also prevent declaring two new independent clusters. No sequences or full qualifiers are committed.

**Single delegated review request:** review [panel v2](../data/manifests/prospective-2026-11-12/panel.v2.json), canonical proposal hash `f9111dd17ec028fdc170552aacc8af7275fccaf46d1e9712e3f7a681bdc0d70b`. Requested scope: E1 source-annotation preservation on the two exact supplemental records; E3 reference-composition sensitivity on the existing four queries with the 18-entry catalogue; E2 same AF298787 supported subset with original primary merge plus the separately proposed round-robin/component-attributed diagnostic. Exact new reference, original exclusions, per-query eligibility, model/weight identity, five-hit budget, seed and resource limits are frozen. This does **not** request fusion retrieval, ESM-C, Renilla or tier 2/3. All changed scored comparisons remain **HOLD pending actual external ChatGPT clearance**; no approval is manufactured. Reproduce preparation with `.venv/bin/python scripts/audit_independent_evidence.py` after restoring the hash-identified private provider originals. This expansion does not solve the independent-query sample-size shortfall.

## Reproduce and maintain

E1 uses the existing extension runner and separately admitted manifests:

```bash
.venv/bin/python scripts/run_bounded_stage.py --output results/cleared-rerun/E1.new.resources.json -- \
  .venv/bin/python scripts/run_evaluation_extension.py \
  --protocol config/execution_2026_11_11/E1.json \
  --inputs data/manifests/execution_2026_11_11/E1.new.inputs.json --output results/cleared-rerun/E1.new
```

Use `E1.control.inputs.json` separately for the control. E3/E2 use the existing parser/scorers through the narrow cleared-scope launcher; natural protein references are read directly from their exact frozen UniProt JSON, avoiding fabricated parent-DNA lengths or extra FASTA copies:

```bash
HF_HOME="$PWD/results/week2/model-cache" HF_HUB_OFFLINE=1 TRANSFORMERS_OFFLINE=1 \
  .venv/bin/python scripts/run_bounded_stage.py --output results/cleared-rerun/E3.resources.json -- \
  .venv-models/bin/python scripts/execute_cleared_retrieval.py --stage E3 --output results/cleared-rerun/E3
# Repeat with --stage E2 and fresh E2 paths for an independent encoder execution.
# Or reuse the original E3 scores, with checksum/model/source/panel validation and no encoder loading:
.venv/bin/python scripts/execute_cleared_retrieval.py --stage E2 \
  --reuse-e3 results/cleared-2026-11-11/E3/run01 --output results/cleared-rerun/E2.reuse
```

Restoring exact ignored caches is required; changed bytes block instead of being substituted. Each retrieval run verifies the original review freeze, source checksums, all 17 reference originals, protein fingerprints and input budgets before scoring. The executor stores the exact protocol/panel before results and refuses incomplete score support. Full diagnostics/test transcripts stay ignored; Git retains the authorization, minimal manifests, source-free stage evidence and this single status summary. Redundant full diagnostics were archived to ignored storage and remain exactly recoverable from reviewed Git commit `3fb4262`; [archive hashes](results/analysis-2026-11-12/archived-artifacts.json) document each move. The launcher retains the complete original authority digest, then verifies the separately pinned [stage dependency inventory](../config/execution_2026_11_11/retrieval-dependencies.v1.json). E2 requires its one cleared parent; E3 requires its three parents. Both require the frozen 17-reference catalogue and runtime/configuration evidence. Inference requires the pinned ESM-2 weights; E2 score replay requires independently anchored cached evidence instead. Archived reports, held ESM-C weights and unrelated iGEM originals are not live stage requirements.

Targeted cleanup added shared score reuse instead of another model implementation, removed duplicated per-input diagnostic details from committed stage summaries, and formatted the small launcher. Cache tampering is regression-tested without expensive models. Earlier clean-checkout checks: **249 run, 244 passed, five optional skips**; dependency-free **249 run, 190 passed, 59 skips**. Historical release verification passed. [Reproducibility evidence](results/analysis-2026-11-12/reproducibility.v1.json) records commands, pinned parsing dependencies, unchanged E2 conditions, byte-identical E4/derived archetypes and bounded replay resources. The fresh clone used a separate parsing environment and restored exact private caches; it is not a fresh model install or new encoder benchmark. Full logs remain ignored in the clean checkout's `results/` directory.

A final bounded biological E3 rerun after cleanup reproduced identical raw
alignment and ESM-2 score matrices. This is reproducibility evidence, not new
independent observations. Canonical score-matrix hashes and repeat resources
are retained in the E3 evidence.


## Scientific report and workflow handoff

The pilot supports **faithful annotation preservation and additional represented components**, but does not yet support a general improvement in family retrieval. Source decomposition retained two eligible documented proteins instead of one in E2; it did not improve the frozen five-hit merged endpoint. Attribution and rank/budget sensitivities identify how the endpoint obscures that distinction. Equal E3 family-hit rates and higher observed ESM-2 work time describe this small run only. None of these outcomes establishes biological function.

Failure analysis separates source limitations (missing iGEM translations, four held GenBank records), unsupported labels (AmpR, independent Renilla positive), representation limits (protein-only model, literal-versus-eligible longest), evaluation limits (off-family parent credit, component score calibration, crowded K=5, uneven catalogue), dependence (copies/shared lineage) and measurement limits (sequential loading, uncontrolled OS caches, absolute sampled RSS). All scored E3 queries have well-annotated multi-CDS parents: source CDS counts are 2, 2 and 3, with one, one and two supported query units. Single-CDS, backbone and genuine-fusion retrieval subgroup effects are unmeasured; strict exclusion tiers remain held. Family types, method differences and cost observations can be listed descriptively, with no powered subgroup comparison.

The [original workbook](../planning/4-week-workflow.xlsx), whose unchanged SHA is captured in [its extraction](week2/WORKFLOW_ACCEPTANCE.json), defines these task criteria. The old frozen Week 2 ledger is historical; this table is current.

| Task | Current assessment against the workbook definition of done |
| --- | --- |
| T15 baseline runner | Complete for cleared protein-compatible pilot: machine-readable longest eligible CDS outputs, pinned ESM-2 and configuration. Whole-plasmid DNA baseline is incompatible; original intact-record goal remains partial. |
| T16 decomposed runner | Complete for pilot: same model, traceable source CDSs, frozen parent merge and outputs. Unsupported components retained in coverage. |
| T17 scoring | Complete for pilot: overall endpoints, coverage, family labels, leakage exclusions and versioned derived archetype report. Broad archetype validation remains partial. |
| T18 paired evaluation | Cleared AF298787 longest/all-CDS comparison and matched E3 methods complete; reproducible run records retained. General decomposition hypothesis remains unresolved. |
| T19 failure/strata analysis | Complete descriptively above and in E2 attribution; small/absent strata explicitly unestimated. |
| T20 bounded improvement | Stable round-robin prototype implemented with focused budget/origin regression; exploratory mechanism benefit observed. Primary prospective evaluation held. |
| T21 regression comparison | Cached before/after attribution saved; hardened original E2 replay unchanged. Confirmatory comparison on a newly frozen panel held pending actual clearance. |
| T22 stable CLI/config | Existing annotation, cleared retrieval and analysis CLIs retained; README contains working one-command demos and reproducible outputs. Packaging is scoped to the completed pilot, without waiving T21's scientific hold. |
| T23 README/demo | Current README replaces stale readiness claims. Clean checkout annotation demo and source-free E4 demo pass. Model inference reproduction requires exact private caches and recorded model environment. |
| T24 clean rerun/environment | Separate parsing environment and fresh clone passed demos, optimized replay and full tests; dependencies/config/model identities recorded. Fresh model installation portability remains untested. |
| T25 final benchmark | Final cleared-pilot metrics and repeat evidence saved; expanded final benchmark waits for T21/panel clearance. No new independent observations claimed. |
| T26 limitations/recommendations | Pilot failure taxonomy and scientific interpretation complete. Next priorities: independently annotated multi-family parents, component-attributed secondary metric, pre-score calibration/aggregation criteria; no post hoc relabeling. |
| T27 demo/report | Scoped research report, evidence links and clean demos ready in this file/README; complete-project acceptance remains partial because coverage and paired population evidence are insufficient. |
| T28 cleanup/handoff | Duplicate full diagnostics archived with Git hashes; source-free results, exact holds and reproduction commands retained. PR left open. Next action is one exact-panel external scientific review, not routine researcher reapproval. |

MVP-01–03 have source-bound parser/unit/preservation evidence. MVP-04 remains partial across scored architectures. MVP-05 has the scientifically cleared protein-compatible proxy, not a raw-DNA/intact-protein-model claim. MVP-06–07 have scoped pilot reporting/tests/CLI/demo and explicit reproducibility limits. No unreported parser/runtime failure remains in the executed pilot; inadequate independence, source truth and unreviewed expansion are the material research gaps.

Handoff: retain the private cache originals and exact model snapshot; review the source-free results and the single panel-v2 request; do not rerun held panels, change primary criteria or merge without authorization. If v2 is cleared, run its exact bounds and preserve a new execution freeze. If an independent-query population study remains infeasible, report the central hypothesis as unresolved rather than adding sibling vectors or treating more references as more samples.


Review follow-up: both P2 findings at `b86b44b` are fixed without changing scientific scope. A minimal workspace with only 50 dependency/evidence files, no model weights and no iGEM originals reproduced the exact original E2 condition outputs under `python -O`, with zero model calls. Adding only E3's additional cleared sources/configuration and ESM-2 weight permitted E3 inference preflight; held ESM-C/iGEM assets remained absent. Focused regressions reject missing/changed required sources, references and configuration, a modified dependency inventory, and missing inference weights.

The alignment sensitivity now verifies the independently anchored raw matrix and each reference's original bytes/sequence fingerprint, recomputes reference-self scores using the frozen BLOSUM62 matrix/gaps/epsilon, rejects any mismatched cached normalized values, and consumes only the recomputation. The normalized-only mutation regression is rejected; regenerated E2 attribution JSON is unchanged. No model inference or primary experimental metric was rerun or altered. Follow-up tests: **251 run, 246 passed, five optional skips**; dependency-free **251 run, 191 passed, 60 skips**. Focused authorization/cache/stage/normalization regressions: **five passed**. Full transcripts remain ignored under `results/review-fixes-2026-11-12/`.
