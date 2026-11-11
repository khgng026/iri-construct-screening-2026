# Cleared experimental results

Authority: [ChatGPT delegated decision](https://github.com/khgng026/iri-construct-screening-2026/pull/11#issuecomment-6800553713), bound to reviewed commit `32e0928` and evidence freeze `c90d063f…47c2e`. Original proposal, sources and frozen criteria are preserved. Raw sources, qualifiers, proteins and detailed diagnostics remain ignored/private.

| Experiment | Actual result |
| --- | --- |
| E1, seven new GenBank sources | Completed: 67/67 observable features and qualifiers preserved; 15/15 source translations validated; all 74 original feature coordinates/strands and parent links agree. No parser errors or invented features. |
| E1, historical pBR322 control repeat | Completed separately: 48/48 observable features/qualifiers and 3/3 translations; all 50 original feature coordinates/strands agree. Not an additional novel source. |
| E3, tier-1 alignment/ESM-2 descriptive partial-method pilot | Cleared; execution next. |
| E2, AF298787 supported URA3/EGFP3 subset | Cleared; execution follows E3. |
| E4 / ESM-C / stricter tiers / Renilla | HOLD as specified in the reviewer decision. |

E1 measures agreement with deposited annotations, not biological correctness. All seven new records are in the well-annotated stratum; uncertainty/justified-abstention denominators are zero and their rates remain null. Duplicated vector proteins do not become independent family observations. [E1 results and measured costs](results/cleared-2026-11-11/E1.json) include source-specific denominators and full local report fingerprints.

Reproduce using new output paths:

```bash
.venv/bin/python scripts/run_bounded_stage.py --output results/cleared-rerun/E1.new.resources.json -- \
  .venv/bin/python scripts/run_evaluation_extension.py \
  --protocol config/execution_2026_11_11/E1.json \
  --inputs data/manifests/execution_2026_11_11/E1.new.inputs.json --output results/cleared-rerun/E1.new
```

Use `E1.control.inputs.json` separately for the historical control. Source hashes are enforced by the admitted manifest on every invocation. The bounded wrapper monitors RSS, kills on a sampled 6-GiB/1,800-second violation, and checks the OS high-water RSS before accepting success. No new raw-source Git copy is authorized.
