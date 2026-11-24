# Standing autonomy and data policy

Updated November 22, 2026. This policy governs work within the active [research scope](SCOPE.md). No additional routine experiment, source-record or separate ChatGPT signoff is required. User instructions take precedence.

## Proceed autonomously

- Inspect, implement, evaluate and document the active research workstreams and existing toolkit; fix prerequisites; create commits and reviewable PRs.
- Process suitable user-supplied and public benign research records; discover a modest number of examples when required for behavior coverage. Review source-specific annotations and notices, record evidence and proceed when the source guidance and record facts support use. Source-level review can be reused unless terms or record-specific facts change.
- Install suitable public dependencies in isolated environments. Select pretrained models appropriate to the active research scope, verify versions and weight hashes, and record use terms. No model training or fine-tuning is required by the scope.
- Version shareable research results in Git through routine reviewable PRs. Include a report, the protocol or configuration, and the commands needed to reproduce it. Include compact CSV or JSON tables when they help reviewers verify counts, comparisons, or uncertainty. Record the input provenance, denominators, and model version where applicable.
- Use the ignored workspace for source records, intermediate outputs, embeddings, and model caches, subject to source and model terms. Keep active held-out answer keys in evaluator-only storage while a benchmark is blind. Release underlying data or model artifacts only when redistribution is suitable and publication is authorized.
- Run regression, exploratory and descriptive evaluation within scope; choose useful examples and version new runs without per-panel approval. Invalid records or unsupported metrics are reported individually. A disputed family label blocks that metric, not extraction of otherwise valid source annotations.

## Resource guidance

Prefer available free compute. Choose CPU or GPU and model concurrency based on the workload and available hardware. If resource needs are unclear, time a representative batch and check memory and disk use before scaling up. Reuse cached artifacts where practical, verify downloads, and respect provider rate limits. For substantial runs, checkpoint progress and share periodic updates. Adjust batch size or concurrency if the run strains the environment. Run size and duration should reflect what the research question requires and what the environment can support. Ask before incurring costs or using paid compute.

## Automated evidence, not approval loops

Record provider/accession/version when available, retrieval date, source hash, relevant terms/notice evidence and disposition. For ordinary public molecular records with no conflicting notice and suitable research context, document the source-level rationale and continue. Public visibility alone is not proof of reuse rights, and an absent conflict is not a legal guarantee. Escalate an actual material uncertainty; do not demand new affirmative permission for every routine record.

For local parser inputs, a source URL or license approval field is not a technical parsing prerequisite. The toolkit records supplied provenance and missing fields. The agent's research data acquisition policy still applies to examples it chooses to obtain.

Compute current hashes automatically. A new provider snapshot becomes a new identified input; validate it and record the change. An exact historical rerun requires its original bytes, but their absence does not block new development or a new labeled run. Keep completed results immutable and distinguish an exploratory adjustment from an untouched confirmatory evaluation. Family/lineage exclusions matter for performance claims, not annotation preservation. Report self-review as self-review, without implying human or independent specialist review.

## Ask only for material exceptions

Spending or using paid compute; new restricted/nonpublic data or account agreements; substantive suitability or rights conflicts; public dataset/source release; production deployment; destructive external changes; a new research direction outside scope; a contested scientific judgment central to a reported conclusion; merging or publishing when not already explicitly authorized. Isolate the affected task and continue independent core work. Do not ask again for authority already granted.

## Legacy transition

Old T02/T07 approval mechanisms, approval-comment identifiers, candidate whitelists and blanket no-download/cache restrictions do not govern new in-scope work. This does not retroactively change a historical clearance or its results. Existing historical scripts enforce their original contract and may fail outside it; build the general toolkit separately rather than weakening those checks. Migrate useful engineering validation, not historical approval gates.
