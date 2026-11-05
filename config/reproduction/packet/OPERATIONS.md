# IRI overnight execution and handoff

## Boundaries

Prepare and run only the locally authorized protein-retrieval research. The user has not yet instructed this preparation chat to launch the overnight session. Do not launch from this packet alone. Do not push, merge, publish, contact collaborators, buy compute, use paid APIs, generate biological sequences, or change persistent system settings. Preserve prior experiments and the original environment. No source PDF or raw sequence enters Git.

Hardware setup is owned by **Prepare Mac for IRI benchmarking**. Use its final handoff, environment paths, probes, access/download approvals, and tested resource guard. The preparation chat's earlier observations are not current readiness proof. Recheck power, free disk, memory pressure, and competing processes at launch. Never stop unrelated jobs or close the user's apps autonomously.

## Ownership and concurrency

Use a lead and at most three concurrent subagents where available. Hardware resources are shared. Assign file ownership explicitly; the shared checkout does not isolate edits.

| Owner | Responsibilities | Writes |
| --- | --- | --- |
| Lead | Freeze configuration, integrate adapters, schedule inference, manage ledger, write final report | Runner, common interfaces, final protocol/report |
| Model agent | Verify official interfaces/licences, implement isolated model adapters, record dependencies | Dedicated new adapter files and adapter-specific tests |
| Benchmark agent | Validate prepared manifest, implement suite construction, catalogue subsets, metrics | Dedicated new benchmark/analysis files and tests |
| Audit agent | Review method fairness, independently recalculate metrics, inspect anomalies | Separate audit artifacts; request corrections from lead |

The lead owns shared configuration and orchestration files. Subagents must not switch branches, run Git mutations, edit another owner's files, modify environments, or launch heavy inference without the lead's resource slot. Reassign the model agent to alignment implementation or the benchmark agent to context controls as tasks finish. Use agents for useful work, not duplicated opinions.

Use actual available collaboration tools for substantive independent work, not simulated role labels in one response. Delegate adapter development, benchmark implementation and independent audit when those tasks can proceed without file or compute conflicts. Record real task IDs and outputs. If this runtime lacks subagent tools, state that limitation and execute a sequential fallback; a prompt cannot create missing tool capabilities. Keep the lead focused on integration, experimental interpretation and next-job selection rather than duplicating delegated work.

Create a new branch from the verified current research branch only after checking worktree state. The lead alone commits reviewed changes. Work remains local. If simultaneous chats modify the repository, reconcile ownership before editing shared files; never reset others' changes.

## Evidence and adaptation loop

Maintain a compact shared `SESSION_STATE.md` in the new session's artifact folder, with goals, current branch/config, completed/active/blocked jobs, accepted findings linked to evidence, unresolved questions, resource status, and exact next actions. Keep a `LESSONS.md` for verified implementation lessons, failed approaches and reproducible recovery steps. Store scientific interpretations separately from confirmed software facts. These files support within-session learning and recovery; they do not change model weights or authorize writing persistent personal memory.

After each substantive completed work package:

1. Check its outputs against its acceptance criteria and have a separate agent inspect material new conclusions or critical code where available.
2. Record what was observed, which explanation remains uncertain, and whether it changes priorities for the remaining research questions.
3. Update the compact state and lessons with evidence paths. Pass the relevant findings and constraints to agents handling dependent tasks; require them to acknowledge consequential corrections.
4. Choose the next feasible experiment by its expected ability to distinguish explanations, remaining compute and verification cost. Preserve the common benchmark, and label designs prompted by outcomes as exploratory.
5. Improve adapters, caching, scheduling and checks when evidence reveals a defect or bottleneck. Verify the change before applying it to later jobs. Invalidate or rerun earlier results if the defect affects them; never mix results across incompatible implementations silently.

This is an iterative research and workflow-improvement loop, not recursive self-modification of Astra. Do not add meta-agents, rewrite prompts repeatedly, or spawn redundant reviewers simply to look busy. Preserve the experimental goals, resource boundaries, raw evidence and final audit. Reuse idle agents on the next independent task, and keep computing only when the work is scientifically or operationally useful.

At context compaction, interruption or agent replacement, restore from compact state, the ledger and required evidence rather than rereading full tool transcripts. Load only context relevant to the current task. The session should accumulate useful findings without accumulating unbounded duplicated instructions.

## Execution order

The user requests approximately 10–12 hours of useful work. At actual launch, record start time and a hard end time 12 hours later, using America/New_York for display and a monotonic deadline for execution. Aim to use the first 9–10 hours for preparation and scientific work and the final 1–2 hours for verification and reporting. Do not enforce the hardware handoff's earlier suggested eight-hour session cap: that scheduling suggestion is superseded by this user request. Retain all of its power, pressure, disk and swap limits. The existing guard accepts finite timeouts up to 12 hours; use shorter measured per-shard bounds as well. Prevent-sleep duration must match the actual finite session/runner rather than the earlier eight-hour example.

Preflight, ordinary reversible software repairs, representative tests, and transition to the main run form one unattended workflow once launch is authorized. Do not ask for another confirmation between stages. Repairs may adjust implementation, isolated dependencies, batch size, scheduling and validated device fallback. Record changed conditions as distinct attempts. Do not silently change examples, model identity, pooling, precision or labels to make tests pass. Refreeze code/configuration after material repairs and rerun affected checks. A failed memory/power gate is not permission to weaken the guard, close apps or restart the Mac.

If the initial resource gate is blocked, lightweight code preparation and review may proceed. Recheck resource state at bounded intervals, keep checkpoints, and run heavy work only after stable normal pressure and power are present. Do not let waiting consume the entire scientific window without reporting the resource blocker. Optional blocked models should not stall other feasible work. Account/tool interruptions can end oversight; preserve resumable local outputs rather than claim guaranteed unattended continuation.

1. Read current project instructions, this protocol, prepared checks, and the final hardware handoff. Record access/download approvals and compute a fresh resource snapshot.
2. Verify source hashes and original outputs. Materialize prepared queries into ignored storage only after checks pass. Copy/version the approved protocol/config in the experiment branch; hash it before scoring. This packet is a preparation artifact, not a substitute for a run freeze.
3. Integrate common adapter/metric interfaces and complete the representative trial below.
4. Reproduce existing 8M results; run local alignment and 8M controls. Complete the core for 650M once admitted.
5. Admit alternative models one at a time, using the same frozen suite and reporting common subsets. Prioritize getting each model through the shared core before extending one to a much larger catalogue.
6. Run nested catalogue tests and additional controls. Prepare labelled holdout or larger-catalogue extensions only if core comparisons and verification are on track.
7. Complete fresh-process repeats and independent audit. Reserve enough of the user-specified session window for final validation and reporting; a suggested reserve is the final 20%, not a reason to discard running evidence.

If the shared core completes early, work through EXTENSION_QUEUE.md rather than ending automatically. Its first priority is automatic segmentation with fair output budgets, followed by credible catalogue/reference tests and representation sensitivity. Stop taking new expensive jobs when the reporting reserve begins. Remain within the defined scientific and resource scope.

Within each model, finish clean and existing pair cases, controls, triples, then quadruples. Save outputs incrementally so interruption leaves interpretable work. Predeclare a deterministic minimum subset in addition to the full suite. If only that subset completes, compare every model on the same subset and label full/partial results separately.

## Representative trial

Choose before inference the shortest eligible pair, longest eligible pair, most balanced pair, and most unequal pair; stable ID order breaks ties. Use both orders, all relevant clean targets, and the four averaging/context arms. Add one eligible triple and quadruple chosen by fixed hash for boundary/metric checks.

Require:

- A single documented command writes versioned run metadata, complete rankings, metrics, and a job state.
- Original source/sequence hashes match and token-to-residue mapping is exact.
- Clean vectors and reference vectors agree as expected for matching configurations.
- Whole vector equals the weighted contextual-vector mean; isolated mixture is invariant to order within tolerance.
- Multi-component recall and all-component-hit denominators are independently checked.
- Alignment spans and scores pass the independent subset check.
- A deliberately interrupted owned job can resume without duplicate rows, skipped records, or stale-cache reuse.
- Runtime and resource observations support a feasible full run. Failed probes lead to smaller batches, CPU fallback where validated, or a recorded skip, not silent changes to model identity or examples.

This trial is pending. Do not claim launch readiness until it passes.

## Job ledger and recovery

### Unattended decisions

The user will be asleep. Do not request routine preferences, confirm reversible implementation choices, or stop after a progress update. Use the protocol and recorded launch authorization to resolve ordinary choices. If a path requires credentials, a new agreement, an unavailable tool, paid services, additional download allowance, publication, or a privileged/destructive action, mark that path blocked, record the exact missing requirement for the morning, and proceed with other independent authorized work. Do not repeatedly ask or wait indefinitely on an optional blocked path. This rule does not bypass platform approval prompts or authorize crossing scope boundaries.

An optional model failure should not halt the shared runner or erase completed jobs. Catch per-job errors and record them. If an essential resource gate fails, stop heavy work, checkpoint, and continue only feasible lightweight analysis/reporting. A resource failure is a reason to pause owned computation, not permission to alter the user's apps or system protections.

### Data durability

Create a versioned session directory with an immutable input manifest, model/configuration metadata, code revision, raw source references, checkpoints, outputs, stdout/stderr logs, and audit artifacts. Preserve each completed unit immediately rather than holding an entire model run in memory until the end. Keep prior runs and original data unchanged. Do not overwrite a successful attempt with a retry.

Use temporary files, flush and fsync where appropriate, then atomic rename on the same filesystem. Record checksums and completion only after data files are durable. For a multi-file unit, publish a final completion descriptor last; resume trusts only units with valid descriptors and checksums. Use an append-only event log or a transactional ledger, and reconcile it against saved files on restart. Disk exhaustion must become a recorded failure, not a partially written result marked complete.

Keep catalogue embeddings and derived data in ignored local storage with source/model/configuration hashes. A cache is replaceable acceleration, not the sole copy of irreplaceable experiment results or source records. Preserve full raw rankings and summary denominators so reports can be regenerated. Store enough hidden-state summaries to reproduce the mechanism analysis without retaining unnecessary massive tensors. Explicitly list any artifact omitted due to storage limits and its regeneration route.

At the end of each model or major suite, produce a checksum index and consistent snapshot/export of the ledger, state and verified outputs in a separate session recovery directory. Such a same-disk snapshot protects against accidental overwrites or partial edits, not drive loss. No cloud upload or external backup is authorized. Do not claim redundant copies on the same disk guarantee data survival.

Before full execution, interrupt the representative trial's owned process deliberately, restart it, and confirm that already completed records survive, incomplete records rerun, cache corruption is detected, and no record is duplicated or silently omitted. Verify that an optional-model exception permits the next job to run. Verify that a simulated denied-access job is recorded as blocked without requesting input. Test these behaviours before marking the session unattended-ready. The planning packet specifies these requirements; it does not itself implement or prove them.

Implement a ledger with immutable job IDs from model revision, code/config/input hashes, suite/catalogue IDs, arm, device, dtype, and seed. States are `planned`, `running`, `complete`, `failed`, `skipped`, or `blocked`. Store start/end times, output paths/checksums, retry count, supported/excluded query counts, exception summary, and resource exit reason. `complete` requires verified output integrity, not merely process exit zero.

Write outputs atomically and checkpoint after each query or bounded batch. Cache catalogue embeddings per model, representation, tokenizer, settings, code, and sequence hash; validate cache bytes and shape before reuse. No sharing vectors across different checkpoints. Keep completed job files intact and append retry evidence rather than overwriting it.

Allow at most two retries for a transient download/network error with bounded backoff. Do not retry deterministic invalid input, access rejection, licence restriction, or repeated out-of-memory errors unchanged. Timeouts should be estimated from the representative probe with recorded slack and an absolute bound from the hardware handoff. Skip a blocked optional model and continue other independent work. Record missing access for the user without trying to bypass it.

One heavy model process at a time is the default. Use the hardware guard for owned jobs; treat unavailable telemetry conservatively. Avoid storing all layers or attention matrices. If a resource gate fails, save state, stop the owned job as the guard specifies, and continue lightweight reporting/review rather than repeatedly loading weights.

## Decisions required at launch

- User start instruction and desired session window or stop time, in America/New_York. No artificial requirement to keep working once useful planned work is complete.
- Final machine handoff and passing representative probe.
- Which large downloads/checkpoints are approved and the aggregate storage limit. Preparation does not override approval requests already pending in the hardware chat.
- Default is local-only, no paid services. Any different spending or external-service scope must be explicitly authorized.
- Model admission table including exact revisions, supported devices/dtypes, environments, access status, and skip reasons.

If the hardware chat establishes bounded numeric limits, reuse them. Do not invent incompatible thresholds or run a second independent watchdog. Update this checklist with concrete values before launch.

## Morning report

The first paragraph states what completed and the strongest supported finding. Then include:

1. A four-question table with answer, evidence/run IDs, limitation, and unresolved items.
2. Comparable model/method results on identical catalogue/query intersections, separated by component count and fraction.
3. Alignment comparison and context/averaging-control interpretation, including counterexamples.
4. Runtime, memory, numerical repeatability, completed/failed/skipped job counts, and provenance.
5. One next justified scientific experiment and a plain-language explanation for Rudra.
6. An autonomy audit: delegation ownership, interventions needed, failures recovered, tests performed, and work left unfinished.

Link the saved report, plots, ledger, and rerun commands. Do not claim all questions answered merely because jobs ran. Do not infer a new unsolved problem or model superiority without comparisons that support it.
