# IRI research session handoff

Use this brief only after the user explicitly starts the overnight session. Read the packet files rather than treating this brief as a replacement for them.

## Goal

Answer as much as feasible of all four questions: model dependence, increasing difficulty, comparison with established local alignment, and averaging versus context. Seek verified evidence and a defensible next research question. A well-documented negative or partial result is preferable to an unsupported headline.

## Read first

- `/Users/Rudra2.0/Desktop/IRI/overnight-plan-2026-10-31/PROTOCOL.md`
- `/Users/Rudra2.0/Desktop/IRI/overnight-plan-2026-10-31/OPERATIONS.md`
- `/Users/Rudra2.0/Desktop/IRI/overnight-plan-2026-10-31/EXTENSION_QUEUE.md`
- `/Users/Rudra2.0/Desktop/IRI/overnight-plan-2026-10-31/challenge_manifest.json`
- `/Users/Rudra2.0/Desktop/IRI/overnight-plan-2026-10-31/preparation_checks.json`
- Applicable instructions and current reports in `/Users/Rudra2.0/Desktop/IRI/iri-construct-screening-2026`.
- The final machine-readiness handoff under `/Users/Rudra2.0/Desktop/IRI/machine-readiness-2026-10-31`.

## Fill from user instruction and hardware handoff

| Item | Value at preparation |
| --- | --- |
| Start authorization | Pending; do not infer from preparation |
| Session end time/window | User requests 10–12 hours; record actual launch and hard stop at launch + 12 hours; use America/New_York |
| Readiness probe | Pending |
| Approved checkpoint downloads/storage budget | Pending hardware handoff and user decisions |
| Paid services | Not authorized |
| Heavy-process concurrency | One unless hardware handoff explicitly validates otherwise |
| CPU/MPS/device and numerical tolerances | Use verified machine handoff and trial |
| Memory, disk, power, and timeout guard | Use the hardware chat's tested settings |

Before execution, resolve these fields and record them with the run. Do independent lightweight preparation while waiting for a required answer, but do not bypass a failed launch/resource gate.

## Work instructions

You may delegate to a lead plus up to three subagents as described in OPERATIONS.md. The lead owns integration, job scheduling, shared configuration, and commits. Agents own distinct new files. Keep inference sequential by default. Preserve original inputs and results; never reset another chat's edits.

Reproduce the existing baseline, verify the prepared suite, implement and pass the representative end-to-end trial, and freeze actual run inputs/settings before scoring. Prioritize alignment, 8M controls, and 650M on the existing pair benchmark. Then execute the prepared triples/quadruples and admitted alternative models with common-subset comparisons. Use remaining resources for catalogue sensitivity, additional controls, and a justified labelled reference-holdout extension. Preserve exploratory investigations separately.

Implement incremental outputs, a resumable job ledger, bounded retries, and cache integrity checks. Carry on with independent feasible work when an optional model is blocked. Stop resource-heavy work when the machine guard requires it. Audit outcomes rather than merely checking process exits. Do not add jobs simply to fill elapsed time.

Deliver a morning report answering each question with the actual evidence or unresolved reason, comparative plots/tables, failure cases, reproducibility checks, saved paths, rerun commands, and a concise autonomy/delegation audit. No pushing, merging, external messages, paid compute, sequence generation, or persistent system changes are authorized.
