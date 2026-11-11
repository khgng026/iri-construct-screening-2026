# T02 v1.2 | Delegated AI-assisted source and experiment admission

**Policy request:** 2026-11-10, researcher asks to revise the original review rules so ChatGPT can approve running experimental analysis on data **after asking the researcher for explicit permission**. This document is the proposed versioned mechanism. **It is not a recorded approval of any scored run.** Effective only after this change is merged/accepted and the separate scoped researcher run authorization is documented.

## Purpose and binding prior rules

Change **who may perform documented source/label reviews** under a researcher-authorized project execution scope, not **what evidence counts**. Retain T01 scope, T02 v1.0 + v1.1 permitted public-benign data/rights/safety/provenance, T03 schema, and T07 v1.0/v1.1 frozen scientific endpoints and leakage guards. Preserve historical records and signatures; don't retroactively treat legacy HOLDs as admitted.

## Approved candidate pool versus executable panel

- **Candidate-selection authority already granted:** eleven versioned NCBI GenBank records (merged PR #7), ten iGEM records with publicly retrievable sequences and source features (merged PR #8). The nine earlier Addgene HOLDs were retired by researcher direction, retained as historical records, not pursued. pBR322 J01749.1 was previously E1-admitted under merged PR #6; its previous E1 result is not an E2/E3 claim.
- The above candidate approvals are **not** data/scoring admissions. A candidate may be selected but excluded after rights, safety, annotation, or independence review. Do not force the entire 21-record candidate pool into the earlier approximately 8–12-record frozen pilot. Prospective changes to scored panel, archetype coverage, hypothesis or budget require a dated versioned plan before scoring.
- Final "run approval" must identify a named versioned dataset/subset, eligible experiment(s), fixed model and reference rules, output/storage policy and material cost/resource limits. It must come from the researcher as an explicit affirmative reply after a clear request. Silence, selection approvals, merging a PR or this amendment alone do **not** confer run authorization.

## Delegation and approval roles

1. **Codex (evidence preparer):** retrieve exact public source under permitted access; log original bytes' hashes and source versions; check annotated CDS translation, source completeness, reuse restrictions and benign research context; propose per-protein family labels and independent references with external evidence; generate machine-readable candidate/review manifests; freeze a proposed protocol *before* model inference. Codex may run software unit tests and non-scored preparatory checks to prepare this evidence, but may not impersonate an independent reviewer or mark missing evidence as approved.
2. **ChatGPT (delegated evidence reviewer):** after a researcher-authorized scope exists, inspect the evidence packet, the applicable frozen protocol and the rights/safety facts. Explicitly approve, reject or HOLD each proposed eligible source, family and reference relation; record reasons, evidence links, uncertainty, exclusions, review date, and a specific panel/hash. ChatGPT can issue documented execution clearance for the eligible subset when **every applicable gate passes**, without another researcher signature on each record. It must identify its decision as **AI-assisted delegated review**, not an independent human or specialist sign-off.
3. **Researcher (ultimate authority):** authorizes the specified experimental scope and can withdraw delegation. New hazardous/dual-use concerns, nonpublic data, unlicensed hosting, public release, material hypothesis/model-budget change, or contested biological judgment must be escalated, not decided by Codex on implicit consent.
4. **Codex (executor):** once a recorded researcher authorization and ChatGPT clearance both exist for the exact frozen panel, run only the allowed experiments, retain raw results, reproduce hashes, and report abstentions, failures and limitations. A result does not retroactively approve a label.

## Nonwaivable per-record gates

- Confirm public provider and immutable accession/part ID; exact downloaded byte checksum/retrieval date; type of record (full plasmid versus cassette/part); source annotations and caveats; ordinary benign suitability; license and **separately** local computation and private cloud/GitHub storage rights. Unclear copyright/attribution or unknown hazards => HOLD and no source upload/scoring.
- Validate source-provided CDS locations/translation where present. Never invent gene/protein components or turn an unannotated fusion-ready vector into a proven fusion. Unsupported noncoding content stays metadata.
- For protein-family scoring, cite independent family evidence for every expected class; freeze positive and competing reference IDs and conservative independence/exact-copy/related-vector exclusions. Do not use the query construct's own components, clones or derivatives as independent references.
- For source-only annotation preservation E1, family/reference panel is not required; source/annotation/permission/benign gates still apply. E2/E3 require reviewed family panel and independent references. E4 needs actual measured E2/E3 data and compatible environment; no fixture scores can qualify.
- Record any scientific shortfall as blocked/descriptive-only/unknown and continue only with defensibly eligible subset. No model-snapshot substitution, post hoc family relabeling or threshold tuning after inspecting results.

## Approval record contract (human-readable and machine-readable)

Every dataset/scored-run clearance requires: (a) explicit researcher authorization quote/date and scoped named experiment/panel, (b) ChatGPT review ID/date/evidence pointers, (c) frozen manifest, input/checksum, labels/exclusions/protocol and config/model hashes, (d) distinct stage statuses `candidate_approved`, `source_admitted`, `cloud_copy_permitted`, `family_panel_approved`, `run_authorized`, (e) reasoned exclusion list and resource/storage bounds. Each field is independent and must not be inferred from a previous one.

An editable audit cannot replace a real approval message. If ChatGPT is unavailable to verify an evidence packet, Codex stops at review-ready status. A prospective change to the frozen panel requires a new review record before scoring.

**Current application (2026-11-10):** 11 GenBank + 10 iGEM identities selected; the earlier pBR322 E1 is recorded; new source files and iGEM annotations have acquisition evidence but no complete current T02 admission and no reviewed independent retrieval panel. **No new experimental execution clearance is given here.**
