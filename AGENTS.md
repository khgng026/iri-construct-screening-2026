# Active project instructions

The researcher authorized the toolkit-first reset on November 16, 2026 (America/New_York), after reviewing the preceding workflow critique and action items. Follow the current conversation first, then these active project documents:

1. [Scope and acceptance](docs/project/SCOPE.md)
2. [Autonomy and data rules](docs/project/AUTONOMY.md)
3. [Toolkit contract](docs/project/TOOLKIT_CONTRACT.md)
4. [Implementation plan](docs/project/PLAN.md) and [current status](docs/project/STATUS.md)

These supersede older T01–T28 planning, the workbook's release dependencies, T02/T07 per-run approval mechanisms, dated mission records and benchmark completion gates for new toolkit work. Those documents and their frozen results remain historical evidence. A historical HOLD is a fact about that historical run; it is not a prohibition on separately versioned work allowed by the active policy. Do not silently modify a frozen experiment to take advantage of this change.

Deliver the preprocessing toolkit before extensions. Keep source-guided CDS extraction, explicit unknowns, validated translations and full provenance. No domain splitting, gene prediction, model training, new model families, public dataset release or production service is part of the current release.

Work autonomously within AUTONOMY.md. Do not request routine panel/source/AI-review signoffs. Checksums establish input identity, not permission. Label self-review honestly; external review is not required for routine in-scope work. Continue independent core tasks when an optional model or record is unavailable.

Do not edit historical approval messages, frozen manifests, original scientific reports or their checksum-bound files to make new code pass. Put new policy, runs and results in new versions. Historical scripts may remain available for reproduction; new toolkit code must not depend on their PR comment IDs, dates, Git objects, research labels or private caches.

Run relevant tests after implementation. Core parsing/export/install tests require no model weights or provider network access. Keep model tests optional and disclose skips. Preserve useful existing regressions; avoid removing tests merely to reduce counts. Verify the built wheel in a clean environment outside the checkout before calling the toolkit release-ready.

Use reviewable commits and pull requests. Merge, public publication, spending and destructive external actions require explicit authorization. Update STATUS.md at meaningful checkpoints; do not add another competing status ledger.
