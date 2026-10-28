# T04 | Git repository and reproducible environment

**Week:** 1  |  **Budget:** 2 hours  |  **Status:** Implemented and hosted privately; GitHub Actions passing  |  **Prerequisite:** T01  |  **Calendar source:** `Weekly To-Dos!A8:K8`.

## What has been set up

- Local Git repository named `iri-construct-screening-2026`, branch `main`, with tracked Week 1 docs, metadata fixtures, tests and package skeleton. A private GitHub repository is live at https://github.com/khgng026/iri-construct-screening-2026.
- `pyproject.toml`: Python 3.10+ metadata; **no runtime dependencies** for the Week 1 checker. Installation is optional for the first-run test.
- `config/default.json`: human-readable, versioned defaults for manifest path and log level. Paths are repository-relative to make the smoke check independent of current working directory.
- `src/iri_constructs/`: metadata manifest validation and a small logging/configuration demo; **not** a sequence parser or hazard detector.
- `scripts/smoke_test.py`: built-in-standard-library smoke command that reads configuration, validates a fictional manifest and schema examples, and exits nonzero on failure.
- `tests/`: Python `unittest` checks. `.github/workflows/ci.yml` runs tests on GitHub push/pull request; the verified run passed (11 unit tests).
- `.gitignore`: excludes downloaded biological sequences, raw/processed datasets, secrets, Python caches, local virtual environments and generated results.

## Fresh-checkout test (Linux/macOS)

From the extracted or cloned repository root:

```bash
python3 -m venv .venv
.venv/bin/python scripts/smoke_test.py
.venv/bin/python -m unittest discover -s tests -v
```

On Windows, use `.venv\\Scripts\\python.exe` instead of `.venv/bin/python`. Neither test downloads data/models nor requires `pip install`; only Python standard-library modules are used. An optional editable install for development is `python -m pip install -e .` when package tooling is available.

**Expected success:** smoke check prints `Week 1 smoke check passed` and tests report `OK`. An intentional invalid example/manifest should generate an explanatory failure and nonzero exit status. Verification is available in GitHub Actions and the local test command.

## Configuration and logging conventions

- Configuration is JSON and checked for known settings. `config/default.json` is committed; secrets and machine-specific overrides live outside the repo.
- Use Python's `logging` module with `LEVEL logger event=...` messages. Default `INFO`; error detail at `ERROR`, diagnostic detail at `DEBUG` only when explicitly configured.
- Log **task name, file/record ID, count, outcome, model/version (when later used), and elapsed time**. Never log sequence strings, customer identifiers, secrets, private paths, or unapproved biological content.
- Commands return exit status `0` on success and nonzero on invalid configuration/metadata; the current checker does not assert the biology behind example labels.

## Git controls and handoff

The local repository has an initial commit and subsequent alignment commit. Confirm with `git log --oneline -2` and `git status --short`. The remote is a private GitHub repository at https://github.com/khgng026/iri-construct-screening-2026; the `main` branch contains the Week 1 scaffold and planning workbook. Changes can be pushed via Git or the connected GitHub integration.

**T04 done when:** The archived Git repository opens with history, the above commands pass in a fresh Python environment, and configuration/logging rules are documented. This task unblocks T05, T10, T13 and T15; it does not imply any model inference or screened data exists.