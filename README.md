# IRI Fall 2026 | Complex-construct preprocessing (Week 1)

**Status (2026-10-28):** T01–T03 drafts aligned to the actual `4-week workflow.xlsx` calendar. T04 reproducible metadata-only environment and private GitHub repository implemented. **No model or biological data processing has been run.**

## Research question

Does preprocessing **public, benign, annotated DNA constructs** into traceable, model-compatible units improve a fixed pretrained model's characterization of their components compared with intact/minimally processed records under a fair paired evaluation?

The IRI proposal frames this as an approximately **80-hour integration/evaluation MVP**, not model training or operational biosecurity screening. The calendar allocates 20 hours/week including 2 hours of parallel technical skills. This repository covers **only T01–T04** of Week 1; T05–T08 and all later tasks remain outside this initial handoff.

## Week 1 official tasks (from `Weekly To-Dos`)

| ID | Calendar task | Artifact | Status |
| --- | --- | --- | --- |
| T01 | Write MVP question, scope, and acceptance criteria | [T01 scope and acceptance](docs/week1/T01_scope_and_acceptance.md) | Draft awaiting agreement |
| T02 | Set public-data and safety rules | [T02 governance checklist](docs/week1/T02_data_governance.md) | Draft awaiting approval |
| T03 | Define construct archetypes and normalized output schema | [T03 schema](docs/week1/T03_archetypes_and_schema.md), [three fictional examples](data/examples/) | Draft awaiting pilot review |
| T04 | Create repository and reproducible environment | [T04 setup](docs/week1/T04_repository_setup.md), config, smoke test and unit tests | Done; hosted privately with passing CI |

[Four-week workflow with progress notes](planning/4-week-workflow.xlsx)

## Reproducible quick start (Python 3.10+, no packages to download)

```bash
python3 -m venv .venv
.venv/bin/python scripts/smoke_test.py
.venv/bin/python -m unittest discover -s tests -v
```

The smoke check reads `config/default.json`, validates a **fictional metadata-only** example manifest and the three fictional normalized schema examples, and writes safe status messages via the standard-library logger. It does not parse sequences or run a model. On Windows, use `.venv\\Scripts\\python.exe` rather than `.venv/bin/python`.

## Project layout

- `docs/week1/` — exact-calendar T01–T04 drafts and task status.
- `planning/4-week-workflow.xlsx` — four-week planning workbook, including dated notes for T01–T28 and current task statuses.
- `config/default.json` — checked-in defaults; no secrets or per-machine credentials.
- `data/manifests/` — fictional example inventory in CSV format, with provenance/review gate columns.
- `data/examples/` — three schematic examples of the normalized output schema, with no biological sequences.
- `src/iri_constructs/` — standard-library metadata and configuration validators.
- `scripts/smoke_test.py` and `tests/` — local and CI checks.
- `.github/workflows/ci.yml` — Python tests on push and pull request; latest verified run passed.

## Boundaries and limitations

Only approved public benign records may be used for the future benchmark. Do not commit sequences, restricted or customer data, or secrets. Protein- and DNA-specific interfaces must be confirmed in later T05; labels and coordinates may be unknown. T04 is hosted in a private GitHub repository with a passing CI run. Early schema examples are fictional and establish structure, not biological correctness or model performance.