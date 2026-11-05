# IRI Fall 2026 | Complex-construct preprocessing (Week 1)

**Status (2026-11-05):** T01-T04 are approved. T05 is in progress: reviewed local-input parsing, ESM-2 fusion/probe and modular ESM-2/ESM-C/ESM-3 feasibility runners are implemented. No reviewed-data pretrained inference or scientific model comparison has been executed. See [current status](docs/week1/STATUS.md).

## Research question

Does preprocessing **public, benign, annotated DNA constructs** into traceable, model-compatible units improve a fixed pretrained model's characterization of their components compared with intact/minimally processed records under a fair paired evaluation?

The IRI proposal frames this as an approximately **80-hour integration/evaluation MVP**, with 20 hours/week including 2 hours of parallel technical skills. This repository retains T01-T04 and adds bounded T05 feasibility modules and reusable skills exercises; later dataset and benchmark tasks remain incomplete.

## Week 1 official tasks (from `Weekly To-Dos`)

| ID | Calendar task | Artifact | Status |
| --- | --- | --- | --- |
| T01 | Write MVP question, scope, and acceptance criteria | [T01 scope and acceptance](docs/week1/T01_scope_and_acceptance.md) | Approved |
| T02 | Set public-data and safety rules | [T02 governance checklist](docs/week1/T02_data_governance.md) | Approved (policy v1.0, 2026-10-28); record-by-record reviews pending T06/T08 |
| T03 | Define construct archetypes and normalized output schema | [T03 schema](docs/week1/T03_archetypes_and_schema.md), [three fictional examples](data/examples/) | Approved baseline v0.1; real-record pilot checks pending |
| T04 | Create repository and reproducible environment | [T04 setup](docs/week1/T04_repository_setup.md), config, smoke test and unit tests | Done; hosted privately with passing CI |
| T05 | Run a pretrained model feasibility spike | [Experiments and commands](experiments/README.md) | In progress; reviewed inputs and actual model-run evidence pending |

[Four-week workflow with progress notes](planning/4-week-workflow.xlsx)

## Reproducible quick start (Python 3.10+, no packages to download)

```bash
python3 -m venv .venv
.venv/bin/python scripts/smoke_test.py
.venv/bin/python -m unittest discover -s tests -v
```

The smoke check reads `config/default.json`, validates a **fictional metadata-only** example manifest and the three fictional normalized schema examples, and writes safe status messages via the standard-library logger. It does not parse sequences or run a model. On Windows, use `.venv\\Scripts\\python.exe` rather than `.venv/bin/python`.

## Project layout

- `docs/week1/` - approved scope, governance, schema and current task status.
- `planning/4-week-workflow.xlsx` — four-week planning workbook, including dated notes for T01–T28 and current task statuses.
- `config/default.json` — checked-in defaults; no secrets or per-machine credentials.
- `data/manifests/` — fictional example inventory in CSV format, with provenance/review gate columns.
- `data/examples/` — three schematic examples of the normalized output schema, with no biological sequences.
- `src/iri_constructs/` — standard-library metadata and configuration validators.
- `scripts/smoke_test.py` and `tests/` — local and CI checks.
- `.github/workflows/ci.yml` — Python tests on push and pull request; latest verified run passed.

## Boundaries and limitations

The T02 public-data and safety policy is approved; approval of individual records, permissions, and dataset releases remains future work. Only approved public benign records may be used for the future benchmark. Do not commit sequences, restricted or customer data, or secrets. T05 model suitability and reviewed-input inference evidence remain pending; labels and coordinates may be unknown. T04 is hosted in a private GitHub repository with a passing CI run. Early schema examples are fictional and establish structure, not biological correctness or model performance.

## T05 feasibility modules

See [experiments/README.md](experiments/README.md) for the annotation parser, reviewed-input contract, separate model dependencies, CLI commands, metrics and limitations. Run `python scripts/t05_environment.py` to record local dependency/device/cache availability without loading sequences or downloading models.

Week 1 checks remain dependency-free. Optional Biopython tests skip when it is absent; install `requirements-parsing.txt` to exercise real parsing on fictional software fixtures. Pretrained tests require independently approved local source records, labels and fusion coordinates. Generated sequences, embeddings, source qualifiers and reports remain in ignored local directories.
