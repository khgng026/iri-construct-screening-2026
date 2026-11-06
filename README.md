# Latest research handoff

See [START_HERE.md](START_HERE.md) for the November 5 local research handoff and shared evidence. The Week 1 status below is historical.

# IRI Fall 2026 | Complex-construct preprocessing (Week 1)

**Status (2026-11-06 UTC):** T01-T04 remain approved and T05 remains in progress. The researcher selected ten candidate identities; zero source files are admitted. T07 v1.0 is approved in principle and remains frozen. The [v1.1 E1-E4 addendum](docs/week1/T07_v1_1_evaluation_addendum.md) now has tested annotation-quality, multi-CDS, independent-reference and cost-reporting implementations. Scientific experiments remain blocked by source/label review. See [current status](docs/week1/STATUS.md) and the [oversight handoff](docs/coordination/CODEX_STATUS.md).

## Research question

Does preprocessing **public, benign, annotated DNA constructs** into traceable, model-compatible units improve a fixed pretrained model's characterization of their components compared with intact/minimally processed records under a fair paired evaluation?

The IRI proposal frames this as an approximately **80-hour integration/evaluation MVP**, with 20 hours/week including 2 hours of parallel technical skills. This repository retains T01-T04 and adds bounded T05 feasibility modules and reusable skills exercises; later dataset and benchmark tasks remain incomplete.

## Week 1 official tasks (from `Weekly To-Dos`)

| ID | Calendar task | Artifact | Status |
| --- | --- | --- | --- |
| T01 | Write MVP question, scope, and acceptance criteria | [T01 scope and acceptance](docs/week1/T01_scope_and_acceptance.md) | **Done — approved 2026-10-29** |
| T02 | Set public-data and safety rules | [T02 governance checklist](docs/week1/T02_data_governance.md) | **Done — approved 2026-10-28; reaffirmed 2026-10-29** (record-level reviews still required in T06/T08) |
| T03 | Define construct archetypes and normalized output schema | [T03 schema](docs/week1/T03_archetypes_and_schema.md), [three fictional examples](data/examples/) | **Done — v0.1 baseline approved 2026-10-29** (real-data pilot validation pending) |
| T04 | Create repository and reproducible environment | [T04 setup](docs/week1/T04_repository_setup.md), config, smoke test and unit tests | **Done — accepted 2026-10-29**; hosted privately with passing CI |
| T05 | Run a pretrained model feasibility spike | [Experiments and commands](experiments/README.md) | In progress; reviewed inputs and actual model-run evidence pending |
| T06 | Build public benign construct candidate inventory | [Candidate inventory](docs/week1/T06_candidate_inventory.md) | Implemented: 20 metadata candidates; researcher admission pending |
| T07 | Define supported reference labels and benchmark metrics | [Frozen v1.0](docs/week1/T07_benchmark_protocol.md), [v1.1 extensions](docs/week1/T07_v1_1_evaluation_addendum.md) | v1.0 approved in principle; E1-E4 implemented/tested; admitted-panel validation pending |
| T08 | Curate and manually spot-check a pilot | [Pilot pipeline and queue](docs/week1/T08_pilot_pipeline.md) | Ten selected identities; zero admitted files; scientific validation blocked |

[Four-week workflow with progress notes](planning/4-week-workflow.xlsx) | [Task status and approval history](docs/week1/STATUS.md)

The [v1.1 implementation report](docs/week1/T07_v1_1_implementation_report.md)
records test evidence and blocked experiments. Reproduction commands are in the
[addendum](docs/week1/T07_v1_1_evaluation_addendum.md); first-run source and policy
requirements are in [DECISIONS_NEEDED.md](docs/coordination/DECISIONS_NEEDED.md).

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

## Real-construct pilot contribution

Start with the [T06 inventory](docs/week1/T06_candidate_inventory.md),
[T07 protocol](docs/week1/T07_benchmark_protocol.md) and
[T08 review/preparation commands](docs/week1/T08_pilot_pipeline.md).
The shared `scripts/run_pilot_experiment.py` freezes designs/panels before scoring
and supports annotation preservation (A), documented fusion comparisons (B),
alignment versus embeddings (C), and sequential model-cost comparisons (D).
All real-data commands currently stop at admission. Software fixtures and metadata
validation are distinct from real biological or pretrained-model evidence.

This contribution is isolated on `codex/real-construct-pilot`; PR #1 was read as
research context and is not modified or merged. Its artificial-combination results
are not treated as real-construct validation. See the
[implementation report](docs/week1/T06_T08_implementation_report.md) for provenance,
test evidence, actual blockers and scope limits.
