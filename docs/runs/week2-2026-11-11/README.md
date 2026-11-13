# Week 2 reproducible evidence

This directory preserves **technical validation and preparatory source QA**.
No new biological model inference, scientific alignment rankings or scored E1–E4
results were produced. [execution-summary.json](execution-summary.json) records
each stage as blocked with null scientific performance/cost fields. Historical
pBR322 E1 artifacts remain unchanged and pass their original verifier.

## Actual checks and source scope

| Evidence | Result | Artifact |
| --- | --- | --- |
| Fresh core Python 3.12 environment | Install and dependency checks passed | [core lock](core-environment.lock.txt), [dependency check](clean-environment-check.txt) |
| Integrated full suite | 246 run, 241 passed, five optional skips, zero failures/errors | [test log](technical/tests.txt) |
| Smoke and historical source/artifact verification | Passed | [smoke](technical/smoke.txt), [historical verifier](technical/historical-release.txt) |
| Main actual-source QA twice | 18 records, 271 features, 18 provided-source proteins; identical reports | [first](technical/curated-1/report.json), [repeat](technical/curated-2/report.json) |
| Separate actual fusion QA twice | One record, 27 features, one 825-aa provided-source protein; identical reports | [first](technical/supplemental-1/report.json), [repeat](technical/supplemental-2/report.json) |
| Clean committed-source archive | Full suite and all nine checks passed; actual-source QA reports exactly match original checkout | [origin](clean-checkout/source-origin.json), [verification](clean-checkout/verification.json), [tests](clean-checkout/tests.txt) |
| Offline ESM-2/ESM-C fictional fixtures | Both exact snapshots passed shape, finite output, float32, pooling, repeat and input-rejection checks | [original report](model/report.json), [fixture log](model/fixture-run.txt) |
| Separate model tests and dependency check | Five runtime tests, four adapter tests and `pip check` passed | [runtime tests](model/runtime-tests.txt), [adapter tests](model/adapter-tests.txt), [dependency check](model/pip-check.txt) |

The dependency-free CI path also passed: 246 tests, 187 passes, 59 optional skips,
zero errors. See [its log](dependency-free-tests.txt). Five parser-environment
core-suite tests skip optional PyTorch/plot dependencies; the model tests
were independently exercised in the dedicated Python 3.11 environment. The model
fixture contains no biological source, family annotation or function expectation.
No fixture can enter E4 as scientific cost/quality evidence.

The original-main-plus-control acquisition contains 22 source records, four of
which remain source-scope HOLD. Only the 18 permitted preparatory inputs enter
the main QA run. A [native Registry export check](registry-export/report.json) additionally
prepared two alternate-format snapshots of existing parts: 33 features and zero
provided translations/protein exports. It does not add biological records.
The genuine fusion is a separately scoped supplement, not an
extra independent scored observation. Original source bytes, full qualifier
sidecars, normalized source objects, protein FASTA, weights and embeddings remain
ignored/private; committed reports contain counts, identities and hashes.

## Costs and reproducibility limits

The integrated nine-command technical validation took **4.61050 s wall** and
**5.85703 child CPU s**, with a child-process high-water RSS of **44,634,112 bytes**.
The clean committed-source archive took **3.56114 s wall**, **4.76680 child CPU s**,
and **46,522,368-byte** child high-water RSS. These are separate technical
verification runs, not independent biological observations or per-model memory
measurements. Acquisition, installation, prior model validation, inspection and
tool overhead are outside that timing scope. GPU time and monetary cost were
not measured and remain null.

The [original model report](model/report.json) separates cached initialization,
first fixture forward and deterministic repeat costs. ESM-2 initialization/first
forward were 3.12688/0.01831 s; ESM-C were 8.16608/0.33019 s. Highest sampled process
RSS was 480,595,968 and 3,115,356,160 bytes respectively. Models were sequential;
weight/OS caches were not experimentally controlled. RSS sampling can miss brief
peaks. These small-fixture timings cannot establish biological-protein latency,
quality-versus-cost tradeoffs or a model winner.

The initial integrated reports identify the pre-publication Git revision and
dirty working tree; their exact code hashes are retained. The later clean
archive is bound to source commit
`a53d77a685e2e30b7f0d35f3ede20dba376b4a59` and its archive SHA-256 in the origin
record. It uses a fresh source archive and isolated core virtual environment;
exact private source-cache bytes were copied locally, not reacquired or presumed
approved. It has no `.git` directory, so its internal Git metadata is null and
the external archive-origin record supplies provenance. Documentation-only
handoff changes follow that code commit.

## Commands

The [validation driver](reproduce.py) never initiates a scored biological
experiment. For a clean checkout with the already committed pBR322 source:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-parsing.txt
.venv/bin/python docs/runs/week2-2026-11-11/reproduce.py \
  --output results/week2/verification-rerun01
```

After the separately documented source acquisition/cache has been restored:

```bash
.venv/bin/python docs/runs/week2-2026-11-11/reproduce.py \
  --output results/week2/verification-rerun02 \
  --curated-inputs data/manifests/week2/inputs.preparatory.json \
  --supplemental-inputs results/week2/supplemental-fusion/inputs.preparatory.json
```

The supplemental local manifest must bind the source and rights/scope evidence
described in [the fusion report](../../week2/SUPPLEMENTAL_FUSION.md); it is not
generated by guessing protein components. Source/annotation response changes
produce different fingerprints and require a new evidence review. Full source
acquisition, reference preparation and model installation/offline commands are
in [the dataset report](../../week2/DATASET.md),
[reference report](../../week2/REFERENCE_EVIDENCE.md), and
[model report](../../week2/MODEL_RUNTIME.md).

Validate the exact scientific-review binding separately:

```bash
.venv/bin/python scripts/freeze_week2_review.py \
  --verify docs/week2/review-freeze.json
```

That verification requires exact private source/reference/model cache bytes.
Missing or changed files block verification. A successful byte check **does not
grant scientific clearance** or authenticate a reviewer message. The canonical
review hash and the recorded real ChatGPT decision must both be checked before
any future scored stage.

[artifact_checksums.json](artifact_checksums.json) binds every committed artifact
in this run directory. [The complete packet](../../coordination/EXPERIMENT_APPROVAL_PACKET.md)
contains the prospective stage scope and remaining external/scientific gates.
