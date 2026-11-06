# T08: reviewed pilot preparation and annotation-preservation audit

The [ten-record queue](../../data/manifests/t08_pilot_review_queue.json) is proposed
metadata, **not an admitted dataset**: two backbone-heavy, three principal-reporter,
two multi-reporter and three tagged/fusion candidates. Every record needs researcher
decisions. All twenty [T06 candidates](T06_candidate_inventory.md) remain excluded
from scoring. No raw biological source was acquired or approved in this assignment.

## Manual admission checklist

For each record, retain reviewer/date and evidence for:

1. Exact public construct/material identity, immutable source-file version and URL,
   distinguishing a whole construct from an insert accession or catalog variant.
2. Ordinary-benign purpose and record-specific safety suitability under T02.
3. Permitted access, local analysis/storage, attribution and redistribution scope.
   Public visibility or a material-transfer agreement alone does not decide this.
4. A permitted local source, retrieval date and SHA256 under ignored `data/raw/`.
5. Original CDS/noncoding annotations, parent length, feature count, qualifiers,
   coordinates, strand, multipart spans and documented translations/conflicts.
6. Separately supported family labels, independence groups and fusion amino-acid
   boundaries. No model-derived labels or guessed missing features are permitted.

Only then create a complete T02 CSV manifest and local `t05-inputs-1` configuration
under ignored `data/processed/`, following the [existing contract](../../experiments/README.md).
Unknown discovery counts must not be replaced with guessed values to pass validation.
Protein-only independent references are admitted separately through T05; the T08
annotation preparer requires parent DNA records. Discovery history stays unapproved.

## Commands and expected pending outcomes

```powershell
python scripts/inventory_candidates.py --inventory data/manifests/t06_candidate_inventory.json --shortlist data/manifests/pilot_shortlist_v1_1.json --check-pilot-queue data/manifests/t08_pilot_review_queue.json --output results/t06_inventory/run01
python scripts/prepare_pilot.py --inputs config/pilot/inputs.pending.json --output results/t08_pilot/pending01
python scripts/run_pilot_experiment.py --protocol config/pilot/experiment_A.json --inputs config/pilot/inputs.pending.json --output results/pilot/A/pending01
python scripts/run_pilot_experiment.py --protocol config/pilot/experiment_B.json --inputs config/pilot/inputs.pending.json --output results/pilot/B/pending01
python scripts/run_pilot_experiment.py --protocol config/pilot/experiment_C.json --inputs config/pilot/inputs.pending.json --output results/pilot/C/pending01
python scripts/run_pilot_experiment.py --protocol config/pilot/experiment_D.json --inputs config/pilot/inputs.pending.json --output results/pilot/D/pending01
```

Use a new output directory each time. Inventory exit zero means **metadata** checking
succeeded. Other commands intentionally exit nonzero with `no_inputs` or
`blocked_no_reviewed_inputs`: the checked-in pending config contains no admitted
records and no fabricated reviewed CSV. Reports record zero parsing/inference and
scientific scores. These are actual admission-blocker checks, not successful science.

After admission, install `requirements-parsing.txt` and substitute a separate
`data/processed/reviewed-inputs.json` in the same commands. A needs no model packages.
See the [frozen protocol](T07_benchmark_protocol.md) and existing separate model
requirements for B-D; inference defaults to offline and runs sequentially.

## Annotation audit and local artifacts

The preparer reuses `reviewed_sources`, `read_source` and `normalize_record`.
All four T02 gates must pass before scored preservation; aliases, incomplete
provenance and changed files are rejected. Shared source reading checks SHA256
immediately before and after parsing, protecting the existing T05 runners too.

T03 remains v0.1. Reverse, joined and circular spans preserve source order.
Fuzzy/remote/missing/mixed-strand annotations retain uncertainty and their original
audit locations. Only documented `/translation` values verified against exact
CDS/frame/table are exported. No ORF inference, completion or source approval occurs.
Original coding/noncoding qualifiers remain recoverable in the local sidecar.

Experiment A independently checks feature/qualifier/location retention, coordinate
and strand agreement, parent-feature-unit links, usable annotated translations,
withheld/missing/ambiguous CDS and per-archetype failures. Rates retain denominators;
agreement is with supplied annotations, not biological truth. Annotation-free FASTA
produces null rates/zero scored records; shared scientific A returns non-success
when it has no observable annotation evidence.

Successful local artifacts are `normalized.json`, `annotation-audit.json`,
`units.fasta`, `experiment-a.json` and `report.json`. The shared A runner puts these
under `annotation/` beside its freeze/report. Qualifiers and exports stay ignored.

## Retrieval safeguards and limitations

B uses complete documented fusions and reviewed target spans; it never concatenates
CDSs. All target/fusion parents, exact target/full-fusion copies and linked groups
are excluded from a common reference set. Query summaries use the fusion parent
and archetype. Shared target/fusion sources form conservative bootstrap clusters.

C compares raw local BLOSUM62 alignment, predeclared reference-normalized scores
and embeddings on identical references. D runs ESM-2/ESM-C sequentially, reporting
dependency, access, cap and memory failures. Larger-model science was not run here.

Each shared run saves `protocol-freeze.json` before admission, `panel-freeze.json`
before scoring, and `model-freeze-<name>.json` with resolved revision/digest before
that model's first forward. Checked-in `config/pilot/design_freeze.json` hashes the
written protocol and designs. Initial symbolic ESM-C `main` is a placeholder:
retain its resolved immutable SHA/weight digest for reproducible reruns.
The checked-in design hashes normalize text to UTF-8/LF so Git line-ending
conversion does not invalidate them. Biological source-file checksums always
identify the original bytes, as T02 requires; they are never normalized.

Missing reference scores cause abstention rather than favorable panel shrinkage.
K is a total hit budget. A query returning every available reference is descriptive
only and excluded from advantage intervals. Cosine/family retrieval cannot establish
function or improved real-construct characterization without admitted independent
evidence. Multi-CDS parent-level score merging was prospective in frozen v1.0;
the separate [v1.1 E2 implementation](E2_multi_cds.md) now provides tested execution
without claiming a real construct benchmark. The original v1.0 protocol/configs
are unchanged.

## Approved selection and v1.1 extensions

The current queue uses the researcher's exact ten selected identities, including
pGreen2x35SLUC. pUC19c is an optional annotation-limited control after its own T02
review. No selected identity implies approval of a file, annotation or label.
The original queue is retained as `t08_pilot_review_queue_v1_0.json`; validate that
historical queue without `--shortlist`. The original 20-record inventory is intact.

[E1](E1_annotation_quality.md) audits annotation quality and appropriate abstention;
[E2](E2_multi_cds.md) compares longest/all eligible CDSs under five total parent hits;
[E3/E4](E3_E4_generalization_cost.md) freeze reference difficulty and report genuinely
comparable measured costs. Current pending commands are in the
[v1.1 addendum](T07_v1_1_evaluation_addendum.md). All real runs remain blocked.
