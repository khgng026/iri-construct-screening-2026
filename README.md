# Construct preprocessing for function-based screening

The research project investigates whether preprocessing unannotated synthetic DNA constructs into model-ready units improves identification of documented protein components, including in the low-similarity regime. See the [research scope](docs/project/SCOPE.md) and [parallel workstreams](docs/project/PLAN.md). The package described below is the existing annotation-aware toolkit; it helps prepare evaluation material but does not solve sequence-only decomposition.

## Existing toolkit: iri-constructs

`iri-constructs` reads annotated GenBank records and writes traceable features,
**validated** protein units and plain diagnostics. It reuses what the source
record documents. A protein is exported only when the record's own
`/translation` agrees with the documented coordinates, frame and translation
table; everything else is kept and explained. Nothing is predicted.

It needs no review manifest, family labels, reference panel, model, Git checkout
or network access.

## Install and run

```bash
python -m pip install .            # from a checkout; Python 3.10+, installs Biopython
```

```bash
iri-constructs preprocess plasmid.gb --out out/
```

```bash
iri-constructs preprocess a.gb b.gbk more_records/ --out batch_out/ --max-model-residues 1022
```

Inputs can be single-record files, multi-record files, nucleotide FASTA
(sequence only; no annotations are invented) and directories (their `.gb`,
`.gbk`, `.gbff`, `.genbank`, `.fasta`, `.fa`, `.fna` files, not recursive).
A fictional example set is in [`examples/synthetic/`](examples/README.md):

```bash
iri-constructs preprocess examples/synthetic --out demo_out/ --max-model-residues 1022
```

That batch deliberately contains one malformed record and one non-sequence
file, so it finishes with exit status 2 and still writes the other seven records.

From Python:

```python
from iri_constructs import PreprocessOptions, preprocess_files

result = preprocess_files(["plasmid.gb"], "out/", options=PreprocessOptions(model_max_residues=1022))
result.exit_code                 # 0, 2 or 1 as below; the library never exits the process
result.proteins                  # {unit_id: validated source protein}
for record in result.failed_records:
    print(record.input_name, record.failure.code, record.failure.message)
```

`preprocess_file`, `process_inputs` (no files written), `write_bundle`, the
result dataclasses and `PreprocessError`/`OutputError` are importable from
`iri_constructs`.

## What you get

| File | Content |
| --- | --- |
| `normalized.json` | Constructs, features and units with stable IDs, parent links, source hashes, original feature type beside the normalized label |
| `annotations.json` | Original qualifiers and every source location part (including fuzzy and remote parts), keyed by the same IDs |
| `proteins.fasta` | Validated source translations only; header carries unit, construct, feature, length and SHA-256 |
| `features.tsv`, `units.tsv` | One row per source feature / per unit: strand, ordered segments, extraction state, reason code |
| `report.json` | Counts, per-record status, warnings, failures, every reason code with its meaning |
| `run.json` | Timestamps, paths, versions. The only file that differs between identical runs |

Conventions: coordinates are zero-based half-open with strand; a joined or
origin-spanning feature keeps its **ordered segment list** (`[1497:1563](+);[0:57](+)`)
and is never flattened to a bounding interval. TSV is UTF-8, tab separated, one
header row; tab, newline and backslash inside a value are written `\t`, `\n`,
`\\`; an unknown value is empty; segments are joined by `;`. IDs are
`<record id>-f-<digest>` / `-u-<digest>` from the feature's own content, so they
do not depend on file order or on what else is in the batch.

Each feature has one unit. `extraction_state` is `exported` (validated protein),
`withheld` (a CDS that is not exported, with a reason) or `metadata_only`
(noncoding feature, kept as metadata).

## What happens when the input is imperfect

| Situation | Outcome |
| --- | --- |
| CDS without `/translation` | `withheld`, `missing_translation`. No translation is synthesized |
| Translation disagrees with the coordinates | `withheld`, `translation_mismatch` |
| Coordinates are not a complete CDS (no stop, partial codon, no start) | `withheld`, `incomplete_cds`; the detail says which defect and whether the translation still agrees in frame |
| Partial/fuzzy boundary (`<1..300`) | `withheld`, `location_fuzzy`; the uncertain location is kept as written |
| `/pseudo`, `/transl_except`, `/exception` | `withheld`, `exceptional_cds_annotation` |
| Repeated `/translation`, `/product`, `/gene`, … | `withheld`, `conflicting_cds_qualifiers` |
| Remote, `order(...)`, between-bases or out-of-range location | feature kept with the matching `location_*` code |
| Identical duplicate features | both kept; the second gets a `-dup2` suffix and a warning |
| Protein longer than a model limit | **still exported**; `model_compatibility` says `exceeds_max_residues` |
| Record without usable ID (tool export) | ID from the file name plus a content digest |
| Same record twice in a batch | processed once, second marked `skipped_duplicate` |
| Different records with the same ID | both kept, IDs disambiguated by content digest |
| Malformed record in a multi-record file | that record fails with a code; the others are processed |
| File that is not GenBank/FASTA, or empty | that input fails; the batch continues |
| Output directory exists and is not empty | refused. `--overwrite` replaces a previous bundle only, never other files |
| Output directory would contain an input | refused |

Exit status: **0** every record processed; **2** some records failed and usable
output was written; **1** invalid invocation, unsafe output location or no
usable record. Withholding a unit is an expected outcome and never changes the
exit status.

## Limits

- It reproduces source annotations. Agreement between a translation and its
  coordinates is not evidence that the annotation is biologically right, that
  the protein is functional, or of anything about safety.
- No gene prediction, no inferred translations, no domain splitting. A fusion
  annotated as one CDS is one protein; separate CDSs are never concatenated.
- Noncoding features are metadata only. On 1,000 public vector records an
  average of one third of the construct length lay inside exported CDSs.
- Records were validated on NCBI GenBank flat files and on fictional fixtures.
  Exports from plasmid editors are covered by fixtures, not by a public sample.
- Qualifier values are preserved as parsed by Biopython (line wrapping is
  joined), not as raw text.

## Develop and verify

```bash
python -m pip install -e . && python -m unittest discover -s tests
```

```bash
python scripts/verify_wheel_install.py
```

The second command is the release gate: it builds the wheel, installs it into a
fresh environment outside the checkout, checks that no model package is pulled
in, runs the CLI and reruns the toolkit tests against the installed copy with
sockets disabled. See the [release assessment](docs/project/RELEASE_ASSESSMENT.md).

## Research use

The toolkit helps build a traceable answer key and an annotation-derived comparison. The current research compares sequence-only approaches with a model-compatible baseline at matched budgets or across budget curves, and requires a separate low-similarity evaluation. See the [scope](docs/project/SCOPE.md) and [plan](docs/project/PLAN.md).

The earlier [component study](docs/results/component-study-v1/README.md) and [extension](docs/results/component-study-v1/extension-1/README.md) remain evidence of their original conditions: greater coverage came with more queries, and an oracle segmentation result motivated further fusion research. They do not establish the outcome of the new study.

## Project documents

[Scope](docs/project/SCOPE.md) · [status](docs/project/STATUS.md) ·
[low-identity evidence](docs/project/Q3_EVIDENCE.md) ·
[contract](docs/project/TOOLKIT_CONTRACT.md) · [plan](docs/project/PLAN.md) ·
[autonomy and data rules](docs/project/AUTONOMY.md) · [agent instructions](AGENTS.md)

Earlier work — the E1–E4 pilot ([report](docs/EXPERIMENT_STATUS.md)), panel v2
([index](docs/results/panel-v2-2026-11-13/README.md)), the T01–T28 workbook and
per-panel approval records — is historical. [HISTORY.md](docs/project/HISTORY.md)
says what each superseded document is; their results and conditions are
unchanged. Historical runners (`scripts/parse_annotations.py`,
`scripts/execute_cleared_retrieval.py`, …) still enforce their original review
manifests and private caches and are kept for reproduction only.

Public records, reference sequences, embeddings and weights stay in ignored
private storage. Only the fictional fixtures and the previously reviewed
`data/approved_public/J01749.1.gb` are committed.
