# T05 local feasibility experiments

> **Historical — superseded.** This page describes the pre-reset T05 feasibility
> experiment interface, from before the toolkit-first reset. The delivered,
> installed preprocessing toolkit is `iri-constructs` (see the top-level
> [README](../README.md) and [docs/project/STATUS.md](../docs/project/STATUS.md));
> the current component-coverage research lives in
> [experiments/component_study_v1](component_study_v1/) with its results under
> [docs/results/component-study-v1](../docs/results/component-study-v1/README.md).
> The commands below are retained for reproducing the historical T05 runs and are
> not the toolkit's interface. Statements such as "T05 remains in progress" or the
> Python/environment notes describe their own time, not the current state.

These modules support the four-week MVP and the parallel skills exercises. They do
not complete the later dataset, paired benchmark, scoring or release tasks. No
model is trained or fine-tuned.

## Environment and commands

Week 1 checks still require only Python's standard library:

```powershell
python scripts/smoke_test.py
python -m unittest discover -s tests -v
```

Parsing is independently installable:

```powershell
python -m pip install -r requirements-parsing.txt
python scripts/parse_annotations.py --inputs data/processed/reviewed-inputs.json --output results/annotation-demo/run01
python -m unittest discover -s tests -v
```

Use a separate **Python 3.11** environment for models. The supplied direct pins
target that interpreter; they are not a tested complete dependency lock. The
existing local environment uses Python 3.14 and was not replaced. On Windows:

```powershell
py -3.11 -m venv .venv-models
.venv-models/Scripts/python.exe -m pip install -r experiments/requirements-esm2.txt
# For ESM-C and ESM-3, use the SDK set instead (includes ESM-2 dependencies):
.venv-models/Scripts/python.exe -m pip install -r experiments/requirements-esm-sdk.txt
.venv-models/Scripts/python.exe -m pip check
.venv-models/Scripts/python.exe -m pip freeze > results/model-environment-lock.txt
.venv-models/Scripts/python.exe scripts/t05_environment.py
# Optional offline encoder/tokenizer load check, without sequence inference:
.venv-models/Scripts/python.exe scripts/t05_environment.py --load-esm2
```

On Linux/macOS, use `python3.11 -m venv .venv-models` and
`.venv-models/bin/python` in place of the Windows commands. CPU is the default.
For CUDA, install the matching platform build of the pinned PyTorch/torchvision
versions using [official PyTorch instructions](https://pytorch.org/get-started/previous-versions/),
then install the requirements. Models run sequentially. ESM-3's memory and CPU
latency requirements have not been measured here; reduce the panel/cap or use
appropriate hardware based on actual observed failures.

```powershell
python experiments/t05_esm2_fusion/compare_fusion.py --inputs data/processed/reviewed-inputs.json --validate-only --output results/fusion-validation/run01
python experiments/t05_esm2_fusion/compare_fusion.py --inputs data/processed/reviewed-inputs.json --allow-download --device cpu --output results/fusion/run01
python experiments/t05_model_comparison/compare_models.py --inputs data/processed/reviewed-inputs.json --models esm2,esmc,esm3 --allow-download --device cpu --output results/model-comparison/run01
python scripts/t05_model_probe.py --inputs data/processed/reviewed-inputs.json --device cpu --output results/probe/run01
```

In these commands, `python` must be the chosen environment's interpreter. All
entry points accept `--help` without model packages. Source-file paths and the
review-manifest path resolve relative to the input JSON, not the working directory.
Use absolute CLI paths when running from another directory. Output must be a new
subdirectory of this repository's ignored `results/` or `data/processed/`.
Existing experiment reports are not overwritten.

Inference defaults to **offline cached weights**. `--allow-download` authorizes
only model-weight downloads; no command downloads biological inputs. Follow model
license/access conditions and authenticate with `hf auth login` if needed. No
credentials are recorded in reports. Omit `--allow-download` for the frozen rerun.

## Reviewed input contract: `t05-inputs-1`

Create the local JSON only after applying the existing T02 record checklist.
This is a separate experiment configuration, not a change to T02's CSV columns
or T03 v0.1. A minimal shape is shown below. Uppercase placeholders are instructions
to supply reviewed values, not approvals or valid example data:

```json
{
  "schema_version": "t05-inputs-1",
  "review_manifest": "REVIEWED_T02_MANIFEST.csv",
  "records": [
    {
      "input_id": "fusion",
      "record_id": "REVIEWED_FUSION_RECORD_ID",
      "path": "REVIEWED_FUSION.fasta",
      "format": "fasta",
      "molecule_type": "protein",
      "sequence_id": "EXACT_FASTA_HEADER_ID"
    },
    {
      "input_id": "target",
      "record_id": "REVIEWED_TARGET_RECORD_ID",
      "path": "REVIEWED_TARGET.fasta",
      "format": "fasta",
      "molecule_type": "protein",
      "sequence_id": "EXACT_FASTA_HEADER_ID"
    }
  ],
  "fusion_cases": [
    {
      "case_id": "REVIEWED_CASE_ID",
      "fusion_id": "fusion",
      "target_id": "target",
      "start0": 0,
      "end0": 1,
      "coordinate_evidence": "REVIEWER, DATE, SOURCE ANNOTATION AND DOCUMENTED AA BOUNDARIES",
      "review_status": "pending"
    }
  ],
  "families": [
    {
      "input_id": "target",
      "family": "REVIEWED_ANNOTATION_FAMILY",
      "annotation_status": "unknown",
      "review_status": "pending",
      "evidence": "SOURCE QUALIFIER OR CURATOR EVIDENCE; REVIEWER AND DATE",
      "independence_group": "REVIEWED_LINKED_CONSTRUCT_OR_HOMOLOGY_GROUP"
    }
  ]
}
```

Add the independent reference records and their family entries for retrieval.
The record, license and safety T02 gates must each be `approved`; the complete
source file must match its reviewed SHA256. Nonempty provenance, source URL,
version, retrieval date and review notes are required. Review notes must state
reviewer, date, benign rationale and source-specific reuse permission as T02
requires; software validates structure and recorded decisions, not their truth.

Inference also requires T02 `annotation_status=approved`. Protein FASTA inputs
are permitted for the feasibility experiments when independently reviewed; T02's
`sequence_length_nt` remains the documented **parent DNA length**, not protein
length. For DNA, use `molecule_type=dna` with `format=genbank` or `fasta`;
the parent nucleotide length is checked. GenBank accession/version and original
feature count must also match the review row. `sequence_id` selects exactly one
record; multiple FASTA entries are never concatenated. Duplicate headers fail.

For reviewed GenBank DNA, eligible CDS translations receive input IDs of the form
`<input_id>/<unit_id>`. Run the annotation parser first to inspect `units.fasta`
headers and normalized/audit files, then use these IDs in `families` or
`fusion_cases`. The model runner repeats validation from the original file;
it never trusts an edited exported protein without a matching source review.
FASTA DNA has no invented CDS/features and produces no protein inference units.

Family entries are scored only when `annotation_status=known` and
`review_status=approved`, with evidence and an independence group. Missing,
ambiguous and pending labels are excluded explicitly. Fusion cases require
reviewed 0-based half-open **amino-acid** coordinates, documented evidence and an
isolated sequence exactly matching that span. Change `pending` only after the
corresponding human review; no script approves records, coordinates or labels.

## Reproducibility and results

Each `report.json` records source/config/manifest digests, source IDs/versions,
model repository and resolved checkpoint revision, code hashes, Git revision and
dirty state, Python/platform and installed-package versions, seed, device, cap,
timings, sampled process RSS, CUDA allocator peaks when applicable, exclusions,
coverage, failures and limitations. No sequences or embedding vectors are written
to experimental reports. The parser's audit qualifiers and annotated protein
exports are local artifacts under ignored output roots.

ESM-2 defaults to the same 8M checkpoint family as the original probe, pinned to
`c731040fcd8d73dceaa04b0a8e6329b345b0f5df`. ESM-C 300M and ESM-3 Small use the
explicitly supported `esm==3.2.1` SDK builders and local snapshot files. Their
initial `main` requests resolve to snapshot SHAs recorded in the report. Freeze
reruns with `--esmc-revision <SHA>` and `--esm3-revision <SHA>` from that report,
retaining the same local inputs/config and complete environment lock. SDK weight
file digests are also recorded. Current SDK APIs have changed, so upgrading the
SDK requires adapter validation rather than silently changing this environment.

All adapters use canonical amino acids, sequence-only input, float32 weights,
deterministic algorithms, no autocast, no silent truncation, and an experimental
shared cap of 1022 residues (adjust downward with `--max-residues`). This cap is
a controlled feasibility restriction, not a claim about every model's maximum.
Each adapter verifies token/residue alignment, removes BOS/EOS/PAD, and exposes
its final output representation. Architectures and representation conventions
still differ. Only within-model similarities and common annotation-level
metrics are comparable.

Timings include the first inference, pooling and host transfer; there is no warm-up.
Download/load costs are separate. Exact-sequence cache hits are flagged and
excluded from per-unique-protein latency medians. RSS is absolute process memory
sampled every 20 ms and may miss brief peaks. CUDA measurements describe the
PyTorch allocator, not total device use. Missing psutil produces null RSS, never
an invented value. A blocked/partial inference exits nonzero and saves a report;
`--validate-only` explicitly records that no inference ran. Zero eligible
retrieval queries yield null agreement metrics and explicit abstentions.

See [fusion](t05_esm2_fusion/README.md),
[model comparison](t05_model_comparison/README.md),
[annotation audit contract](../docs/week1/T05_annotation_audit_extension.md), and
[current status](../docs/week1/STATUS.md).

## Official implementation references

- [ESM-2 pinned model configuration](https://huggingface.co/facebook/esm2_t6_8M_UR50D/blob/c731040fcd8d73dceaa04b0a8e6329b345b0f5df/config.json)
- [ESM SDK 3.2.1 builders](https://github.com/evolutionaryscale/esm/blob/v3.2.1/esm/pretrained.py)
- [ESM-C 3.2.1 output/embedding semantics](https://github.com/evolutionaryscale/esm/blob/v3.2.1/esm/models/esmc.py)
- [ESM-3 3.2.1 output/embedding semantics](https://github.com/evolutionaryscale/esm/blob/v3.2.1/esm/models/esm3.py)
- [Biopython feature locations and extraction](https://github.com/biopython/biopython/blob/master/Doc/Tutorial/chapter_seq_annot.rst)
