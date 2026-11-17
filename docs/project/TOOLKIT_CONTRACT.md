# Preprocessing toolkit v0.1 contract

This is the accepted target contract for the implementation stretch, not a claim that the new CLI already exists. Keep the biological extraction approach unchanged. Reuse `normalize_record` and strengthen its surrounding interface rather than adding a second normalizer.

## Public interface

- Python API: `preprocess_file(input_path, output_dir, *, options=None)` returns a structured batch result. Public return types and exceptions are documented and importable from the installed package.
- CLI: `iri-constructs preprocess INPUT [INPUT ...] --out DIRECTORY`; accept multiple files and multiple GenBank records per file. Inputs may be anywhere the user can read; output may be any writable user-selected directory. Validate paths, prohibit accidental source overwrite and do not overwrite existing result bundles by default.
- GenBank is primary. Existing FASTA support may remain as sequence-only metadata with zero invented annotations; it must not imply gene prediction. No network call or download occurs in core preprocessing.
- Minimal options cover input/output and explicit overwrite policy if implemented. Add flags only for demonstrated needs. Help, API import and preprocessing do not import model libraries.

## Output bundle

- `normalized.json`: versioned constructs, source features, units, source identity/hash and coordinate conventions; original feature type plus normalized label; full-parent fusion retained.
- `annotations.json`: original qualifiers and locations linked to normalized IDs, or an equally explicit documented lossless representation. Do not put the original biological payload into summary logs.
- `proteins.fasta`: only source-provided translations validated against documented coordinates/frame/table; no silently synthesized replacement.
- `features.tsv` and `units.tsv`: stable IDs, parent links, strand/ordered segment representation, extraction state and reason codes. Specify escaping and how compound values are encoded.
- `report.json`: schema/tool/configuration version, counts, per-record/per-unit outcomes, warnings/failures and run metadata. Mark annotation agreement separately from biological truth. Optional model compatibility is reported separately from extraction correctness.

Generate an input manifest from the supplied files; no external T02 approval CSV, family labels, reference panel or PR decision is required. Preserve optional provenance, distinguish missing from known fields and never manufacture approval values.

## Correctness and partial results

Use zero-based half-open coordinates with strand and ordered source segments. Do not flatten joined/origin-spanning locations into misleading intervals. Preserve uncertain locations as uncertain. Distinguish a raw byte hash from a sequence hash and from scientific provenance assertions.

IDs must be deterministic for the same input/configuration; duplicate records/features must be detected and handled without overwriting one another. Distinguish expected withholding from a parser failure. Keep valid noncoding metadata and other valid proteins when one unit is unsupported. Record-level malformed input may fail that record while a batch continues when safely parseable; do not promise recovery from an unreadable file.

A valid source protein longer than 1022 residues remains an extracted protein; flag incompatibility with a chosen model separately. Missing translations, conflicting qualifiers and exceptional annotations retain the current conservative withholding behavior. Do not introduce gene calling or domain splitting to make coverage look better.

CLI exit codes: 0 when every input record processes without record-level failures; 2 for a completed partial batch with at least one record-level failure and at least one usable normalized record; 1 for invalid invocation, fatal I/O/configuration failure or no processable records. Expected unit withholding is visible in the report and does not alone make a batch fail. Library callers receive the structured results rather than having the process exit.

## Release gate

Build a wheel, install it in a new virtual environment, then invoke help/API/preprocessing on fixtures from a working directory outside the repository with no `.git`, historical manifests, model stack, private data cache or network access required by preprocessing. Verify exports/provenance and partial-batch behavior. Declare parsing dependencies through packaging. Keep model dependencies in optional extras or a separate environment. Preserve existing useful regression checks and add focused tests for the new user-facing behavior.

## Research bridge

The outputs must support joining units to their original construct, selecting the longest eligible protein or all eligible proteins, adding separately reviewed family labels and loading the same exported proteins into a fixed downstream model. Freeze any comparison's actual inputs/configuration and score only supported endpoints. A small optional integration demonstration establishes interface compatibility; evidence of general retrieval improvement requires a separate appropriate study.
