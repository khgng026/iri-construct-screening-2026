# Toolkit v0.1.0 release assessment

Assessed 2026-11-17 against [SCOPE.md](SCOPE.md) R1–R7 and the
[contract](TOOLKIT_CONTRACT.md). Agent self-review; no independent review is
claimed. Commands were run on Linux x86-64 with Python 3.10, 3.11, 3.12 and 3.13.

| Criterion | Result | Evidence |
| --- | --- | --- |
| R1 wheel installs and runs outside the checkout | Pass | `scripts/verify_wheel_install.py`: wheel built, installed in a fresh venv in a temporary directory, imported from `site-packages`, CLI and API exercised there. Passed on 3.10.20, 3.11.18, 3.12.3, 3.13.18 with Biopython 1.88 resolved from the index (development used 1.86) |
| R2 `preprocess` single and batch forms | Pass | Installed CLI: single file exit 0; fixture directory exit 2 with seven records written |
| R3 semantics, qualifiers, hashes, parent links, ordered locations | Pass | `test_single_file_bundle_is_complete_and_traceable`, `test_edge_cases_…`: every table row resolves to its original annotation; origin-spanning and reverse joins keep source order; original feature type kept beside the label |
| R4 validated translations, stable reason codes, independence from labels and model limits | Pass | Reason codes declared in `EXPECTED.json` by construction and matched; 1,100-residue protein exported and only flagged |
| R5 coverage of awkward inputs; partial batches keep successes | Pass | Reverse, joined, fuzzy, remote, ordered, zero-length, out-of-range, duplicate, missing/mismatching/non-canonical/incomplete/pseudo, multi-record, malformed, truncated, non-UTF-8, empty, non-sequence |
| R6 determinism and no silent overwrite | Pass | Two runs byte-identical except `run.json`; IDs independent of batch composition; duplicate and colliding records handled; existing output, non-bundle directory and input-inside-output refused |
| R7 examples, provenance notes, offline tests, clean-install demo, limitations | Pass | [examples](../../examples/README.md), README limits, the wheel gate reruns the tests with sockets disabled |

Tests: 271 in the repository suite (266 pass, 5 optional historical skips);
dependency-free mode 271 run with 74 skips (Biopython-dependent tests,
including all toolkit tests, skip without it). Installed-wheel run: 14 toolkit
tests, 13 pass, 1 skip (the pBR322 record is only in the checkout). The 251
pre-existing tests pass unchanged, and `scripts/verify_pilot_release.py` passes.

Real records: 1,000 of 1,000 sampled NCBI vector records (two disjoint seeded
samples) processed with no record-level failure and no parser warning; 2,338 of
2,491 CDSs exported.

## Not verified, or known gaps

- **Continuous integration** has a new `wheel` job (Python 3.10 and 3.12). It was
  still queued on GitHub when this assessment was written, so its hosted result
  is not recorded here; the same gate passed locally on all four versions.
- **Windows and macOS** were not run. Paths use `pathlib`; `resource` (peak
  memory in `run.json`) is absent on Windows and is recorded as null there.
- **Plasmid-editor exports** (missing or inconsistent translations, nonstandard
  LOCUS lines) are covered only by fictional fixtures. No public sample of such
  files was obtained; the real-record evidence is NCBI flat files, which were
  uniformly well formed.
- The network-free claim is tested by disabling sockets in a subprocess around
  the API call, not by an OS-level sandbox.
- A CDS whose coordinates omit the stop codon is withheld even when its
  translation matches in frame (`incomplete_cds`, with
  `in_frame_translation_agrees`). None occurred in the 1,000 real records, so no
  relaxing option was added; if tool exports make this common, it is the first
  behavior to revisit.
- Qualifier text is preserved as parsed (wrapped lines joined), not byte for byte.
- Very large inputs are read into memory per file; no streaming.

## Interface notes

The historical `normalize_record(record, review, sha)` interface and its
outputs are unchanged. The toolkit uses the same core
(`annotations.normalized_features`) in a non-strict mode in which an
out-of-range location or identical duplicate feature becomes a reported outcome
instead of failing the record. `biopython` is now a declared dependency of the
package; the repository's dependency-free checks still run without installing it.
