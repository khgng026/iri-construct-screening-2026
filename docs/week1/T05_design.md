# T05 frozen design — 2026-10-29

Recorded before any model inference or retrieval. This is an independently written,
protein-only feasibility experiment, inspired by the supplied confidential
presentation (methods pp. 11–14). No code or datasets from its author are available.
His ESM-2 650M results are reported claims, not results reproduced here. The PDF
and its contents are not copied into this repository.

## Scope and data review

The user's T05 instruction explicitly authorizes a small reviewed benign protein
catalogue and artificial examples. This narrow protein-interface test extends
T02's construct-only input scope for T05; it does not amend T02 or approve a
T06/T08 DNA benchmark. Apply T02's provenance, licensing, exclusion, local storage,
and frozen-manifest rules individually. Review is performed by Codex against the
public records under this instruction, not represented as Khang's sign-off or
independent institutional biosafety certification.

Five candidates fixed in advance: UniProt P42212 (GFP), P00560 (yeast PGK1),
P00924 (yeast ENO1), P00942 (yeast TPI1), P00330 (yeast ADH1). Use full canonical
sequences, without removing initiator residues, mutating, completing, or designing
sequences. Confirm reviewed status, identity, ordinary reporter/metabolic purpose,
versions, and reuse terms before admission. Any unclear candidate remains outside
the scored set; report a changed denominator rather than silently replacing it.

## Inputs fixed before scoring

Expected target is GFP (P42212) for every query:

| ID | Ordered full components | Purpose |
| --- | --- | --- |
| clean_gfp | P42212 | Clean self-retrieval sanity check |
| gfp_pgk | P42212, P00560 | Artificial protein concatenation |
| pgk_gfp | P00560, P42212 | Same components with target position reversed |
| gfp_pgk_adh | P42212, P00560, P00330 | Artificial three-component concatenation |

Join amino-acid strings directly, with no linker, separator, translation, random
filler, or biological claim about folding/function. All are digital test inputs,
not observed fusion proteins or proposed wet-lab constructs. The oracle arm slices
the exact GFP interval from the constructed query and **re-encodes that substring**.
It does not pool GFP rows from the whole-query forward pass. Boundaries and target
identity are supplied by construction, not discovered. No segment selection or
max-score aggregation. Clean whole/oracle inputs are intentionally identical.

## Fixed model and scoring

Use `facebook/esm2_t6_8M_UR50D`, immutable revision in `config/t05_model.json`.
CPU float32, four threads, seed zero, eval mode and inference mode. Mean of last
layer embeddings over actual residue tokens only. Rank all five reference vectors
by cosine similarity descending; exact ties use accession ascending and report
tie counts. No training, tuning, identity holdout, clustering, or significance test.
Encode catalogue and each arm separately; do not cache the oracle vector as a result.

Operational input cap: **1,022 residues**, selected conservatively for this pilot,
not claimed to be a hard architectural limit of rotary ESM-2. Reject longer,
empty, noncanonical, or multi-record inputs; never truncate. Record tokenizer and
model configuration limits separately. No raw sequence strings in Git or logs.

## Required observations

Save complete rankings, target ranks/scores, catalogue denominator, sequence and
source checksums, embedding arrays locally, shapes, model/package versions,
per-arm forward/pooling/ranking times, model load/download and total times,
failures/exclusions. Rerun in a fresh process with local-only model loading and
compare all scores (absolute tolerance 1e-6), ranks, vectors, and input/config hashes.
Timing variation is expected. Check pooling using masked sentinel values and a real
single-versus-padded-batch pass. Test invalid input, provenance gates, checksums,
foreign keys and protein coordinates; preserve the original 11 tests.

Self-retrieval is a wiring test. The oracle returns an exact catalogue entry, so
its success is expected by construction. Three related artificial queries with
one target and five references cannot establish general performance, generalize
to unseen families, prove function, or demonstrate automatic segmentation. Model
training-set overlap is unknown. No hazardous or proprietary material is used.
