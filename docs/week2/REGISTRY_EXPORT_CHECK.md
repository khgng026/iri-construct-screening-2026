# Bounded check of native Registry GenBank exports

Native annotated GenBank exports are anonymously available for both inspected
selected parts. They improve the available source-format evidence, but **do not
repair missing provided protein translations**. No existing source manifest,
schema freeze, model freeze or reference proposal was changed by this check.

The provider's [OpenAPI documentation](https://api.registry.igem.org/docs)
documents `GET /v1/parts/{partUUID}.gb` as
`PartsController_getGenbank_v1`. This is a documented export endpoint, not a
guessed URL. The original Swagger initialization response, its exact digest,
the two export URLs, retrieval timestamps and original-byte digests are recorded
in [registry_export_check.json](../../data/manifests/week2/registry_export_check.json).
Seven public provider requests were used, below the twenty-request investigation
bound; the two export requests were separated by thirteen seconds. No login,
restricted endpoint, inference or scientific scoring was used.

| Existing selected part | Native GenBank features | Original REST feature rows | Formal CDSs | Provided translations | DNA matches original part API |
| --- | ---: | ---: | ---: | ---: | --- |
| BBa_26NKZH8J | 7 | 5 | 1 | 0 | Yes, exact equality |
| BBa_25FDWALN | 26 | 24 | 3 | 0 | Yes, exact equality |

The additional two features in each provider export are source/composite-part
bookkeeping. They are not silently substituted into the original REST feature
denominators. These are two format snapshots of **existing** selected parts, not
two independent additional biological records. The exact parent DNA fingerprints
match the original cached part API responses.

The existing preparatory pipeline parsed both exact GenBank files successfully:
two records, 33 source features, zero blocked records and zero exported proteins.
All four formal CDSs lack `/translation` qualifiers, so the parser retains
explicit unsupported model-input states. No ORF was inferred and no protein was
translated to manufacture coverage. The first part's deposited tag/fusion
description does not provide an independently validated complete protein or
reviewed amino-acid boundaries.

The OpenAPI `SequenceFeatureLocationDTO` specifies a minimum of one for `start`
and `end`, but does not explicitly describe end inclusivity. In these two native
exports, all 29 matched feature locations follow the observed one-based-inclusive
mapping: for example, REST `start=1,end=19` becomes parsed GenBank `[0:19)`.
The complete multisets of label/strand/location tuples agree after this
conversion and removal of the two provider bookkeeping features; interpreting
the raw tuples as already zero-based does not agree. Native exports do not
retain REST feature UUIDs, so this verifies location multisets rather than
claiming per-UUID feature identity or cross-matching repeated label text.
This is exact provider-export evidence for these snapshots, **not** a universal
REST indexing contract or independent biological validation. Native GenBank
locations can be preserved using that format's documented semantics in a future
separately frozen input; the original REST QA continues to retain unknown
normalized coordinates rather than changing its earlier assumptions silently.

Original responses, extracted OpenAPI, export files, local input manifest and
sequence-bearing normalized/qualifier artifacts remain ignored under
`results/week2/registry-export-check/`. The source-free manifest binds their
fingerprints and the QA report. It retains the existing part-level CC-BY-4.0
attribution/rights proposals without granting source or experiment admission.

To repeat the technical QA using the exact cached files and a fresh output path:

```bash
.venv/bin/python scripts/validate_week2.py \
  --inputs results/week2/registry-export-check/inputs.native-exports.preparatory.json \
  --output results/week2/registry-export-check/qa-rerun
```

Only these two exports were inspected. The result does not establish that all
Registry formats or parts lack translations. Any future replacement of the
original REST snapshot with a native export requires an explicit new input
freeze and the applicable separate scientific review before scoring.
