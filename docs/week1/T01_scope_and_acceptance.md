# T01 | MVP question, scope and acceptance criteria

**Week 1 / 1 hour / Approved by researcher (2026-10-29), scope baseline v1.0**  |  Calendar: `Weekly To-Dos!A5:K5`, `MVP Criteria!A4:D13`.

## Research question

**Does decomposing public, benign, annotated synthetic DNA constructs into traceable, model-compatible units improve a fixed pretrained model's characterization of their components compared with a valid intact/minimally processed baseline?**

The research proposal identifies an **input-format mismatch**: multi-part constructs combine coding sequences, tags, regulatory elements and vector/backbone features, whereas many existing models expect defined biological units. The smallest useful demonstration is a reproducible **preprocessing-and-evaluation comparison**, not a new model, a safety classifier, or proof of hazardous-sequence detection.

## Boundaries and design

- **Inputs:** public, benign, versioned, annotated records (prefer GenBank if its feature annotations are sufficient). No customer/proprietary data.
- **Four archetypes:** single CDS; multi-gene cassette; tagged fusion; backbone-heavy construct. Multiple archetype tags are allowed where features overlap.
- **Output:** deterministic parent-to-feature-to-unit mappings with source IDs, coordinates/strand when known, annotation origin, and explicit `unknown`/`ambiguous` states. Unsupported noncoding units are retained as metadata rather than passed to a protein-only model.
- **Comparison:** the same selected pretrained model and settings on matched raw/minimally processed versus decomposed inputs; record raw-input incompatibilities separately. No training/fine-tuning.
- **Benchmark:** a reviewed benign pilot followed by a **150–250-construct** target in Week 2 *if quality permits*, with per-archetype coverage and exclusion counts. Favor a smaller defensible reference set over invented labels.
- **Time:** four weeks / ~80 total hours, including 2 hours weekly parallel skills practice.

**Excluded:** hazardous sequence collection, proprietary orders, wet-lab work, fragment completion, cross-order reconstruction, red-teaming, sub-segment indexing, production deployment and autonomous safety decisions. These are workbook scope guardrails; extensions are future work.

## Measurable MVP acceptance criteria

| ID | Completion evidence (all required unless feasibility limits are explicitly reported) |
| --- | --- |
| **MVP-01** | Documented benign input format; valid records load and malformed ones fail cleanly. |
| **MVP-02** | Stable parent-child IDs, normalized model-ready units and source provenance. |
| **MVP-03** | Original public annotations preserved; absent/conflicting labels marked unknown/ambiguous. |
| **MVP-04** | All four archetypes represented where annotations allow; counts/exclusions reported. |
| **MVP-05** | Paired raw-versus-preprocessed inference with **one fixed model/version/config**. |
| **MVP-06** | Coverage, supported agreement/classification metrics, per-archetype results, runtime and failures reported. |
| **MVP-07** | Reproducible tests, environment/config, CLI, README and benign end-to-end demo. |

**Week 1 sign-off — completed 2026-10-29:** The researcher explicitly approved this MVP question, scope, and measurable criteria as the planning baseline (v1.0). This is researcher approval; no separate advisor signature is asserted. T02 policy v1.0 is approved, and T03 schema v0.1 is an approved preliminary baseline. T05 must establish model/interface and viable label metrics; pilot schema validation and actual performance evidence remain future work. Do **not** assume any improvement before paired tests.

**Source basis:** `4-week workflow.xlsx` (Overview; Weekly To-Dos T01; MVP Criteria) and the IRI proposal, pp. 1–2.