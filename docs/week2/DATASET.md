# Week 2 curated source evidence

The bounded curation inspected **22 exact source records**: the eleven selected
NCBI candidates, ten selected iGEM entries, and independently reacquired historical
pBR322. Eighteen records qualify for preparatory source QA, and four remain on
source-scope HOLD. **Zero new scientific admissions or execution clearances are
issued by this curation.** The earlier pBR322 E1 release remains historical and
unchanged; inspecting it again does not add an independent experiment.

The versioned [source review](../../data/manifests/week2/source_review.json) stores
source identities, accession/annotation versions, original response digests,
retrieval dates, feature/CDS evidence, rights proposals, attribution, exclusions,
and related-lineage proposals. The [preparatory input manifest](../../data/manifests/week2/inputs.preparatory.json)
binds exactly the eighteen permitted preparatory inputs and their annotation
responses. It grants no scoring permission. Original sequence-bearing responses
remain in ignored `results/week2/curation/raw_sources/`.

## Coverage and exclusions

| Preparatory architecture proposal | Records | Important qualification |
| --- | ---: | --- |
| Multi-CDS | 11 | Six NCBI records and five iGEM entries; an annotated CDS role does not imply an exportable translation |
| Single-CDS | 4 | iGEM entries; none provides a source translation suitable for protein export |
| Backbone-heavy | 3 | pBR322, pSEVA111, and a circular iGEM plasmid proposal |
| Confirmed full-fusion example | 0 in this pool | Fusion-ready vectors and tag descriptions do not establish reviewed full-fusion boundaries |

The separate [MK448012.1 supplemental proposal](../../data/manifests/week2/supplemental_fusion_candidates.json)
adds a genuine source-annotated fusion for technical parser QA only. It is not
silently added to the original candidate pool or a scored panel. Its complete
825-aa source translation validates, but no isolated fusion component or new
retrieval case has been admitted.

Four NCBI records are excluded from the preparatory QA input set:

- **KJ541668.1 and KJ541670.1:** original annotations include a toxin-selection
  component, conflicting with the current T02 toxin-exclusion rule. Their
  reporter translation discrepancies are also retained in the source evidence.
- **KJ170897.1 and PP098727.1:** original cited publications raise source-context
  questions under T02. Their reporter/vector names alone do not settle those
  questions. These are contextual HOLDs pending specialist/scope review, not a
  claim that the exact reporter vectors are inherently hazardous.

The nine historical Addgene HOLDs remain archival records and were not pursued.
No authenticated, restricted, customer, or paid source was used.

## Rights and source semantics

The exact NCBI policy response was captured from the
[NCBI policies page](https://www.ncbi.nlm.nih.gov/home/about/policies/).
The review retains its molecular-data distribution policy together with the
submitter intellectual-property caveat and exact-file notice search. The
proposed local-analysis/private-copy scope is not a declaration of CC0, patent
clearance, biological truth, or authorization for public project release.

Each iGEM response identifies the same explicit license UUID, resolved from the
[Registry license API](https://api.registry.igem.org/v1/licenses/d6c69ca7-8be4-4bc0-b4a8-d3ae1d428aa6)
to CC-BY-4.0. Per-entry provider/source credits, title, public URL, license URL,
and response digest are retained. Attribution and modification notices must
accompany any permitted derived redistribution. Source-copy and scientific
admission remain distinct decisions; original source bytes were not committed.

NCBI records here explicitly describe complete constructs. iGEM supplies
published parts/cassettes as well as circular vector descriptions, so their
source-type classifications remain preparer proposals. Circular topology alone
does not establish complete-plasmid annotation or independent provenance.

The iGEM feature responses are complete single-page responses, but their raw
coordinate indexing convention has not been independently established. Their
original locations, feature UUIDs, role labels and qualifiers are preserved in
the local audit; normalized locations remain unknown. No annotated protein
translation is provided, so no protein is synthesized from DNA to fill that gap.

## Dependence and dataset size

Exact parent and component fingerprints expose duplicate observations without
scoring them. Related vector/cassette groups are recorded independently of
protein-family labels. The pGL3/pGL4 lineages and the GfasPurple, TyrBm, and
YF1/FixJ construct/cassette pairs survive into the QA reports. These links are
conservative proposals for later review, not proof that all other sources are
independent. Multiple vector records can contain identical reporter proteins.

The 150–250 target is deferred: expanding an already redundant pool would not
repair absent translations, unclear fusion boundaries, unsupported family labels,
or missing independent references. The workbook explicitly allows a smaller
well-labeled set. T09 now has a versioned source ledger, coverage and exclusion
counts, while admission to a scientific dataset remains partially complete.

## Reproduction

```bash
# Anonymous acquisition, provider throttling, original-byte cache, source-free manifests:
.venv/bin/python scripts/curate_week2_dataset.py \
  --output results/week2/curation --publish-manifests data/manifests/week2

# Verify the existing exact originals without any network request:
.venv/bin/python scripts/curate_week2_dataset.py --offline \
  --output results/week2/curation --publish-manifests data/manifests/week2

.venv/bin/python scripts/validate_week2.py \
  --inputs data/manifests/week2/inputs.preparatory.json \
  --output results/week2/new-qa-run
```

Curation exits **1** while the four HOLDs remain, even though every original was
acquired and all eighteen permitted preparatory inputs are available. This is an
intentional incomplete-admission signal, not a download failure. Offline mode
requires the original response/cache-metadata pairs and refuses missing,
tampered, or mismatched-URL bytes. A clean checkout must reacquire originals;
committed metadata does not stand in for source bytes. A changed file must be
reviewed and frozen again before scientific use.

`SourceClient` respects provider throttling and `Retry-After`; waits above sixty
seconds become resumable HOLDs rather than early retries. Acquisition is bounded
to the selected identities and response sizes. `--only-id` and `--source` permit
small targeted repairs without repeatedly downloading the entire pool.

Focused validation: sixteen curation tests pass, including exact identities,
provided translations, source-context exclusions, license/pagination checks,
corrupt-cache refusal, offline operation, provider retry handling, and safe
source-free manifest publication.
