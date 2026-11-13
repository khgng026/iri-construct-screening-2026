# Reference evidence and prospective query panel

The completed reference preparation contains **17 exact reviewed UniProt
entries**, four proposed molecular homolog classes and five source-bound query
proposals. Four query proposals across three parent constructs have proposed
positive and competing references after conservative exclusions. The fifth,
Renilla hRluc, has no independent positive and remains HOLD.

These are **unscored proposals requiring delegated scientific review**, not a
new T02 admission or approved benchmark. No biological model inference,
alignment score, retrieval ranking or E1–E4 result was produced by this work.

## Sources and family evidence

[The catalogue](../../data/manifests/week2/references/catalogue.json) records
explicit accession/taxon/lineage selections and provider-rights evidence.
[The source-free reference proposal](../../data/manifests/week2/references/reference_proposal.json)
records original response SHA-256, retrieval date, entry and sequence versions,
protein fingerprint, curated family annotations, Pfam/InterPro evidence and
original literature citations for each entry. Original sequence-bearing JSON
and the original policy responses remain in ignored private results.

| Proposed class | Reviewed natural reference entries | Scope and caveat |
| --- | --- | --- |
| GFP-like homologs | P42212, Q9U6Y8, Q8ISF8, Q9U6Y4, Q9U6Y3, Q9U6Y6, P83690 | Molecular GFP family, PF01353. Fluorescence colors and nonfluorescent chromoproteins do not become separate families. |
| Coleopteran luciferase homologs | P08659, Q26304, P13129, Q01158, Q27757 | Reviewed luciferase identity, Coleoptera lineage and curated AMP-binding-family/Pfam evidence are jointly required. AMP-binding domains alone do not establish luciferase identity. |
| OMP-decarboxylase homologs | P03962, P08244, P25971, O26232 | Reviewed standalone OMP-decarboxylase/PF00215 family. A molecular family assignment is proposed for the documented URA3 component; growth-selection phenotype alone is insufficient. |
| Renilla luciferase lineage | P27652 | Source-lineage anchor only. The catalogue contains no independent positive for the Renilla-derived query. It is not merged with firefly luciferase solely because both emit light. |

The [curated GFP family description](https://www.ebi.ac.uk/interpro/entry/InterPro/IPR000786/)
and its original literature, including
[PMID 12325128](https://pubmed.ncbi.nlm.nih.gov/12325128/) and
[PMID 10852900](https://pubmed.ncbi.nlm.nih.gov/10852900/), explicitly support a
homology class spanning different colors and nonfluorescent members. The source
U55762.1 cites the GFP mutant study
[PMID 8707053](https://pubmed.ncbi.nlm.nih.gov/8707053/). The explicit
[query specification](../../data/manifests/week2/references/query_specification.json)
retains the factual basis and links for every proposed assignment; the software
does not guess a class from an unreviewed product string.

UniProt's [database license](https://www.uniprot.org/help/license) permits
copyrightable database material under CC-BY-4.0 with attribution. The exact
policy-response digest, retrieval date, modification date and attribution to
the UniProt Consortium are in the catalogue. This is a documented rights
proposal, not a claim about patents, biological correctness or authorization to
publish project source data. New source and family/scoring decisions remain
separate fields and are false/pending in every entry.

## Query identities and prospective eligibility

Queries are bound to original source-provided CDS translations and deterministic
normalized units. The builder verifies source file digest, source protein ID,
protein fingerprint, source segments and protein length against the preparatory
unit manifest. It does not infer ORFs or translate missing iGEM annotations.

| Parent / source protein | Proposed class | Positive / competing / total eligible references | Current evidence status |
| --- | --- | --- | --- |
| U47295.2 / AAA89082.1 | Coleopteran luciferase homologs | 4 / 12 / 16 | Proposed supported subset, pending review |
| U55762.1 / AAB02574.1 | GFP-like homologs | 6 / 10 / 16 | Proposed supported subset, pending review |
| AF298787.1 / AAG34528.1 | OMP-decarboxylase homologs | 3 / 12 / 15 | Proposed supported subset, pending review |
| AF298787.1 / AAG34529.1 | GFP-like homologs | 6 / 9 / 15 | Proposed supported subset, pending review |
| AY864931.1 / AAW66991.1 | Renilla luciferase lineage | 0 / 16 / 16 | **HOLD**, no independent positive |

All proposed catalogues exceed K=5, so none of these prospective rows is
budget-saturated. That does not establish scientific independence or predict
model performance. There are at most **three parent constructs** among the four
proposed-supported query rows. Shared vector components may reduce independent
groups further. The panel cannot meet the eight-independent-cluster threshold
for a comparative bootstrap interval.

Natural source ancestors are excluded even when their exact sequence differs
from an engineered query: P42212 for Aequorea-derived GFP, P08659 for
Photinus-derived luciferase, P03962 for yeast URA3 and P27652 for Renilla. Every
AF298787 component uses the same parent-wide ancestor exclusions, including
both P42212 and P03962. The common parent component fingerprint set includes
the unsupported AmpR component, so a copy cannot leak into any condition.
Same-parent, same-lineage, exact-copy and known-ancestor removals are recorded in
[exclusion_evidence.json](../../data/manifests/week2/references/exclusion_evidence.json).
No reviewed strict sequence-relatedness policy was invented; tier 3 remains HOLD.

## Full denominators and dependence findings

Unsupported source components remain visible in
[parent_component_coverage](../../data/manifests/week2/references/query_proposals.json).
The proposed supported subsets do not claim complete-family coverage:

- U47295 has two validated source CDSs; only its luciferase has a reference-family
  proposal. AmpR is retained as an unsupported component.
- U55762 has two validated source CDSs; only its GFP has a reference-family
  proposal. The 264-aa neomycin marker is longer than the 239-aa reporter.
- AF298787 has three validated source CDSs: 267-aa URA3, 238-aa GFP and 286-aa
  AmpR. The **longest eligible supported** unit would be URA3, while the literal
  source-longest translation is unsupported AmpR. A future supported-subset E2
  comparison must disclose that difference and retain the unresolved full
  documented-family denominator.
- AY864931 has three validated source CDSs; even its proposed Renilla class lacks
  an independent positive, so it is not a retrieval-ready parent.

Exact fingerprint inspection, without alignment or model scores, found that
U47295.2 AAA89082.1, U47296.2 AAA89084.1 and AY738222.1 AAV52869.1 encode the same
550-aa reporter fingerprint
`9023bdf7532e7e7ce555faa66e0ebac09e7c2d8c2e997b51cf4dd572b45a6f43`.
Their 286-aa AmpR fingerprints also match. Only U47295 is retained as a
prospective reporter representative; three vector identities are not three
independent protein observations. The core dataset manifest preserves all
original source identities and their distinct parent-DNA/annotation evidence.

The four T02 source/context HOLD records are excluded. iGEM records without
source-provided translations contribute annotation QA, not invented protein
queries. The supplemental MK448012.1 fusion remains a separate preparatory
schema example and is not silently added to this panel.

## Reproduction and validation

After source curation and preparatory QA, use the cached originals:

```bash
.venv/bin/python scripts/prepare_reference_dataset.py \
  --output results/week2/references \
  --query-specification data/manifests/week2/references/query_specification.json \
  --unit-provenance docs/week2/qa/unit-provenance.json \
  --source-manifest data/manifests/week2/source_review.json
.venv/bin/python -m unittest discover -s tests -p test_reference_dataset.py -v
```

For an initial clean acquisition, add `--download` and use a new ignored output
folder. The downloader retrieves the same accession identities, validates each
actual returned entry/sequence version and records its exact bytes. Provider
updates legitimately produce a different evidence fingerprint and require a
new review; a new download is not presumed identical to an approved snapshot.
Existing cached source bytes are hash-verified before reuse. An incomplete cache
pair is blocked rather than silently overwritten. The committed proposal binds
the source manifest, query specification, catalogue and unit-provenance digests.

Eleven focused software tests pass, covering source/entry identity, versions,
taxon, alphabet/cap, independently citable family evidence, all-parent-component
exclusions, ancestors, duplicates, missing positives, source/coordinate tampering,
unsupported full denominators and deterministic freeze content. Their fixtures
are explicitly fictional and do not constitute scientific family approvals.
