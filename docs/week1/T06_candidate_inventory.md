# T06 candidate inventory and manual admission queue

This version contains **20 publicly documented candidate records**, with **10
prioritized for pilot review** and **zero admitted to parsing or scoring**.
Primary candidate categories cover 4 backbone-heavy, 6 single-principal-reporter,
3 multi-reporter and 7 tagged/fusion constructs. These are discovery labels from
public purpose/product descriptions, not verified feature-table assignments.
No biological source record, family label, license or safety decision is approved.

The versioned [metadata inventory](../../data/manifests/t06_candidate_inventory.json)
is a separate `t06-candidates-1` discovery schema. It does not change T02's approved
CSV fields or T03 v0.1. Unknown source-file versions, nucleotide lengths, feature
counts and checksums remain explicit; no plausible numbers were substituted.

## Candidate records

| ID / primary source | Construct | Proposed primary archetype | Pilot review priority |
| --- | --- | --- | --- |
| [L09137.2](https://www.ncbi.nlm.nih.gov/nuccore/L09137.2) | pUC19c | backbone_heavy | 1 |
| [J01749.1](https://www.ncbi.nlm.nih.gov/nuccore/J01749.1) | pBR322 | backbone_heavy | 2 |
| [Addgene 50004](https://www.addgene.org/50004/) | pUC18 | backbone_heavy | reserve |
| [Addgene 50005](https://www.addgene.org/50005/) | pUC19 | backbone_heavy | reserve |
| [Addgene 212936](https://www.addgene.org/212936/) | pGL3 Basic Vector | single_cds | 3 |
| [Addgene 18964](https://www.addgene.org/18964/) | Luciferase-pcDNA3 | single_cds | 4 |
| [Addgene 87121](https://www.addgene.org/87121/) | pcDNA-RLuc8 | single_cds | 5 |
| [Addgene 12178](https://www.addgene.org/12178/) | pIS0 | single_cds | reserve |
| [Addgene 12179](https://www.addgene.org/12179/) | pIS1 | single_cds | reserve |
| [Addgene 16539](https://www.addgene.org/16539/) | pBV-Luc | single_cds | reserve |
| [Addgene 18673](https://www.addgene.org/18673/) | pYIC | multi_cds | 6 |
| [Addgene 64034](https://www.addgene.org/64034/) | pGL4Luc-RLuc | multi_cds | 7 |
| [Addgene 108228](https://www.addgene.org/108228/) | pGreen2x35SLUC | multi_cds | reserve |
| [Addgene 21073](https://www.addgene.org/21073/) | pEGFP-LC3 | tagged_fusion | 8 |
| [Addgene 21074](https://www.addgene.org/21074/) | ptfLC3 | tagged_fusion | reserve |
| [Addgene 24920](https://www.addgene.org/24920/) | pEGFP-LC3 (human) | tagged_fusion | reserve |
| [Addgene 51970](https://www.addgene.org/51970/) | Nano-lantern/pcDNA3 | tagged_fusion | 9 |
| [Addgene 51819](https://www.addgene.org/51819/) | Syn-ATP | tagged_fusion | reserve |
| [Addgene 113442](https://www.addgene.org/113442/) | pcDNA3.1 NL | tagged_fusion | reserve |
| [Addgene 255006](https://www.addgene.org/255006/) | pcDNA3-Renilla-FLAG-firefly | tagged_fusion | 10 |

Priorities identify ten **review requests**, not an approved 10-record dataset.
Every inventory row records proposed inclusion merit, an explicit exclusion from
scoring, original source identity/link, access date, annotation evidence/limits,
terms links, unresolved issues and a conservative relatedness group. Reserves
provide alternatives if permissions or annotations disqualify prioritized records.

## What was actually verified

Public Addgene catalog identity, purpose/product/tag metadata and source-specific
terms summaries were inspected on **2026-11-05**. Download links request login;
no account authentication or biological sequence download was attempted. Material
IDs are persistent, but catalog pages expose multiple sequence entries without
an admitted immutable source-file version. `accession_version="unknown"` records
that gap rather than treating a material ID or insertion accession as a complete
versioned construct.

Direct NCBI nucleotide pages presented a browser challenge. The pUC19c identity
and version are visible in [indexed NCBI metadata](https://www.ncbi.nlm.nih.gov/nuccore/L09137),
with [official NEB pUC19 map/product support](https://www.neb.com/products/n3041-puc19-vector).
The [official NCBI Entrez Direct examples](https://www.ncbi.nlm.nih.gov/books/NBK565821/)
identify J01749.1/pBR322 and its cloning-marker context. The inventory explicitly
distinguishes these metadata observations from direct source feature-table review.
No original GenBank qualifiers or coordinates were fabricated or copied into Git.

Addgene material terms, including [UBMTA](https://www.addgene.org/agreement/1/) and
record-specific limited-use labels, require a source-specific researcher decision;
they do not independently establish permission to redistribute sequence data.
The [website terms](https://www.addgene.org/terms-of-use/) are an additional reference.
For NCBI, its [usage policy](https://www.ncbi.nlm.nih.gov/home/about/policies/)
distinguishes government-created information from contributed material. These
observations leave licensing `uncertain` for all records. No legal conclusion or
automatic source approval is asserted.

## Recorded ambiguities and independence limits

The inventory retains pYIC's disagreement between tag summaries and depositor
comments, pEGFP-LC3's missing insert start codon note, and ptfLC3's reported
variant/stop-extension issue. Syn-ATP includes a modified luciferase, and several
catalogs provide multiple sequence entries. These require exact-source manual
checks before selecting complete CDS products, mature products, or fusion spans.

`single_cds` means one principal reporter as a discovery hypothesis, not a claim
that the entire plasmid has only one CDS. A vector's capacity for a later ORF
insert does not establish the presence of that insert. A complete fusion can be
one CDS encoding several components; it must not be counted as several CDS solely
because its name contains several proteins. The dual-reporter candidates describe
separate products; their exact source annotations remain pending.

pUC records, LC3/EGFP constructs, related firefly reporters and Rluc8-containing
Nano-lantern are conservatively linked for review. These groups are not validated
homology clusters, and the catalog count is not the number of independent examples.
Reviewers must resolve exact duplicates, shared backbones and family relatedness
before reference-panel freezing. No alignment/embedding advantage is inferred.

## Reproduce validation and the queue

The metadata command needs only the standard library, does not access the network,
and cannot mark any source approved:

```powershell
python scripts/inventory_candidates.py --inventory data/manifests/t06_candidate_inventory.json --output results/t06_inventory/run01
python -m unittest discover -s tests -p test_candidate_inventory.py -v
```

`review_queue.json` lists all 20 records as excluded from scoring, their specific
missing reviews/version/count/checksum fields, the ten proposed pilot IDs and the
manual checklist. `report.json` records inventory SHA256 and execution metadata.
Both are ignored local outputs. Exit 0 means metadata/queue generation succeeded;
it does not mean any biological record was admitted or experiment succeeded.

The loader rejects duplicate IDs/source identities, incomplete provenance,
credential-bearing URLs, invalid dates/counts, raw payload fields and **every
approved gate**. Human decisions belong in a distinct complete T02 admission
manifest, retaining the original candidate inventory as discovery history.

## Researcher checklist for every record

1. Select and identify the exact public construct/source-file version. Resolve
   accession-versus-insert and catalogue variant differences.
2. Record reviewer/date and the ordinary-benign rationale under T02; reject or
   hold any uncertain source.
3. Establish permitted access, analysis, local storage, attribution and any
   redistribution limits for that exact source, not merely the physical material.
4. Obtain only the permitted local source and record its retrieval date/SHA256.
5. Manually inspect original coding/noncoding annotations, parent length, feature
   count, strands, multipart coordinates, translations and any documented conflicts.
6. Resolve fusion component evidence/bounds and keep unknown/ambiguous states when
   evidence is insufficient. Do not infer or complete unsupported features.
7. Separately approve supported family labels and independent query/reference
   groups under the frozen T07 protocol.
8. Create the complete T02 admission row and reviewed local-input mapping only
   after the applicable decisions. The pilot preparer and existing T05 runners
   remain responsible for enforcing those independent gates/checksums.

T06's metadata inventory and review tooling are implemented. Record admission,
biological validation and scientific T08 pilot evidence remain pending researcher
review; zero scored results or successful model inferences are implied.

The metadata CLI also writes a derived `pilot_review_queue.json`. The checked-in
[T08 pending queue](../../data/manifests/t08_pilot_review_queue.json) contains the
ten prioritized records and binds a canonical-JSON inventory SHA256, independent
of Git line endings. The inventory run report separately records the local file-byte digest. Add
`--shortlist data/manifests/pilot_shortlist_v1_1.json --check-pilot-queue data/manifests/t08_pilot_review_queue.json` to reject stale or
altered queue metadata. Editing counters/statuses cannot promote it to admission;
source approvals remain separate T02 decisions.

The v1.1 shortlist reflects the researcher's approved ten identities and includes
`addgene-108228` (pGreen2x35SLUC). pUC19c is optional after its own review.
`t08_pilot_review_queue_v1_0.json` preserves the original proposed queue;
the twenty-record discovery inventory remains unchanged. The approved selection
manifest is metadata, with zero admitted biological records.
