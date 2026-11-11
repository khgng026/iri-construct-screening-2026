# Supplemental real fusion example

The original candidate pool did not establish a complete, translated fusion:
fusion-ready vectors are not completed fusions, and the inspected iGEM feature
responses omit provided protein translations. A separate public record was
therefore acquired for **preparatory parser QA**, using the mission's authority
to investigate suitable alternatives. It is outside the original 21-candidate
scored authorization and is not silently added to an experimental panel.

[GenBank MK448012.1](https://www.ncbi.nlm.nih.gov/nuccore/MK448012.1) is the
complete 7,632-nt pCHX50 cloning vector. Its exact source has 27 features and one
formal CDS, protein QCO92761.1, annotated `3xUbVR-tTA-Dendra2`. The source-provided
825-residue translation agrees with the source coordinates, frame and translation
table. Its deposited reference is Chassin et al., *A modular degron library for
synthetic circuits in mammalian cells*, Nature Communications (2019),
[PMID 31043592](https://pubmed.ncbi.nlm.nih.gov/31043592/),
[DOI 10.1038/s41467-019-09974-5](https://doi.org/10.1038/s41467-019-09974-5).

This is an ordinary synthetic-circuit reporter proposal. TetR/VP16 transactivator
and hCMV/SV40-derived regulatory annotations are retained for case-specific source
review. No source admission, specialist clearance or biological-function claim
is issued by this inspection. The exact source contained no conflicting rights
notice; the NCBI molecular-data policy and submitter-IP caveat remain applicable.
Original bytes and annotation sidecars stay in ignored private project storage.

The source distinguishes the complete fusion CDS from annotated subregions.
The Dendra2 `misc_feature` spans nucleotide coordinates [2454,3192) in normalized
zero-based notation, inside the full CDS [717,3195). These are retained annotation
boundaries, not an independently reviewed isolated-target/fusion retrieval case.
The normalizer exports the source's complete protein only. It neither invents
separate CDSs for the fusion regions nor concatenates proteins.

The notes `NeoR_CDS` and `AmpR_CDS` occur on `misc_feature` annotations, not formal
CDSs. They correctly remain metadata. Turning their names into additional
proteins would invent coding annotations. Repeated preparation preserved all
27 features, one validated protein, stable parent-child provenance and the
unchanged T03 v0.1 schema. See [parser QA](PARSER_QA.md) for final content hashes.

The source-free identity, byte digest, retrieval time, original references,
scope/rights proposals and exact boundary evidence are in
[supplemental_fusion_candidates.json](../../data/manifests/week2/supplemental_fusion_candidates.json).
The original-file SHA-256 is
`038fe3d5fb47562aca2003bc20a06a91809ccab8380221d4573976214b952c07`.

To reacquire the exact version using the shared bounded cache client:

```bash
.venv/bin/python - <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, 'src')
from iri_constructs.dataset_curation import SourceClient, parse_genbank
client = SourceClient(Path('results/week2/supplemental-fusion/raw_sources'))
url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
       'db=nuccore&id=MK448012.1&rettype=gbwithparts&retmode=text')
data, metadata = client.get(url, 'MK448012.1.gb')
record, evidence = parse_genbank(data, 'MK448012.1')
assert metadata['sha256'] == '038fe3d5fb47562aca2003bc20a06a91809ccab8380221d4573976214b952c07'
print(record.id, evidence['feature_count'], evidence['validated_source_translation_count'])
PY
```


The exact source-free preparatory input snapshot is committed separately. After
acquisition, restore it at its original relative location and use a new output:

```bash
cp docs/week2/supplemental-inputs.snapshot.json \
  results/week2/supplemental-fusion/inputs.preparatory.json
.venv/bin/python scripts/validate_week2.py \
  --inputs results/week2/supplemental-fusion/inputs.preparatory.json \
  --output results/week2/supplemental-fusion/qa-new
```

The snapshot contains only identities, checksums and rights/scope metadata, not
biological source bytes. Its relative source path is intentionally resolved after
copying it to the ignored local source directory.

A changed annotated-file snapshot requires new evidence review; an accession
version alone does not freeze every future annotation response. This one source
improves preparatory architecture coverage. It does not establish retrieval
performance, independent family diversity, or unseen-protein generalization.
