# Related enzyme family comparison

The first housekeeping pilot had eight low-identity queries and only one closely competing label. This experiment tests whether the characterization method can distinguish documented enzyme labels within the same broad protein group. It uses natural housekeeping proteins from reviewed UniProtKB with a single EC label and a common InterPro entry for the group. The five groups cover dehydrogenases, aminotransferases, and related lyases. [First pilot](README.md)

## Comparison fixed before scoring

- Choose one reference for each EC label by the taxon priority and accession order in the runner. References come from model bacteria and archaea. Queries come from a separate set of plants, fungi, and animals. Keep one query per taxon and label, deduplicate exact sequences, and exclude exact query accessions and sequence hashes seen in the first pilot.
- Rank only references in the query's group. Every method sees the same competitors. The main metric is the fraction of queries whose top reference has the correct documented EC label. Rank ties break by accession.
- As in the first pilot, call a query low-identity if its same-label alignment has less than 40% identity and covers at least 70% of both proteins. Report intermediate, high, and insufficient-coverage cases separately. A strong competitor has an adequately covered alignment with identity at least as high as the same-label alignment's identity.
- Compare raw local BLOSUM62, reference-self-normalized BLOSUM62, cosine on mean-pooled ESM-2 8M embeddings, and a local alignment on residue embeddings. The two embedding arms share the same encoder pass. For local embedding alignment, subtract the row and column means of the cosine matrix and restore its overall mean. Use affine gaps of -0.15 and -0.02 in those score units, then divide by the shorter protein length. This is a simple prototype, with settings fixed before scores were computed.

The catalog uses curated labels as the answer key. Some source labels may rely on homology. Query taxa and enzyme labels can also be related across cases. The reported counts are descriptive. This study addresses characterization of known proteins.

## Reproduce

Use the project in a Python environment with the dependencies and cached model described in the [pilot handoff](../../docs/project/HANDOFF.md). Acquisition records the provider release, versions, and hashes. Keep the source snapshots and sequence payloads under the ignored results directory. [UniProt terms](https://www.uniprot.org/help/license)

```sh
run=results/related-housekeeping-v1
python -m experiments.benign_protein_pilot.related_study acquire "$run"
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python -m experiments.benign_protein_pilot.related_study \
  score "$run" --out "$run/report.json"
```
