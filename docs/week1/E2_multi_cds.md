# E2: parent-level documented-CDS comparison, T07 addendum v1.1

**Implemented; real-data validation blocked by source/label admission.** This
extension implements the prospective multi-CDS comparison without modifying
the frozen [T07 v1.0 protocol](T07_benchmark_protocol.md) or T03 v0.1 schema.
No biological record or family is approved by this implementation. It creates
no protein concatenations, inferred translations or DNA inputs for protein models.

```powershell
python scripts/run_multi_cds.py --protocol config/extensions/E2.json --inputs config/pilot/inputs.pending.json --output results/extensions/E2/pending01
```

The checked-in pending input config intentionally produces
`blocked_no_reviewed_inputs`, a nonzero CLI exit, and zero scientific scores.
After exact-record admission, replace it with a local `t05-inputs-1` config under
ignored `data/processed/`; use a new output directory for every execution. The
shared extension runner calls the same `run_multi_cds` API.

## Inputs and frozen comparison

The existing T02 source CSV, checksums and all four review gates apply before
parsing. Parent sources are reviewed DNA FASTA/GenBank records. Only their
documented `/translation` qualifiers passing existing coordinate/strand/frame
validation become proteins. FASTA with no CDS cannot supply a parent protein
representation. Independent protein references use the existing reviewed FASTA
interface. Source bytes, annotations, proteins and model weights stay local.

Use explicit `reference_ids` in the input JSON and, preferably, `parent_ids`
listing T02 construct IDs; the latter otherwise selects all reviewed DNA sources.
For a GenBank unit, its existing input ID is `<source input_id>/<T03 unit_id>`.
The `families` entries independently declare `known`/`approved` family annotations,
evidence and reviewed independence groups. No label is derived from model output.

**Scored eligibility includes a supported annotated translation and a known,
approved family label.** The baseline is the longest such eligible CDS within
the common residue cap, with ascending stable unit ID breaking length ties.
A longer translation with an unreviewed family is explicitly withheld from the
scored baseline, rather than silently selecting an easier protein. Source CDS,
parsed-translation and model/label-eligible denominators are separately retained.

The panel freezes every original CDS ID, source count, usable units, baseline ID,
source/sequence hashes, reviewed groups, full known family set and missing-family
unit IDs. Before any method score, it removes references sharing the parent,
**any** documented parent-component sequence, or any reviewed component group.
This includes hashes of existing but withheld source translations. All conditions
and models use exactly the remaining reference IDs for that parent.

For each independent method space:

1. Baseline uses one longest eligible CDS score row.
2. Decomposed condition takes the maximum score for each reference across all
   eligible source CDSs.
3. References are unique, ties use ascending reference ID, and only the final
   **five total parent hits** are returned. With fewer available references, return
   every available reference and flag the panel as descriptive/budget-saturated.

Every required CDS/reference score must exist and be finite. A failed extra CDS
abstains the decomposed parent while preserving a valid baseline result; missing
references cannot favorably shrink a ranking. No eligible source CDS or no valid
independent positive/competing family produces abstention/null, never success.

## Endpoints and attribution

Reports retain distinct expected/recovered family counts, family recall@5,
all-family recovery, parent coverage and unit exclusions. Duplicate CDS families
contribute one family to the denominator. `full` uses all independently reviewed
documented family expectations and blocks when unknown family annotations prevent
a complete definition. `matched` freezes the subset represented by eligible CDSs
and independently supported references; excluded families remain explicit.
These subsets never change after seeing model failures or scores.

The same longest-CDS query/reference scores are reported as a matched per-unit
control in both conditions. This is not extra independent evidence: improvement
from representing additional CDSs is **coverage**, not proof of improved per-unit
retrieval accuracy or biological-function prediction. Parents sharing component
sequences or reviewed groups form conservative connected independence clusters.
Paired all-CDS minus longest-CDS intervals use those clusters, the predeclared
seed and replicate count, and the minimum of eight independent clusters.
Parents returning the entire eligible catalogue are excluded from gain intervals;
insufficient clusters leave interval bounds null.

## Resources and artifacts

Protocol and parent-panel JSON freezes precede scoring; actual model revision,
weight digest and configuration are saved before the first forward. Models run
sequentially, using the existing adapters. Each run measures cold loading,
reference encoding, baseline encoding and additional-CDS encoding separately.
Exact sequence cache hits, successful/attempted calls and processed residues are
recorded. The all-CDS run reuses the measured baseline/reference prefix; it is not
an independent cold run. Baseline-only total wall/peak costs remain null where
not independently observed, and CPU/GPU usage seconds remain null if unavailable.

`report.json` contains conditions, abstentions, per-unit controls and actual
resource events. `cost-quality-records.json` contains a `records` array using
`cost-quality-record-1`, with full panel/protocol/source/config identity, hardware,
cache context and actual cost/quality fields for E4. An undeclared `evidence_kind`
defaults to `unclassified`; fictional tests are `fictional_software_fixture` and
cannot be promoted to admitted research measurements or scientific plots.

Low-level APIs `freeze_multi_cds_panel`, `score_multi_cds`, and
`produce_multi_cds_scores` are reusable with already admitted inputs. The public
execution API performs T02 gating; low-level numerical tests do not approve data.
Tests cover multi-CDS max merging, duplicate families, longest/tie choices,
missing labels/CDSs, empty parents, all-component leakage, strict score budgets,
incomplete/failing inputs, frozen hashes, stage ordering and source-admission
integration on generated, explicitly fictional software fixtures. No real
pretrained result is claimed.
