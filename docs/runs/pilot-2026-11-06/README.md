# Verified-subset pilot execution

One exact source was admitted for E1; nine remain HOLD. The
[release ledger](../../../data/manifests/pilot_execution_2026_11_06/release.json)
accounts for all ten individually. Researcher approval remains recorded in
merged PR #5; technical evidence determines executable scope.

| Attempt | Report status | Scope |
| --- | --- | --- |
| [E1](E1/report.json) | `ok` | 48/48 substantive features and qualifiers; 3/3 validated source translations; zero inventions; pBR322 only |
| [E2](E2/report.json) | `blocked_insufficient_parent_panel` | Three parsed CDSs; zero reviewed family-eligible units/references; no model calls |
| [E3](E3/report.json) | `blocked_insufficient_panel` | No positive/competing-family panel; no retrieval/inference scores |
| [E4](E4/report.json) | `blocked_no_measurements` | No scientific cost/quality measurement or plot |
| [Software suite](baseline-tests.txt) | 187 tests; 182 passed; five optional skips | Software fixtures only |
| [Model preflight](environment-global/report.json) | Cached ESM-2 loaded offline | No forward pass; ESM-C not validated |

E1's `scientific_scores_produced` flag means source-annotation agreement only,
as its `evidence_class` states. It does not mean family retrieval, biological
function prediction or pretrained inference. Two source bookkeeping features
are outside the 48-feature observable denominator. Zero-denominator abstention/
uncertainty rates remain null. One source cannot support a multi-archetype claim.

Sequence-free reports/freezes were copied byte-for-byte from ignored local
results. [Artifact digests](artifact_checksums.json) preserve those bytes;
[attempts.json](attempts.json) documents commands, outcomes, source/access blockers
and environment limitations. Reports identify main `4352009703dcd1825e181d203edf6088c882182e`
and a dirty working tree from new data/configuration. Existing parser/evaluation
modules and protocols were unchanged; source hashes are recorded. Original
qualifier sidecars, normalized objects and protein exports remain ignored locally.

## Reproduce

Create an isolated environment with `requirements-parsing.txt`, then run:

```powershell
.venv/Scripts/python.exe scripts/verify_pilot_release.py
.venv/Scripts/python.exe scripts/run_evaluation_extension.py --protocol config/extensions/E1.json --inputs config/pilot/inputs.verified_e1_2026_11_06.json --output results/my-new-pilot/E1
```

Use a new output directory for each attempt. E1 should exit 0. Changing E1 to E2
or E3 in protocol/output blocks with these inputs. No family labels or references
were fabricated to force inference. E4 consumes E3's empty measurement file and
blocks. Original commands are in the attempt ledger.

See [private-copy rights](../../../data/approved_public/README.md),
[independence review](../../../data/manifests/pilot_execution_2026_11_06/independence_review.json)
and [remaining decisions](../../coordination/DECISIONS_NEEDED.md).
PR #1 stays read-only; its artificial-query results/Drive bundle are not source
admissions or outcomes of this pilot.
