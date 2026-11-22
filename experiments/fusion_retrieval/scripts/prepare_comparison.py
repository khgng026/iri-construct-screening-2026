"""Freeze a benign, bounded adaptation of the published benchmark before scoring."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
UP = ROOT / 'sources/nishit-working'
sys.path.insert(0, str(UP / 'src'))
from pilot.reference import iter_fasta
from budget_compare import atomic_json, digest

PAYLOADS = ['B2GDS5', 'P11542', 'Q30U33', 'Q6F1K1', 'Q6N9W0', 'Q92QI7']
EXCLUDED_FILLERS = {'P20531': 'viral host-range/anti-apoptosis factor',
                    'Q02HI0': 'DING protein; conservative unresolved-scope exclusion',
                    'Q9FWW5': 'uncertain function'}


def csvrows(path):
    return list(csv.DictReader(path.open(encoding='utf-8-sig', newline='')))


def fasta(text):
    out, acc, parts = {}, None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if acc:
                out[acc] = ''.join(parts)
            acc, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if acc:
        out[acc] = ''.join(parts)
    return out


def main():
    if (ROOT / 'configs/frozen.json').exists():
        raise SystemExit('Freeze already exists; preserve it rather than silently overwrite')
    seqs = {}
    # Preserve first observed sequence per accession, and reject discrepancies for
    # any reference actually selected below.
    conflicts = set()
    paths = [UP / x for x in ['results/homologs_all173/reference_twilight.fasta',
             'data/reference/reference_twilight.fasta', 'data/reference_handcheck/reference_twilight.fasta',
             'data/reference/reference.fasta', 'data/reference/payloads.fasta']]
    origins = {}
    for path in paths:
        for acc, seq in fasta(path.read_text()).items():
            if acc in seqs and seqs[acc] != seq:
                conflicts.add(acc)
            else:
                seqs[acc] = seq
                origins.setdefault(acc, str(path.relative_to(ROOT)))
    manifest = csvrows(UP / 'data/synthetic_real/synthetic_manifest.csv')
    source_by_id = {x['construct_id']: x for x in manifest}
    queries = fasta((UP / 'benchmark/v1/queries.faa').read_text())
    items, excluded = [], []
    for row in csvrows(UP / 'benchmark/v1/items.csv'):
        if row['tier'] != 'synth_real' or row['payload_accession'] not in PAYLOADS:
            continue
        sid = row['construct_id'].split(':', 1)[1]
        source = source_by_id[sid]
        fillers = [a for a in source['filler_accessions'].split('|') if a]
        reason = [EXCLUDED_FILLERS[a] for a in fillers if a in EXCLUDED_FILLERS]
        seq = queries[row['construct_id']]
        if len(seq) > 2046:
            reason.append('above common 2046 residue cap; no silent model truncation')
        if reason:
            excluded.append({'id': row['construct_id'], 'reasons': reason})
            continue
        # Published sequence is the input; no rebuilding or sequence generation.
        if len(seq) != int(row['orf_len']):
            raise ValueError('Published sequence and item length differ')
        items.append({**row, 'sequence': seq, 'sequence_sha256': hashlib.sha256(seq.encode()).hexdigest(),
                      'source_manifest': source, 'filler_accessions': fillers})
    hs = csvrows(UP / 'benchmark/v1/homologs.csv')
    selected, ref_reason = set(PAYLOADS), {a: ['payload_self_control'] for a in PAYLOADS}
    for x in items:
        for a in x['filler_accessions']:
            selected.add(a)
            ref_reason.setdefault(a, []).append('documented_existing_filler')
    # Fixed accession hash subsampling within identity bands, no retrieval outcome.
    bands = [(0, .25), (.25, .35), (.35, .5), (.5, .7), (.7, .9), (.9, 1.00001)]
    for acc in PAYLOADS:
        for lo, hi in bands:
            candidates = {h['homolog_accession'] for h in hs if h['payload_accession'] == acc
                          and lo < float(h['identity']) <= hi and h['homolog_accession'] in seqs
                          and 50 <= len(seqs[h['homolog_accession']]) <= 1020
                          and h['homolog_accession'] not in conflicts}
            pick = sorted(candidates, key=lambda a: hashlib.sha256(('iri-v1:' + a).encode()).hexdigest())[:8]
            for a in pick:
                selected.add(a)
                ref_reason.setdefault(a, []).append(f'published_homolog:{acc}:{lo}-{hi}')
    # Reuse exact reviewed October source bytes as unrelated catalogue competitors.
    old = ROOT.parents[1] / 'iri-construct-screening-2026/data/manifests/t05_benign_catalogue.json'
    old_records = json.loads(old.read_text())['records']
    archive = ROOT.parents[1] / 'iri-reviewed-core-inputs-2026-10-31.zip'
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
        for r in old_records:
            suffix = r['sequence_ref']['path']
            hits = [n for n in names if n.endswith(suffix)]
            if len(hits) != 1:
                raise ValueError('Missing/ambiguous old reference in reviewed archive: ' + suffix)
            ss = list(fasta(z.read(hits[0]).decode()).values())[0]
            if hashlib.sha256(ss.encode()).hexdigest() != r['sequence_sha256']:
                raise ValueError('Old sequence hash mismatch')
            a = r['accession']
            if a in seqs and seqs[a] != ss:
                raise ValueError('Cross-source selected reference mismatch')
            if len(ss) <= 1020:
                seqs[a], origins[a] = ss, 'reviewed-core-inputs-2026-10-31.zip:' + hits[0]
                selected.add(a)
                ref_reason.setdefault(a, []).append('reviewed_october_competitor')
    absent_fillers = sorted(selected - set(seqs))
    if any(a in PAYLOADS for a in absent_fillers) or selected & conflicts:
        raise ValueError({'conflicts': sorted(selected & conflicts), 'missing_payload': sorted(selected - set(seqs))})
    # Some full filler sources were not shipped. Keep every query unchanged and
    # explicitly record missing filler references; this is a single-target task.
    selected -= set(absent_fillers)
    refs = [{'accession': a, 'sequence': seqs[a], 'sequence_sha256': hashlib.sha256(seqs[a].encode()).hexdigest(),
             'origin': origins[a], 'selection_reasons': sorted(set(ref_reason[a]))} for a in sorted(selected)]
    ordered_payloads = sorted(PAYLOADS, key=lambda a: hashlib.sha256(('iri-split-v1:' + a).encode()).hexdigest())
    dev = set(ordered_payloads[:3])
    for item in items:
        item['split'] = 'development' if item['payload_accession'] in dev else 'evaluation'
    config = {'schema': 'iri-nishit-adaptation-1', 'created_before_scoring': True,
              'claim_scope': 'Published alignment-derived target retrieval; not function, screening or AI-designed construct validation',
              'canonical_replication': False, 'canonical_gap': 'Exact gate_all173 FASTA and canonical maskfix item-rank files absent; release assets empty',
              'source_commit': '5dfef4734f07b981fa6305383fee170d6af28a6d',
              'model': 'facebook/esm2_t33_650M_UR50D', 'model_revision': '08e4846e537177426273712802403f7ba8261b6c',
              'max_length': 2048, 'max_query_residues': 2046, 'device': 'cpu', 'batch_size': 1, 'threads': 4,
              'model_pooling': 'upstream attention-mask mean includes CLS and EOS',
              'thresholds': [1.0, .7, .25], 'primary_budget': 10, 'secondary_budget': 5,
              'window': 300, 'step': 150, 'min_segment': 50, 'max_evalue': .001,
              'alignment': {'matrix': 'BLOSUM62', 'open_gap_score': -11, 'extend_gap_score': -1,
                            'backend': 'Biopython PairwiseAligner exact local; Parasail build unavailable'},
              'segmentation': 'Unmodified Nishit best_hit_segments; its span locator uses upstream Bio blastp defaults (-12/-1)',
              'arms': ['whole', 'windows', 'bhs', 'oracle'],
              'aggregations': ['max_score', 'round_robin', 'per_segment_diagnostic'],
              'improvement_hypothesis': 'Coordinate-ordered round-robin reduces family crowding at a fixed total review budget without target labels or cross-segment score calibration',
              'decision_rule': 'No parameter tuning on either split; report paired evaluation cells by f stratum and payload; three evaluation payloads preclude generalization',
              'split': {'development': sorted(dev), 'evaluation': sorted(set(PAYLOADS)-dev)},
              'references': refs, 'items': items, 'exclusions': excluded,
              'documented_fillers_without_shipped_reference': absent_fillers,
              'homologs': [h for h in hs if h['payload_accession'] in PAYLOADS and h['homolog_accession'] in selected],
              'source_hashes': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths + [UP/'benchmark/v1/items.csv', UP/'benchmark/v1/queries.faa', UP/'benchmark/v1/homologs.csv', old, archive] if p.is_relative_to(ROOT)},
              'negative_label_limit': 'Non-target-labelled returned hits are not biological false positives; other real components may be present',
              'complete_recovery_limit': 'Only one designated target per item; cannot claim recovery of all biological components',
              'abstention': 'No fitted similarity threshold; only missing/nonfinite computation abstains; absent-target references remain misses'}
    config['freeze_sha256'] = digest(config)
    atomic_json(ROOT/'configs/frozen.json', config)
    atomic_json(ROOT/'evidence/cohort_summary.json', {'n_items':len(items), 'n_references':len(refs), 'split':config['split'], 'exclusions':excluded, 'freeze_sha256':config['freeze_sha256']})
    print(json.dumps({'n_items':len(items), 'n_references':len(refs), 'split':config['split'], 'excluded':len(excluded), 'sha256':config['freeze_sha256']}))


if __name__ == '__main__':
    main()
