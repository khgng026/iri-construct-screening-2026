"""Versioned E4 and exploratory E2 attribution; no new inference or primary criteria."""
import argparse, hashlib, json, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.runtime import InputError, write_json, digest_file
E3_SHA = 'f13e363599c4173ea7973301a696600d7dcf15c742e3d572e8c2b82ea87f1a53'

def require(value, reason):
    if not value:
        raise InputError(reason)

def attributed_hits(scores, families, reference_families, k=5, balanced=False):
    ordered = {q: sorted(row, key=lambda r: (-row[r], r)) for q, row in scores.items()}
    if balanced:
        hits = []
        origins = {}
        for rank in range(max(map(len, ordered.values()))):
            for q in sorted(ordered):
                if rank < len(ordered[q]) and ordered[q][rank] not in hits:
                    hits.append(ordered[q][rank])
                    origins[ordered[q][rank]] = q
                if len(hits) == k:
                    break
            if len(hits) == k:
                break
    else:
        hits = sorted(next(iter(scores.values())), key=lambda r: (-max((row[r] for row in scores.values())), r))[:k]
    rows = []
    for r in hits:
        winner = min(scores, key=lambda q: (-scores[q][r], q)) if not balanced else origins[r]
        rows.append({'reference_id': r, 'winning_component': winner, 'reference_family': reference_families[r], 'component_family': families[winner], 'component_attributed_correct': reference_families[r] == families[winner]})
    return rows

def analyze_e4(out):
    e3path = ROOT / 'docs/results/cleared-2026-11-11/E3.json'
    require(digest_file(e3path) == E3_SHA, 'Changed approved E3 evidence')
    protocol = ROOT / 'config/week2/E4.proposed.json'
    require(digest_file(protocol) == '41587c43d414a53cbe9d39dcba0df4fcfd7361d5719e5b1eec732557b5d6655e', 'Changed E4 protocol')
    e3 = json.loads(e3path.read_text())
    source = json.loads((ROOT / 'data/manifests/week2/source_review.json').read_text())
    archetypes = {r['source_id']: r.get('proposed_archetype') for r in source['records']}
    derived = {'schema_version': 'derived-archetype-report-1', 'original_E3_sha256': E3_SHA, 'score_matrices_unchanged': True, 'rows': []}
    rows = []
    for name, m in e3['methods'].items():
        for r in m['result']['rows']:
            derived['rows'].append({'method': name, 'query_id': r['query_id'], 'archetype': archetypes[r['construct_id']], 'family_recall_at_5': r['documented_family_recall_at_k']})
        c = m['costs']
        rows.append({'method': name, 'family_recall_at_5': m['result']['overall']['documented_family_recall_at_k'], 'coverage': m['result']['coverage'], 'wrong_family_hits': m['result']['false_family_hit_count'], 'returned_hits': 20, 'work_seconds': c['wall_seconds'], 'load_seconds': c.get('cold_load_seconds'), 'sampled_work_RSS_bytes': c['sampled_peak_process_rss_bytes'], 'calls': c.get('model_calls', c.get('alignment_calls')), 'processed_residues': c.get('processed_residues')})
    write_json(out / 'archetypes.v1.json', derived)
    manifest = {'clearance': 'https://github.com/khgng026/iri-construct-screening-2026/pull/11#issuecomment-6815404726', 'decision': 'CLEAR_FOR_E4_DESCRIPTIVE_MATCHED_E3_ONLY', 'E3_sha256': E3_SHA, 'E4_protocol_sha256': digest_file(protocol), 'panel_sha256': e3['panel']['panel_sha256'], 'methods': ['alignment_raw', 'esm2'], 'excluded': ['E2', 'fixture measurements', 'ESM-C'], 'confidence_interval': None}
    write_json(out / 'E4.execution.json', manifest)
    result = {'execution': manifest, 'rows': rows, 'inference_vs_alignment_work_ratio': rows[1]['work_seconds'] / rows[0]['work_seconds'], 'conclusion': 'Equal observed family-recall endpoints in four queries; ESM-2 inference/pooling used more method-work time in this run, plus separate initialization. No general winner/equivalence claim.', 'limitations': ['Different implementation initialization is not symmetrically measured; alignment initialization remains null.', 'RSS is absolute process memory, including prior work/imports; not incremental model memory.', 'OS cache uncontrolled; two dependence clusters; no significance or hardware-independent efficiency.']}
    write_json(out / 'E4.json', result)

def analyze_e2(out):
    proposal = json.loads((ROOT / 'data/manifests/week2/references/reference_proposal.json').read_text())
    rf = {r['reference_id']: r['proposed_family'] for r in proposal['references']}
    qs = [q for q in proposal['query_proposals'] if q['parent_id'] == 'ncbi-AF298787.1']
    families = {q['query_id']: q['proposed_family'] for q in qs}
    scores = json.loads((ROOT / 'results/cleared-2026-11-11/E2/run01/unit-scores.json').read_text())
    baseline = next((q for q in scores if families[q] == 'omp_decarboxylase_homologs'))
    design = {'status': 'PROSPECTIVE_PRIMARY_EVALUATION_PENDING_CLEARANCE', 'method': 'round-robin per-component ranked references, deduplicate, total K=5; stable unit/reference ordering', 'primary_endpoint': 'Unchanged historical E2; any future replacement requires separate clearance', 'scores_sha256': digest_file(ROOT / 'results/cleared-2026-11-11/E2/run01/unit-scores.json')}
    write_json(out / 'balanced_merge.prospective.json', design)
    diagnostic = {'evidence_kind': 'exploratory_cached_score_attribution_not_replacement_primary', 'conditions': {}}
    for name, subset in [('longest', {baseline: scores[baseline]}), ('all', scores)]:
        hits = attributed_hits(subset, families, rf)
        diagnostic['conditions'][name] = {'hits': hits, 'parent_families_returned': sorted({h['reference_family'] for h in hits} & set(families.values())), 'families_recovered_through_matching_component': sorted({h['reference_family'] for h in hits if h['component_attributed_correct']})}
    diagnostic['exploratory_balanced_hits'] = attributed_hits(scores, families, rf, balanced=True)
    diagnostic['score_scale_by_component'] = {q: {'minimum': min(row.values()), 'maximum': max(row.values())} for q, row in scores.items()}
    diagnostic['reference_composition'] = {f: sum((rf[r] == f for r in next(iter(scores.values())))) for f in set(rf.values())}
    write_json(out / 'E2.attribution.json', diagnostic)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage', choices=['E4', 'E2'], required=True)
    a = p.parse_args()
    out = ROOT / 'docs/results/analysis-2026-11-12'
    out.mkdir(parents=True, exist_ok=True)
    (analyze_e4 if a.stage == 'E4' else analyze_e2)(out)
if __name__ == '__main__':
    main()
