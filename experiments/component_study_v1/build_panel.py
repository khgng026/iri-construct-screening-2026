"""Label units, select references and freeze the panel — before any retrieval scoring.

Reads a toolkit bundle (units.tsv, proteins.fasta, normalized.json) and
families.v1.json, fetches the declared anchors and reviewed family members from
UniProt (cached privately) and writes

* a private panel with reference sequences (ignored storage), and
* a source-free freeze (IDs, hashes, labels, exclusions) for Git.

    python experiments/component_study_v1/build_panel.py

Follows SPEC.md; nothing here looks at a retrieval score.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
SEED = '20261117'
IDENTITY, COVERAGE, PER_FAMILY, MAX_RESIDUES = 0.90, 0.80, 5, 1022
CANONICAL = set('ACDEFGHIKLMNPQRSTVWY')


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def file_sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fetch(url: str, cache: Path) -> bytes:
    if not cache.exists():
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=120) as response:
                    cache.write_bytes(response.read())
                time.sleep(0.3)
                break
            except Exception:
                time.sleep(3 * (attempt + 1))
        else:
            raise SystemExit('UniProt request failed after bounded retries')
    return cache.read_bytes()


def evidence(entry: dict) -> dict:
    return {
        'accession': entry['primaryAccession'], 'reviewed': 'Swiss-Prot' in entry.get('entryType', ''),
        'pfam': sorted(x['id'] for x in entry.get('uniProtKBCrossReferences', []) if x['database'] == 'Pfam'),
        'similarity': [t['value'] for c in entry.get('comments', []) if c['commentType'] == 'SIMILARITY' for t in c['texts']],
        'name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
        'fragment': 'Fragment' in str(entry.get('proteinDescription', {}).get('flag', '')), 'sequence': entry['sequence']['value'],
        'entry_version': entry.get('entryAudit', {}).get('entryVersion'),
        'sequence_version': entry.get('entryAudit', {}).get('sequenceVersion')}


def satisfies(item: dict, rule: dict) -> bool:
    return (item['reviewed'] and set(rule.get('pfam_all', [])) <= set(item['pfam'])
            and (not rule.get('similarity_contains') or any(rule['similarity_contains'] in s for s in item['similarity']))
            and (not rule.get('recommended_name') or item['name'] == rule['recommended_name']))


class Aligner:
    """Local BLOSUM62 alignment statistics with the project's historical gap settings."""
    def __init__(self):
        from Bio.Align import PairwiseAligner, substitution_matrices
        self.aligner = PairwiseAligner(mode='local')
        self.aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
        self.aligner.open_gap_score, self.aligner.extend_gap_score = -10.0, -0.5
        self.self_scores, self.calls = {}, 0

    def self_score(self, s: str) -> float:
        if s not in self.self_scores:
            self.self_scores[s] = self.aligner.score(s, s)
        return self.self_scores[s]

    def stats(self, a: str, b: str) -> dict | None:
        """Identity over the aligned region and the fraction of each sequence it covers."""
        self.calls += 1
        if self.aligner.score(a, b) < 0.3 * min(self.self_score(a), self.self_score(b)):
            return None
        alignment = self.aligner.align(a, b)[0]
        same = columns = 0
        for (a0, a1), (b0, b1) in zip(*alignment.aligned):
            columns += a1 - a0
            same += sum(x == y for x, y in zip(a[a0:a1], b[b0:b1]))
        blocks_a, blocks_b = alignment.aligned
        span_a = blocks_a[-1][1] - blocks_a[0][0]
        span_b = blocks_b[-1][1] - blocks_b[0][0]
        columns_with_gaps = max(span_a, span_b)
        return {'identity': same / columns_with_gaps, 'coverage_a': span_a / len(a), 'coverage_b': span_b / len(b)}


def near_copy(stats: dict | None) -> bool:
    return bool(stats) and stats['identity'] >= IDENTITY and max(stats['coverage_a'], stats['coverage_b']) >= COVERAGE


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--bundle', type=Path, default=ROOT / 'results' / 'component_study_v1' / 'bundle')
    parser.add_argument('--cache', type=Path, default=ROOT / 'data' / 'raw' / 'component_study_v1' / 'uniprot')
    parser.add_argument('--private', type=Path, default=ROOT / 'results' / 'component_study_v1' / 'panel.private.json')
    parser.add_argument('--freeze', type=Path, default=ROOT / 'docs' / 'results' / 'component-study-v1' / 'panel.freeze.json')
    args = parser.parse_args(argv)
    args.cache.mkdir(parents=True, exist_ok=True)
    table = json.loads((HERE / 'families.v1.json').read_text())
    families = table['families']
    aligner = Aligner()

    # ---- units from the toolkit bundle (the research bridge: plain TSV/FASTA/JSON) ----
    with (args.bundle / 'units.tsv').open(newline='') as handle:
        rows = [r for r in csv.DictReader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)]
    with (args.bundle / 'features.tsv').open(newline='') as handle:
        names = {r['unit_id']: r for r in csv.DictReader(handle, delimiter='\t', quoting=csv.QUOTE_NONE)}
    sequences, current = {}, None
    for line in (args.bundle / 'proteins.fasta').read_text().splitlines():
        if line.startswith('>'):
            current = line[1:].split()[0]
            sequences[current] = ''
        else:
            sequences[current] += line
    normalized = json.loads((args.bundle / 'normalized.json').read_text())
    assert all(sha(sequences[r['unit_id']]) == r['protein_sha256'] for r in rows if r['extraction_state'] == 'exported')

    # ---- anchors: fetched individually and checked against the declared rule ----
    anchors = {}
    for family, rule in families.items():
        for accession in rule['anchors']:
            entry = json.loads(fetch(f'https://rest.uniprot.org/uniprotkb/{accession}.json', args.cache / f'{accession}.json'))
            item = evidence(entry)
            if not satisfies(item, rule):
                raise SystemExit(f'Declared anchor {accession} does not satisfy the rule of {family}')
            anchors[accession] = {**item, 'family': family}

    # ---- label every distinct exported protein by near-identity to an anchor ----
    labels = {}
    for sequence in sorted(set(sequences.values())):
        contained, fragments, best = {}, {}, {}
        for accession, anchor in anchors.items():
            stats = aligner.stats(sequence, anchor['sequence'])
            if not stats or stats['identity'] < IDENTITY:
                continue
            target = contained if stats['coverage_b'] >= COVERAGE else fragments
            if anchor['family'] not in target or stats['identity'] > target[anchor['family']]['identity']:
                target[anchor['family']] = {'anchor': accession, 'identity': round(stats['identity'], 4),
                                            'unit_coverage': round(stats['coverage_a'], 4),
                                            'anchor_coverage': round(stats['coverage_b'], 4)}
        if len(contained) > 1:
            label_class = 'multi_family'
        elif contained:
            label_class = 'full' if next(iter(contained.values()))['unit_coverage'] >= COVERAGE else 'embedded'
        else:
            label_class = 'fragment' if fragments else 'unknown'
        labels[sha(sequence)] = {'label_class': label_class, 'families': sorted(contained),
                                 'evidence': contained or fragments}

    # ---- references: reviewed members, seeded order, no lineage copies ----
    panel, family_report = {}, {}
    for family, rule in families.items():
        clauses = ['reviewed:true'] + [f'xref:pfam-{p}' for p in rule.get('pfam_all', [])]
        if rule.get('similarity_contains'):
            clauses.append('cc_similarity:"' + rule['similarity_contains'].replace('"', '') + '"')
        if rule.get('recommended_name'):
            clauses.append('protein_name:"' + rule['recommended_name'] + '"')
        query = ' AND '.join(f'({c})' for c in clauses)
        url = 'https://rest.uniprot.org/uniprotkb/search?' + urllib.parse.urlencode(
            {'query': query, 'format': 'json', 'size': 500, 'sort': 'accession asc'})
        raw = fetch(url, args.cache / f'members.{family}.json')
        members = [evidence(e) for e in json.loads(raw)['results']]
        lineage = [a['sequence'] for a in anchors.values() if a['family'] == family]
        lineage += sorted({s for s in set(sequences.values()) if family in labels[sha(s)]['families']})
        reasons = {'not_rule_member': 0, 'anchor': 0, 'fragment_or_length_or_residues': 0,
                   'lineage_copy_of_anchor_or_unit': 0, 'near_copy_of_chosen_reference': 0, 'not_needed': 0}
        chosen = []
        for item in sorted(members, key=lambda m: sha(f"{SEED}:{m['accession']}")):
            s = item['sequence']
            if len(chosen) == PER_FAMILY:
                reasons['not_needed'] += 1
            elif not satisfies(item, rule):
                reasons['not_rule_member'] += 1
            elif item['accession'] in anchors:
                reasons['anchor'] += 1
            elif item['fragment'] or not 50 <= len(s) <= MAX_RESIDUES or set(s) - CANONICAL:
                reasons['fragment_or_length_or_residues'] += 1
            elif any(near_copy(aligner.stats(s, other)) for other in lineage):
                reasons['lineage_copy_of_anchor_or_unit'] += 1
            elif any(near_copy(aligner.stats(s, c['sequence'])) for c in chosen):
                reasons['near_copy_of_chosen_reference'] += 1
            else:
                chosen.append(item)
        family_report[family] = {
            'rule': {k: v for k, v in rule.items() if k != 'anchors'}, 'query': query,
            'members_response_sha256': hashlib.sha256(raw).hexdigest(), 'members_returned': len(members),
            'anchors': [{'accession': a, 'sequence_sha256': sha(anchors[a]['sequence']),
                         'entry_version': anchors[a]['entry_version']} for a in rule['anchors']],
            'candidate_dispositions': reasons, 'scorable': len(chosen) >= 2,
            'references': [{'reference_id': 'uniprot-' + c['accession'], 'length_aa': len(c['sequence']),
                            'sequence_sha256': sha(c['sequence']), 'entry_version': c['entry_version'],
                            'sequence_version': c['sequence_version']} for c in chosen]}
        if len(chosen) >= 2:
            for c in chosen:
                panel['uniprot-' + c['accession']] = {'family': family, 'sequence': c['sequence']}
    scorable = {f for f, r in family_report.items() if r['scorable']}

    # ---- units and parents ----
    units, by_parent = [], {}
    for r in rows:
        if r['molecule_type'] != 'protein':
            continue
        exported = r['extraction_state'] == 'exported'
        label = labels[r['protein_sha256']] if exported else {'label_class': None, 'families': [], 'evidence': {}}
        eligible = exported and int(r['length_aa']) <= MAX_RESIDUES
        supported_families = sorted(set(label['families']) & scorable)
        unit = {'unit_id': r['unit_id'], 'construct_id': r['construct_id'], 'extraction_state': r['extraction_state'],
                'reason_code': r['reason_code'], 'length_aa': int(r['length_aa']) if exported else None,
                'sequence_sha256': r['protein_sha256'] or None, 'eligible': eligible,
                'label_class': label['label_class'], 'families': label['families'],
                'label_evidence': label['evidence'], 'supported': eligible and bool(supported_families),
                'supported_families': supported_families if eligible else []}
        units.append(unit)
        by_parent.setdefault(r['construct_id'], []).append(unit)
    parents = []
    for construct in normalized['constructs']:
        cid = construct['construct_id']
        own = by_parent.get(cid, [])
        eligible = sorted((u for u in own if u['eligible']), key=lambda u: (-u['length_aa'], u['unit_id']))
        supported = [u for u in eligible if u['supported']]
        covered = set()
        exported_ids = {u['unit_id'] for u in own if u['extraction_state'] == 'exported'}
        for feature in construct['features']:
            if feature['unit_id'] in exported_ids:
                for segment in feature['segments']:
                    covered.update(range(segment['start0'], segment['end0']))
        parents.append({
            'construct_id': cid, 'sequence_length_nt': construct['sequence_length_nt'], 'topology': construct['topology'],
            'source_features': len(construct['features']),
            'noncoding_features': sum(f['source_feature_type'] != 'CDS' for f in construct['features']),
            'source_cds': len(own), 'exported_cds': len(exported_ids),
            'withheld_by_reason': {code: sum(u['reason_code'] == code for u in own if u['extraction_state'] == 'withheld')
                                   for code in sorted({u['reason_code'] for u in own if u['extraction_state'] == 'withheld'})},
            'over_model_limit': sum(u['extraction_state'] == 'exported' and not u['eligible'] for u in own),
            'eligible_unit_ids': [u['unit_id'] for u in eligible],
            'round_robin_order': [u['unit_id'] for u in eligible],
            'baseline_unit_id': eligible[0]['unit_id'] if eligible else None,
            'baseline_supported_unit_id': supported[0]['unit_id'] if supported else None,
            'supported_unit_ids': [u['unit_id'] for u in supported],
            'exported_cds_nt_fraction': round(len(covered) / construct['sequence_length_nt'], 4),
            'analysis_set': ('multi_unit' if len(eligible) >= 2 else 'single_unit') if supported else 'no_supported_component'})

    freeze = {
        'schema_version': 'component-study-panel-1', 'spec': 'experiments/component_study_v1/SPEC.md',
        'spec_sha256': file_sha(HERE / 'SPEC.md'), 'families_sha256': file_sha(HERE / 'families.v1.json'),
        'bundle': {name: file_sha(args.bundle / name) for name in ('units.tsv', 'features.tsv', 'proteins.fasta', 'normalized.json')},
        'thresholds': {'identity': IDENTITY, 'coverage': COVERAGE, 'references_per_family': PER_FAMILY,
                       'max_residues': MAX_RESIDUES, 'seed': SEED, 'budget_k': 5},
        'families': family_report, 'reference_ids': sorted(panel),
        'reference_family': {r: panel[r]['family'] for r in sorted(panel)},
        'units': units, 'parents': parents,
        'counts': {
            'parents': len(parents), 'cds_units': len(units),
            'exported_units': sum(u['extraction_state'] == 'exported' for u in units),
            'eligible_units': sum(u['eligible'] for u in units),
            'distinct_eligible_sequences': len({u['sequence_sha256'] for u in units if u['eligible']}),
            'label_class': {c: sum(u['label_class'] == c for u in units) for c in ('full', 'embedded', 'multi_family', 'fragment', 'unknown')},
            'supported_units': sum(u['supported'] for u in units),
            'distinct_supported_sequences': len({u['sequence_sha256'] for u in units if u['supported']}),
            'scorable_families': len(scorable), 'references': len(panel),
            'analysis_sets': {s: sum(p['analysis_set'] == s for p in parents) for s in ('multi_unit', 'single_unit', 'no_supported_component')},
            'alignment_calls_for_labels_and_exclusions': aligner.calls},
        'status': 'frozen_before_scoring'}
    freeze['panel_sha256'] = sha(json.dumps(freeze, sort_keys=True, separators=(',', ':')))
    args.freeze.parent.mkdir(parents=True, exist_ok=True)
    args.freeze.write_text(json.dumps(freeze, indent=1, sort_keys=True) + '\n')
    args.private.parent.mkdir(parents=True, exist_ok=True)
    args.private.write_text(json.dumps({'panel_sha256': freeze['panel_sha256'], 'references': panel}, sort_keys=True) + '\n')
    # Private review aid: source names beside the labels (qualifier text stays out of Git).
    with (args.private.parent / 'labels.private.tsv').open('w') as handle:
        for u in units:
            n = names[u['unit_id']]
            handle.write('\t'.join([u['unit_id'], str(u['length_aa']), str(u['label_class']), ','.join(u['families']),
                                    n['gene'], n['product'], n['display_name']]) + '\n')
    print(json.dumps(freeze['counts'], indent=1))
    print({f: len(r['references']) for f, r in family_report.items()})


if __name__ == '__main__':
    main()
