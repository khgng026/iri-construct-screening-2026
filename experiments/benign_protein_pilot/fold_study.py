"""Evaluate curated natural protein families with same-fold alternatives."""
from __future__ import annotations

import argparse
import json
import random
import re
import time
from collections import Counter, defaultdict
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from statistics import mean
from urllib.error import HTTPError, URLError
from zoneinfo import ZoneInfo

from .reaction_study import _encode_150m, _prior_sequence_hashes
from .related_study import LOCAL_EXTEND, LOCAL_OPEN, _encode, _residue_local_score
from .structural_study import (CLASSIFICATION, DESCRIPTIONS, RELEASE, TAXA, _domain_sequence, _global_identity,
                               _profile_scores, _top_label, family_folds, metadata_candidates, superfamily)
from .study import ALPHABET, MIN_COVERAGE, _alignment_comparison, sha, write_json

SCHEMA = 'curated-homologous-domains-1'
FREEZE = Path(__file__).with_name('fold_panel.freeze.json')
SEED = 20261123
PRIMARY_IDENTITY = 0.25
MIN_LENGTH = 55
MAX_PER_FINE_FAMILY = 6

# Curated natural domains involved in trafficking, folding, metabolism, translation and chromatin.
FOLD_GROUPS = {
    'a.118': {1: (1, 2, 3, 4, 9, 10, 14, 19, 20, 21, 22, 28, 29), 8: (1, 7, 9, 10), 9: (1, 2, 3, 4)},
    'b.34': {5: (1, 2, 3, 4, 5, 6, 7), 9: (1, 2, 3, 4), 13: (1, 2, 3)},
    'b.43': {3: (1, 2, 3, 4, 5, 6), 4: (1, 2, 3), 5: (1, 2)},
    'd.15': {1: (1, 2, 3, 4, 7, 8), 3: (1, 2), 4: (1, 2)},
    'd.58': {1: (1, 2, 3, 4, 5, 6), 7: (1, 3, 4, 5), 16: (1, 2),
             18: (1, 2, 3, 5, 6, 10), 26: (1, 2, 3, 4, 5, 7)},
    'd.79': {1: (1, 2), 3: (1, 2)},
}
FAMILIES = {f'{fold}.{group}.{family}' for fold, groups in FOLD_GROUPS.items()
            for group, families in groups.items() for family in families}
COMMON_NAME_FILTERS = {
    'a.118.8.1': ('cyclophilin', 'fkbp', 'hop', 'mitochondria', 'peroxin', 'pre-mrna splicing', 'sec17',
                  'phosphatase 5', 'prolyl 4-hydroxylase', 'stip1', 'smg-7'),
    'd.15.1.1': ('ubiquitin$', 'sumo', 'rub1', 'nedd8$', 'dsk2', 'tubulin', 'rad23', 'elongin', 'hub1'),
    'd.58.7.1': ('cbp20', 'cold-inducible', 'translation initiation', 'nuclear ribonucleoprotein',
                  r'poly\(a\)-binding', 'splic', 'nucleolin', 'polypyrimidine', 'nucleoporin'),
}
CLASSICAL = ('raw_alignment', 'normalized_alignment', 'fine_profiles', 'superfamily_profiles')
MODELS = ('esm2_8m_mean', 'esm2_8m_local', 'esm2_150m_mean', 'esm2_150m_local')
MAIN_COMPARISONS = ('raw_alignment', 'fine_profiles', 'superfamily_profiles')
REFERENCE_POLICIES = ('three-way', 'leave-one-family-out')


def arch(family):
    return family.rsplit('.', 2)[0]


def _eligible(row):
    if row['family'] not in FAMILIES:
        return False
    filters = COMMON_NAME_FILTERS.get(row['family'])
    return filters is None or re.search('|'.join(filters), row['protein'], re.IGNORECASE) is not None


def _worklist(buckets):
    by_family = defaultdict(dict)
    for (family, taxon), rows in buckets.items():
        for row in rows:
            if _eligible(row):
                key = row['protein'], taxon
                by_family[family].setdefault(key, []).append(row)
    work = []
    for family, candidates in sorted(by_family.items()):
        name_counts = Counter()
        ranked = []
        for (name, taxon), rows in sorted(candidates.items()):
            ranked.append((name_counts[name], name, taxon, sorted(rows, key=lambda row: row['id'])[:2]))
            name_counts[name] += 1
        # Reserve a few source candidates for records with nonstandard residues or exact duplicates.
        work.extend(rows for _, _, _, rows in sorted(ranked)[:MAX_PER_FINE_FAMILY + 3])
    return work


def _download_domain(rows, cache):
    failed = []
    for row in rows:
        cached = cache / f"{row['id']}.json"
        if cached.exists():
            result = json.loads(cached.read_text())
            if result is not None:
                return result, failed
            failed.append((row['id'], 'length_or_alphabet'))
            continue
        for attempt in range(2):
            try:
                time.sleep(0.85)
                result = _domain_sequence(row, minimum_length=MIN_LENGTH)
                write_json(cached, result)
                if result is not None:
                    return result, failed
                failed.append((row['id'], 'length_or_alphabet'))
                break
            except HTTPError as error:
                if error.code == 429:
                    wait = error.headers.get('Retry-After', 'provider rate window')
                    raise RuntimeError(f'SCOPe requested a pause ({wait}). Successful entries are cached.') from error
                if error.code >= 500 and attempt == 0:
                    time.sleep(1)
                    continue
                failed.append((row['id'], f'http_{error.code}'))
                break
            except (URLError, TimeoutError):
                if attempt == 1:
                    failed.append((row['id'], 'unavailable'))
                time.sleep(1)
    return None, failed


def _admitted_folds(rows):
    subdivisions = family_folds(rows)
    labels = defaultdict(set)
    for family in subdivisions:
        labels[arch(family)].add(superfamily(family))
    return {family: fold for family, fold in subdivisions.items() if len(labels[arch(family)]) >= 2}


def acquire(folder, descriptions, classification):
    desc_source, class_source = descriptions.read_bytes(), classification.read_bytes()
    if (sha(desc_source), sha(class_source)) != (DESCRIPTIONS, CLASSIFICATION):
        raise ValueError('SCOPe source metadata does not match the pinned release')
    work = _worklist(metadata_candidates(desc_source.decode(), class_source.decode(), FAMILIES))
    seen = _prior_sequence_hashes() | {row['sha256'] for row in json.loads(
        Path(__file__).with_name('reaction_panel.freeze.json').read_text())['entries']}
    selections, rejected = [], []
    counts = Counter()
    cache = folder / 'domain-cache'
    cache.mkdir(parents=True, exist_ok=True)
    print(f'Collecting up to {len(work)} curated domain entries', flush=True)
    for index, candidates in enumerate(work, 1):
        result, failures = _download_domain(candidates, cache)
        rejected.extend(failures)
        if result is not None:
            row, sequence = result
            if row['sha256'] not in seen and counts[row['family']] < MAX_PER_FINE_FAMILY:
                counts[row['family']] += 1
                seen.add(row['sha256'])
                selections.append((row, sequence))
        if index % 30 == 0:
            print(f'Checked {index}/{len(work)} source entries', flush=True)
    splits = _admitted_folds([row for row, _ in selections])
    entries = [dict(row, split=splits[row['family']]) for row, _ in selections if row['family'] in splits]
    sequences = {row['id']: sequence for row, sequence in selections if row['family'] in splits}
    if len({arch(row['family']) for row in entries}) < 2:
        raise ValueError('Need at least two folds with multiple represented superfamilies')
    write_json(folder / 'panel.freeze.json', {
        'schema': SCHEMA, 'release': RELEASE, 'descriptions_sha256': DESCRIPTIONS, 'classification_sha256': CLASSIFICATION,
        'source': 'https://scop.berkeley.edu/astral/ver=2.08', 'terms': 'https://scop.berkeley.edu/about',
        'retrieved_on': datetime.now(ZoneInfo('America/New_York')).date().isoformat(), 'entries': entries,
        'rejected_sources': rejected, 'candidate_fine_families': sorted(counts),
    })
    write_json(folder / 'sequences.private.json', sequences)
    retained_families = {entry['family'] for entry in entries}
    print(f'Frozen {len(entries)} unique domains in {len(retained_families)} fine families')
    by_fold = Counter(arch(row['family']) for row in entries)
    print(dict(sorted(by_fold.items())))


def restore(folder):
    frozen = FREEZE.read_bytes()
    entries = json.loads(frozen)['entries']
    cache = folder / 'domain-cache'
    cache.mkdir(parents=True, exist_ok=True)
    sequences = {}
    for index, row in enumerate(entries, 1):
        result, failures = _download_domain([row], cache)
        if result is None:
            raise ValueError(f"{row['id']}: pinned ASTRAL snapshot unavailable ({failures})")
        observed, sequence = result
        if any(observed[key] != row[key] for key in ('id', 'family', 'length', 'sha256', 'source_sha256')):
            raise ValueError(f"{row['id']}: downloaded data differs from the frozen panel")
        sequences[row['id']] = sequence
        if index % 30 == 0:
            print(f'Verified {index}/{len(entries)} frozen source entries', flush=True)
    with (folder / 'panel.freeze.json').open('xb') as handle:
        handle.write(frozen)
    write_json(folder / 'sequences.private.json', sequences)
    _panel(folder)
    print(f'Verified all {len(entries)} frozen source entries')


def _panel(folder):
    frozen = (folder / 'panel.freeze.json').read_bytes()
    panel = json.loads(frozen)
    entries = {row['id']: row for row in panel['entries']}
    sequences = json.loads((folder / 'sequences.private.json').read_text())
    splits = _admitted_folds(entries.values())
    if (panel['schema'] != SCHEMA or not entries or entries.keys() != sequences.keys()
            or len(entries) != len(panel['entries']) or len({row['sha256'] for row in entries.values()}) != len(entries)
            or any(not _eligible(row) or row['taxon'] not in TAXA
                   or row['family'] not in splits or row['split'] != splits[row['family']]
                   or row['length'] != len(sequences[key]) or set(sequences[key]) - ALPHABET
                   or sha(sequences[key].encode('ascii')) != row['sha256'] for key, row in entries.items())):
        raise ValueError('Sequence contents or cross-validation groups differ from the frozen panel')
    return entries, sequences, sha(frozen)


def comparisons(entries, sequences, reference_policy):
    group_key = 'family' if reference_policy == 'leave-one-family-out' else 'split'
    for fold in sorted({arch(row['family']) for row in entries.values()}):
        fold_ids = [key for key, row in entries.items() if arch(row['family']) == fold]
        expected = {superfamily(entries[key]['family']) for key in fold_ids}
        for group in sorted({entries[key][group_key] for key in fold_ids}):
            query = {key: sequences[key] for key in fold_ids if entries[key][group_key] == group}
            refs = {key: sequences[key] for key in fold_ids if entries[key][group_key] != group}
            if {superfamily(entries[key]['family']) for key in refs} != expected:
                raise ValueError(f'{fold}: a split lost a reference superfamily')
            yield fold, query, refs


def _classical(entries, sequences, reference_policy):
    from Bio.Align import PairwiseAligner, substitution_matrices

    fine = {key: row['family'] for key, row in entries.items()}
    broad = {key: superfamily(family) for key, family in fine.items()}
    fine_labels = {family: superfamily(family) for family in set(fine.values())}
    broad_labels = {label: label for label in set(broad.values())}
    aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)
    rows, pairs, profile_count = {}, 0, 0
    for fold, query, refs in comparisons(entries, sequences, reference_policy):
        seq_scores, local_identity = _alignment_comparison(query, refs)
        seq_profile, fine_count = _profile_scores(query, refs, fine)
        family_profile, broad_count = _profile_scores(query, refs, broad)
        profile_count += fine_count + broad_count
        pairs += len(query) * len(refs)
        for key, sequence in query.items():
            label = broad[key]
            best_global = max(_global_identity(aligner, sequence, ref) for ref in refs.values())
            covered = [hit['identity'] for ref, hit in local_identity[key].items()
                       if broad[ref] == label and hit['coverage'] >= MIN_COVERAGE]
            predictions = {
                'raw_alignment': _top_label(seq_scores['alignment'][key], broad),
                'normalized_alignment': _top_label(seq_scores['alignment_ref_normalized'][key], broad),
                'fine_profiles': _top_label(seq_profile[key], fine_labels),
                'superfamily_profiles': _top_label(family_profile[key], broad_labels),
            }
            rows[key] = {'query': key, 'fold': fold, 'fine_family': fine[key], 'label': label,
                         'split': entries[key]['split'],
                         'max_global_identity': best_global, 'covered_positive_local_identity': max(covered, default=None),
                         'correct': {method: pred == label for method, pred in predictions.items()}}
    return rows, {'pairs': pairs, 'profiles': profile_count}


def _embeddings(entries, sequences, size, reference_policy):
    from Bio.Align import PairwiseAligner

    encode = _encode if size == '8m' else _encode_150m
    start = time.perf_counter()
    pooled, local, model = encode(sequences)
    embedded = time.perf_counter()
    aligner = PairwiseAligner(mode='local', open_gap_score=LOCAL_OPEN, extend_gap_score=LOCAL_EXTEND)
    labels = {key: superfamily(row['family']) for key, row in entries.items()}
    rows = {}
    for _, query, refs in comparisons(entries, sequences, reference_policy):
        for key in query:
            scores = {ref: float(pooled[key] @ pooled[ref]) for ref in refs}
            residue_scores = {ref: _residue_local_score(local[key], local[ref], aligner) for ref in refs}
            rows[key] = {f'esm2_{size}_mean': _top_label(scores, labels) == labels[key],
                         f'esm2_{size}_local': _top_label(residue_scores, labels) == labels[key]}
    return rows, {'model': model, 'encodings': len(sequences), 'encoding_seconds': round(embedded - start, 3),
                  'comparison_seconds': round(time.perf_counter() - embedded, 3)}


def evaluate(folder, output, arm, reference_policy):
    entries, sequences, fingerprint = _panel(folder)
    if arm == 'conventional':
        cases, work = _classical(entries, sequences, reference_policy)
    else:
        cases, work = _embeddings(entries, sequences, arm, reference_policy)
    write_json(output, {'schema': SCHEMA, 'panel_sha256': fingerprint, 'arm': arm,
                        'reference_policy': reference_policy, 'cases': cases, 'work': work,
                        'environment': {pkg: version(pkg) for pkg in ('biopython', 'numpy', 'pyfamsa', 'pyhmmer',
                                                                        'torch', 'transformers')}})
    print(f'{arm}: scored {len(cases)} domains')


def _fold_scores(rows, methods):
    families = defaultdict(list)
    for row in rows:
        families[row['fold'], row['label'], row['fine_family']].append(row)
    broad_groups = defaultdict(list)
    for (fold, label, _), cases in families.items():
        broad_groups[fold, label].append({method: mean(case['correct'][method] for case in cases) for method in methods})
    folds = defaultdict(list)
    for (fold, _), fine_scores in broad_groups.items():
        folds[fold].append({method: mean(row[method] for row in fine_scores) for method in methods})
    return {fold: {method: mean(row[method] for row in broad_scores) for method in methods}
            for fold, broad_scores in folds.items()}


def _cluster_summary(rows, methods):
    fold_scores = _fold_scores(rows, methods)
    return {'domains': len(rows), 'fine_families': len({row['fine_family'] for row in rows}),
            'folds': len(fold_scores), 'hits': {method: sum(row['correct'][method] for row in rows) for method in methods},
            'equal_fold_accuracy': {method: mean(row[method] for row in fold_scores.values()) if fold_scores else None
                                    for method in methods}}


def _paired(rows, comparison):
    model = 'esm2_150m_local'
    estimates = _fold_scores(rows, (model, comparison))
    differences = [values[model] - values[comparison] for values in estimates.values()]
    randomizer = random.Random(SEED)
    samples = sorted(mean(randomizer.choices(differences, k=len(differences))) for _ in range(4000)) if differences else []
    omit_one = [mean(differences[:i] + differences[i + 1:]) for i in range(len(differences))] if len(differences) > 1 else []
    return {'equal_fold_gain': mean(differences) if differences else None,
            'fold_bootstrap95': [samples[99], samples[3899]] if len(differences) >= 4 else None,
            'omit_one_fold_range': [min(omit_one), max(omit_one)] if omit_one else None,
            'model_only': sum(row['correct'][model] and not row['correct'][comparison] for row in rows),
            'comparison_only': sum(row['correct'][comparison] and not row['correct'][model] for row in rows)}


def summarize(conventional, small, large, output):
    sources = [json.loads(path.read_text()) for path in (conventional, small, large)]
    fingerprint = sources[0]['panel_sha256']
    reference_policy = sources[0].get('reference_policy', 'three-way')
    cases = sources[0]['cases']
    if (tuple(source['arm'] for source in sources) != ('conventional', '8m', '150m')
            or any(source['panel_sha256'] != fingerprint or source['cases'].keys() != cases.keys()
                   or source.get('reference_policy', 'three-way') != reference_policy for source in sources)):
        raise ValueError('The three reports must use the same panel, case IDs, and reference policy')
    for source in sources[1:]:
        for key, correct in source['cases'].items():
            cases[key]['correct'].update(correct)
    methods = (*CLASSICAL, *MODELS)
    rows = list(cases.values())
    primary = [row for row in rows if row['max_global_identity'] < PRIMARY_IDENTITY]
    subgroups = {'under_25': primary, '25_to_under_30': [r for r in rows if .25 <= r['max_global_identity'] < .30],
                 '30_and_over': [r for r in rows if r['max_global_identity'] >= .30]}
    conventional_misses = [r for r in primary if not any(r['correct'][m] for m in CLASSICAL)]
    report = {'schema': SCHEMA, 'panel_sha256': fingerprint, 'reference_policy': reference_policy,
              'primary_max_global_identity_below': PRIMARY_IDENTITY,
              'overall': _cluster_summary(rows, methods),
              'by_tier': {key: _cluster_summary(entries, methods) for key, entries in subgroups.items()},
              'primary_by_fold': {fold: _cluster_summary([r for r in primary if r['fold'] == fold], methods)
                                  for fold in sorted({r['fold'] for r in primary})},
              'primary_paired_150m_local': {method: _paired(primary, method) for method in MAIN_COMPARISONS},
              'primary_with_covered_positive': _cluster_summary(
                  [r for r in primary if r['covered_positive_local_identity'] is not None], methods),
              'primary_all_conventional_misses': _cluster_summary(conventional_misses, methods),
              'work': {source['arm']: source['work'] for source in sources},
              'environment': {source['arm']: source['environment'] for source in sources}, 'cases': rows}
    write_json(output, report)
    print(json.dumps({key: report[key] for key in ('overall', 'by_tier', 'primary_paired_150m_local',
                                                   'primary_all_conventional_misses')}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest='stage', required=True)
    enroll = stages.add_parser('acquire')
    enroll.add_argument('folder', type=Path)
    enroll.add_argument('--descriptions', required=True, type=Path)
    enroll.add_argument('--classification', required=True, type=Path)
    stages.add_parser('restore').add_argument('folder', type=Path)
    score = stages.add_parser('score')
    score.add_argument('folder', type=Path)
    score.add_argument('--arm', required=True, choices=('conventional', '8m', '150m'))
    score.add_argument('--reference-policy', choices=REFERENCE_POLICIES, default='three-way')
    score.add_argument('--out', required=True, type=Path)
    group = stages.add_parser('summarize')
    group.add_argument('--conventional', required=True, type=Path)
    group.add_argument('--small', required=True, type=Path)
    group.add_argument('--large', required=True, type=Path)
    group.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    if args.stage == 'acquire':
        acquire(args.folder, args.descriptions, args.classification)
    elif args.stage == 'restore':
        restore(args.folder)
    elif args.stage == 'score':
        evaluate(args.folder, args.out, args.arm, args.reference_policy)
    else:
        summarize(args.conventional, args.small, args.large, args.out)


if __name__ == '__main__':
    main()
