"""Evaluate held-out families of curated, naturally occurring metabolic protein domains."""
from __future__ import annotations

import argparse
import json
import re
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from urllib.request import Request, urlopen
from zoneinfo import ZoneInfo

from .related_study import LOCAL_EXTEND, LOCAL_OPEN, _encode, _residue_local_score
from .study import ALPHABET, MIN_COVERAGE, _alignment_comparison, sha, write_json

SCHEMA = 'natural-metabolic-structural-families-1'
RELEASE = '2.08'
DESCRIPTIONS = '41aad433fda2d30eb05fb5a4d03692345e0cce39134a8f7cddb2ec140b5c8af8'
CLASSIFICATION = '979eceafbc2cb8c50a6549d156d888478a43003e8f9e21e1940d11bcf3d728ca'
TAXA = frozenset({1423, 1582, 1590, 2190, 2238, 2287, 2336, 3055, 3562, 3702, 3847, 3880, 4530,
                  4577, 4932, 7227, 9606, 10090, 10116, 83333, 100226, 187420, 224308, 224325,
                  243232, 243274, 262724, 272944, 273057, 273116, 309800})
# Curated carbohydrate, nucleotide, amino acid and cofactor metabolism in the TIM-barrel fold.
FAMILIES = frozenset(f'c.1.{superfamily}.{family}' for superfamily, families in {
    2: (1, 2, 3, 4, 6), 8: (1, 3, 13), 9: (1, 6, 7, 8, 10), 10: (1, 2, 3, 4, 5, 8),
    11: (1, 2), 12: (1, 2, 3, 8), 15: (2, 3, 6), 22: (1, 2), 23: (1, 2),
}.items() for family in families)
MAX_DOMAINS_PER_FAMILY = 4
PRIMARY_IDENTITY = 0.30
METHODS = ('pairwise_raw', 'pairwise', 'profile', 'mean_cosine', 'residue_local')


def superfamily(family):
    return family.rsplit('.', 1)[0]


def metadata_candidates(descriptions, classification):
    annotations = {}
    for line in descriptions.splitlines():
        if not line.startswith('#'):
            key, level, _, _, name = line.split('\t', 4)
            if level in ('sf', 'fa', 'dm', 'sp'):
                annotations[key] = name

    candidates = defaultdict(list)
    for line in classification.splitlines():
        if line.startswith('#'):
            continue
        domain, pdb, location, family, node, hierarchy = line.split('\t')
        if family not in FAMILIES or ',' in location:
            continue
        nodes = dict(level.split('=') for level in hierarchy.split(','))
        protein, species = (annotations[nodes[key]] for key in ('dm', 'sp'))
        taxon = re.search(r'\[TaxId: (\d+)\]', species)
        if (not taxon or int(taxon[1]) not in TAXA or
                protein.lower().startswith(('automated', 'hypothetical', 'putative', 'probable')) or
                re.search(r'engineered|\bmutant\b|\bchimera\b|synthetic', species, re.IGNORECASE)):
            continue
        candidates[family, int(taxon[1])].append({
            'id': domain, 'node': int(node), 'pdb': pdb, 'family': family, 'taxon': int(taxon[1]),
            'protein': protein, 'species': species, 'family_name': annotations[nodes['fa']],
            'superfamily_name': annotations[nodes['sf']],
        })
    return candidates


def _source_file(url):
    request = Request(url, headers={'User-Agent': 'iri-constructs/natural-structural-domain-study'})
    with urlopen(request, timeout=30) as response:
        return response.read()


def domain_source(node):
    return f'https://scop.berkeley.edu/astral/seq/ver={RELEASE}&id={node}&seqOption=0&output=text'


def _domain_sequence(row):
    payload = _source_file(domain_source(row['node']))
    header, *lines = payload.decode('ascii').splitlines()
    if header.split()[:2] != [f">{row['id']}", row['family']]:
        raise ValueError(f"ASTRAL response did not identify {row['id']}")
    sequence = ''.join(lines).upper()
    if not 100 <= len(sequence) <= 620 or set(sequence) - ALPHABET:
        return None
    return dict(row, length=len(sequence), sha256=sha(sequence.encode('ascii')), source_sha256=sha(payload)), sequence


def _prior_sequence_hashes():
    here = Path(__file__).parent
    return {row['sha256'] for name in ('panel.freeze.json', 'related_panel.freeze.json')
            for row in json.loads((here / name).read_text())['entries']}


def _choose_taxon(candidates):
    for row in sorted(candidates, key=lambda item: item['id'])[:3]:
        result = _domain_sequence(row)
        if result is not None:
            return result
    return None


def family_folds(entries):
    """Assign all records of each fine family to one fold before scoring."""
    families = defaultdict(set)
    for row in entries:
        families[superfamily(row['family'])].add(row['family'])
    return {family: index % 3 for related in families.values() if len(related) >= 2
            for index, family in enumerate(sorted(related, key=lambda item: int(item.rsplit('.', 1)[1])))}


def acquire(folder, descriptions, classification):
    description_bytes, class_bytes = descriptions.read_bytes(), classification.read_bytes()
    if (sha(description_bytes), sha(class_bytes)) != (DESCRIPTIONS, CLASSIFICATION):
        raise ValueError('SCOPe classification metadata does not match the pinned release')
    candidates = metadata_candidates(description_bytes.decode(), class_bytes.decode())
    seen = _prior_sequence_hashes()
    counts = defaultdict(int)
    selected = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        for result in pool.map(_choose_taxon, [rows for _, rows in sorted(candidates.items())]):
            if result is None:
                continue
            row, sequence = result
            if row['sha256'] not in seen and counts[row['family']] < MAX_DOMAINS_PER_FAMILY:
                selected.append((row, sequence))
                seen.add(row['sha256'])
                counts[row['family']] += 1

    folds = family_folds([row for row, _ in selected])
    entries = [dict(row, fold=folds[row['family']]) for row, _ in selected if row['family'] in folds]
    sequences = {row['id']: sequence for row, sequence in selected if row['family'] in folds}
    folder.mkdir(parents=True, exist_ok=True)
    write_json(folder / 'panel.freeze.json', {
        'schema': SCHEMA, 'source': f'SCOPe and ASTRAL {RELEASE}', 'descriptions_sha256': DESCRIPTIONS,
        'classification_sha256': CLASSIFICATION, 'source_documentation': 'https://scop.berkeley.edu/astral/ver=2.08',
        'terms': 'https://scop.berkeley.edu/about',
        'retrieved_on': datetime.now(ZoneInfo('America/New_York')).date().isoformat(), 'entries': entries,
        'primary_max_global_identity_below': PRIMARY_IDENTITY,
    })
    write_json(folder / 'sequences.private.json', sequences)
    print(f'Froze {len(entries)} domains, {len(folds)} families and '
          f'{len({superfamily(family) for family in folds})} superfamilies')


def _global_identity(aligner, query, reference):
    hit = next(iter(aligner.align(query, reference)))
    counts = hit.counts()
    return counts.identities / hit.length


def _identity_tier(identity):
    if identity < 0.25:
        return 'under_25'
    if identity < PRIMARY_IDENTITY:
        return '25_to_under_30'
    if identity < 0.40:
        return '30_to_under_40'
    return '40_or_more'


def _profile_scores(queries, references, entries):
    from pyfamsa import Aligner, Sequence
    from pyhmmer import easel, plan7

    alphabet = easel.Alphabet.amino()
    aligner = Aligner(threads=1)
    builder = plan7.Builder(alphabet, seed=42)
    background = plan7.Background(alphabet)
    pipeline = plan7.Pipeline(alphabet, background=background, F1=1.0, F2=1.0, F3=1.0,
                              T=-1000.0, domT=-1000.0)
    grouped = defaultdict(list)
    for key, sequence in references.items():
        grouped[entries[key]['family']].append((key, sequence))
    block = easel.DigitalSequenceBlock(alphabet, [
        easel.TextSequence(name=key.encode(), sequence=sequence).digitize(alphabet)
        for key, sequence in queries.items()
    ])
    hits = {key: {} for key in queries}
    for family, members in sorted(grouped.items()):
        if len(members) == 1:
            _, sequence = members[0]
            digital = easel.TextSequence(name=family.encode(), sequence=sequence).digitize(alphabet)
            hmm, _, _ = builder.build(digital, background)
        else:
            msa = aligner.align(Sequence(key.encode(), sequence.encode()) for key, sequence in members)
            aligned = [easel.TextSequence(name=item.id, sequence=item.sequence.decode()) for item in msa]
            digital = easel.TextMSA(name=family.encode(), sequences=aligned).digitize(alphabet)
            hmm, _, _ = builder.build_msa(digital, background)
        for hit in pipeline.search_hmm(hmm, block):
            name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            hits[name][family] = hit.score
    return hits, len(grouped)


def _top_label(scores, labels):
    best = {}
    for key, score in scores.items():
        label = labels[key]
        best[label] = max(best.get(label, float('-inf')), score)
    if not best:
        return None
    highest = max(best.values())
    winners = [label for label, score in best.items() if score == highest]
    return winners[0] if len(winners) == 1 else None


def _summary(rows):
    by_family = defaultdict(list)
    for row in rows:
        by_family[row['family']].append(row)
    return {
        'queries': len(rows), 'fine_families': len(by_family),
        'superfamilies': len({row['label'] for row in rows}),
        'correct': {method: sum(row['correct'][method] for row in rows) for method in METHODS},
        'equal_fine_family_accuracy': {
            method: sum(sum(row['correct'][method] for row in group) / len(group) for group in by_family.values())
            / len(by_family) if by_family else None for method in METHODS
        },
        'residue_local_disagreements': {method: {
            'residue_only': sum(row['correct']['residue_local'] and not row['correct'][method] for row in rows),
            'other_only': sum(row['correct'][method] and not row['correct']['residue_local'] for row in rows),
        } for method in METHODS if method != 'residue_local'},
    }


def _panel(folder):
    freeze_bytes = (folder / 'panel.freeze.json').read_bytes()
    freeze = json.loads(freeze_bytes)
    sequences = json.loads((folder / 'sequences.private.json').read_text())
    entries = {row['id']: row for row in freeze['entries']}
    folds = family_folds(entries.values())
    if (freeze['schema'] != SCHEMA or len({superfamily(family) for family in folds}) < 2
            or len(entries) != len(freeze['entries']) or entries.keys() != sequences.keys()
            or freeze['primary_max_global_identity_below'] != PRIMARY_IDENTITY
            or len({row['sha256'] for row in entries.values()}) != len(entries)
            or any(row['family'] not in FAMILIES or row['taxon'] not in TAXA or row['fold'] != folds.get(row['family'])
                   or len(sequences[key]) != row['length'] or set(sequences[key]) - ALPHABET
                   or sha(sequences[key].encode('ascii')) != row['sha256'] for key, row in entries.items())):
        raise ValueError('Private domains or split do not match the frozen panel')
    return entries, sequences, sha(freeze_bytes)


def score(folder, output):
    from Bio.Align import PairwiseAligner, substitution_matrices

    entries, sequences, freeze_hash = _panel(folder)
    began = time.perf_counter()
    pooled, local, model = _encode(sequences)
    encoded = time.perf_counter()
    residue_aligner = PairwiseAligner(mode='local', open_gap_score=LOCAL_OPEN, extend_gap_score=LOCAL_EXTEND)
    global_aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                                     open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)
    labels = {key: superfamily(row['family']) for key, row in entries.items()}
    profile_labels = {row['family']: superfamily(row['family']) for row in entries.values()}
    cases, pairs, models = [], 0, 0
    for fold in sorted({row['fold'] for row in entries.values()}):
        queries = {key: sequences[key] for key, row in entries.items() if row['fold'] == fold}
        refs = {key: sequences[key] for key, row in entries.items() if row['fold'] != fold}
        if {labels[key] for key in refs} != set(labels.values()):
            raise ValueError(f'Fold {fold} is missing a reference superfamily')
        pairs += len(queries) * len(refs)
        sequence_scores, distances = _alignment_comparison(queries, refs)
        profile_scores, count = _profile_scores(queries, refs, entries)
        models += count
        for query, sequence in queries.items():
            label = labels[query]
            global_identity = max(_global_identity(global_aligner, sequence, ref) for ref in refs.values())
            covered_identities = [row['identity'] for key, row in distances[query].items()
                                  if labels[key] == label and row['coverage'] >= MIN_COVERAGE]
            predictions = {
                'pairwise_raw': _top_label(sequence_scores['alignment'][query], labels),
                'pairwise': _top_label(sequence_scores['alignment_ref_normalized'][query], labels),
                'profile': _top_label(profile_scores[query], profile_labels),
                'mean_cosine': _top_label({key: float(pooled[query] @ pooled[key]) for key in refs}, labels),
                'residue_local': _top_label({key: _residue_local_score(local[query], local[key], residue_aligner)
                                             for key in refs}, labels),
            }
            cases.append({'query': query, 'family': entries[query]['family'], 'label': label, 'fold': fold,
                          'max_global_identity': global_identity, 'tier': _identity_tier(global_identity),
                          'max_covered_positive_local_identity': max(covered_identities, default=None),
                          'predicted': predictions,
                          'correct': {method: value == label for method, value in predictions.items()}})

    primary = [row for row in cases if row['max_global_identity'] < PRIMARY_IDENTITY]
    report = {
        'schema': SCHEMA, 'panel_sha256': freeze_hash, 'model': model,
        'environment': {name: version(name) for name in ('biopython', 'numpy', 'pyfamsa', 'pyhmmer',
                                                         'torch', 'transformers')},
        'overall': _summary(cases), 'primary': _summary(primary),
        'by_tier': {tier: _summary([row for row in cases if row['tier'] == tier])
                    for tier in ('under_25', '25_to_under_30', '30_to_under_40', '40_or_more')},
        'primary_by_superfamily': {label: _summary([row for row in primary if row['label'] == label])
                                   for label in sorted(set(labels.values()))},
        'cost': {'model_encodings': len(sequences), 'pairwise_comparisons': pairs, 'profiles_built': models,
                 'encoding_seconds': round(encoded - began, 3),
                 'scoring_seconds': round(time.perf_counter() - encoded, 3)},
        'cases': cases,
    }
    write_json(output, report)
    print(json.dumps({key: report[key] for key in ('overall', 'primary', 'by_tier', 'cost')}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest='stage', required=True)
    download = stages.add_parser('acquire')
    download.add_argument('folder', type=Path)
    download.add_argument('--descriptions', required=True, type=Path)
    download.add_argument('--classification', required=True, type=Path)
    compare = stages.add_parser('score')
    compare.add_argument('folder', type=Path)
    compare.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    if args.stage == 'acquire':
        acquire(args.folder, args.descriptions, args.classification)
    else:
        score(args.folder, args.out)


if __name__ == '__main__':
    main()
