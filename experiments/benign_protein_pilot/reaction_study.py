"""Test exact reaction labels on reviewed natural enzymes with experimental evidence."""
from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from zoneinfo import ZoneInfo

from .related_study import LOCAL_EXTEND, LOCAL_OPEN, _encode, _representations, _residue_local_score
from .structural_study import METHODS, PRIMARY_IDENTITY, _global_identity, _profile_scores, _top_label
from .study import ALPHABET, MIN_COVERAGE, _alignment_comparison, sha, write_json

SCHEMA = 'natural-experimental-enzyme-reactions-1'
TAXON_GROUPS = {
    'microbial': (1423, 1582, 1590, 2190, 2238, 2287, 2336, 83333, 100226, 187420, 224308, 224325,
                  243232, 243274, 262724, 272944, 273057, 273116, 309800),
    'plant_and_algal': (3055, 3218, 3562, 3702, 3847, 3880, 4113, 4530, 4577, 39947),
    'animal_and_fungal': (4896, 4932, 559292, 284812, 6239, 7227, 7955, 9031, 9606, 9913, 10090, 10116),
}
GROUP_BY_TAXON = {taxon: group for group, taxa in TAXON_GROUPS.items() for taxon in taxa}
REACTIONS = {
    'carbohydrate_hydrolysis': ('3.2.1.1', '3.2.1.2', '3.2.1.3', '3.2.1.4', '3.2.1.8', '3.2.1.10',
                                '3.2.1.20', '3.2.1.22', '3.2.1.23'),
    'sugar_phosphate_metabolism': ('2.2.1.2', '4.1.2.4', '4.1.2.13', '5.1.3.1', '5.3.1.1', '5.3.1.6',
                                   '5.3.1.8', '5.3.1.9'),
    'nucleobase_and_nucleoside_deamination': ('3.5.4.2', '3.5.4.3', '3.5.4.4', '3.5.4.5', '3.5.4.6',
                                            '3.5.4.9', '3.5.4.12'),
}
GROUP_BY_REACTION = {label: group for group, labels in REACTIONS.items() for label in labels}
MAX_PER_TAXON_GROUP = 6
ESM2_150M = {
    'repository': 'facebook/esm2_t30_150M_UR50D',
    'revision': 'a695f6045e2e32885fa60af20c13cb35398ce30c',
    'weight_sha256': 'c3f1da8aea53bddd32c246c86168c23b9fd72341fb9db9a94436f855f5053566',
}


def _fetch(label):
    taxa = ' OR '.join(f'organism_id:{taxon}' for taxon in GROUP_BY_TAXON)
    url = 'https://rest.uniprot.org/uniprotkb/search?' + urlencode({
        'query': f'reviewed:true AND ec:{label} AND ({taxa})', 'format': 'json', 'size': 250,
        'fields': 'accession,protein_name,organism_id,sequence,ec,cc_catalytic_activity,version,sequence_version',
    })
    request = Request(url, headers={'User-Agent': 'iri-constructs/experimental-metabolic-enzymes'})
    with urlopen(request, timeout=25) as response:
        payload = response.read(2_000_001)
        release = response.headers.get('X-UniProt-Release')
        if not release or response.headers.get('Link') or len(payload) > 2_000_000:
            raise ValueError(f'{label}: incomplete or unexpected UniProt response')
    return json.loads(payload)['results'], {'sha256': sha(payload), 'release': release, 'url': url}


def _evidence(row, label):
    reactions = [comment['reaction'] for comment in row.get('comments', ())
                 if comment['commentType'] == 'CATALYTIC ACTIVITY']
    if {reaction['ecNumber'] for reaction in reactions if 'ecNumber' in reaction} != {label}:
        return []
    return sorted({citation['id'] for reaction in reactions if reaction.get('ecNumber') == label
                   for citation in reaction.get('evidences', ())
                   if citation['evidenceCode'] == 'ECO:0000269' and citation.get('source') == 'PubMed'
                   and citation.get('id')})


def _prior_sequence_hashes():
    here = Path(__file__).parent
    return {row['sha256'] for name in ('panel.freeze.json', 'related_panel.freeze.json', 'structural_panel.freeze.json',
                                     'structural_panel.c.23.freeze.json', 'structural_panel.c.26.freeze.json')
            for row in json.loads((here / name).read_text())['entries']}


def _candidates(rows, label):
    for row in rows:
        sequence = row['sequence']['value']
        citations = _evidence(row, label)
        taxon = row['organism']['taxonId']
        if not citations or not 100 <= len(sequence) <= 850 or set(sequence) - ALPHABET:
            continue
        if taxon not in GROUP_BY_TAXON or not row['entryType'].startswith('UniProtKB reviewed'):
            raise ValueError(f'{label}: source returned an unexpected record')
        audit = row['entryAudit']
        yield {'id': row['primaryAccession'], 'label': label, 'taxon': taxon, 'group': GROUP_BY_TAXON[taxon],
               'entry_version': audit['entryVersion'], 'sequence_version': audit['sequenceVersion'],
               'length': len(sequence), 'sha256': sha(sequence.encode('ascii')), 'reaction_pubmed_ids': citations}, sequence


def acquire(folder):
    labels = sorted(GROUP_BY_REACTION)
    sources, selected = {}, []
    seen = _prior_sequence_hashes()
    with ThreadPoolExecutor(max_workers=3) as workers:
        for label, (rows, source) in zip(labels, workers.map(_fetch, labels), strict=True):
            sources[label] = source
            taxa, sizes = set(), defaultdict(int)
            candidates = sorted(_candidates(rows, label), key=lambda pair: (pair[0]['taxon'], pair[0]['id']))
            for row, sequence in candidates:
                if (row['taxon'] not in taxa and row['sha256'] not in seen
                        and sizes[row['group']] < MAX_PER_TAXON_GROUP):
                    selected.append((row, sequence))
                    sizes[row['group']] += 1
                    taxa.add(row['taxon'])
                    seen.add(row['sha256'])

    groups_by_label = defaultdict(set)
    for row, _ in selected:
        groups_by_label[row['label']].add(row['group'])
    shared_labels = {label for label, groups in groups_by_label.items() if len(groups) >= 2}
    eligible_labels = {label for labels in REACTIONS.values() if len(shared_labels.intersection(labels)) >= 2
                       for label in shared_labels.intersection(labels)}
    entries = [row for row, _ in selected if row['label'] in eligible_labels]
    releases = {source['release'] for source in sources.values()}
    if len(releases) != 1 or not entries:
        raise ValueError('UniProt release changed or no comparable, experimentally documented reactions were found')
    folder.mkdir(parents=True, exist_ok=True)
    write_json(folder / 'panel.freeze.json', {
        'schema': SCHEMA, 'retrieved_on': datetime.now(ZoneInfo('America/New_York')).date().isoformat(),
        'release': releases.pop(), 'terms': 'https://www.uniprot.org/help/license',
        'entries': entries, 'sources': sources,
    })
    write_json(folder / 'sequences.private.json', {row['id']: sequence for row, sequence in selected
                                                  if row['label'] in eligible_labels})
    print(f'Froze {len(entries)} experimentally supported proteins across {len(eligible_labels)} reactions')


def _summarize(rows):
    labels = defaultdict(list)
    for row in rows:
        labels[row['label']].append(row)
    return {
        'queries': len(rows), 'reactions': len(labels),
        'correct': {method: sum(row['correct'][method] for row in rows) for method in METHODS},
        'equal_reaction_accuracy': {
            method: sum(sum(row['correct'][method] for row in group) / len(group) for group in labels.values())
            / len(labels) if labels else None for method in METHODS
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
    if (freeze['schema'] != SCHEMA or not entries or len(entries) != len(freeze['entries'])
            or sequences.keys() != entries.keys() or len({row['sha256'] for row in entries.values()}) != len(entries)
            or any(row['label'] not in GROUP_BY_REACTION or not row['reaction_pubmed_ids']
                   or GROUP_BY_TAXON.get(row['taxon']) != row['group'] or len(sequences[key]) != row['length']
                   or set(sequences[key]) - ALPHABET or sha(sequences[key].encode('ascii')) != row['sha256']
                   for key, row in entries.items())):
        raise ValueError('Private protein sequences or reaction evidence do not match the frozen panel')
    return entries, sequences, sha(freeze_bytes)


def _encode_150m(sequences):
    import torch
    from huggingface_hub import snapshot_download
    from transformers import AutoModel, AutoTokenizer

    from iri_constructs.model_adapters import prepare_device
    from iri_constructs.runtime import digest_file

    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    prepare_device('cpu', 20261122)
    spec = ESM2_150M
    snapshot = Path(snapshot_download(spec['repository'], revision=spec['revision'], local_files_only=True,
        allow_patterns=['config.json', 'tokenizer_config.json', 'special_tokens_map.json', 'vocab.txt',
                        'model.safetensors']))
    if snapshot.name != spec['revision'] or digest_file(snapshot / 'model.safetensors') != spec['weight_sha256']:
        raise ValueError('ESM-2 150M cache does not match the pinned checkpoint')
    tokenizer = AutoTokenizer.from_pretrained(snapshot, local_files_only=True)
    encoder, loading = AutoModel.from_pretrained(snapshot, local_files_only=True, add_pooling_layer=False,
                                                 output_loading_info=True)
    if (loading.get('missing_keys') or loading.get('mismatched_keys')
            or any(not key.startswith('lm_head.') for key in loading.get('unexpected_keys', ()))):
        raise ValueError('Checkpoint parameters do not match the ESM-2 encoder')
    encoder = encoder.float().eval()
    pooled, residues = {}, {}
    with torch.inference_mode():
        for key, sequence in sequences.items():
            inputs = tokenizer(sequence, return_tensors='pt', truncation=False, return_special_tokens_mask=True)
            if tokenizer.unk_token_id in inputs['input_ids'][0].tolist():
                raise ValueError('Tokenizer returned an unknown residue token')
            keep = ~inputs.pop('special_tokens_mask')[0].bool() & inputs['attention_mask'][0].bool()
            states = encoder(**inputs).last_hidden_state[0][keep]
            if states.shape[0] != len(sequence) or not torch.isfinite(states).all():
                raise ValueError('Invalid ESM-2 residue representation')
            pooled[key], residues[key] = _representations(states.float().numpy())
    return pooled, residues, {**spec, 'device': 'cpu', 'representation': 'Final layer; special tokens excluded'}


def score(folder, output, esm2_size='8m'):
    from Bio.Align import PairwiseAligner, substitution_matrices

    entries, sequences, freeze_hash = _panel(folder)
    began = time.perf_counter()
    encoder = _encode if esm2_size == '8m' else _encode_150m
    pooled, local, model = encoder(sequences)
    encoded = time.perf_counter()
    residue_aligner = PairwiseAligner(mode='local', open_gap_score=LOCAL_OPEN, extend_gap_score=LOCAL_EXTEND)
    global_aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                                     open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)
    labels = {key: row['label'] for key, row in entries.items()}
    cases, pairs, profiles = [], 0, 0
    for comparison, comparison_labels in REACTIONS.items():
        keys = [key for key, label in labels.items() if label in comparison_labels]
        for group in TAXON_GROUPS:
            queries = {key: sequences[key] for key in keys if entries[key]['group'] == group}
            if not queries:
                continue
            refs = {key: sequences[key] for key in keys if entries[key]['group'] != group}
            if {labels[key] for key in refs} != {labels[key] for key in keys}:
                raise ValueError(f'{comparison}/{group}: a reaction lacks reference proteins')
            pairs += len(queries) * len(refs)
            seq_scores, distances = _alignment_comparison(queries, refs)
            profile_scores, count = _profile_scores(queries, refs, labels)
            profiles += count
            for query, sequence in queries.items():
                label = labels[query]
                positives = [key for key in refs if labels[key] == label]
                best = max(_global_identity(global_aligner, sequence, refs[key]) for key in positives)
                local_identities = [distances[query][key]['identity'] for key in positives
                                    if distances[query][key]['coverage'] >= MIN_COVERAGE]
                predictions = {
                    'pairwise_raw': _top_label(seq_scores['alignment'][query], labels),
                    'pairwise': _top_label(seq_scores['alignment_ref_normalized'][query], labels),
                    'profile': _top_label(profile_scores[query], {ec: ec for ec in comparison_labels}),
                    'mean_cosine': _top_label({key: float(pooled[query] @ pooled[key]) for key in refs}, labels),
                    'residue_local': _top_label({key: _residue_local_score(local[query], local[key], residue_aligner)
                                                 for key in refs}, labels),
                }
                cases.append({'query': query, 'comparison': comparison, 'group': group, 'label': label,
                              'max_positive_global_identity': best,
                              'max_covered_positive_local_identity': max(local_identities, default=None),
                              'correct': {method: result == label for method, result in predictions.items()}})

    low = [row for row in cases if row['max_positive_global_identity'] < PRIMARY_IDENTITY]
    report = {'schema': SCHEMA, 'panel_sha256': freeze_hash, 'model': model,
              'environment': {name: version(name) for name in ('biopython', 'numpy', 'pyfamsa', 'pyhmmer', 'torch',
                                                               'transformers')},
              'overall': _summarize(cases), 'primary': _summarize(low),
              'under_25': _summarize([row for row in low if row['max_positive_global_identity'] < 0.25]),
              'low_with_covered_positive': _summarize([row for row in low
                                                      if row['max_covered_positive_local_identity'] is not None]),
              'primary_by_comparison': {label: _summarize([row for row in low if row['comparison'] == label])
                                        for label in REACTIONS},
              'cost': {'encoder_calls': len(sequences), 'pairwise_comparisons': pairs, 'profiles_built': profiles,
                       'encoding_seconds': round(encoded - began, 3),
                       'scoring_seconds': round(time.perf_counter() - encoded, 3)}, 'cases': cases}
    write_json(output, report)
    print(json.dumps({key: report[key] for key in ('overall', 'primary', 'low_with_covered_positive', 'cost')}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stages = parser.add_subparsers(dest='stage', required=True)
    stages.add_parser('acquire').add_argument('folder', type=Path)
    compare = stages.add_parser('score')
    compare.add_argument('folder', type=Path)
    compare.add_argument('--out', required=True, type=Path)
    compare.add_argument('--esm2-size', choices=('8m', '150m'), default='8m')
    args = parser.parse_args()
    if args.stage == 'acquire':
        acquire(args.folder)
    else:
        score(args.folder, args.out, esm2_size=args.esm2_size)


if __name__ == '__main__':
    main()
