"""Compare simple input segmentation on four documented natural housekeeping proteins."""
from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from urllib.request import Request, urlopen
from zoneinfo import ZoneInfo

from .study import ALPHABET, MODEL_HASH, sha, write_json

SCHEMA = 'native-housekeeping-fusions-1'
QUERIES = {
    'P11172': (9606, 'ump_synthase', {'OPRTase': '2.4.2.10', 'OMPdecase': '4.1.1.23'}),
    'Q01637': (7227, 'ump_synthase', {'OPRTase': '2.4.2.10', 'OMPdecase': '4.1.1.23'}),
    'P07244': (559292, 'ade5_7', {'GARS': '6.3.4.13', 'AIRS': '6.3.3.1'}),
    'P00909': (83333, 'trpcf', {'Indole-3-glycerol phosphate synthase': '4.1.1.48',
                               "N-(5'-phosphoribosyl)anthranilate isomerase": '5.3.1.24'}),
}
REFERENCES = {
    'P0A7E3': (83333, '2.4.2.10'), 'P08244': (83333, '4.1.1.23'),
    'P08178': (83333, '6.3.3.1'), 'P15640': (83333, '6.3.4.13'),
    'P03964': (224308, '4.1.1.48'), 'P20167': (224308, '5.3.1.24'),
    'P0A799': (83333, '2.7.2.3'), 'P0A858': (83333, '5.3.1.1'), 'P0A6P9': (83333, '4.2.1.11'),
}


def _labels(description):
    sections = [description, *description.get('includes', []), *description.get('contains', [])]
    return {value['value'] for section in sections
            for value in section.get('recommendedName', {}).get('ecNumbers', [])}


def _regions(record, names, length):
    regions = []
    for feature in record.get('features', []):
        if feature['type'] != 'Region' or feature['description'] not in names:
            continue
        start, end = (feature['location'][key] for key in ('start', 'end'))
        left, right = start['value'] - 1, end['value']
        if not (start['modifier'] == end['modifier'] == 'EXACT' and 0 <= left < right <= length):
            raise ValueError('Expected an exact region inside the protein')
        regions.append({'name': feature['description'], 'label': names[feature['description']], 'span': [left, right]})
    regions.sort(key=lambda region: region['span'])
    if (len(regions) != len(names) or set(names) != {region['name'] for region in regions} or
            any(a['span'][1] > b['span'][0] for a, b in zip(regions, regions[1:]))):
        raise ValueError('Expected two distinct documented regions')
    return regions


def acquire(folder):
    folder.mkdir(parents=True)
    entries, sequences, releases = {}, {}, set()
    for accession, spec in (QUERIES | REFERENCES).items():
        url = f'https://rest.uniprot.org/uniprotkb/{accession}.json'
        request = Request(url, headers={'User-Agent': 'iri-constructs/native-housekeeping'})
        with urlopen(request, timeout=20) as response:
            raw = response.read(2_000_001)
            releases.add(response.headers.get('X-UniProt-Release'))
        time.sleep(0.35)
        if len(raw) > 2_000_000:
            raise ValueError(f'{accession}: source record exceeds the size bound')
        record = json.loads(raw)
        sequence = record['sequence']['value']
        if (record['entryType'] != 'UniProtKB reviewed (Swiss-Prot)' or record['primaryAccession'] != accession or
                record['organism']['taxonId'] != spec[0] or not 1 <= len(sequence) <= 1022 or set(sequence) - ALPHABET):
            raise ValueError(f'{accession}: source differs from the selected natural protein')
        label_set = set(spec[2].values()) if accession in QUERIES else {spec[1]}
        if _labels(record['proteinDescription']) != label_set:
            raise ValueError(f'{accession}: source does not have the expected documented EC labels')
        entry = {'source': url, 'source_sha256': sha(raw), 'sha256': sha(sequence.encode('ascii')),
                 'entry_version': record['entryAudit']['entryVersion'],
                 'sequence_version': record['entryAudit']['sequenceVersion'], 'taxon': spec[0], 'length': len(sequence)}
        if accession in QUERIES:
            entry.update(role='query', architecture=spec[1], regions=_regions(record, spec[2], len(sequence)))
        else:
            entry.update(role='reference', label=spec[1])
        entries[accession], sequences[accession] = entry, sequence
        (folder / f'{accession}.private.json').write_bytes(raw)
    if len(releases) != 1 or None in releases:
        raise ValueError('UniProt release could not be fixed for the acquisition')
    write_json(folder / 'panel.freeze.json', {'schema': SCHEMA, 'entries': entries, 'release': releases.pop(),
               'retrieved_on': datetime.now(ZoneInfo('America/New_York')).date().isoformat(),
               'terms': 'https://www.uniprot.org/help/license'})
    write_json(folder / 'sequences.private.json', sequences)


def _top_labels(spans, scores):
    best = {}
    for span, row in zip(spans, scores, strict=True):
        for label, score in row.items():
            if label not in best or score > best[label][0]:
                best[label] = (score, span)
    return {label: best[label][1] for label in sorted(best, key=lambda key: (-best[key][0], key))[:2]}


def _recovery(calls, regions):
    hits, supported = 0, 0
    for region in regions:
        if region['label'] in calls:
            hits += 1
            left, right = region['span']
            start, end = calls[region['label']]
            supported += max(0, min(end, right) - max(start, left)) >= 0.80 * (right - left)
    return {'label_hits': hits, 'region_supported_hits': supported}


def score(folder, output):
    import torch
    from Bio.Align import PairwiseAligner, substitution_matrices
    from iri_constructs.metrics import cosine, mean_pool
    from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device

    freeze_bytes = (folder / 'panel.freeze.json').read_bytes()
    freeze = json.loads(freeze_bytes)
    seqs = json.loads((folder / 'sequences.private.json').read_text())
    entries = freeze['entries']
    if (freeze['schema'] != SCHEMA or seqs.keys() != entries.keys() or seqs.keys() != (QUERIES | REFERENCES).keys()
            or any(sha(sequence.encode('ascii')) != entries[key]['sha256'] for key, sequence in seqs.items())):
        raise ValueError('Sequence snapshot differs from the frozen natural protein panel')
    aligner = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)
    targets = {entries[r]['label']: (seqs[r], aligner.score(seqs[r], seqs[r])) for r in REFERENCES}
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    model = LocalAdapter('esm2', prepare_device('cpu', 20261122), ESM2_REVISION, offline=True)
    cases = []
    try:
        if model.metadata['weight_sha256'] != MODEL_HASH:
            raise ValueError('Pinned ESM-2 weight hash does not match')
        vectors = {entries[key]['label']: mean_pool(model.residues(seqs[key])) for key in REFERENCES}

        def embedding_scores(states):
            vector = mean_pool(states)
            return {family: cosine(vector, ref) for family, ref in vectors.items()}

        for accession in QUERIES:
            protein, record = seqs[accession], entries[accession]
            length, mid = len(protein), len(protein) // 2
            regions = record['regions']
            plans = {'whole': [(0, length)], 'midpoint': [(0, mid), (mid, length)],
                     'annotated_regions': [tuple(region['span']) for region in regions]}
            scores = {}
            states = model.residues(protein)
            for method, spans in plans.items():
                words = [protein[start:end] for start, end in spans]
                aligned = [{label: aligner.score(word, ref) / self_score for label, (ref, self_score) in targets.items()}
                           for word in words]
                scores[f'blosum62_{method}'] = _recovery(_top_labels(spans, aligned), regions)
                encoded = [states] if method == 'whole' else [model.residues(word) for word in words]
                scores[f'esm2_{method}'] = _recovery(_top_labels(spans, map(embedding_scores, encoded)), regions)
            spans = plans['midpoint']
            pooled = [embedding_scores(states[start:end]) for start, end in spans]
            scores['esm2_midpoint_pooled_from_whole'] = _recovery(_top_labels(spans, pooled), regions)
            cases.append({'query': accession, 'architecture': record['architecture'], 'scores': scores})
        model_metadata = model.metadata
    finally:
        model.close()
    totals = defaultdict(lambda: defaultdict(int))
    for case in cases:
        for method, result in case['scores'].items():
            for metric, value in result.items():
                totals[method][metric] += value
    report = {'schema': SCHEMA, 'panel_sha256': sha(freeze_bytes), 'model': model_metadata, 'queries': len(cases),
              'eligible_regions': 2 * len(cases), 'total': totals, 'cases': cases,
              'encoder_calls': {'references': len(REFERENCES), 'whole_shared_with_midpoint_pooling': len(cases),
                                'separately_encoded_midpoint_halves': 2 * len(cases), 'annotated_region_inputs': 2 * len(cases)}}
    write_json(output, report)
    print(json.dumps(totals, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    sub.add_parser('acquire').add_argument('folder', type=Path)
    run = sub.add_parser('score')
    run.add_argument('folder', type=Path)
    run.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.stage == 'acquire':
        acquire(args.folder)
    else:
        score(args.folder, args.out)


if __name__ == '__main__':
    main()
