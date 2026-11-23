"""Compare characterization of natural enzymes against related competing labels."""
from __future__ import annotations

import argparse
import csv
import io
import json
import time
from collections import Counter
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from zoneinfo import ZoneInfo

from .study import (ALPHABET, HIGH_IDENTITY, LOW_IDENTITY, MIN_COVERAGE, MODEL_HASH,
                    _alignment_comparison, _tier, sha, write_json)

SCHEMA = 'related-housekeeping-families-1'
REFERENCE_TAXA = (83333, 224308, 416870, 196627, 220668, 100226, 300852, 262724, 309800, 243232)
QUERY_TAXA = (3702, 39947, 4577, 3847, 4113, 3880, 3218, 3055, 559292, 284812, 6239, 7227,
              7955, 9913, 9031, 9606, 10090, 10116)
COHORTS = {
    'hydroxyacid_dehydrogenases': ('IPR001236', ('1.1.1.27', '1.1.1.37')),
    'isocitrate_like_dehydrogenases': ('IPR019818', ('1.1.1.42', '1.1.1.85')),
    'aminotransferases': ('IPR015424', ('2.6.1.1', '2.6.1.2', '2.6.1.9', '2.6.1.11', '2.6.1.13', '2.6.1.57')),
    'aspartase_fumarase_lyases': ('IPR000362', ('4.2.1.2', '4.3.2.1', '4.3.2.2')),
    'medium_chain_alcohol_dehydrogenases': ('IPR002328', ('1.1.1.1', '1.1.1.103')),
}
LOCAL_OPEN, LOCAL_EXTEND = -0.15, -0.02


def _fetch(label):
    taxa = ' OR '.join(f'organism_id:{taxon}' for taxon in (*REFERENCE_TAXA, *QUERY_TAXA))
    query = {'query': f'reviewed:true AND ec:{label} AND ({taxa})', 'format': 'tsv', 'size': 250,
             'fields': 'accession,protein_name,organism_id,length,xref_interpro,ec,version,sequence_version,sequence'}
    url = 'https://rest.uniprot.org/uniprotkb/search?' + urlencode(query)
    request = Request(url, headers={'User-Agent': 'iri-constructs/related-housekeeping-families'})
    with urlopen(request, timeout=20) as response:
        data = response.read(1_000_001)
        release = response.headers.get('X-UniProt-Release')
        if response.headers.get('Link') or not release or len(data) > 1_000_000:
            raise ValueError(f'{label}: incomplete or unexpected UniProt response')
    return data, {'url': url, 'release': release, 'sha256': sha(data)}


def _candidates(data, label, interpro, excluded):
    for row in csv.DictReader(io.StringIO(data.decode()), delimiter='\t'):
        sequence = row['Sequence']
        if (row['EC number'] != label or interpro not in row['InterPro'].split(';') or
                row['Protein names'].lower().startswith(('putative', 'probable', 'uncharacterized')) or
                not 100 <= len(sequence) <= 750 or len(sequence) != int(row['Length']) or set(sequence) - ALPHABET):
            excluded['annotation_or_sequence_rules'] += 1
            continue
        taxon = int(row['Organism (ID)'])
        if taxon not in (*REFERENCE_TAXA, *QUERY_TAXA):
            raise ValueError(f'{label}: unexpected taxon')
        yield {'accession': row['Entry'], 'entry_version': int(row['Entry version']),
               'sequence_version': int(row['Sequence version']), 'sha256': sha(sequence.encode('ascii')),
               'taxon': taxon, 'length': len(sequence), 'label': label}, sequence


def acquire(folder):
    folder.mkdir(parents=True)
    previous = json.loads(Path(__file__).with_name('panel.freeze.json').read_text())
    previous_queries = [row for row in previous['entries'] if row['role'] == 'query']
    prior_ids = {row['accession'] for row in previous_queries}
    prior_hashes = {row['sha256'] for row in previous_queries}
    refs, queries, sources, excluded = [], [], {}, Counter()
    for cohort, (interpro, labels) in COHORTS.items():
        for label in labels:
            data, sources[label] = _fetch(label)
            (folder / f'{label}.tsv').write_bytes(data)
            time.sleep(0.35)
            usable = [(dict(entry, cohort=cohort, interpro=interpro), sequence)
                      for entry, sequence in _candidates(data, label, interpro, excluded)]
            pool = [pair for pair in usable if pair[0]['taxon'] in REFERENCE_TAXA]
            if not pool:
                raise ValueError(f'{label}: missing a reference for the frozen cohort')
            entry, sequence = min(pool, key=lambda pair: (REFERENCE_TAXA.index(pair[0]['taxon']), pair[0]['accession']))
            refs.append((dict(entry, role='reference'), sequence))
            queries.extend((dict(entry, role='query'), seq) for entry, seq in usable if entry['taxon'] in QUERY_TAXA)

    seen_hashes = {row['sha256'] for row, _ in refs}
    if len(seen_hashes) != len(refs):
        raise ValueError('Reference labels share an identical sequence')
    seen_taxa = set()
    selected = list(refs)
    for row, sequence in sorted(queries, key=lambda pair: (pair[0]['label'], pair[0]['taxon'], pair[0]['accession'])):
        taxon_label = row['taxon'], row['label']
        if row['accession'] in prior_ids or row['sha256'] in prior_hashes:
            excluded['v0_query_overlap'] += 1
        elif row['sha256'] in seen_hashes or taxon_label in seen_taxa:
            excluded['repeated_sequence_or_taxon'] += 1
        else:
            selected.append((row, sequence))
            seen_hashes.add(row['sha256'])
            seen_taxa.add(taxon_label)
    if len({source['release'] for source in sources.values()}) != 1:
        raise ValueError('UniProt release changed during acquisition')
    freeze = {'schema': SCHEMA, 'retrieved_on': datetime.now(ZoneInfo('America/New_York')).date().isoformat(),
              'source': 'Reviewed UniProtKB/Swiss-Prot', 'terms': 'https://www.uniprot.org/help/license',
              'sources_by_label': sources, 'entries': [entry for entry, _ in selected], 'excluded': dict(excluded)}
    write_json(folder / 'panel.freeze.json', freeze)
    write_json(folder / 'sequences.private.json', {entry['accession']: seq for entry, seq in selected})
    print(f'Frozen {len(selected) - len(refs)} queries and {len(refs)} references across {len(COHORTS)} cohorts')


def _residue_local_score(query, reference, aligner):
    """Local alignment over centred cosine scores, per shorter-protein residue."""
    import numpy as np

    cosine = query @ reference.T
    centred = cosine - cosine.mean(axis=0) - cosine.mean(axis=1, keepdims=True) + cosine.mean()
    n, m = centred.shape
    # C aligner operates on integer residue positions with a position-specific matrix.
    matrix = np.zeros((n + m, n + m), dtype='float64')
    matrix[:n, n:] = centred
    aligner.substitution_matrix = matrix
    aligned = aligner.score(np.arange(n, dtype='int32'), np.arange(n, n + m, dtype='int32'))
    return aligned / min(n, m)


def _encode(sequences):
    import numpy as np
    import torch
    from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device

    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    adapter = LocalAdapter('esm2', prepare_device('cpu', 20261122), ESM2_REVISION, offline=True)
    pooled, residues = {}, {}
    try:
        if adapter.metadata['weight_sha256'] != MODEL_HASH:
            raise ValueError('Pinned ESM-2 weight hash does not match')
        for accession, sequence in sequences.items():
            states = np.asarray(adapter.residues(sequence), dtype=np.float32)
            vector = states.mean(axis=0)
            pooled[accession] = vector / np.linalg.norm(vector)
            residues[accession] = states / np.linalg.norm(states, axis=1, keepdims=True)
        return pooled, residues, adapter.metadata
    finally:
        adapter.close()


def _summarize(rows):
    methods = ('alignment', 'alignment_ref_normalized', 'mean_cosine', 'residue_local')
    hit_counts = {method: sum(row['matches'][method] for row in rows) for method in methods}
    return {'queries': len(rows), 'labels': len({row['label'] for row in rows}),
            'strong_competing_label': sum(row['strong_competing_label'] for row in rows),
            'hit_at_1': hit_counts,
            'paired_vs_alignment': {method: {
                'method_only': sum(row['matches'][method] and not row['matches']['alignment'] for row in rows),
                'alignment_only': sum(row['matches']['alignment'] and not row['matches'][method] for row in rows),
            } for method in methods[1:]}}


def score(folder, output):
    from Bio.Align import PairwiseAligner

    frozen_bytes = (folder / 'panel.freeze.json').read_bytes()
    freeze = json.loads(frozen_bytes)
    sequences = json.loads((folder / 'sequences.private.json').read_text())
    entries = {entry['accession']: entry for entry in freeze['entries']}
    if (freeze['schema'] != SCHEMA or len(entries) != len(freeze['entries']) or sequences.keys() != entries.keys()
            or any(sha(sequence.encode('ascii')) != entries[key]['sha256'] for key, sequence in sequences.items())):
        raise ValueError('Private sequences do not match the frozen panel')

    began = time.perf_counter()
    pooled, residues, model = _encode(sequences)
    encoded = time.perf_counter()
    aligner = PairwiseAligner(mode='local', open_gap_score=LOCAL_OPEN, extend_gap_score=LOCAL_EXTEND)
    cases, alignment_seconds, local_seconds, pairs = [], 0.0, 0.0, 0
    for cohort in COHORTS:
        queries = {key: sequences[key] for key, row in entries.items() if row['cohort'] == cohort and row['role'] == 'query'}
        refs = {key: sequences[key] for key, row in entries.items() if row['cohort'] == cohort and row['role'] == 'reference'}
        labels = {entries[key]['label']: key for key in refs}
        if len(refs) != len(COHORTS[cohort][1]) or set(labels) != set(COHORTS[cohort][1]):
            raise ValueError(f'{cohort}: requires the frozen one-reference-per-label panel')
        pairs += len(queries) * len(refs)
        started = time.perf_counter()
        scores, distances = _alignment_comparison(queries, refs)
        alignment_seconds += time.perf_counter() - started
        for qid in queries:
            label = entries[qid]['label']
            target = labels[label]
            positive = distances[qid][target]
            identity = positive['identity'] if positive['coverage'] >= MIN_COVERAGE else None
            competing = any(key != target and value['coverage'] >= MIN_COVERAGE and value['identity'] >= identity
                            for key, value in distances[qid].items()) if identity is not None else False
            started = time.perf_counter()
            local = {key: _residue_local_score(residues[qid], residues[key], aligner) for key in refs}
            local_seconds += time.perf_counter() - started
            rows = {**{method: values[qid] for method, values in scores.items()},
                    'mean_cosine': {key: float(pooled[qid] @ pooled[key]) for key in refs}, 'residue_local': local}
            cases.append({'query_id': qid, 'cohort': cohort, 'label': label, 'tier': _tier(identity),
                          'strong_competing_label': competing,
                          'matches': {method: min(row, key=lambda key: (-row[key], key)) == target for method, row in rows.items()}})

    report = {'schema': SCHEMA, 'panel_sha256': sha(frozen_bytes), 'model': model,
              'environment': {package: version(package) for package in ('biopython', 'numpy', 'torch', 'transformers')},
              'criteria': {'low_identity_below': LOW_IDENTITY, 'high_identity_at_least': HIGH_IDENTITY,
                           'min_coverage_each': MIN_COVERAGE, 'local_gap_open': LOCAL_OPEN, 'local_gap_extend': LOCAL_EXTEND},
              'overall': _summarize(cases),
              'by_tier': {tier: _summarize([row for row in cases if row['tier'] == tier])
                          for tier in ('low', 'intermediate', 'high', 'insufficient_coverage')},
              'by_cohort': {cohort: _summarize([row for row in cases if row['cohort'] == cohort]) for cohort in COHORTS},
              'cost': {'model_encodings': len(sequences), 'reference_pairs': pairs,
                       'encoding_seconds': round(encoded - began, 3), 'alignment_seconds': round(alignment_seconds, 3),
                       'residue_local_seconds': round(local_seconds, 3)}, 'cases': cases}
    write_json(output, report)
    print(json.dumps({key: report[key] for key in ('overall', 'by_tier', 'cost')}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    sub.add_parser('acquire').add_argument('folder', type=Path)
    compare = sub.add_parser('score')
    compare.add_argument('folder', type=Path)
    compare.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.stage == 'acquire':
        acquire(args.folder)
    else:
        score(args.folder, args.out)


if __name__ == '__main__':
    main()
