"""Compare known housekeeping protein labels against a separate reference panel."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import platform
import time
from collections import Counter
from datetime import datetime
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from zoneinfo import ZoneInfo

SCHEMA = 'benign-housekeeping-protein-pilot-1'
REFERENCE_TAXA = (83333, 224308)  # E. coli K-12 and B. subtilis 168
QUERY_TAXA = (3702, 39947, 559292, 6239, 7227, 7955, 9606)
GROUPS = {
    '1.1.1.27': 'IPR001236',  # lactate dehydrogenase
    '1.1.1.37': 'IPR001236',  # malate dehydrogenase
    '1.1.1.42': 'IPR019818',  # isocitrate dehydrogenase
    '1.1.1.85': 'IPR004429',  # isopropylmalate dehydrogenase
    '1.8.1.7': 'IPR001100',   # glutathione reductase
    '2.2.1.1': 'IPR005475',   # transketolase
    '2.6.1.1': 'IPR004839',   # aspartate aminotransferase
    '2.6.1.2': 'IPR004839',   # alanine aminotransferase
    '2.6.1.9': 'IPR001917',   # histidinol-phosphate aminotransferase
    '2.7.2.3': 'IPR001576',   # phosphoglycerate kinase
    '2.7.4.3': 'IPR000850',   # adenylate kinase
    '4.2.1.11': 'IPR000941',  # enolase
    '5.3.1.1': 'IPR000652',   # triosephosphate isomerase
    '5.3.1.6': 'IPR004788',   # ribose-5-phosphate isomerase A
    '5.3.1.9': 'IPR001672',   # glucose-6-phosphate isomerase
}
ALPHABET = frozenset('ACDEFGHIKLMNPQRSTVWY')
MODEL_HASH = '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189'
MIN_COVERAGE = 0.70
LOW_IDENTITY = 0.40
HIGH_IDENTITY = 0.60


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, document):
    with path.open('x', encoding='utf-8') as handle:
        json.dump(document, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def _get_group(ec):
    taxa = ' OR '.join(f'organism_id:{taxon}' for taxon in (*REFERENCE_TAXA, *QUERY_TAXA))
    url = 'https://rest.uniprot.org/uniprotkb/search?' + urlencode({
        'query': f'reviewed:true AND ec:{ec} AND ({taxa})', 'format': 'tsv', 'size': 100,
        'fields': 'accession,protein_name,organism_id,length,xref_interpro,ec,version,sequence_version,sequence',
    })
    with urlopen(Request(url, headers={'User-Agent': 'iri-constructs/benign-protein-pilot'}), timeout=20) as response:
        if response.headers.get('Link'):
            raise ValueError(f'{ec}: source returned more than 100 records')
        data = response.read(1_000_001)
        release = response.headers.get('X-UniProt-Release')
    if not release or len(data) > 1_000_000:
        raise ValueError(f'{ec}: unexpected UniProt response')
    return data, {'url': url, 'release': release, 'sha256': sha(data)}


def acquire(out):
    out.mkdir(parents=True)
    entries, sequences, sources, excluded = [], {}, {}, Counter()
    for ec, interpro in GROUPS.items():
        data, sources[ec] = _get_group(ec)
        (out / f'{ec}.tsv').write_bytes(data)
        rows = csv.DictReader(io.StringIO(data.decode()), delimiter='\t')
        taxa_seen, signatures = set(), set()
        # Reference first, then query. Selection precedes calculation of either retrieval score.
        order = sorted(rows, key=lambda row: (int(row['Organism (ID)']) not in REFERENCE_TAXA, row['Entry']))
        for row in order:
            seq = row['Sequence']
            taxon = int(row['Organism (ID)'])
            if (row['EC number'] != ec or interpro not in row['InterPro'].split(';') or
                    row['Protein names'].lower().startswith(('putative', 'probable', 'uncharacterized')) or
                    not (100 <= len(seq) <= 750) or set(seq) - ALPHABET or len(seq) != int(row['Length'])):
                excluded['annotation_or_sequence_rules'] += 1
                continue
            if taxon not in (*REFERENCE_TAXA, *QUERY_TAXA):
                raise ValueError(f'{ec}: unexpected taxon from UniProt')
            digest = sha(seq.encode('ascii'))
            if taxon in taxa_seen or digest in signatures:
                excluded['repeated_taxon_or_sequence_within_label'] += 1
                continue
            accession = row['Entry']
            if accession in sequences:
                raise ValueError(f'{accession}: more than one EC label')
            taxa_seen.add(taxon)
            signatures.add(digest)
            entries.append({'accession': accession, 'entry_version': int(row['Entry version']),
                            'sequence_version': int(row['Sequence version']), 'sha256': digest, 'taxon': taxon,
                            'role': 'reference' if taxon in REFERENCE_TAXA else 'query', 'length': len(seq),
                            'label': ec, 'interpro': interpro})
            sequences[accession] = seq
        time.sleep(0.25)
    if len({source['release'] for source in sources.values()}) != 1:
        raise ValueError('UniProt release changed during acquisition')
    freeze = {'schema': SCHEMA, 'retrieved_on': datetime.now(ZoneInfo('America/New_York')).date().isoformat(),
              'source': 'Reviewed UniProtKB/Swiss-Prot', 'terms': 'https://www.uniprot.org/help/license',
              'sources_by_label': sources, 'entries': entries, 'excluded': dict(excluded)}
    write_json(out / 'sequences.private.json', sequences)
    write_json(out / 'panel.freeze.json', freeze)
    print(f'Frozen {len(entries)} proteins across {len(GROUPS)} labels in {out}')


def _alignment_comparison(queries, references):
    from Bio.Align import PairwiseAligner, substitution_matrices

    aligner = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)
    raw, normalized, distances = {}, {}, {}
    self_scores = {rid: aligner.score(ref, ref) for rid, ref in references.items()}
    for qid, query in queries.items():
        raw[qid], normalized[qid], distances[qid] = {}, {}, {}
        for rid, reference in references.items():
            alignments = iter(aligner.align(query, reference))
            hit = next(alignments, None)
            raw[qid][rid] = hit.score if hit is not None else 0.0
            normalized[qid][rid] = raw[qid][rid] / self_scores[rid]
            pairs, identical = 0, 0
            if hit is not None:
                for (a, b), (c, d) in zip(*hit.aligned, strict=True):
                    pairs += int(b - a)
                    identical += sum(x == y for x, y in zip(query[a:b], reference[c:d], strict=True))
            distances[qid][rid] = {'identity': identical / pairs if pairs else 0.0,
                                   'coverage': min(pairs / len(query), pairs / len(reference))}
    return {'alignment': raw, 'alignment_ref_normalized': normalized}, distances


def _embeddings(queries, references):
    import torch
    from iri_constructs.metrics import cosine, mean_pool
    from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device

    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    model = LocalAdapter('esm2', prepare_device('cpu', 20261122), ESM2_REVISION, offline=True)
    try:
        if model.metadata['weight_sha256'] != MODEL_HASH:
            raise ValueError('Pinned ESM-2 weight hash does not match')
        query_vecs = {qid: mean_pool(model.residues(query)) for qid, query in queries.items()}
        ref_vecs = {rid: mean_pool(model.residues(ref)) for rid, ref in references.items()}
        rows = {qid: {rid: cosine(query, ref) for rid, ref in ref_vecs.items()}
                for qid, query in query_vecs.items()}
        return rows, model.metadata
    finally:
        model.close()


def _tier(best_identity):
    if best_identity is None:
        return 'insufficient_coverage'
    if best_identity < LOW_IDENTITY:
        return 'low'
    if best_identity < HIGH_IDENTITY:
        return 'intermediate'
    return 'high'


def _summary(cases, methods):
    return {
        'queries': len(cases), 'labels': len({case['label'] for case in cases}),
        'strong_competing_label': sum(case['strong_competing_label'] for case in cases),
        'hits': {method: {'hit_at_1': sum(case['ranks'][method] == 1 for case in cases),
                          'hit_at_5': sum(case['ranks'][method] <= 5 for case in cases)} for method in methods},
        'paired_hit_at_1': {
            'esm2_only': sum(c['ranks']['esm2'] == 1 and c['ranks']['alignment'] != 1 for c in cases),
            'alignment_only': sum(c['ranks']['alignment'] == 1 and c['ranks']['esm2'] != 1 for c in cases),
        },
    }


def score(folder, output):
    from importlib.metadata import version
    from iri_constructs.component_eval import unit_retrieval

    metadata = (folder / 'panel.freeze.json').read_bytes()
    freeze = json.loads(metadata)
    if freeze['schema'] != SCHEMA:
        raise ValueError('Invalid panel schema')
    seqs = json.loads((folder / 'sequences.private.json').read_text())
    entries = {row['accession']: row for row in freeze['entries']}
    if seqs.keys() != entries.keys() or len(freeze['entries']) != len(entries) or any(
            sha(seq.encode('ascii')) != entries[key]['sha256'] for key, seq in seqs.items()):
        raise ValueError('Sequence snapshot differs from its freeze')
    queries = {key: seqs[key] for key, e in entries.items() if e['role'] == 'query'}
    refs = {key: seqs[key] for key, e in entries.items() if e['role'] == 'reference'}
    labels = {key: e['label'] for key, e in entries.items()}
    ref_labels = {key: labels[key] for key in refs}
    if (not queries or set(labels[q] for q in queries) - set(ref_labels.values()) or
            set(queries.values()) & set(refs.values())):
        raise ValueError('Panel lacks matching labels or has exact query/reference copies')

    started = time.perf_counter()
    scores, distances = _alignment_comparison(queries, refs)
    alignment_seconds = time.perf_counter() - started
    started = time.perf_counter()
    scores['esm2'], model = _embeddings(queries, refs)
    esm_seconds = time.perf_counter() - started
    cases = []
    for qid in queries:
        label = labels[qid]
        valid = {rid: d['identity'] for rid, d in distances[qid].items() if d['coverage'] >= MIN_COVERAGE}
        positive = max((v for rid, v in valid.items() if labels[rid] == label), default=None)
        competing = max((v for rid, v in valid.items() if labels[rid] != label), default=None)
        cases.append({
            'query_id': qid, 'label': label, 'tier': _tier(positive), 'best_same_label_identity': positive,
            'strong_competing_label': positive is not None and competing is not None and competing >= positive,
            'ranks': {method: unit_retrieval(rows[qid], {label}, ref_labels, 5)['first_same_family_rank']
                      for method, rows in scores.items()},
        })
    report = {
        'schema': SCHEMA, 'panel_sha256': sha(metadata), 'model': model, 'references': len(refs),
        'settings': {'minimum_coverage_each': MIN_COVERAGE, 'low_identity_below': LOW_IDENTITY,
                     'high_identity_at_least': HIGH_IDENTITY, 'alignment': 'local BLOSUM62; open -10, extend -0.5'},
        'environment': {**{p: version(p) for p in ('biopython', 'torch', 'transformers', 'huggingface-hub', 'numpy')},
                        'python': platform.python_version()},
        'cost': {'alignment_seconds': round(alignment_seconds, 3), 'esm2_seconds': round(esm_seconds, 3),
                 'pairwise_alignments': len(queries) * len(refs), 'esm2_calls': len(queries) + len(refs)},
        'overall': _summary(cases, scores),
        'by_tier': {tier: _summary([c for c in cases if c['tier'] == tier], scores)
                    for tier in ('low', 'intermediate', 'high', 'insufficient_coverage')},
        'cases': cases,
    }
    write_json(output, report)
    print(json.dumps({key: report[key] for key in ('overall', 'by_tier', 'cost')}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    fetch = sub.add_parser('acquire')
    fetch.add_argument('folder', type=Path, help='new local folder for source snapshots and panel freeze')
    compare = sub.add_parser('score')
    compare.add_argument('folder', type=Path, help='previously acquired panel folder')
    compare.add_argument('--out', type=Path, required=True, help='new local report path')
    args = parser.parse_args()
    if args.stage == 'acquire':
        acquire(args.folder)
    else:
        score(args.folder, args.out)


if __name__ == '__main__':
    main()
