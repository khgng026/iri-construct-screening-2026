"""Run a paired validation pilot on the documented reporter constructs in reporter_panel.json."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .verify import SCHEMA

MODEL_HASH = '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189'
MAX_RESIDUES = 1022
MIN_RESIDUES = 30
AMINO_ACIDS = frozenset('ACDEFGHIKLMNPQRSTVWY')


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, value):
    with path.open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(panel, output):
    """Reacquire the exact deposited records and reviewed protein references."""
    def fetch(url, expected):
        with urlopen(Request(url, headers={'User-Agent': 'iri-constructs/reporter-pilot'}), timeout=20) as response:
            data = response.read(200_001)
        if len(data) > 200_000 or hashlib.sha256(data).hexdigest() != expected:
            raise ValueError('Source snapshot differs from the frozen panel')
        return data

    root = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
    records = {}
    for row in panel['constructs']:
        query = urlencode(dict(db='nuccore', id=row['accession'], rettype='gb', retmode='text', tool='iri_constructs'))
        records[row['accession']] = fetch(root + query, row['sha256'])
        time.sleep(0.4)
    reference_source = panel['reference_source']
    raw = fetch(reference_source['url'], reference_source['response_sha256'])
    sequences = {r['Entry']: r['Sequence'] for r in csv.DictReader(io.StringIO(raw.decode()), delimiter='\t')}
    references = [{**row, 'sequence': sequences[row['accession']]} for row in panel['references']]
    out_bytes = (json.dumps(references, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
    if hashlib.sha256(out_bytes).hexdigest() != reference_source['prepared_sha256']:
        raise ValueError('Reference panel differs from the frozen input')
    (output / 'sources').mkdir(parents=True)
    for accession, data in records.items():
        (output / 'sources' / f'{accession}.gb').write_bytes(data)
    (output / 'sources' / 'reference-snapshot.tsv').write_bytes(raw)
    (output / 'references.private.json').write_bytes(out_bytes)


def prepare(sources, panel, output):
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from iri_constructs import process_inputs

    fasta, constructs, annotation_units = [], [], {}
    for item in panel['constructs']:
        accession, key = item['accession'], item['id']
        source = sources / f'{accession}.gb'
        if digest(source) != item['sha256']:
            raise ValueError(f'{accession}: GenBank snapshot differs from the panel')
        sequence = SeqIO.read(source, 'genbank')
        result = process_inputs([source])
        if sequence.id != accession or len(result.ok_records) != 1:
            raise ValueError(f'{accession}: source cannot be parsed as one expected record')
        validated = result.ok_records[0]
        features = {f.source_feature_index: f for f in validated.features}
        units = {u.feature_id: u for u in validated.units}
        components, proteins = [], []
        for label in item['labels']:
            feature = features[label['feature_index']]
            unit = units[feature.feature_id]
            if unit.extraction_state != 'exported' or feature.names['product'] != label['source_product']:
                raise ValueError(f'{accession}: documented component no longer matches the answer key')
            locus = {'id': f'{key}-{label["feature_index"]}',
                     'strand': {'+': 1, '-': -1}[feature.strand],
                     'spans': [[s['start0'], s['end0']] for s in feature.segments]}
            components.append({**locus, 'family': label['family']})
            proteins.append({**locus, 'protein': unit.sequence})
        documented_cds = {f.source_feature_index for f in validated.features if f.source_feature_type == 'CDS'}
        if documented_cds != {label['feature_index'] for label in item['labels']}:
            raise ValueError(f'{accession}: panel must account for every documented CDS in this pilot')
        fasta.append(SeqRecord(sequence.seq.upper(), id=key, name=key, description=''))
        constructs.append({'id': key, 'archetype': item['archetype'], 'length': len(sequence), 'components': components})
        annotation_units[key] = proteins

    output.mkdir(parents=True)
    (output / 'method_inputs').mkdir()
    (output / 'evaluator').mkdir()
    SeqIO.write(fasta, output / 'method_inputs' / 'constructs.fasta', 'fasta')
    write(output / 'evaluator' / 'answer-key.json', {'schema': SCHEMA, 'constructs': constructs})
    write(output / 'evaluator' / 'annotation-units.json', annotation_units)


def _candidates(sequence, start_aware):
    from Bio.Seq import Seq

    units = []
    length = len(sequence)
    for strand, dna in ((1, Seq(sequence)), (-1, Seq(sequence).reverse_complement())):
        for offset in range(3):
            translated = str(dna[offset:length - (length - offset) % 3].translate())
            start = 0
            for stop_free in translated.split('*'):
                end = start + len(stop_free)
                # Both arms use the same translation and model length handling.
                # The alternative starts at the first methionine of each block.
                index = stop_free.find('M') if start_aware else 0
                if index >= 0 and not set(stop_free) - AMINO_ACIDS:
                    for left in range(start + index, end, MAX_RESIDUES):
                        right = min(left + MAX_RESIDUES, end)
                        if right - left < MIN_RESIDUES:
                            continue
                        a, b = offset + left * 3, offset + right * 3
                        span = [a, b] if strand == 1 else [length - b, length - a]
                        units.append({'id': f'q{len(units) + 1}', 'strand': strand, 'spans': [span],
                                      'protein': translated[left:right]})
                start = end + 1
    return sorted(units, key=lambda unit: (-len(unit['protein']), unit['spans'], -unit['strand'], unit['id']))


def _predictions(selected, model, references, limit):
    from iri_constructs.metrics import cosine, mean_pool

    rows = []
    for query in selected:
        vector = mean_pool(model.residues(query['protein']))
        scores = {}
        for family, reference in references:
            score = cosine(vector, reference)
            scores[family] = max(scores.get(family, -1.0), score)
        rows.append((query['id'], sorted(scores, key=lambda label: (-scores[label], label))))
    # Each model input can contribute one family per round. A repeated family
    # may refer to a different coordinate and is scored one to one by the verifier.
    predictions = []
    for rank in range(max((len(labels) for _, labels in rows), default=0)):
        for qid, labels in rows:
            if rank < len(labels):
                predictions.append({'query_id': qid, 'family': labels[rank]})
                if len(predictions) == limit:
                    return predictions
    return predictions


def run(method, inputs, references, protocol, output):
    import torch
    from Bio import SeqIO
    from iri_constructs.metrics import mean_pool
    from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device

    if digest(references) != protocol['reference_panel_sha256']:
        raise ValueError('Reference panel differs from the fixed characterization protocol')
    input_records = (read(inputs) if method == 'annotation_reference' else
                     {r.id: str(r.seq) for r in SeqIO.parse(inputs, 'fasta')})
    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    model = LocalAdapter('esm2', prepare_device('cpu', 20261122), ESM2_REVISION, MAX_RESIDUES, offline=True)
    try:
        if model.metadata['weight_sha256'] != MODEL_HASH:
            raise ValueError('Pinned model weights differ')
        vectors = [(ref['family'], mean_pool(model.residues(ref['sequence']))) for ref in read(references)]
        run_data = {'schema': SCHEMA, 'method': method, 'characterization_id': protocol['characterization_id'],
                    'tiers': [], 'input_sha256': digest(inputs), 'model': model.metadata}
        for cap in protocol['query_caps']:
            traces = []
            for key, original in sorted(input_records.items()):
                started = time.perf_counter()
                units = (sorted(original, key=lambda u: (-len(u['protein']), u['id']))
                         if method == 'annotation_reference' else _candidates(original, method == 'start_aware'))
                usable = [unit for unit in units if MIN_RESIDUES <= len(unit['protein']) <= MAX_RESIDUES]
                selected = usable[:cap]
                prepared = time.perf_counter()
                predictions = _predictions(selected, model, vectors, protocol['prediction_limit'])
                traces.append({'id': key, 'queries': [{**{k: q[k] for k in ('id', 'strand', 'spans')},
                                                       'input_residues': len(q['protein'])} for q in selected],
                               'predictions': predictions, 'candidate_count': len(usable), 'status': 'ok',
                               'preprocess_seconds': prepared - started,
                               'characterize_seconds': time.perf_counter() - prepared})
            run_data['tiers'].append({'query_cap': cap, 'constructs': traces})
        write(output, run_data)
    finally:
        model.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    for stage in ('acquire', 'prepare'):
        command = sub.add_parser(stage)
        if stage == 'prepare':
            command.add_argument('--sources', type=Path, required=True)
        command.add_argument('--panel', type=Path, default=Path(__file__).with_name('reporter_panel.json'))
        command.add_argument('--out', type=Path, required=True)
    execute = sub.add_parser('run')
    execute.add_argument('--method', choices=('minimal', 'start_aware', 'annotation_reference'), required=True)
    execute.add_argument('--inputs', type=Path, required=True)
    execute.add_argument('--references', type=Path, required=True)
    execute.add_argument('--protocol', type=Path, required=True)
    execute.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.stage == 'acquire':
        acquire(read(args.panel), args.out)
    elif args.stage == 'prepare':
        prepare(args.sources, read(args.panel), args.out)
    else:
        run(args.method, args.inputs, args.references, read(args.protocol), args.out)


if __name__ == '__main__':
    main()
