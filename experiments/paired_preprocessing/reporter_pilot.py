"""Run a paired validation pilot on the documented reporter constructs in reporter_panel.json."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import time
from contextlib import contextmanager
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from .verify import SCHEMA

MODEL_HASH = '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189'
MAX_RESIDUES = 1022
MIN_RESIDUES = 30
AMINO_ACIDS = frozenset('ACDEFGHIKLMNPQRSTVWY')
PYRODIGAL_VERSION = '3.6.3'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path, value):
    with path.open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def acquire(panel, output, *, canonical_record_end=False):
    """Reacquire frozen deposited records and protein references."""
    def fetch(url, expected, canonical_end=False):
        with urlopen(Request(url, headers={'User-Agent': 'iri-constructs/reporter-pilot'}), timeout=20) as response:
            data = response.read(200_001)
        if canonical_end and len(data) <= 200_000:
            data = data.rstrip(b'\n') + b'\n'
        if len(data) > 200_000 or hashlib.sha256(data).hexdigest() != expected:
            raise ValueError('Source snapshot differs from the frozen panel')
        return data

    root = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
    records = {}
    for row in panel['constructs']:
        query = urlencode(dict(db='nuccore', id=row['accession'], rettype='gb', retmode='text', tool='iri_constructs'))
        records[row['accession']] = fetch(root + query, row['sha256'], canonical_record_end)
        time.sleep(0.4)
    reference_source = panel['reference_source']
    if reference_source.get('provider') == 'SCOPe ASTRAL':
        snapshots, sequences = {}, {}
        for row in panel['references']:
            accession = row['accession']
            raw = fetch(row['url'], row['source_sha256'])
            description, *lines = raw.decode('ascii').splitlines()
            sequence = ''.join(lines).upper()
            if (description.split()[:1] != [f'>{accession}']
                    or hashlib.sha256(sequence.encode('ascii')).hexdigest() != row['sequence_sha256']):
                raise ValueError(f'{accession}: ASTRAL record differs from the frozen panel')
            snapshots[f'{accession}.fasta'] = raw
            sequences[accession] = sequence
            time.sleep(0.85)
    else:
        raw = fetch(reference_source['url'], reference_source['response_sha256'])
        snapshots = {'reference-snapshot.tsv': raw}
        sequences = {r['Entry']: r['Sequence'] for r in csv.DictReader(io.StringIO(raw.decode()), delimiter='\t')}
    references = [{**row, 'sequence': sequences[row['accession']]} for row in panel['references']]
    out_bytes = (json.dumps(references, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
    if hashlib.sha256(out_bytes).hexdigest() != reference_source['prepared_sha256']:
        raise ValueError('Reference panel differs from the frozen input')
    (output / 'sources').mkdir(parents=True)
    for accession, data in records.items():
        (output / 'sources' / f'{accession}.gb').write_bytes(data)
    for filename, data in snapshots.items():
        (output / 'sources' / filename).write_bytes(data)
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


def _append_chunks(units, protein, coding_start, strand):
    """Use the translated residues' span, including when a candidate is split for the model."""
    for offset in range(0, len(protein), MAX_RESIDUES):
        chunk = protein[offset:offset + MAX_RESIDUES]
        if len(chunk) < MIN_RESIDUES:
            continue
        start = coding_start + strand * offset * 3
        end = start + strand * len(chunk) * 3
        units.append({'id': f'q{len(units) + 1}', 'strand': strand, 'spans': [sorted((start, end))],
                      'protein': chunk})


def _length_order(unit):
    return -len(unit['protein']), unit['spans'], -unit['strand'], unit['id']


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
                    anchor = offset + (start + index) * 3
                    _append_chunks(units, stop_free[index:], anchor if strand == 1 else length - anchor, strand)
                start = end + 1
    return sorted(units, key=_length_order)


def _rank_by_start(units):
    """Prefer an early start without trimming or dropping any candidate."""
    def score(unit):
        protein = unit['protein']
        first_start = protein.find('M')
        penalty = MIN_RESIDUES if first_start < 0 else min(first_start, MIN_RESIDUES)
        return len(protein) - penalty

    # Stable ties keep the original length and coordinate order.
    return sorted(units, key=score, reverse=True)


def _gene_caller_candidates(sequence, finder):
    units = []
    for gene in finder.find_genes(sequence.encode('ascii')):
        protein = gene.translate(include_stop=False)
        if not protein or set(protein) - AMINO_ACIDS:
            continue
        # Pyrodigal gives forward coordinates that are one-based and inclusive.
        start, end, strand = gene.begin - 1, gene.end, gene.strand
        remainder = end - start - len(protein) * 3
        if not (0 <= start < end <= len(sequence) and strand in (-1, 1) and 0 <= remainder <= 3):
            raise ValueError('Pyrodigal translation and source coordinates disagree')
        _append_chunks(units, protein, start if strand == 1 else end, strand)
    return sorted(units, key=_length_order)


def _rank_families(reference_scores):
    scores = {}
    for family, score in reference_scores:
        scores[family] = max(scores.get(family, -1.0), score)
    return sorted(scores, key=lambda label: (-scores[label], label))


@contextmanager
def _characterizer(name, references):
    if name == 'blosum62':
        from Bio.Align import PairwiseAligner, substitution_matrices

        aligner = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                                  open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)
        targets = [(ref['family'], ref['sequence'], aligner.score(ref['sequence'], ref['sequence']))
                   for ref in references]
        if any(self_score <= 0 for _, _, self_score in targets):
            raise ValueError('Alignment reference self-score must be positive')

        def align(protein):
            return _rank_families((family, aligner.score(protein, sequence) / self_score)
                                  for family, sequence, self_score in targets)

        yield align, {'alignment': {'matrix': 'BLOSUM62', 'mode': 'local', 'open_gap_score': -10.0,
                                    'extend_gap_score': -0.5, 'normalization': 'reference self-score'}}
        return

    if name != 'esm2':
        raise ValueError('Unknown characterization method')
    import torch
    from iri_constructs.metrics import cosine, mean_pool
    from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device

    torch.set_num_threads(2)
    torch.set_num_interop_threads(1)
    model = LocalAdapter('esm2', prepare_device('cpu', 20261122), ESM2_REVISION, MAX_RESIDUES, offline=True)
    try:
        if model.metadata['weight_sha256'] != MODEL_HASH:
            raise ValueError('Pinned model weights differ')
        vectors = [(ref['family'], mean_pool(model.residues(ref['sequence']))) for ref in references]

        def embed(protein):
            vector = mean_pool(model.residues(protein))
            return _rank_families((family, cosine(vector, reference)) for family, reference in vectors)

        yield embed, {'model': model.metadata}
    finally:
        model.close()


def _predictions(selected, rank_families, limit):
    rows = [(query['id'], rank_families(query['protein'])) for query in selected]
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
    input_sha256 = digest(inputs)
    if (method != 'annotation_reference' and 'sequence_input_sha256' in protocol and
            input_sha256 != protocol['sequence_input_sha256']):
        raise ValueError('Sequence input differs from the fixed characterization protocol')
    reference_sha256 = digest(references)
    if reference_sha256 != protocol['reference_panel_sha256']:
        raise ValueError('Reference panel differs from the fixed characterization protocol')

    finder = None
    if method == 'gene_caller':
        try:
            import pyrodigal
        except ImportError as error:
            raise ValueError(f'The gene_caller arm requires pyrodigal=={PYRODIGAL_VERSION}') from error
        if pyrodigal.__version__ != PYRODIGAL_VERSION:
            raise ValueError(f'The gene_caller arm requires pyrodigal=={PYRODIGAL_VERSION}')
        finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)

    from Bio import SeqIO

    input_records = (read(inputs) if method == 'annotation_reference' else
                     {r.id: str(r.seq) for r in SeqIO.parse(inputs, 'fasta')})
    characterizer = protocol.get('characterizer', 'esm2')
    reference_rows = read(references)
    selected_inputs = hashlib.sha256()
    with _characterizer(characterizer, reference_rows) as (rank_families, metadata):
        run_data = {'schema': SCHEMA, 'method': method, 'characterization_id': protocol['characterization_id'],
                    'characterizer': characterizer, 'reference_panel_sha256': reference_sha256,
                    'reference_count': len(reference_rows), 'tiers': [], 'input_sha256': input_sha256, **metadata}
        if finder is not None:
            run_data['gene_caller'] = {'engine': 'pyrodigal', 'version': PYRODIGAL_VERSION,
                                       'meta': True, 'closed': False, 'mask': False}
        for cap in protocol['query_caps']:
            traces = []
            for key, original in sorted(input_records.items()):
                started = time.perf_counter()
                if method == 'annotation_reference':
                    units = sorted(original, key=lambda u: (-len(u['protein']), u['id']))
                elif finder is not None:
                    units = _gene_caller_candidates(original, finder)
                else:
                    units = _candidates(original, method == 'start_aware')
                    if method == 'start_ranked':
                        units = _rank_by_start(units)
                usable = [unit for unit in units if MIN_RESIDUES <= len(unit['protein']) <= MAX_RESIDUES]
                selected = usable[:cap]
                prepared = time.perf_counter()
                predictions = _predictions(selected, rank_families, protocol['prediction_limit'])
                scored = time.perf_counter()
                # Includes the protein, order and source spans in a single private digest.
                encoded = json.dumps((cap, key, selected), sort_keys=True, separators=(',', ':')).encode('utf-8')
                selected_inputs.update(encoded + b'\n')
                traces.append({'id': key, 'queries': [{**{k: q[k] for k in ('id', 'strand', 'spans')},
                                                       'input_residues': len(q['protein'])} for q in selected],
                               'predictions': predictions, 'candidate_count': len(usable), 'status': 'ok',
                               'preprocess_seconds': prepared - started,
                               'characterize_seconds': scored - prepared})
            run_data['tiers'].append({'query_cap': cap, 'constructs': traces})
        run_data['selected_inputs_sha256'] = selected_inputs.hexdigest()
        write(output, run_data)


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
    execute.add_argument('--method', choices=('minimal', 'start_aware', 'start_ranked',
                                              'gene_caller', 'annotation_reference'), required=True)
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
