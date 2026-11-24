"""Compare protein scoring on the fixed, documented public reporter panel."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from importlib.metadata import version
from pathlib import Path
from statistics import mean

from experiments.benign_protein_pilot.study import MIN_COVERAGE, _alignment_comparison

from .cohort_scorers import Encoder, _normalized_rows, _residue_scorer
from .reporter_pilot import PYRODIGAL_VERSION, _candidates, _gene_caller_candidates, acquire, digest, read, write
from .verify import _covered, _match

PANEL = Path(__file__).with_name('distant_reporter_panel.json')
ALTERNATE = PANEL.with_name('distant_reporter_alternate_references.json')
RENILLA = PANEL.with_name('distant_renilla_panel.json')
GUS_CONTEXT = PANEL.with_name('gus_context_panel.json')
CAPS = (1, 2, 4)
PRIMARY_CAP = 2
MIN_OVERLAP = 0.80
REMOTE_IDENTITY = 0.30
MODELS = ('esm2_8m', 'esm2_150m')


def _signature(sequence):
    return hashlib.sha256(sequence.encode('ascii')).hexdigest()


def _domain_locus(spans, strand, start, end):
    if len(spans) != 1 or not 0 <= start < end or 3 * end > spans[0][1] - spans[0][0]:
        raise ValueError('The archived domain needs a contiguous in-frame source interval')
    left, right = spans[0]
    return [[left + 3 * start, left + 3 * end]] if strand == 1 else [[right - 3 * end, right - 3 * start]]


def load(folder, panel):
    from Bio import SeqIO

    from iri_constructs import process_inputs

    ref_path = folder / 'references.private.json'
    if digest(ref_path) != panel['reference_source']['prepared_sha256']:
        raise ValueError('The reference snapshot differs from the frozen reporter panel')
    refs = read(ref_path)
    records = {}
    for row in panel['constructs']:
        accession = row['accession']
        path = folder / 'sources' / f'{accession}.gb'
        if digest(path) != row['sha256']:
            raise ValueError(f'{accession}: source snapshot differs from the frozen reporter panel')
        record = SeqIO.read(path, 'genbank')
        validated_records = process_inputs([path]).ok_records
        if record.id != accession or len(record) != row['length'] or len(validated_records) != 1:
            raise ValueError(f'{accession}: expected one validated deposited record')
        validated, = validated_records
        features = {feature.source_feature_index: feature for feature in validated.features}
        units = {unit.feature_id: unit for unit in validated.units}
        targets = []
        for target in row['components']:
            feature = features[target['feature_index']]
            unit = units[feature.feature_id]
            sequence = unit.sequence if unit.extraction_state == 'exported' else None
            complete = sequence
            source_feature = record.features[target['feature_index']]
            qualifiers = source_feature.qualifiers
            translation_table = int(qualifiers.get('transl_table', ['1'])[0])
            if 'context_feature_index' in target:
                # The documented reporter segment is internal to a validated fusion CDS.
                context = features[target['context_feature_index']]
                context_unit = units[context.feature_id]
                enclosing = record.features[target['context_feature_index']].location
                location = source_feature.location
                sequence = str(source_feature.extract(record.seq).translate(table=translation_table)).removesuffix('*')
                complete = context_unit.sequence
                offset = (location.start - enclosing.start if location.strand == 1
                          else enclosing.end - location.end)
                if (context_unit.extraction_state != 'exported' or sequence not in complete
                        or qualifiers.get('translation') != [sequence] or qualifiers.get('codon_start', ['1']) != ['1']
                        or context.names['product'] != target['context_product']
                        or _signature(complete) != target['context_sha256'] or location.strand != enclosing.strand
                        or location.strand not in (-1, 1) or len(location.parts) != 1 or len(enclosing.parts) != 1
                        or offset % 3 or complete[offset // 3:offset // 3 + len(sequence)] != sequence
                        or not enclosing.start <= location.start < location.end <= enclosing.end):
                    raise ValueError(f'{accession}: the documented embedded reporter did not validate')
                strand = location.strand
                spans = [[int(location.start), int(location.end)]]
            else:
                strand = {'+': 1, '-': -1}[feature.strand]
                spans = [[s['start0'], s['end0']] for s in feature.segments]
            if (feature.source_feature_type != 'CDS' or not sequence or feature.names['product'] != target['source_product']
                    or _signature(sequence) != target['translation_sha256']):
                raise ValueError(f'{accession}: documented target protein could not be validated')
            case = {**target, 'id': f"{accession}:{target['feature_index']}", 'strand': strand,
                    'spans': spans, 'protein': sequence}
            if 'domain_residue_start0' in target:
                start, end = target['domain_residue_start0'], target['domain_residue_end0']
                region = sequence[start:end]
                if qualifiers.get('codon_start', ['1']) != ['1'] or _signature(region) != target['domain_sha256']:
                    raise ValueError(f'{accession}: the documented domain no longer matches its record')
                domain_spans = _domain_locus(spans, strand, start, end)
                (left, right), = domain_spans
                domain_dna = record.seq[left:right]
                if strand == -1:
                    domain_dna = domain_dna.reverse_complement()
                if str(domain_dna.translate(table=translation_table)) != region:
                    raise ValueError(f'{accession}: the domain locus does not reproduce its protein')
                case.update(protein=region, spans=domain_spans, complete_product=complete,
                            layout=row['layout'], submission_group=row['submission_group'])
            targets.append(case)
        records[accession] = {'dna': str(record.seq).upper(), 'targets': targets}
    return records, refs


def _identities(targets, refs):
    from Bio.Align import PairwiseAligner, substitution_matrices

    positive = defaultdict(list)
    for ref in refs:
        positive[ref['family']].append(ref)
    aligner = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)
    _, local = _alignment_comparison({t['id']: t['protein'] for t in targets},
                                    {r['accession']: r['sequence'] for r in refs})
    identities = {}
    for target in targets:
        matches = positive[target['family']]
        global_hits = [next(iter(aligner.align(target['protein'], ref['sequence']))) for ref in matches]
        global_identity = max(hit.counts().identities / hit.length for hit in global_hits)
        local_hits = [local[target['id']][ref['accession']] for ref in matches]
        identities[target['id']] = {
            'max_same_family_global_identity': global_identity,
            'max_covered_same_family_local_identity': max(
                (hit['identity'] for hit in local_hits if hit['coverage'] >= MIN_COVERAGE), default=None),
            'low_identity': global_identity < REMOTE_IDENTITY,
        }
    return identities


def _top(scores):
    if not scores:
        return None
    best = max(scores.values())
    winners = [family for family, score in scores.items() if score == best]
    return winners[0] if len(winners) == 1 else None


def _rank_reference_scores(scores, families):
    best = {}
    for key, score in scores.items():
        family = families[key]
        best[family] = max(best.get(family, float('-inf')), score)
    return _top(best)


def _profiles(queries, refs):
    from pyfamsa import Aligner, Sequence
    from pyhmmer import easel, plan7

    alphabet = easel.Alphabet.amino()
    background = plan7.Background(alphabet)
    builder = plan7.Builder(alphabet, seed=42)
    aligner = Aligner(threads=1)
    pipeline = plan7.Pipeline(alphabet, background=background, F1=1.0, F2=1.0, F3=1.0,
                              T=-1000.0, domT=-1000.0)
    db = easel.DigitalSequenceBlock(alphabet, [
        easel.TextSequence(name=key.encode(), sequence=sequence).digitize(alphabet)
        for key, sequence in queries.items()
    ])
    groups = defaultdict(list)
    for ref in refs:
        groups[ref['family']].append(Sequence(ref['accession'].encode(), ref['sequence'].encode()))
    hits = {key: {} for key in queries}
    for family, group in groups.items():
        aligned = [easel.TextSequence(name=row.id, sequence=row.sequence.decode()) for row in aligner.align(group)]
        msa = easel.TextMSA(name=family.encode(), sequences=aligned).digitize(alphabet)
        hmm, _, _ = builder.build_msa(msa, background)
        for hit in pipeline.search_hmm(hmm, db):
            key = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
            hits[key][family] = hit.score
    return {key: _top(scores) for key, scores in hits.items()}


def _predictions(queries, refs):
    import numpy as np

    families = {ref['accession']: ref['family'] for ref in refs}
    references = {ref['accession']: ref['sequence'] for ref in refs}
    raw, _ = _alignment_comparison(queries, references)
    predictions = {method: {key: _rank_reference_scores(scores, families) for key, scores in hits.items()}
                   for method, hits in raw.items()}
    predictions['profile'] = _profiles(queries, refs)
    model_metadata = {}
    for name in MODELS:
        encoder = Encoder(name, 2)
        model_metadata[name] = encoder.metadata
        pooled, residues = {}, {}
        for key, sequence in {**references, **queries}.items():
            states = encoder.states(sequence)
            vector = states.mean(axis=0)
            pooled[key] = vector / np.linalg.norm(vector)
            residues[key] = _normalized_rows(states)
        del encoder
        scorer = _residue_scorer()
        predictions[f'{name}_mean'] = {key: _rank_reference_scores(
            {ref: float(pooled[key] @ pooled[ref]) for ref in references}, families) for key in queries}
        predictions[f'{name}_local'] = {key: _rank_reference_scores(
            {ref: scorer(residues[key], residues[ref]) for ref in references}, families) for key in queries}
    return predictions, model_metadata


def _compatible(target, query):
    return _covered((target['strand'], target['spans']), (query['strand'], query['spans'])) >= MIN_OVERLAP


def _score_targets(targets, pool, selected, rankings):
    allowed = [{t['id'] for t in targets if _compatible(t, q)} for q in selected]
    statuses = {}
    for method, calls in rankings.items():
        edges = [[t['id'] for t in targets if t['id'] in locations and t['family'] == calls[_signature(q['protein'])]]
                 for q, locations in zip(selected, allowed, strict=True)]
        hits = _match(edges)
        statuses[method] = {t['id']: t['id'] in hits for t in targets}
    covered = set().union(*allowed)
    details = {}
    for target in targets:
        key = target['id']
        details[key] = {
            'pool_covering': any(_compatible(target, q) for q in pool),
            'selected_covering': key in covered,
            'selected_exact_protein': any(_compatible(target, q) and q['protein'] == target['protein'] for q in selected),
            'correct': {method: hits[key] for method, hits in statuses.items()},
        }
    return details


def _summarize(cases, methods):
    groups = defaultdict(list)
    for case in cases:
        groups[case['ancestry']].append(case)
    return {'targets': len(cases), 'ancestries': len(groups),
            'pool_covering': sum(case['pool_covering'] for case in cases),
            'selected_covering': sum(case['selected_covering'] for case in cases),
            'selected_exact_protein': sum(case['selected_exact_protein'] for case in cases),
            'correct': {method: sum(case['correct'][method] for case in cases) for method in methods},
            'equal_ancestry': {method: mean(mean(case['correct'][method] for case in group) for group in groups.values())
                              if groups else None for method in methods}}


def score(records, refs, panel_hash):
    import pyrodigal

    if version('pyrodigal') != PYRODIGAL_VERSION:
        raise ValueError(f'Expected Pyrodigal {PYRODIGAL_VERSION}')
    finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    pools = {'minimal': {key: _candidates(row['dna'], start_aware=False) for key, row in records.items()},
             'gene_caller': {key: _gene_caller_candidates(row['dna'], finder) for key, row in records.items()}}
    targets = [target for entry in records.values() for target in entry['targets']]
    distances = _identities(targets, refs)
    scored = [q['protein'] for arm in pools.values() for pool in arm.values() for q in pool[:max(CAPS)]]
    scored.extend(target['protein'] for target in targets)
    complete = {target['id']: target['complete_product'] for target in targets if 'complete_product' in target}
    scored.extend(complete.values())
    ranked, models = _predictions({_signature(sequence): sequence for sequence in scored}, refs)
    metadata_fields = ('source_product', 'family', 'ancestry', 'layout', 'submission_group')
    metadata = {t['id']: {**{key: t[key] for key in metadata_fields if key in t},
                           **distances[t['id']], 'joined_cds': len(t['spans']) > 1} for t in targets}
    diagnostics = {t['id']: {method: hits[_signature(t['protein'])] == t['family']
                              for method, hits in ranked.items()} for t in targets}
    results = {}
    for arm, by_record in pools.items():
        results[arm] = {}
        for cap in CAPS:
            cases = {key: {**metadata[key], **value} for record_id, pool in by_record.items()
                     for key, value in _score_targets(records[record_id]['targets'], pool, pool[:cap], ranked).items()}
            results[arm][cap] = {
                'all': _summarize(list(cases.values()), ranked),
                'low': _summarize([row for row in cases.values() if row['low_identity']], ranked),
                'under_25': _summarize([row for row in cases.values() if row['max_same_family_global_identity'] < 0.25],
                                       ranked),
                'selected_inputs_sha256': _signature(json.dumps(
                    {key: [(_signature(q['protein']), q['strand'], q['spans']) for q in pool[:cap]]
                     for key, pool in by_record.items()}, sort_keys=True)),
                'actual_queries': sum(min(cap, len(pool)) for pool in by_record.values()),
                'query_residues': sum(len(q['protein']) for pool in by_record.values() for q in pool[:cap]),
                'cases': cases,
            }
    report = {'schema': 'distant-public-reporter-results-1', 'panel_sha256': panel_hash, 'primary_cap': PRIMARY_CAP,
            'low_identity_below': REMOTE_IDENTITY, 'locus_overlap_at_least': MIN_OVERLAP,
            'unique_candidate_and_diagnostic_proteins': len({_signature(sequence) for sequence in scored}),
            'reference_proteins': len(refs), 'models': models, 'targets': metadata,
            'environment': {name: version(name) for name in ('biopython', 'numpy', 'pyfamsa', 'pyhmmer', 'pyrodigal',
                                                               'torch', 'transformers')},
            'known_boundaries_correct': diagnostics, 'construct_results': results}
    if complete:
        report['known_full_products_correct'] = {
            key: {method: hits[_signature(protein)] == metadata[key]['family'] for method, hits in ranked.items()}
            for key, protein in complete.items()
        }
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel', choices=('initial', 'alternate', 'renilla', 'gus-context'), default='initial')
    stages = parser.add_subparsers(dest='stage', required=True)
    stages.add_parser('acquire').add_argument('folder', type=Path)
    stages.add_parser('audit').add_argument('folder', type=Path)
    runner = stages.add_parser('score')
    runner.add_argument('folder', type=Path)
    runner.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    panel_path = {'renilla': RENILLA, 'gus-context': GUS_CONTEXT}.get(args.panel, PANEL)
    panel = read(panel_path)
    panel_hash = digest(panel_path)
    if args.panel == 'alternate':
        panel.update(read(ALTERNATE))
        panel_hash = _signature(panel_hash + digest(ALTERNATE))
    if args.stage == 'acquire':
        acquire(panel, args.folder, canonical_record_end=True)
        return
    records, refs = load(args.folder, panel)
    if args.stage == 'audit':
        targets = [target for row in records.values() for target in row['targets']]
        distances = _identities(targets, refs)
        for target in targets:
            print(target['id'], target['ancestry'], distances[target['id']])
        return
    write(args.out, score(records, refs, panel_hash))


if __name__ == '__main__':
    main()
