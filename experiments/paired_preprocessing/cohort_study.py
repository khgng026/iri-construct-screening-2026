"""Freeze, prepare and run the fresh grouped construct cohort (cohort v2).

Stages, each reading local bytes except where noted:

    freeze    evaluator only: eligible components and their evidence, lineage and
              shared-protein groups, the development/holdout partition and the
              low-similarity tier, from the PR #19 audit and its source snapshots
    prepare   sequence-only method inputs with opaque IDs, evaluator answer keys,
              annotation-derived units, group manifests, the fusion subregion key
              and the reference panel, all under the ignored run directory
    run       one preprocessing arm with one scorer, writing a verifier trace

Methods under test receive only the opaque-ID FASTA. Annotations, accessions,
groups and tiers stay in the evaluator files.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
import time
from collections import Counter, defaultdict
from pathlib import Path

from . import corpus_batch as cb
from .reporter_pilot import (AMINO_ACIDS, MAX_RESIDUES, MIN_RESIDUES, PYRODIGAL_VERSION, _candidates,
                             _characterizer, _gene_caller_candidates, _length_order, _predictions, _rank_by_start)
from .verify import SCHEMA

ROOT = Path(__file__).resolve().parents[2]
FREEZE_SCHEMA = 'paired-cohort-freeze-1'
GROUP_SCHEMA = 'paired-recall-groups-1'

# Fixed before any method output on this cohort exists (see COHORT_STUDY.md).
RULES = {
    'version': 'cohort-v2',
    'seed': '20261123',
    'evidence_tiers': {
        'corroborated': 'reviewed UniProt anchor (identity >= 0.90, anchor coverage >= 0.80) and the deposited '
                        'CDS name or note names the same family',
        'anchor_only': 'reviewed UniProt anchor only; deposited text is silent or names only a partner',
    },
    'lineage_links': [
        'same direct submission (authors and submission text, including date and address)',
        'same published reference title',
        'same unpublished or catalog reference title',
        'near-copy constructs: at least 90% of the sampled canonical 25-mers of each record are in the other',
    ],
    'conservative_links': 'lineage links plus any identical eligible protein shared by two records',
    'partition': 'connected components of the conservative graph; the largest goes to holdout; of the rest, '
                 'those ranking in the lowest 20% by sha256(seed:first accession) go to development',
    'development_fraction': 0.2,
    'kmer': {'k': 25, 'sample_modulus': 8, 'near_copy_fraction': 0.90},
    'low_similarity': {'low_identity_below': 0.40, 'high_identity_at_least': 0.60, 'min_coverage_each': 0.70,
                       'query': 'whole protein; anchor region residues for fusion-class components'},
    'segmentation': {'window_residues': 300, 'stride_residues': 150,
                     'rule': 'Pyrodigal candidates in length order; a candidate longer than the window is replaced '
                             'by windows at the stride, the last ending at the candidate end'},
}


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def read(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path: Path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def file_digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------ grouping

class UnionFind:
    def __init__(self, keys):
        self.parent = {k: k for k in keys}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            # Deterministic root: the lexicographically smaller accession.
            if rb < ra:
                ra, rb = rb, ra
            self.parent[rb] = ra

    def groups(self):
        out = defaultdict(list)
        for key in sorted(self.parent):
            out[self.find(key)].append(key)
        return dict(out)


def _norm(text: str) -> str:
    return re.sub(r'\W+', ' ', text.lower()).strip()


def lineage_keys(references):
    """Evaluator-side source lineage keys from GenBank REFERENCE blocks."""
    keys = set()
    for ref in references:
        journal, title = ref.journal or '', ref.title or ''
        if title == 'Direct Submission' or journal.startswith('Submitted'):
            keys.add(('submission', _norm(ref.authors or ''), _norm(journal)))
        elif title and 'unpublished' in journal.lower():
            keys.add(('unpublished_title', _norm(title)))
        elif title:
            keys.add(('publication_title', _norm(title)))
    return keys


_COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


def sampled_kmers(sequence: str, k: int, modulus: int):
    """Deterministic sample of canonical k-mers, treating the record as circular."""
    wrapped = sequence + sequence[:k - 1]
    reverse = wrapped.translate(_COMPLEMENT)[::-1]
    total = len(wrapped)
    out = set()
    for i in range(len(sequence)):
        forward = wrapped[i:i + k]
        backward = reverse[total - i - k:total - i]
        value = int.from_bytes(hashlib.blake2b(min(forward, backward).encode(), digest_size=8).digest(), 'little')
        if value % modulus == 0:
            out.add(value)
    return out


def near_copy_pairs(kmer_sets, fraction):
    index = defaultdict(list)
    for accession in sorted(kmer_sets):
        for kmer in kmer_sets[accession]:
            index[kmer].append(accession)
    shared = defaultdict(Counter)
    for accessions in index.values():
        for i, a in enumerate(accessions):
            for b in accessions[i + 1:]:
                shared[a][b] += 1
    pairs = []
    for a in sorted(shared):
        for b, count in sorted(shared[a].items()):
            if count >= fraction * max(len(kmer_sets[a]), len(kmer_sets[b]), 1):
                pairs.append((a, b))
    return pairs


def group_stats(groups, weights):
    sizes = sorted(((len(members), sum(weights[m] for m in members)) for members in groups.values()), reverse=True)
    total = sum(w for _, w in sizes)
    return {'groups': len(groups), 'groups_with_weight': sum(1 for _, w in sizes if w),
            'largest_records': sizes[0][0] if sizes else 0,
            'largest_share': sizes[0][1] / total if total else None,
            'effective_groups': (1 / sum((w / total) ** 2 for _, w in sizes)) if total else None}


# ------------------------------------------------------------------- freeze

def _toolkit_units(paths):
    """Validated protein sequences by (accession, feature index) from the annotation-aware toolkit."""
    from iri_constructs import process_inputs

    out = {}
    for path in paths:
        result = process_inputs([path])
        for record in result.records:
            if record.status != 'ok':
                continue
            units = {u.feature_id: u for u in record.units}
            for feature in record.features:
                if feature.source_feature_type == 'CDS':
                    unit = units[feature.feature_id]
                    if unit.extraction_state == 'exported':
                        out[(record.source_record_id, feature.source_feature_index)] = unit.sequence
    return out


def _reference_panel(corpus_run: Path, audit):
    rows = {}
    for path in sorted((corpus_run / 'sources' / 'uniprot').glob('chunk_*.tsv')):
        for row in csv.DictReader(io.StringIO(path.read_text()), delimiter='\t'):
            rows[row['Entry']] = row
    freeze = read(cb.STUDY / 'panel.freeze.json')
    panel = []
    for family, ids in sorted(audit['scorable_families'].items()):
        frozen = {r['reference_id']: r for r in freeze['families'][family]['references']}
        for ref_id in ids:
            sequence = rows[ref_id.removeprefix('uniprot-')]['Sequence']
            if sha(sequence) != frozen[ref_id]['sequence_sha256']:
                raise ValueError(f'{ref_id}: reference differs from the frozen component-study panel')
            panel.append({'accession': ref_id.removeprefix('uniprot-'), 'family': family, 'sequence': sequence})
    return panel


def _aligner():
    from Bio.Align import PairwiseAligner, substitution_matrices

    return PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                           open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)


def identity_coverage(aligner, query, reference):
    hit = next(iter(aligner.align(query, reference)), None)
    pairs = identical = 0
    if hit is not None:
        for (a, b), (c, d) in zip(*hit.aligned, strict=True):
            pairs += int(b - a)
            identical += sum(x == y for x, y in zip(query[a:b], reference[c:d], strict=True))
    return {'identity': identical / pairs if pairs else 0.0,
            'coverage': min(pairs / len(query), pairs / len(reference))}


def similarity_tier(aligner, protein, family, references):
    rule = RULES['low_similarity']
    rows = [{'accession': r['accession'], 'family': r['family'], **identity_coverage(aligner, protein, r['sequence'])}
            for r in references]
    same = [r for r in rows if r['family'] == family and r['coverage'] >= rule['min_coverage_each']]
    best = max((r['identity'] for r in same), default=None)
    if best is None:
        tier = 'insufficient_coverage'
    elif best < rule['low_identity_below']:
        tier = 'low'
    elif best < rule['high_identity_at_least']:
        tier = 'intermediate'
    else:
        tier = 'high'
    competitor = best is not None and any(r['family'] != family and r['coverage'] >= rule['min_coverage_each']
                                          and r['identity'] >= best for r in rows)
    best_any = max(rows, key=lambda r: (r['coverage'] >= rule['min_coverage_each'], r['identity']))
    return {'tier': tier, 'best_same_family_identity': best, 'strong_competitor': competitor,
            'max_identity_any_reference': max(r['identity'] for r in rows),
            'closest_adequate_family': best_any['family'] if best_any['coverage'] >= rule['min_coverage_each'] else None}


def freeze(corpus_run: Path, out: Path):
    from Bio import SeqIO

    audit = cb.read(corpus_run / 'audit.private.json')
    acquisition = cb.read(corpus_run / 'acquisition.private.json')
    if not all(f['byte_identical'] for s in acquisition['samples'].values() for f in s['files']):
        raise ValueError('Corpus snapshot differs from the frozen component-study batches')
    scorable = set(audit['scorable_families'])
    reporter = cb.read(Path(__file__).with_name('reporter_panel.json'))
    pilot = {row['accession'] for row in reporter['constructs']}

    records = {}
    for construct in audit['constructs']:
        if construct['status'] != 'ok' or construct['accession'] in pilot:
            continue
        components = cb.eligible_components(construct, scorable)
        if components:
            records[construct['accession']] = (construct, components)

    batch_paths, references, sequences = {}, {}, {}
    for path in sorted(corpus_run.glob('sources/*/batch_*.gb')):
        for record in SeqIO.parse(path, 'genbank'):
            if record.id in records:
                batch_paths[record.id] = path
                references[record.id] = record.annotations.get('references', [])
                sequences[record.id] = str(record.seq).upper()
    if batch_paths.keys() != records.keys():
        raise ValueError('Eligible records are missing from the source snapshots')
    proteins = _toolkit_units(sorted(set(batch_paths.values())))
    panel = _reference_panel(corpus_run, audit)

    # Lineage, conservative and partition graphs.
    kmer_rule = RULES['kmer']
    kmers = {a: sampled_kmers(s, kmer_rule['k'], kmer_rule['sample_modulus']) for a, s in sequences.items()}
    copies = near_copy_pairs(kmers, kmer_rule['near_copy_fraction'])
    lineage = UnionFind(records)
    owner = {}
    for accession in sorted(records):
        for key in sorted(lineage_keys(references[accession])):
            lineage.union(accession, owner.setdefault(key, accession))
    for a, b in copies:
        lineage.union(a, b)
    conservative = UnionFind(records)
    for a, b in ((a, lineage.find(a)) for a in records):
        conservative.union(a, b)
    owner = {}
    for accession in sorted(records):
        for component in records[accession][1]:
            conservative.union(accession, owner.setdefault(component['protein_sha256'], accession))
    lineage_groups, conservative_groups = lineage.groups(), conservative.groups()
    lineage_of = {a: f'L{root}' for root, members in lineage_groups.items() for a in members}
    conservative_of = {a: f'C{root}' for root, members in conservative_groups.items() for a in members}

    largest = max(conservative_groups, key=lambda r: (len(conservative_groups[r]), r))
    rest = sorted((r for r in conservative_groups if r != largest), key=lambda r: sha(f'{RULES["seed"]}:{r}'))
    n_dev = round(RULES['development_fraction'] * len(rest))
    partition_of_group = {r: ('development' if i < n_dev else 'holdout') for i, r in enumerate(rest)}
    partition_of_group[largest] = 'holdout'
    partition = {a: partition_of_group[conservative.find(a)] for a in records}

    aligner = _aligner()
    tier_cache = {}
    reporter_proteins = set()
    rows = []
    for accession in sorted(records):
        construct, components = records[accession]
        out_components = []
        for component in components:
            protein = proteins[(accession, component['feature_index'])]
            if sha(protein) != component['protein_sha256']:
                raise ValueError(f'{accession}: toolkit translation differs from the audit')
            query = protein
            if component['fusion_class']:
                start, end = component['anchor_region']['residues']
                query = protein[start:end]
            key = (query, component['family'])
            if key not in tier_cache:
                tier_cache[key] = similarity_tier(aligner, query, component['family'], panel)
            evidence_tier = ('corroborated' if component['name_evidence'] in ('deposited_name_agrees',
                                                                                'deposited_note_agrees')
                             else 'anchor_only')
            row = {'id': f'{component["feature_index"]}-{component["family"]}', 'family': component['family'],
                   'feature_index': component['feature_index'], 'strand': component['strand'],
                   'spans': component['spans'], 'length_aa': component['length_aa'],
                   'protein_sha256': component['protein_sha256'],
                   'translation_evidence': 'toolkit validated the deposited translation against the DNA',
                   'family_evidence': {'anchor': component['anchor'], 'label_class': component['label_class'],
                                       'deposited_name': component['name_evidence']},
                   'evidence_tier': evidence_tier, 'fusion_class': component['fusion_class'],
                   'over_model_limit': component['over_model_limit'],
                   'similarity': tier_cache[key]}
            if component['fusion_class']:
                row['anchor_region'] = component['anchor_region']
                boundary = component['deposited_subregion']
                row['documented_subregion'] = (boundary if boundary and boundary['usable']
                                               and boundary['overlaps_anchor_region'] else None)
            out_components.append(row)
        rows.append({'accession': accession, 'sample': construct['sample'],
                     'batch_file': f'{construct["sample"]}/{construct["batch_file"]}',
                     'batch_sha256': file_digest(batch_paths[accession]),
                     'sequence_sha256': construct['sequence_sha256'], 'length': construct['length'],
                     'topology': construct['topology'], 'archetype': cb.archetype_of(construct, components),
                     'lineage_group': lineage_of[accession], 'conservative_group': conservative_of[accession],
                     'partition': partition[accession], 'components': out_components})

    # Repeated proteins across partitions and against the reporter development panel.
    by_protein = defaultdict(set)
    for row in rows:
        for component in row['components']:
            by_protein[component['protein_sha256']].add(row['partition'])
    crossing = sorted(p for p, parts in by_protein.items() if len(parts) > 1)
    if crossing:
        raise ValueError('An eligible protein crosses the development/holdout partition')
    reporter_run = ROOT / 'results' / 'paired-reporter-pilot-v0' / 'prepared' / 'evaluator' / 'annotation-units.json'
    if reporter_run.exists():
        reporter_proteins = {sha(u['protein']) for units in read(reporter_run).values() for u in units}
    holdout_reporter_overlap = sorted(p for p in reporter_proteins if 'holdout' in by_protein.get(p, ()))

    weights = {row['accession']: sum(c['evidence_tier'] == 'corroborated' for c in row['components']) for row in rows}
    summary = {}
    for part in ('development', 'holdout'):
        members = {r['accession'] for r in rows if r['partition'] == part}
        sub = lambda groups: {k: [m for m in v if m in members] for k, v in groups.items()  # noqa: E731
                              if any(m in members for m in v)}
        summary[part] = {'records': len(members),
                         'lineage': group_stats(sub(lineage_groups), weights),
                         'conservative': group_stats(sub(conservative_groups), weights)}
    document = {
        'schema': FREEZE_SCHEMA, 'rules': RULES,
        'corpus_audit_sha256': file_digest(corpus_run / 'audit.private.json'),
        'corpus_acquisition_sha256': file_digest(corpus_run / 'acquisition.private.json'),
        'reference_panel': [{k: r[k] for k in ('accession', 'family')} | {'sequence_sha256': sha(r['sequence'])}
                            for r in panel],
        'near_copy_construct_pairs': len(copies),
        'reporter_development_proteins_in_holdout': holdout_reporter_overlap,
        'summary': summary, 'records': rows,
    }
    write(out / 'freeze.private.json', document)
    print(json.dumps(summary, indent=1))


# ------------------------------------------------------------------ prepare

def _opaque_ids(accessions, prefix):
    ordered = sorted(accessions, key=lambda a: sha(f'{RULES["seed"]}:opaque:{a}'))
    return {a: f'{prefix}{i + 1:04d}' for i, a in enumerate(ordered)}


def _subregion_component(component):
    sub = component['documented_subregion']
    return {'id': f'{component["id"]}-region{sub["feature_index"]}', 'family': component['family'],
            'strand': component['strand'], 'spans': sub['spans'], 'residues': sub['residues'],
            'parent_component': component['id']}


def prepare(corpus_run: Path, run: Path):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    document = read(run / 'freeze.private.json')
    records = document['records']
    sequences, paths = {}, sorted({corpus_run / 'sources' / r['batch_file'] for r in records})
    wanted = {r['accession'] for r in records}
    for path in paths:
        for record in SeqIO.parse(path, 'genbank'):
            if record.id in wanted:
                sequences[record.id] = str(record.seq).upper()
    proteins = _toolkit_units(paths)
    panel = _reference_panel(corpus_run, cb.read(corpus_run / 'audit.private.json'))
    if [{k: r[k] for k in ('accession', 'family')} | {'sequence_sha256': sha(r['sequence'])} for r in panel] \
            != document['reference_panel']:
        raise ValueError('Reference panel differs from the freeze')
    prepared = run / 'prepared'
    reference_bytes = (json.dumps(panel, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
    prepared.mkdir(parents=True)
    (prepared / 'references.private.json').write_bytes(reference_bytes)

    hashes = {'references': hashlib.sha256(reference_bytes).hexdigest()}
    for part, prefix in (('development', 'd'), ('holdout', 'h')):
        rows = [r for r in records if r['partition'] == part]
        ids = _opaque_ids([r['accession'] for r in rows], prefix)
        folder = prepared / part
        (folder / 'method_inputs').mkdir(parents=True)
        (folder / 'evaluator').mkdir()
        fasta, keys, units, groups, map_rows = [], defaultdict(list), {}, defaultdict(list), []
        fusion_fasta, fusion_key, region_units = [], [], {}
        for row in sorted(rows, key=lambda r: ids[r['accession']]):
            key, sequence = ids[row['accession']], sequences[row['accession']]
            if sha(sequence) != row['sequence_sha256'] or len(sequence) != row['length']:
                raise ValueError(f'{row["accession"]}: construct sequence differs from the freeze')
            fasta.append(SeqRecord(Seq(sequence), id=key, name=key, description=''))
            for tier in ('corroborated', 'anchor_only', 'all'):
                chosen = [c for c in row['components'] if tier == 'all' or c['evidence_tier'] == tier]
                keys[tier].append({'id': key, 'archetype': row['archetype'], 'length': row['length'],
                                   'components': [{k: c[k] for k in ('id', 'family', 'strand', 'spans')}
                                                  for c in chosen]})
            seen = {}
            for component in row['components']:
                seen.setdefault(component['feature_index'], component)
            units[key] = [{'id': f'{key}-{c["feature_index"]}', 'strand': c['strand'], 'spans': c['spans'],
                           'protein': proteins[(row['accession'], c['feature_index'])]} for c in seen.values()]
            for kind, field in (('lineage', 'lineage_group'), ('conservative', 'conservative_group')):
                groups[kind].append({'id': key, 'group': row[field], 'partition': part})
            map_rows.append({'id': key, 'accession': row['accession'], 'lineage_group': row['lineage_group'],
                             'conservative_group': row['conservative_group'], 'archetype': row['archetype'],
                             'components': [{k: c[k] for k in ('id', 'family', 'evidence_tier', 'fusion_class',
                                                               'protein_sha256', 'over_model_limit', 'similarity')}
                                            | {'documented_subregion': bool(c.get('documented_subregion'))}
                                            for c in row['components']]})
            regions = [c for c in row['components'] if c.get('documented_subregion')]
            if regions:
                fusion_fasta.append(SeqRecord(Seq(sequence), id=key, name=key, description=''))
                fusion_key.append({'id': key, 'archetype': row['archetype'], 'length': row['length'],
                                   'components': [{k: v for k, v in _subregion_component(c).items()
                                                   if k in ('id', 'family', 'strand', 'spans')} for c in regions]})
                region_units[key] = []
                for c in regions:
                    sub = _subregion_component(c)
                    protein = proteins[(row['accession'], c['feature_index'])]
                    start, end = sub['residues']
                    region_units[key].append({'id': sub['id'], 'strand': c['strand'], 'spans': sub['spans'],
                                              'protein': protein[start:end]})
        SeqIO.write(fasta, folder / 'method_inputs' / 'constructs.fasta', 'fasta')
        hashes[f'{part}_inputs'] = file_digest(folder / 'method_inputs' / 'constructs.fasta')
        for tier, constructs in keys.items():
            write(folder / 'evaluator' / f'answer-key.{tier}.json', {'schema': SCHEMA, 'constructs': constructs})
            hashes[f'{part}_answer_key_{tier}'] = file_digest(folder / 'evaluator' / f'answer-key.{tier}.json')
        write(folder / 'evaluator' / 'annotation-units.json', units)
        for kind, rows_ in groups.items():
            write(folder / 'evaluator' / f'groups.{kind}.json', {'schema': GROUP_SCHEMA, 'constructs': rows_})
        write(folder / 'evaluator' / 'construct-map.private.json', map_rows)
        hashes[f'{part}_annotation_units'] = file_digest(folder / 'evaluator' / 'annotation-units.json')
        if fusion_fasta:
            (folder / 'fusion').mkdir()
            SeqIO.write(fusion_fasta, folder / 'fusion' / 'constructs.fasta', 'fasta')
            write(folder / 'fusion' / 'subregion-key.json', {'schema': SCHEMA, 'constructs': fusion_key})
            write(folder / 'fusion' / 'documented-regions.json', region_units)
            hashes[f'{part}_fusion_inputs'] = file_digest(folder / 'fusion' / 'constructs.fasta')
            hashes[f'{part}_fusion_key'] = file_digest(folder / 'fusion' / 'subregion-key.json')
    hashes['freeze'] = file_digest(run / 'freeze.private.json')
    write(prepared / 'hashes.json', hashes)
    print(json.dumps(hashes, indent=1))


def fusion_units(partition_dir: Path):
    """Whole documented units for the constructs in the fusion subset (annotation-derived diagnostic)."""
    from Bio import SeqIO

    keys = {r.id for r in SeqIO.parse(partition_dir / 'fusion' / 'constructs.fasta', 'fasta')}
    units = read(partition_dir / 'evaluator' / 'annotation-units.json')
    write(partition_dir / 'fusion' / 'annotation-units.json', {k: v for k, v in units.items() if k in keys})


PROTOCOL = Path(__file__).with_name('cohort_protocol.json')
PROTOCOL_FIELDS = ('schema', 'cohort', 'query_caps', 'primary_query_cap', 'prediction_limit',
                   'minimum_component_overlap', 'minimum_interesting_gain', 'reference_panel_sha256')


def scorer_protocol(inputs_key: str, scorer: str, template=None):
    """The verifier protocol for one partition and scorer, from the committed template."""
    template = template or read(PROTOCOL)
    protocol = {k: template[k] for k in PROTOCOL_FIELDS}
    protocol.update(template['scorers'][scorer])
    protocol['sequence_input_sha256'] = template['inputs'][inputs_key]
    if 'fusion' in inputs_key:
        protocol.update({k: template['fusion'][k] for k in ('query_caps', 'primary_query_cap')})
    return protocol


def public_freeze(run: Path, out: Path):
    """Source identities, partitions and groups, bound to the private evaluator files by hash.

    Family labels, loci and tiers stay in the private files; the committed freeze holds no
    sequence, qualifier text or label.
    """
    document = read(run / 'freeze.private.json')
    hashes = read(run / 'prepared' / 'hashes.json')
    records = [{k: r[k] for k in ('accession', 'sample', 'batch_file', 'batch_sha256', 'sequence_sha256',
                                  'partition', 'lineage_group', 'conservative_group')}
               | {'eligible_components': len(r['components'])} for r in document['records']]
    write(out, {'schema': 'paired-cohort-public-freeze-1', 'rules': document['rules'],
                'corpus_audit_sha256': document['corpus_audit_sha256'],
                'corpus_acquisition_sha256': document['corpus_acquisition_sha256'],
                'reference_panel': document['reference_panel'], 'private_file_sha256': hashes,
                'summary': document['summary'], 'near_copy_construct_pairs': document['near_copy_construct_pairs'],
                'reporter_development_proteins_in_holdout': len(document['reporter_development_proteins_in_holdout']),
                'records': records})


# ---------------------------------------------------------------------- run

def _windows(units):
    """Fixed windows over each candidate, keeping the parent length order (segmentation arm)."""
    size, stride = RULES['segmentation']['window_residues'], RULES['segmentation']['stride_residues']
    out = []
    for unit in units:
        protein, strand = unit['protein'], unit['strand']
        if len(protein) <= size:
            out.append(unit)
            continue
        starts = list(range(0, len(protein) - size, stride)) + [len(protein) - size]
        (lo, hi), = unit['spans']
        anchor = lo if strand == 1 else hi
        for offset in starts:
            start = anchor + strand * offset * 3
            end = start + strand * size * 3
            out.append({'id': f'{unit["id"]}w{offset}', 'strand': strand, 'spans': [sorted((start, end))],
                        'protein': protein[offset:offset + size]})
    return out


ARMS = ('minimal', 'start_ranked', 'gene_caller', 'gene_caller_windows', 'annotation_reference',
        'documented_regions')
UNIT_ARMS = ('annotation_reference', 'documented_regions')


def chunk_unit(unit):
    """Apply the sequence arms' model length rule to a documented unit, including joined spans.

    Fixed after the development check (COHORT_STUDY.md, Checkpoint C): the earlier diagnostic
    dropped documented units longer than the model limit while the sequence arms split them.
    Coordinates assume translation from the first base of the first segment.
    """
    protein = unit['protein']
    if len(protein) <= MAX_RESIDUES:
        return [unit]
    segments = sorted(map(tuple, unit['spans']), reverse=unit['strand'] == -1)
    positions = [p for a, b in segments for p in (range(a, b) if unit['strand'] == 1 else reversed(range(a, b)))]
    chunks = []
    for index, offset in enumerate(range(0, len(protein), MAX_RESIDUES)):
        piece = protein[offset:offset + MAX_RESIDUES]
        if len(piece) < MIN_RESIDUES:
            continue
        covered = sorted(positions[3 * offset:3 * (offset + len(piece))])
        spans = []
        for p in covered:
            if spans and spans[-1][1] == p:
                spans[-1][1] = p + 1
            else:
                spans.append([p, p + 1])
        chunks.append({'id': f'{unit["id"]}c{index + 1}', 'strand': unit['strand'], 'spans': spans, 'protein': piece})
    return chunks


def candidate_units(method, original, finder):
    if method in UNIT_ARMS:
        units = [chunk for unit in original for chunk in chunk_unit(unit)]
        return sorted(units, key=lambda u: (-len(u['protein']), u['id']))
    if method in ('gene_caller', 'gene_caller_windows'):
        units = _gene_caller_candidates(original, finder)
        return _windows(units) if method == 'gene_caller_windows' else units
    units = _candidates(original, False)
    return _rank_by_start(units) if method == 'start_ranked' else units


def run(method, inputs, references, protocol, output, characterizer_factory=None):
    """Run one arm. ``characterizer_factory`` substitutes another scorer with the reporter interface."""
    input_sha256 = file_digest(inputs)
    if method not in UNIT_ARMS and input_sha256 != protocol['sequence_input_sha256']:
        raise ValueError('Sequence input differs from the fixed protocol')
    reference_sha256 = file_digest(references)
    if reference_sha256 != protocol['reference_panel_sha256']:
        raise ValueError('Reference panel differs from the fixed protocol')
    finder = None
    if method.startswith('gene_caller'):
        import pyrodigal
        if pyrodigal.__version__ != PYRODIGAL_VERSION:
            raise ValueError(f'The gene caller arms require pyrodigal=={PYRODIGAL_VERSION}')
        finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    from Bio import SeqIO

    records = (read(inputs) if method in UNIT_ARMS else {r.id: str(r.seq) for r in SeqIO.parse(inputs, 'fasta')})
    reference_rows = read(references)
    reference_residues = sum(len(r['sequence']) for r in reference_rows)
    characterizer = protocol.get('characterizer', 'esm2')
    selected_inputs = hashlib.sha256()
    cache, unique_seconds = {}, 0.0
    started_run = time.perf_counter()
    factory = characterizer_factory or _characterizer
    with factory(characterizer, reference_rows) as (rank_families, metadata):
        def memo(protein):
            # Identical inputs give identical rankings; caching changes cost, not results.
            nonlocal unique_seconds
            if protein not in cache:
                t = time.perf_counter()
                cache[protein] = rank_families(protein)
                unique_seconds += time.perf_counter() - t
            return cache[protein]

        trace = {'schema': SCHEMA, 'method': method, 'characterization_id': protocol['characterization_id'],
                 'characterizer': characterizer, 'reference_panel_sha256': reference_sha256,
                 'reference_count': len(reference_rows), 'tiers': [], 'input_sha256': input_sha256, **metadata}
        if finder is not None:
            trace['gene_caller'] = {'engine': 'pyrodigal', 'version': PYRODIGAL_VERSION,
                                    'meta': True, 'closed': False, 'mask': False}
        if method == 'gene_caller_windows':
            trace['segmentation'] = RULES['segmentation']
        pools = {}
        for cap in protocol['query_caps']:
            rows = []
            for key, original in sorted(records.items()):
                started = time.perf_counter()
                if key not in pools:
                    units = candidate_units(method, original, finder)
                    pools[key] = [u for u in units if MIN_RESIDUES <= len(u['protein']) <= MAX_RESIDUES
                                  and not set(u['protein']) - AMINO_ACIDS]
                usable = pools[key]
                selected = usable[:cap]
                prepared = time.perf_counter()
                predictions = _predictions(selected, memo, protocol['prediction_limit'])
                scored = time.perf_counter()
                encoded = json.dumps((cap, key, selected), sort_keys=True, separators=(',', ':')).encode('utf-8')
                selected_inputs.update(encoded + b'\n')
                residues = sum(len(q['protein']) for q in selected)
                rows.append({'id': key, 'status': 'ok', 'candidate_count': len(usable),
                             'candidate_residues': sum(len(u['protein']) for u in usable),
                             'queries': [{**{k: q[k] for k in ('id', 'strand', 'spans')},
                                          'input_residues': len(q['protein'])} for q in selected],
                             'predictions': predictions,
                             # Full family order per query, for diagnosing covered loci with no correct call.
                             'query_family_rankings': {q['id']: memo(q['protein']) for q in selected},
                             'preprocess_seconds': prepared - started,
                             'characterize_seconds': scored - prepared,
                             'work': {'scorer_inputs': len(selected),
                                      'reference_comparisons': len(selected) * len(reference_rows),
                                      'alignment_cells': residues * reference_residues}})
            trace['tiers'].append({'query_cap': cap, 'constructs': rows})
        trace['selected_inputs_sha256'] = selected_inputs.hexdigest()
        # Loci of every usable candidate before the cap, for separating pool and selection misses.
        trace['candidate_pools'] = {key: [{'strand': u['strand'], 'spans': u['spans'], 'residues': len(u['protein'])}
                                          for u in pool] for key, pool in sorted(pools.items())}
        trace['unique_scorer_inputs'] = len(cache)
        trace['unique_scoring_seconds'] = unique_seconds
        trace['wall_seconds'] = time.perf_counter() - started_run
        write(output, trace)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    proto = sub.add_parser('protocol')
    proto.add_argument('--inputs-key', required=True)
    proto.add_argument('--scorer', choices=('esm2', 'blosum62'), required=True)
    proto.add_argument('--out', type=Path, required=True)
    public = sub.add_parser('public-freeze')
    public.add_argument('--run', type=Path, default=ROOT / 'results' / 'paired-cohort-v2')
    public.add_argument('--out', type=Path, required=True)
    units = sub.add_parser('fusion-units')
    units.add_argument('--partition-dir', type=Path, required=True)
    for stage in ('freeze', 'prepare'):
        command = sub.add_parser(stage)
        command.add_argument('--corpus-run', type=Path, default=ROOT / 'results' / 'paired-corpus-batch-v1')
        command.add_argument('--run', type=Path, default=ROOT / 'results' / 'paired-cohort-v2')
    execute = sub.add_parser('run')
    execute.add_argument('--method', choices=ARMS, required=True)
    execute.add_argument('--inputs', type=Path, required=True)
    execute.add_argument('--references', type=Path, required=True)
    execute.add_argument('--protocol', type=Path, required=True)
    execute.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.stage == 'freeze':
        freeze(args.corpus_run, args.run)
    elif args.stage == 'prepare':
        prepare(args.corpus_run, args.run)
    elif args.stage == 'protocol':
        write(args.out, scorer_protocol(args.inputs_key, args.scorer))
    elif args.stage == 'public-freeze':
        public_freeze(args.run, args.out)
    elif args.stage == 'fusion-units':
        fusion_units(args.partition_dir)
    else:
        run(args.method, args.inputs, args.references, read(args.protocol), args.out)


if __name__ == '__main__':
    main()
