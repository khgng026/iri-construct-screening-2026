"""Audit the component-study corpus and prepare a verified paired-preprocessing evaluation batch.

Stream 1 (evaluation and verification) only. Nothing here implements or scores a
preprocessing method. The stages are separate so the network step can be rerun on
its own and every later stage works from local bytes:

    acquire   reacquire the 1,000 component-study records with the exact GI batches
              of the frozen acquisition manifests, compare hashes, and fetch the
              declared UniProt family anchors and panel references
    audit     run the annotation-aware toolkit on every record, compare its
              validated proteins with the frozen unit hashes, recompute the anchor
              alignments with coordinates, read deposited names and subregion
              features, and group related constructs and repeated components
    select    apply the pre-declared selection rule (SELECTION below)
    prepare   write sequence-only method inputs with opaque IDs, the evaluator-only
              answer key in the verifier format, a separate subregion key and a
              source-free summary

Every biological payload, the accession to opaque-ID map and the answer key stay
under the ignored run directory. Only code and the summary are committed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import hmac
import io
import json
import re
import secrets
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .verify import SCHEMA

ROOT = Path(__file__).resolve().parents[2]
STUDY = ROOT / 'docs' / 'results' / 'component-study-v1'
SAMPLES = (('v1', STUDY), ('ext1', STUDY / 'extension-1'))
FAMILY_TABLE = ROOT / 'experiments' / 'component_study_v1' / 'families.v1.json'
EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'
UNIPROT = 'https://rest.uniprot.org/uniprotkb/'
AGENT = {'User-Agent': 'iri-constructs/corpus-batch'}
AUDIT_SCHEMA = 'paired-corpus-audit-1'
SUBREGION_SCHEMA = 'paired-subregion-key-1'

# Same label thresholds as the frozen component study (SPEC.md); not retuned here.
IDENTITY, COVERAGE, MAX_RESIDUES, MIN_RESIDUES = 0.90, 0.80, 1022, 30

# Pre-declared before any method output exists. Archetypes are assigned from the
# evaluator-side evidence in this priority order; quotas bound the first batch.
SELECTION = {
    'version': 'corpus-batch-v1',
    'seed': '20261122',
    'archetype_priority': ['tagged_fusion', 'multi_cds', 'backbone_heavy', 'single_cds'],
    'quotas': {'tagged_fusion': 24, 'multi_cds': 24, 'backbone_heavy': 8, 'single_cds': 8},
    'rules': {
        'tagged_fusion': 'at least one eligible component whose anchor alignment covers under 80% of its CDS '
                         '(label class embedded or multi_family)',
        'multi_cds': 'eligible components on at least two distinct CDSs and no fusion-class component',
        'backbone_heavy': 'one eligible component, eligible component bases under 15% of the construct, and at '
                          'least one further exported or withheld CDS without an eligible label',
        'single_cds': 'exactly one eligible component and no other CDS in the record',
    },
    'diversity': 'within an archetype, take at most one construct per identical eligible-component set and '
                 'round-robin over dependence groups (seeded order) so that one common marker group cannot '
                 'fill the quota; tagged_fusion additionally round-robins over the fusion family set first',
    'backbone_fraction': 0.15,
}


# --------------------------------------------------------------------- utilities

def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def write(path: Path, value, *, exclusive=True):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x' if exclusive else 'w', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def fetch(url: str, limit: int = 20_000_000) -> bytes:
    """Bounded retries, a size cap, and a pause that keeps NCBI/UniProt requests well under their limits."""
    last = None
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=AGENT), timeout=120) as response:
                data = response.read(limit + 1)
            if len(data) > limit:
                raise ValueError('response exceeds size cap')
            time.sleep(0.5)
            return data
        except Exception as error:  # noqa: BLE001 - reported after bounded retries
            last = error
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f'request failed after bounded retries: {last}')


# ----------------------------------------------------------------------- acquire

def acquire(out: Path):
    """Reacquire each frozen batch with its exact GI list and compare it with the manifest."""
    report = {'schema': 'paired-corpus-acquisition-1', 'retrieved_on': datetime.now(timezone.utc).date().isoformat(),
              'samples': {}, 'uniprot': {}}
    for name, folder in SAMPLES:
        manifest = read(folder / 'acquisition.json')
        target = out / 'sources' / name
        target.mkdir(parents=True, exist_ok=True)
        files = []
        for item in manifest['files']:
            path = target / item['file']
            if not path.exists():
                query = urllib.parse.urlencode({'db': 'nuccore', 'id': ','.join(item['requested_gi']),
                                                'rettype': manifest['rettype'], 'retmode': 'text',
                                                'tool': 'iri_constructs'})
                data = fetch(EUTILS + query)
                if not data.lstrip().startswith(b'LOCUS'):
                    raise RuntimeError(f'{name}/{item["file"]}: unexpected EFetch response; nothing cached')
                path.write_bytes(data)
            data = path.read_bytes()
            versions = [line.split()[1] for line in data.decode('utf-8', 'replace').splitlines()
                        if line.startswith('VERSION')]
            files.append({'file': item['file'], 'frozen_sha256': item['sha256'], 'sha256': hashlib.sha256(data).hexdigest(),
                          'byte_identical': hashlib.sha256(data).hexdigest() == item['sha256'],
                          'frozen_versions': item['accession_versions'], 'returned_versions': versions,
                          'missing_versions': sorted(set(item['accession_versions']) - set(versions)),
                          'new_versions': sorted(set(versions) - set(item['accession_versions']))})
        report['samples'][name] = {'manifest_sha256': digest(folder / 'acquisition.json'), 'files': files}

    # Family anchors (independent label source) and the frozen reference panel (label scorability).
    table = read(FAMILY_TABLE)
    anchors = sorted({a for rule in table['families'].values() for a in rule['anchors']})
    freeze = read(STUDY / 'panel.freeze.json')
    references = {r['reference_id'].removeprefix('uniprot-'): {**r, 'family': family}
                  for family, row in freeze['families'].items() for r in row.get('references', [])}
    accessions = sorted(set(anchors) | set(references))
    uniprot_dir = out / 'sources' / 'uniprot'
    uniprot_dir.mkdir(parents=True, exist_ok=True)
    for offset in range(0, len(accessions), 50):
        chunk = accessions[offset:offset + 50]
        path = uniprot_dir / f'chunk_{offset // 50 + 1:02d}.tsv'
        if not path.exists():
            query = ' OR '.join(f'accession:{a}' for a in chunk)
            url = UNIPROT + 'search?' + urllib.parse.urlencode(
                {'query': f'({query})', 'format': 'tsv', 'size': 500,
                 'fields': 'accession,reviewed,protein_name,version,sequence_version,length,sequence'})
            path.write_bytes(fetch(url))
        report['uniprot'][path.name] = {'sha256': digest(path), 'requested': chunk}
    report['uniprot_terms'] = 'https://www.uniprot.org/help/license (CC BY 4.0)'
    report['ncbi_terms'] = read(STUDY / 'acquisition.json')['terms_note']
    write(out / 'acquisition.private.json', report, exclusive=False)
    for name, sample in report['samples'].items():
        same = sum(f['byte_identical'] for f in sample['files'])
        missing = sum(len(f['missing_versions']) for f in sample['files'])
        print(f'{name}: {same}/{len(sample["files"])} batch files byte-identical; {missing} frozen records missing')


# ------------------------------------------------------------------------- audit

# Deposited-name vocabulary. It records whether the submitter's own CDS or
# subregion names agree with the independently supported family; it does not
# assign labels. Case-insensitive; matched against product, gene, label and
# standard_name first, then note.
NAME_PATTERNS = {
    'class_a_beta_lactamase': r'beta[- ]?lactamase|\bbla\b|\bbla[- ]|ampicillin|\bamp[r]?\b|\btem[- ]?1\b',
    'aminoglycoside_phosphotransferase': r'phosphotransferase|\bneo[r]?\b|\bnpt|neomycin|kanamycin|\bkan[r]?\b|\baph|g418|geneticin',
    'chloramphenicol_acetyltransferase': r'chloramphenicol|\bcat\b|\bcm[r]?\b|\bcam[r]?\b',
    'gfp_like': r'fluorescen|\b[a-z]{0,3}[gyrcb]fp|venus|citrine|cerulean|cherry|dsred|tomato|mkate|staygold|azami|kaede|eos\b|gfp',
    'omp_decarboxylase': r'orotidine|\bura3\b|\bura4\b|\bpyr[fg4]\b|omp decarboxylase',
    'laci_galr_repressor': r'\blac ?i|lac repressor|lactose operon repressor|\blac ?iq',
    'tetr_repressor': r'\btet ?r\b|tetr\b|tetracycline repressor|\br?t[ts]a\b|transactivator|transrepressor',
    'tcr_tet_efflux': r'\btet ?[abcl]?\b|tetracycline (resistance|efflux)',
    'coleopteran_luciferase': r'(?<!renilla )(?<!gaussia )luciferase|\bluc\b|\bluc[+2]|\bfluc\b',
    'glycosyl_hydrolase_2': r'galactosidase|\blac ?z|glucuronidase|\bgus[a]?\b|\buid ?a\b',
    'ant3_adenylyltransferase': r'adenylyltransferase|nucleotidyltransferase|\baad|spectinomycin|streptomycin|\bspc|\bsm[r]?\b|\bspec',
    'phage_integrase': r'integrase|\bcre\b|recombinase|\bint\b',
    'bleomycin_resistance': r'bleomycin|zeocin|phleomycin|\bble\b|sh ?ble',
    'pat_bar_acetyltransferase': r'phosphinothricin|\bbar\b|\bpat\b|bialaphos|basta',
    'igp_dehydratase': r'\bhis3\b|imidazoleglycerol',
    'arac_regulator': r'\bara ?c\b|arabinose',
    'isopropylmalate_dehydrogenase': r'\bleu2\b|isopropylmalate',
    'trpf': r'\btrp1\b|phosphoribosylanthranilate|\btrp ?f\b',
    'phage_rna_polymerase': r'rna polymerase|\bt7 ?rnap|\bt[37] ?pol|\bsp6 ?pol',
    'gst_mu': r'glutathione|\bgst\b',
    'solute_binding_protein_1': r'maltose|\bmal ?e\b|\bmbp\b',
    'thioredoxin': r'thioredoxin|\btrx',
    'reca': r'\brec ?a\b',
    'biotin_protein_ligase': r'\bbir ?a\b|biotin.{0,20}ligase',
    'luxr_autoinducer_regulator': r'\blux ?r\b|autoinducer|\btra ?r\b|\blas ?r\b',
    'glycosyl_hydrolase_68': r'levansucrase|\bsac ?b\b',
    'erm_rrna_methyltransferase': r'\berm|erythromycin|macrolide|rrna.{0,20}methyl',
    'initiator_repb': r'\brep ?b\b|replication (initiat|protein)|\brep\b|\brep ?a\b|initiator',
    'trfa_initiator': r'\btrf ?a\b',
    'dhfr_type_2': r'dihydrofolate|\bdhfr|trimethoprim|\bdfr',
}
_NAME_RE = {family: re.compile(pattern, re.I) for family, pattern in NAME_PATTERNS.items()}
PRIMARY_NAMES = ('product', 'gene', 'label', 'standard_name', 'allele')
TAG_RE = re.compile(r'fusion|fused|tag(ged)?\b|\bhis[- ]?6|6x?his|hexahis|polyhis|\bflag\b|\bha\b|\bmyc\b|\bv5\b|\bstrep', re.I)
SUBREGION_TYPES = {'misc_feature', 'mat_peptide', 'sig_peptide', 'transit_peptide', 'proprotein', 'regulatory'}


def _names(qualifiers, fields):
    return ' | '.join(v for k in fields for v in qualifiers.get(k, []))


def name_support(qualifiers):
    """Families named by the deposited text, separately for primary name fields and notes."""
    primary = _names(qualifiers, PRIMARY_NAMES)
    note = _names(qualifiers, ('note', 'function'))
    return ({f for f, rx in _NAME_RE.items() if rx.search(primary)},
            {f for f, rx in _NAME_RE.items() if rx.search(note)}, bool(TAG_RE.search(primary + ' ' + note)))


def _load_aligner():
    import importlib.util
    spec = importlib.util.spec_from_file_location('_component_panel', ROOT / 'experiments' / 'component_study_v1' / 'build_panel.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.Aligner()   # the frozen study's exact alignment settings


def _uniprot(run: Path):
    rows = {}
    for path in sorted((run / 'sources' / 'uniprot').glob('chunk_*.tsv')):
        for row in csv.DictReader(io.StringIO(path.read_text()), delimiter='\t'):
            rows[row['Entry']] = row
    return rows


def transcript_positions(segments, codon_start):
    """Forward-sequence coordinate of every coding nucleotide, in 5'->3' order of the CDS."""
    positions = []
    for part in segments:
        span = range(part['start0'], part['end0'])
        positions.extend(span if part['strand'] == 1 else reversed(span))
    return positions[codon_start - 1:]


def protein_to_spans(positions, start_aa, end_aa):
    """Zero-based half-open forward spans for residues [start_aa, end_aa), merged where contiguous."""
    covered = sorted(positions[3 * start_aa:3 * end_aa])
    spans = []
    for p in covered:
        if spans and spans[-1][1] == p:
            spans[-1][1] = p + 1
        else:
            spans.append([p, p + 1])
    return spans


def coding_dna(sequence: str, segments, codon_start):
    parts = []
    for part in segments:
        piece = sequence[part['start0']:part['end0']]
        if part['strand'] == -1:
            piece = piece[::-1].translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))
        parts.append(piece)
    return ''.join(parts)[codon_start - 1:]


def occurrence_anchor_evidence(cached_families, protein, positions, coding, table_id):
    """Per-occurrence anchor evidence for one CDS, without mutating the shared cache.

    ``unit_region_aa`` is the anchor hit in protein coordinates and is a property of the
    protein itself, so it is cached and shared by every construct that carries this exact
    protein. ``dna_spans`` and ``dna_region_retranslates`` are properties of *this* CDS
    occurrence's coordinates. They must be written into a fresh dict, never back into the
    cached entry: otherwise a later construct that repeats the same protein would overwrite
    the location saved for this one (the counts can still pass because they are tallied as
    each occurrence is computed, while the saved row keeps whichever write happened last).
    """
    from Bio.Seq import Seq

    evidence = {}
    for family, cached in cached_families.items():
        start_aa, end_aa = cached['unit_region_aa']
        peptide = str(Seq(coding[3 * start_aa:3 * end_aa]).translate(table=table_id))
        retranslates = (peptide[1:] == protein[start_aa + 1:end_aa] if start_aa == 0
                        else peptide == protein[start_aa:end_aa])
        evidence[family] = {**cached, 'unit_region_aa': list(cached['unit_region_aa']),
                            'dna_spans': protein_to_spans(positions, start_aa, end_aa),
                            'dna_region_retranslates': retranslates}
    return evidence


def audit(run: Path):
    from Bio import SeqIO
    from iri_constructs import process_inputs

    acquisition = read(run / 'acquisition.private.json')
    if not all(f['byte_identical'] for s in acquisition['samples'].values() for f in s['files']):
        print('note: at least one batch differs from the frozen bytes; it is audited as a new snapshot')
    table = read(FAMILY_TABLE)['families']
    uniprot = _uniprot(run)
    anchors = {a: {'family': family, 'sequence': uniprot[a]['Sequence']}
               for family, rule in table.items() for a in rule['anchors']}
    freeze_v1 = read(STUDY / 'panel.freeze.json')
    scorable = {f: row for f, row in freeze_v1['families'].items() if row['scorable']}
    references = {r['reference_id']: {'family': f, 'sequence': uniprot[r['reference_id'].removeprefix('uniprot-')]['Sequence']}
                  for f, row in scorable.items() for r in row['references']}
    for ref_id, row in references.items():
        expected = next(r['sequence_sha256'] for r in scorable[row['family']]['references'] if r['reference_id'] == ref_id)
        if sha(row['sequence']) != expected:
            raise ValueError(f'{ref_id}: reference sequence differs from the frozen panel')
    frozen_units = {}
    for _, folder in SAMPLES:
        for unit in read(folder / 'panel.freeze.json')['units']:
            frozen_units[unit['unit_id']] = unit
    aligner = _load_aligner()
    label_cache = {}

    def label(protein):
        if protein in label_cache:
            return label_cache[protein]
        contained, partial = {}, {}
        for accession, anchor in anchors.items():
            stats = aligner.stats(protein, anchor['sequence'])
            if not stats or stats['identity'] < IDENTITY:
                continue
            target = contained if stats['coverage_b'] >= COVERAGE else partial
            if anchor['family'] not in target or stats['identity'] > target[anchor['family']]['identity']:
                blocks = aligner.aligner.align(protein, anchor['sequence'])[0].aligned[0]
                target[anchor['family']] = {'anchor': accession, 'identity': round(stats['identity'], 4),
                                            'unit_coverage': round(stats['coverage_a'], 4),
                                            'anchor_coverage': round(stats['coverage_b'], 4),
                                            'unit_region_aa': [int(blocks[0][0]), int(blocks[-1][1])]}
        if len(contained) > 1:
            label_class = 'multi_family'
        elif contained:
            label_class = 'full' if next(iter(contained.values()))['unit_coverage'] >= COVERAGE else 'embedded'
        else:
            label_class = 'fragment' if partial else 'unknown'
        def near_copy(reference):
            stats = aligner.stats(protein, reference)
            return bool(stats) and stats['identity'] >= IDENTITY and max(stats['coverage_a'], stats['coverage_b']) >= COVERAGE

        # Same lineage rule as the frozen extension: a near copy of a panel reference is not scored.
        near = {f: any(near_copy(r['sequence']) for r in references.values() if r['family'] == f) for f in contained}
        label_cache[protein] = {'label_class': label_class, 'families': contained, 'partial': partial,
                                'near_copy_of_reference': near}
        return label_cache[protein]

    constructs, checks = [], Counter()
    for sample, _ in SAMPLES:
        for path in sorted((run / 'sources' / sample).glob('batch_*.gb')):
            result = process_inputs([path])
            sequences = {record.id: str(record.seq).upper() for record in SeqIO.parse(path, 'genbank')}
            for record in result.records:
                checks['records'] += 1
                if record.status != 'ok':
                    constructs.append({'accession': record.source_record_id, 'sample': sample, 'status': record.status,
                                       'failure': getattr(record.failure, 'code', None)})
                    continue
                meta = record.construct
                dna = sequences[record.source_record_id]
                checks['construct_sequence_hash_match'] += sha(dna) == meta['sequence_sha256']
                annotations = record.annotations['features']
                units = {u.feature_id: u for u in record.units}
                features = list(record.features)
                cds_rows, subregions = [], []
                for feature in features:
                    qualifiers = annotations[feature.feature_id]['qualifiers']
                    if feature.source_feature_type in SUBREGION_TYPES or (
                            feature.source_feature_type not in ('CDS', 'source', 'gene', 'mRNA', 'rep_origin',
                                                                'promoter', 'terminator', 'primer_bind')
                            and len(feature.segments) == 1):
                        primary, note, _ = name_support(qualifiers)
                        if primary or note:
                            subregions.append({'feature_index': feature.source_feature_index,
                                               'type': feature.source_feature_type, 'strand': feature.strand,
                                               'segments': [[s['start0'], s['end0']] for s in feature.segments],
                                               'location_status': feature.location_status,
                                               'named_families': sorted(primary | note)})
                    if feature.source_feature_type != 'CDS':
                        continue
                    unit = units[feature.feature_id]
                    primary, note, tagged = name_support(qualifiers)
                    row = {'feature_index': feature.source_feature_index, 'unit_id': unit.unit_id,
                           'strand': {'+': 1, '-': -1}.get(feature.strand),
                           'segments': [[s['start0'], s['end0']] for s in feature.segments],
                           'segment_strands': sorted({s['strand'] for s in feature.segments}),
                           'extraction_state': unit.extraction_state, 'reason_code': unit.reason_code,
                           'length_aa': unit.length_aa, 'codon_start': int(qualifiers.get('codon_start', ['1'])[0]),
                           'transl_table': int(qualifiers.get('transl_table', ['11'])[0]),
                           'deposited_names_families': sorted(primary),
                           'deposited_note_families': sorted(note), 'deposited_fusion_or_tag_term': tagged,
                           'product_text_sha256': sha(_names(qualifiers, PRIMARY_NAMES))}
                    frozen = frozen_units.get(unit.unit_id)
                    if unit.extraction_state == 'exported':
                        protein = unit.sequence
                        row['protein_sha256'] = sha(protein)
                        row['frozen_unit_match'] = bool(frozen) and frozen['sequence_sha256'] == row['protein_sha256']
                        checks['exported_units'] += 1
                        checks['frozen_unit_match'] += row['frozen_unit_match']
                        evidence = label(protein)
                        row['label_class'] = evidence['label_class']
                        row['partial_anchor_evidence'] = evidence['partial']
                        row['near_copy_of_reference'] = evidence['near_copy_of_reference']
                        if frozen:
                            same = (frozen['label_class'] == evidence['label_class'] and
                                    sorted(frozen['families']) == sorted(evidence['families']))
                            row['frozen_label_match'] = same
                            checks['frozen_label_match'] += same
                        # The anchor region in protein coordinates (unit_region_aa) is a property of the
                        # protein and is cached. The DNA spans and the re-translation depend on THIS CDS
                        # occurrence's coordinates, so build a fresh per-occurrence dict and never write
                        # back into the shared cache: a repeated protein in another construct must not
                        # overwrite this record's saved location.
                        positions = transcript_positions(feature.segments, row['codon_start'])
                        coding = coding_dna(dna, feature.segments, row['codon_start'])
                        anchor_evidence = occurrence_anchor_evidence(evidence['families'], protein, positions, coding, row['transl_table'])
                        row['anchor_evidence'] = anchor_evidence
                        for ev in anchor_evidence.values():
                            checks['anchor_regions'] += 1
                            checks['anchor_regions_retranslate'] += ev['dna_region_retranslates']
                    cds_rows.append(row)
                constructs.append({'accession': record.source_record_id, 'sample': sample, 'status': 'ok',
                                   'batch_file': path.name, 'length': len(dna), 'topology': meta['topology'],
                                   'sequence_sha256': sha(dna), 'description_sha256': sha(meta['description'] or ''),
                                   'cds': cds_rows, 'named_subregions': subregions})
    audit_doc = {'schema': AUDIT_SCHEMA, 'thresholds': {'identity': IDENTITY, 'anchor_coverage': COVERAGE},
                 'scorable_families': {f: [r['reference_id'] for r in row['references']] for f, row in scorable.items()},
                 'checks': dict(checks), 'constructs': constructs}
    write(run / 'audit.private.json', audit_doc, exclusive=False)
    print(json.dumps(dict(checks), indent=1))


# ---------------------------------------------------------------------- shortlist

def _overlap(a, b):
    return sum(max(0, min(y, d) - max(x, c)) for x, y in a for c, d in b)


def _bases(spans):
    return sum(end - start for start, end in spans)


def _cds_transcript_index(cds):
    """Forward coordinate -> position in the CDS coding sequence (5'->3', after codon_start)."""
    strand = cds['strand']
    ranges = (range(s, e) for s, e in cds['segments'])
    coords = [p for r in ranges for p in (r if strand == 1 else reversed(r))]
    return {p: i for i, p in enumerate(coords[cds['codon_start'] - 1:])}


def deposited_subregion_boundary(construct, cds, family, anchor_region):
    """A deposited sub-feature inside this CDS that names the family and gives a usable location.

    This boundary comes from the record itself, independent of the anchor alignment. Returns the
    best such feature (by overlap with the anchor region) with whether it is codon-aligned, long
    enough to score, and consistent with the anchor region, or None.
    """
    index = _cds_transcript_index(cds)
    strand_text = {1: '+', -1: '-'}[cds['strand']]
    best = None
    for sub in construct['named_subregions']:
        if family not in sub['named_families'] or sub['strand'] not in (strand_text, '') or sub['location_status'] != 'exact':
            continue
        covered = [p for a, b in sub['segments'] for p in range(a, b)]
        if not covered or any(p not in index for p in covered):
            continue  # a named feature that is not wholly inside this CDS is not a subregion boundary for it
        first, last = min(index[p] for p in covered), max(index[p] for p in covered)
        residues = [first // 3, (last + 1) // 3]
        spans = sorted([list(s) for s in sub['segments']])
        candidate = {
            'feature_index': sub['feature_index'], 'type': sub['type'], 'spans': spans, 'residues': residues,
            'codon_aligned': first % 3 == 0 and (last + 1) % 3 == 0,
            'long_enough': residues[1] - residues[0] >= MIN_RESIDUES,
            'overlaps_anchor_region': _overlap(spans, anchor_region) >= COVERAGE * min(_bases(spans), _bases(anchor_region))}
        candidate['usable'] = candidate['codon_aligned'] and candidate['long_enough']
        if best is None or _overlap(spans, anchor_region) > _overlap(best['spans'], anchor_region):
            best = candidate
    return best


def eligible_components(construct, scorable):
    """Whole-CDS components with a validated translation, an independent family label and no name conflict."""
    components = []
    for cds in construct['cds']:
        if cds['extraction_state'] != 'exported' or cds['strand'] is None or len(cds['segment_strands']) != 1:
            continue
        deposited = set(cds['deposited_names_families']) | set(cds['deposited_note_families'])
        for family, ev in sorted(cds.get('anchor_evidence', {}).items()):
            if family not in scorable or cds['near_copy_of_reference'].get(family):
                continue
            if family in cds['deposited_names_families']:
                name = 'deposited_name_agrees'
            elif family in cds['deposited_note_families']:
                name = 'deposited_note_agrees'
            elif deposited and cds['label_class'] != 'multi_family':
                continue  # deposited name names a different family for this single-label CDS
            else:
                name = 'deposited_name_partner' if deposited else 'deposited_name_silent'
            fusion = cds['label_class'] in ('embedded', 'multi_family')
            component = {'feature_index': cds['feature_index'], 'family': family, 'protein_sha256': cds['protein_sha256'],
                         'spans': cds['segments'], 'strand': cds['strand'], 'length_aa': cds['length_aa'],
                         'label_class': cds['label_class'], 'name_evidence': name, 'fusion_class': fusion,
                         'over_model_limit': cds['length_aa'] > MAX_RESIDUES,
                         'anchor': {k: ev[k] for k in ('anchor', 'identity', 'unit_coverage', 'anchor_coverage')}}
            if fusion:
                region = ev['dna_spans']
                component['anchor_region'] = {'spans': region, 'residues': ev['unit_region_aa'],
                                              'retranslates': ev['dna_region_retranslates'],
                                              'usable': ev['unit_region_aa'][1] - ev['unit_region_aa'][0] >= MIN_RESIDUES
                                              and ev['dna_region_retranslates']}
                component['deposited_subregion'] = deposited_subregion_boundary(construct, cds, family, region)
            components.append(component)
    return components


def archetype_of(construct, components):
    cds_with_label = {c['feature_index'] for c in components}
    fusion = [c for c in components if c['fusion_class']]
    labelled_bases = sum(_bases(cds['segments']) for cds in construct['cds'] if cds['feature_index'] in cds_with_label)
    fraction = labelled_bases / construct['length']
    if not components:
        return None
    if fusion:
        return 'tagged_fusion'
    if len(cds_with_label) >= 2:
        return 'multi_cds'
    if fraction < 0.15 and len(construct['cds']) - len(cds_with_label) >= 1:
        return 'backbone_heavy'
    if len(construct['cds']) == 1:
        return 'single_cds'
    return 'single_label_other'


def _groups(accessions, protein_sets):
    """Connected components of records that share any eligible-component protein hash."""
    parent = {a: a for a in accessions}

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    owner = {}
    for accession in accessions:
        for protein in protein_sets[accession]:
            other = owner.setdefault(protein, accession)
            parent[root(accession)] = root(other)
    groups = defaultdict(list)
    for accession in accessions:
        groups[root(accession)].append(accession)
    sizes = sorted((len(g) for g in groups.values()), reverse=True)
    return {'records': len(accessions), 'related_groups': len(groups),
            'largest_group': sizes[0] if sizes else 0, 'singleton_groups': sum(1 for s in sizes if s == 1)}


def shortlist(run: Path):
    """Metadata-only shortlist: independently reviewed benign records, with dependence-aware group counts."""
    audited = read(run / 'audit.private.json')
    acquisition = read(run / 'acquisition.private.json')
    scorable = set(audited['scorable_families'])
    pilot = {row['accession'] for row in read(Path(__file__).with_name('reporter_panel.json'))['constructs']}

    protein_sets, rows, by_archetype = {}, [], defaultdict(set)
    documented_fusion_records, anchor_only_fusion_records = set(), set()
    documented_fusion_components = anchor_only_fusion_components = 0
    whole_component_records = set()
    for construct in audited['constructs']:
        if construct['status'] != 'ok' or construct['accession'] in pilot:
            continue
        components = eligible_components(construct, scorable)
        if not components:
            continue
        archetype = archetype_of(construct, components)
        accession = construct['accession']
        protein_sets[accession] = {c['protein_sha256'] for c in components}
        whole_component_records.add(accession)
        by_archetype[archetype].add(accession)
        fusion_rows = []
        for component in components:
            if not component['fusion_class']:
                continue
            boundary = component['deposited_subregion']
            documented = bool(boundary and boundary['usable'] and boundary['overlaps_anchor_region'])
            if documented:
                documented_fusion_components += 1
                documented_fusion_records.add(accession)
            else:
                anchor_only_fusion_components += 1
                anchor_only_fusion_records.add(accession)
            fusion_rows.append({'family': component['family'], 'documented_boundary': documented,
                                'anchor_region_usable': component['anchor_region']['usable']})
        rows.append({'accession': accession, 'sample': construct['sample'], 'archetype': archetype,
                     'eligible_components': len(components), 'fusion_components': fusion_rows})

    # A record can hold both a documented-boundary fusion and an anchor-only one; report components and records.
    categories = {
        'independently_reviewed_records': _groups(sorted(whole_component_records), protein_sets),
        'documented_fusion_boundary': _groups(sorted(documented_fusion_records), protein_sets),
        'anchor_only_fusion': _groups(sorted(anchor_only_fusion_records), protein_sets),
    }
    by_archetype_summary = {a: _groups(sorted(accs), protein_sets) for a, accs in sorted(by_archetype.items())}
    whole_corpus_protein_sets = {}
    for construct in audited['constructs']:
        if construct['status'] != 'ok':
            continue
        whole_corpus_protein_sets[construct['accession']] = {
            cds['protein_sha256'] for cds in construct['cds'] if cds['extraction_state'] == 'exported'}
    report = {
        'schema': 'paired-corpus-shortlist-1',
        'note': 'Metadata-only shortlist. Benign public vectors as deposited; no sequence or qualifier text here.',
        'terms': {'ncbi': acquisition.get('ncbi_terms'), 'uniprot': acquisition.get('uniprot_terms')},
        'independent_label_source': 'UniProtKB/Swiss-Prot reviewed family anchors (identity >= 0.90, anchor '
                                    'coverage >= 0.80); deposited CDS names recorded separately.',
        'audit_sha256': digest(run / 'audit.private.json'),
        'excluded_reporter_pilot_records': sorted(pilot),
        'component_counts': {'documented_fusion_boundary': documented_fusion_components,
                             'anchor_only_fusion': anchor_only_fusion_components},
        'categories': categories,
        'by_archetype': by_archetype_summary,
        'dependence_context': {
            'note': 'Records sharing any identical exported protein are connected; common markers link most '
                    'vectors, so record counts overstate independent cases.',
            'whole_corpus': _groups(sorted(whole_corpus_protein_sets), whole_corpus_protein_sets),
            'shortlist_pool': categories['independently_reviewed_records']},
        'records': sorted(rows, key=lambda r: r['accession']),
    }
    write(run / 'shortlist.private.json', report, exclusive=False)
    summary = {k: {'records': v['records'], 'related_groups': v['related_groups'], 'largest_group': v['largest_group']}
               for k, v in categories.items()}
    print(json.dumps({'component_counts': report['component_counts'], 'categories': summary,
                      'by_archetype': {a: g['records'] for a, g in by_archetype_summary.items()},
                      'whole_corpus_groups': report['dependence_context']['whole_corpus']}, indent=1))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    for stage in ('acquire', 'audit', 'shortlist'):
        command = sub.add_parser(stage)
        command.add_argument('--run', type=Path, default=ROOT / 'results' / 'paired-corpus-batch-v1')
    args = parser.parse_args()
    if args.stage == 'acquire':
        acquire(args.run)
    elif args.stage == 'audit':
        audit(args.run)
    else:
        shortlist(args.run)


if __name__ == '__main__':
    main()
