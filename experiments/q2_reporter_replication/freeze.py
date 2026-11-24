"""Freeze the Q2 reporter-replication panel, reference set, and grouping before scoring.

Reference sequences come from reviewed UniProtKB/Swiss-Prot entries, one or more per
family, chosen as canonical family representatives independent of the deposited
constructs. The panel admits an audited construct only when every documented CDS maps
to a family that has a reference here, so the answer key never contains a component that
no reference could recover. Families without a suitable independent reviewed reference
(notably the engineered NanoLuc) are therefore excluded, and the excluded constructs are
reported as a shortfall.

Reference sequences are written to the ignored results tree; a sequence-free manifest
(accession, family, length, SHA-256) is committed.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

TOOL = 'iri_constructs_q2_reporter'
UNIPROT = 'https://rest.uniprot.org/uniprotkb/search'

# Reviewed Swiss-Prot references per family (verified present and canonical). gfp_like
# carries both Aequorea green and Anthozoa red representatives because the deposited
# "fluorescent protein" family spans both GFP-fold subtypes at low mutual identity.
REFERENCE_ACCESSIONS = {
    'gfp_like': ['P42212', 'Q9U6Y8', 'Q8ISF8'],
    'firefly_luciferase': ['P08659'],
    'renilla_luciferase': ['P27652'],
    'beta_galactosidase': ['P00722'],
    'beta_lactamase': ['P62593'],
    'aminoglycoside_phosphotransferase': ['P00552', 'P00551'],
    'hygromycin_phosphotransferase': ['P00557'],
    'aminoglycoside_nucleotidyltransferase': ['P0AG05'],
    'puromycin_n_acetyltransferase': ['P13249'],
}
REFERENCED_FAMILIES = set(REFERENCE_ACCESSIONS)


def read_json(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def write_json(path: Path, value, *, exist_ok=False):
    with path.open('w' if exist_ok else 'x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def _get(url: str, timeout=30, retries=6) -> bytes:
    from urllib.error import HTTPError, URLError
    delay = 1.0
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': TOOL}), timeout=timeout) as r:
                return r.read()
        except (HTTPError, URLError):
            if attempt == retries - 1:
                raise
            time.sleep(delay)
            delay = min(delay * 2, 30.0)
    raise RuntimeError('unreachable')


def build_references(out_private: Path, out_manifest: Path):
    accession_to_family = {acc: fam for fam, accs in REFERENCE_ACCESSIONS.items() for acc in accs}
    accessions = sorted(accession_to_family)
    query = 'reviewed:true AND (' + ' OR '.join(f'accession:{a}' for a in accessions) + ')'
    params = urllib.parse.urlencode(dict(query=query, format='tsv', size=500,
                                         fields='accession,protein_name,organism_name,length,sequence_version,sequence'))
    raw = _get(f'{UNIPROT}?{params}')
    import csv
    rows = list(csv.DictReader(io.StringIO(raw.decode()), delimiter='\t'))
    by_acc = {r['Entry']: r for r in rows}
    missing = [a for a in accessions if a not in by_acc]
    if missing:
        raise ValueError(f'UniProt did not return references: {missing}')
    references, manifest = [], []
    for acc in accessions:
        row = by_acc[acc]
        sequence = row['Sequence']
        references.append({'accession': acc, 'family': accession_to_family[acc],
                           'protein_name': row['Protein names'], 'organism': row['Organism'],
                           'length': int(row['Length']), 'sequence_version': int(row['Sequence version']),
                           'sequence': sequence,
                           'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest()})
        manifest.append({'accession': acc, 'family': accession_to_family[acc],
                         'length': int(row['Length']), 'sequence_version': int(row['Sequence version']),
                         'sequence_sha256': hashlib.sha256(sequence.encode()).hexdigest()})
    out_bytes = (json.dumps(references, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
    out_private.write_bytes(out_bytes)
    write_json(out_manifest, {'tool': TOOL, 'provider': 'UniProtKB/Swiss-Prot',
                              'query': query, 'response_sha256': hashlib.sha256(raw).hexdigest(),
                              'prepared_sha256': hashlib.sha256(out_bytes).hexdigest(),
                              'families': sorted(REFERENCED_FAMILIES), 'references': manifest}, exist_ok=True)
    print(f'wrote {len(references)} references across {len(REFERENCED_FAMILIES)} families; '
          f'prepared_sha256={hashlib.sha256(out_bytes).hexdigest()}')


def build_panel(audit: Path, fetch_manifest: Path, out: Path):
    records = read_json(audit)['records']
    shas = read_json(fetch_manifest)['fetched']
    admitted, excluded = [], []
    for r in records:
        if not r.get('admit'):
            continue
        families = {c['family'] for c in r['cds']}
        if not families <= REFERENCED_FAMILIES:
            excluded.append({'accession': r['accession'],
                             'reason': 'family without independent reference: '
                                       + ','.join(sorted(families - REFERENCED_FAMILIES))})
            continue
        admitted.append(r)

    admitted.sort(key=lambda r: r['accession'])
    constructs = []
    for i, r in enumerate(admitted, 1):
        acc = r['accession']
        labels = [{'family': c['family'], 'feature_index': c['feature_index'],
                   'source_product': c['product']} for c in sorted(r['cds'], key=lambda c: c['feature_index'])]
        constructs.append({'accession': acc, 'id': f'rc-{i:03d}', 'archetype': r['archetype'],
                           'sha256': shas[acc]['sha256'], 'labels': labels})

    panel = {
        'study': 'q2-reporter-replication',
        'construct_source': 'NCBI GenBank nuccore EFetch, rettype=gb, retmode=text',
        'construct_terms_url': 'https://www.ncbi.nlm.nih.gov/home/about/policies/',
        'reference_manifest': 'references.manifest.json',
        'referenced_families': sorted(REFERENCED_FAMILIES),
        'constructs': constructs,
        'excluded_admitted_for_missing_reference': excluded,
    }
    write_json(out, panel, exist_ok=True)
    reporters = sum(1 for c in constructs for l in c['labels']
                    if l['family'] in {'gfp_like', 'firefly_luciferase', 'renilla_luciferase', 'beta_galactosidase'})
    components = sum(len(c['labels']) for c in constructs)
    print(f'panel: {len(constructs)} constructs, {components} components ({reporters} reporter), '
          f'{len(excluded)} admitted constructs excluded for missing reference')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)

    r = sub.add_parser('references')
    r.add_argument('--out-private', type=Path, required=True)
    r.add_argument('--out-manifest', type=Path, required=True)

    p = sub.add_parser('panel')
    p.add_argument('--audit', type=Path, required=True)
    p.add_argument('--fetch-manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)

    args = parser.parse_args()
    if args.stage == 'references':
        build_references(args.out_private, args.out_manifest)
    else:
        build_panel(args.audit, args.fetch_manifest, args.out)


if __name__ == '__main__':
    main()
