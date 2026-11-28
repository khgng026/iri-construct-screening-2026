"""Discover and audit fresh public benign reporter constructs for the Q2 replication.

This module does NOT assign a frozen panel by itself. It fetches candidate GenBank
records from NCBI, validates their coding translations with the existing toolkit, maps
each documented CDS to a family label from its deposited ``/product`` qualifier, and
writes a sequence-free audit. A human then reviews the audit and freezes the panel.

Family labels come from the depositor's annotation (documentation), never from the
alignment/model scorer under test. Translations are validated against the deposited DNA
by ``iri_constructs.process_inputs``. Source sequences are written only to the ignored
results tree; the audit JSON is sequence-free (lengths and SHA-256 only).

Stages
------
discover : run the committed esearch queries, resolve accessions, drop excluded ones.
fetch    : efetch GenBank for candidate accessions into the ignored sources directory.
audit    : parse + validate each record, map families, decide admissibility, write audit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
TOOL = 'iri_constructs_q2_reporter'
MIN_RESIDUES = 30
MAX_RESIDUES = 1022

# Curated esearch queries over NCBI nuccore. Restricted to plausible plasmid/vector
# lengths. The synthetic division filter (gbdiv_syn) is applied where it does not
# collapse the family to nothing. Each query targets a distinct reporter ancestry so the
# panel spans more than one lineage.
SEARCH_QUERIES = {
    'gfp_syn': '"green fluorescent protein"[All Fields] AND plasmid[Title] AND 1800:9000[SLEN] AND gbdiv_syn[PROP]',
    'fp_vector': '"fluorescent protein"[Title] AND vector[Title] AND 1800:12000[SLEN]',
    'mcherry_rfp': 'mCherry[Title] AND vector[Title] AND 1800:12000[SLEN]',
    'rfp_vector': '"red fluorescent protein"[Title] AND vector[Title] AND 1800:12000[SLEN]',
    'luciferase': 'luciferase[Title] AND (vector[Title] OR plasmid[Title]) AND 1800:12000[SLEN]',
    'renilla_luc': '"Renilla luciferase"[All Fields] AND plasmid[Title] AND 1800:12000[SLEN]',
    'lacz_gus': '(lacZ OR beta-galactosidase OR beta-glucuronidase)[Title] AND vector[Title] AND 1800:14000[SLEN]',
    'reporter_cassette': 'reporter[Title] AND (plasmid[Title] OR vector[Title]) AND 1800:12000[SLEN]',
}

# Ordered (compiled substring/regex, family, role). First match wins. Matching is
# case-insensitive over the ``/product`` string, with the ``/gene`` qualifier as a
# fallback. Roles: reporter | marker | regulator | backbone.
_LEXICON_SPEC = [
    # --- reporters: fluorescent protein superfamily (shared GFP fold) ---
    (r'green fluorescent|superfolder|\bsfgfp\b|\begfp\b|\bgfp\b', 'gfp_like', 'reporter'),
    (r'yellow fluorescent|\byfp\b|\bvenus\b|citrine', 'gfp_like', 'reporter'),
    (r'cyan fluorescent|\bcfp\b|cerulean|turquoise', 'gfp_like', 'reporter'),
    (r'blue fluorescent|\bbfp\b|azurite|tagbfp', 'gfp_like', 'reporter'),
    (r'red fluorescent|\brfp\b|mcherry|mscarlet|dsred|tdtomato|td-tomato|mruby|mkate|tagrfp|morange|mplum', 'gfp_like', 'reporter'),
    (r'mneongreen|neongreen|staygold|mstaygold|emerald|clover|mgreenlantern', 'gfp_like', 'reporter'),
    (r'fluorescent protein', 'gfp_like', 'reporter'),
    # --- reporters: bioluminescent (distinct families) ---
    (r'firefly luciferase|photinus|\bluc2\b|\bffluc\b', 'firefly_luciferase', 'reporter'),
    (r'renilla luciferase|\brluc\b|renilla reniformis luciferase', 'renilla_luciferase', 'reporter'),
    (r'nanoluc|nano-?luc|\bnluc\b|nanoluciferase', 'nanoluc', 'reporter'),
    # --- reporters: enzymatic ---
    (r'beta-?galactosidase|β-?galactosidase|\blacz\b', 'beta_galactosidase', 'reporter'),
    (r'beta-?glucuronidase|β-?glucuronidase|glucuronidase|\bgusa\b|\buida\b', 'beta_glucuronidase', 'reporter'),
    # --- selection markers (common co-occurring context / distractors) ---
    (r'beta-?lactamase|β-?lactamase|\bbla\b|ampicillin|\bampr\b|tem-?1|carbenicillin', 'beta_lactamase', 'marker'),
    (r'aminoglycoside.*phosphotransferase|neomycin|kanamycin|\bkanr\b|\bneo\b|\bnptii?\b|aph\(', 'aminoglycoside_phosphotransferase', 'marker'),
    (r'hygromycin|\bhph\b|\bhygr?\b', 'hygromycin_phosphotransferase', 'marker'),
    (r'puromycin|\bpac\b', 'puromycin_n_acetyltransferase', 'marker'),
    (r'spectinomycin|streptomycin|\baada\b', 'aminoglycoside_nucleotidyltransferase', 'marker'),
    (r'bleomycin|zeocin|sh ble|\bble\b', 'bleomycin_resistance', 'marker'),
    (r'chloramphenicol acetyltransferase|\bcat\b', 'chloramphenicol_acetyltransferase', 'marker'),
    # --- regulators / expression machinery ---
    (r'arabinose operon regulatory|arabinose-binding|\barac\b', 'araC', 'regulator'),
    (r'lac repressor|lactose operon repressor|\blaci\b', 'lacI', 'regulator'),
    (r'tetracycline repressor|\btetr\b', 'tetR', 'regulator'),
    (r't7 rna polymerase', 't7_rna_polymerase', 'regulator'),
    # --- backbone maintenance ---
    (r'repressor of primer|\brop\b|\brom\b', 'rop_rom', 'backbone'),
]
LEXICON = [(re.compile(p, re.IGNORECASE), fam, role) for p, fam, role in _LEXICON_SPEC]


def read_json(path: Path):
    return json.loads(path.read_text(encoding='utf-8'))


def write_json(path: Path, value, *, exist_ok=False):
    mode = 'w' if exist_ok else 'x'
    with path.open(mode, encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def map_family(product: str | None, gene: str | None):
    for text in (product or '', gene or ''):
        for pattern, family, role in LEXICON:
            if pattern.search(text):
                return family, role
    return None, None


def _get(url: str, timeout=30, retries=6) -> bytes:
    from urllib.error import HTTPError, URLError
    delay = 1.0
    for attempt in range(retries):
        try:
            with urlopen(Request(url, headers={'User-Agent': TOOL}), timeout=timeout) as response:
                return response.read()
        except HTTPError as error:
            if error.code in (429, 500, 502, 503) and attempt < retries - 1:
                time.sleep(delay)
                delay = min(delay * 2, 30.0)
                continue
            raise
        except URLError:
            if attempt < retries - 1:
                time.sleep(delay)
                delay = min(delay * 2, 30.0)
                continue
            raise
    raise RuntimeError('unreachable')


def esearch(term: str, retmax=200) -> tuple[int, list[str]]:
    q = urlencode(dict(db='nuccore', term=term, retmax=retmax, retmode='json', tool=TOOL))
    data = json.loads(_get(EUTILS + 'esearch.fcgi?' + q))['esearchresult']
    return int(data['count']), data.get('idlist', [])


def esummary(uids: list[str]) -> dict[str, str]:
    """Map NCBI UID -> accession.version for a batch of UIDs."""
    out = {}
    for i in range(0, len(uids), 200):
        batch = uids[i:i + 200]
        q = urlencode(dict(db='nuccore', id=','.join(batch), retmode='json', tool=TOOL))
        result = json.loads(_get(EUTILS + 'esummary.fcgi?' + q))['result']
        for uid in result.get('uids', []):
            out[uid] = result[uid].get('accessionversion') or result[uid].get('caption')
        time.sleep(0.4)
    return out


def load_exclusions(paths: list[Path]) -> set[str]:
    excluded = set()
    for path in paths:
        if path.suffix == '.txt':
            excluded.update(line.strip() for line in path.read_text().splitlines() if line.strip())
        elif path.name == 'reporter_panel.json':
            excluded.update(c['accession'] for c in read_json(path)['constructs'])
    # normalise to base accession (drop version) as well, for safety
    excluded.update({a.split('.')[0] for a in list(excluded)})
    return excluded


def discover(exclusions: set[str], out: Path):
    found = {}
    report = {}
    for name, term in SEARCH_QUERIES.items():
        count, uids = esearch(term)
        acc = esummary(uids)
        kept = {}
        for uid, accession in acc.items():
            if not accession:
                continue
            base = accession.split('.')[0]
            if accession in exclusions or base in exclusions:
                continue
            kept[accession] = name
        found.update(kept)
        report[name] = {'total_hits': count, 'resolved': len(acc), 'kept_after_exclusion': len(kept)}
        time.sleep(0.4)
    payload = {'tool': TOOL, 'queries': SEARCH_QUERIES, 'per_query': report,
               'candidates': sorted(found), 'candidate_source_query': found}
    write_json(out, payload, exist_ok=True)
    return payload


def fetch(accessions: list[str], sources: Path):
    sources.mkdir(parents=True, exist_ok=True)
    manifest = {}
    for accession in accessions:
        target = sources / f'{accession}.gb'
        if target.exists():
            data = target.read_bytes()
        else:
            q = urlencode(dict(db='nuccore', id=accession, rettype='gb', retmode='text', tool=TOOL))
            data = _get(EUTILS + 'efetch.fcgi?' + q)
            target.write_bytes(data)
            time.sleep(0.5)
        manifest[accession] = {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}
    return manifest


def _infer_archetype(cds_rows, length) -> str:
    reporters = [r for r in cds_rows if r['role'] == 'reporter']
    n = len(cds_rows)
    coding = sum(r['protein_len'] * 3 for r in cds_rows)
    if n == 1:
        return 'single_cds'
    if any(r['protein_len'] > 600 and r['role'] == 'reporter' for r in reporters) and n <= 2:
        return 'tagged_fusion'
    if length and coding / length < 0.45:
        return 'backbone_heavy'
    return 'multi_cds'


def audit_record(gb_path: Path, exclusions: set[str]):
    from Bio import SeqIO
    from iri_constructs import process_inputs

    record = SeqIO.read(gb_path, 'genbank')
    accession = record.id
    result = process_inputs([gb_path])
    reasons = []
    if len(result.ok_records) != 1:
        return {'accession': accession, 'admit': False, 'reasons': ['record did not parse as one ok record'],
                'length': len(record.seq), 'cds': []}
    validated = result.ok_records[0]
    units = {u.feature_id: u for u in validated.units}
    cds_rows = []
    for feature in validated.features:
        if feature.source_feature_type != 'CDS':
            continue
        unit = units.get(feature.feature_id)
        product = feature.names.get('product')
        gene = feature.names.get('gene')
        family, role = map_family(product, gene)
        exported = unit is not None and unit.extraction_state == 'exported'
        protein_len = unit.length_aa if (unit and unit.length_aa) else 0
        cds_rows.append({
            'feature_index': feature.source_feature_index,
            'product': product, 'gene': gene, 'strand': feature.strand,
            'family': family, 'role': role,
            'exported': exported,
            'extraction_state': unit.extraction_state if unit else None,
            'protein_len': protein_len,
            'protein_sha256': unit.sequence_sha256 if (unit and exported) else None,
            'model_eligible': exported and MIN_RESIDUES <= protein_len <= MAX_RESIDUES,
        })

    base = accession.split('.')[0]
    n_cds = len(cds_rows)
    exported_rows = [r for r in cds_rows if r['exported']]
    labelled = [r for r in cds_rows if r['family']]
    reporters = [r for r in cds_rows if r['role'] == 'reporter' and r['exported']]
    unlabelled_cds = [r for r in cds_rows if not r['family']]
    withheld_cds = [r for r in cds_rows if not r['exported']]

    # Strict admissibility: reuse the tested prepare path, which requires every documented
    # CDS to be exported and accounted for. Admit only when every CDS is exported and
    # family-labelled, at least one is a reporter, and the record is fresh.
    if accession in exclusions or base in exclusions:
        reasons.append('accession in exclusion set')
    if n_cds == 0:
        reasons.append('no documented CDS')
    if withheld_cds:
        reasons.append(f'{len(withheld_cds)} CDS not exported by toolkit')
    if unlabelled_cds:
        reasons.append(f'{len(unlabelled_cds)} CDS without a family label')
    if not reporters:
        reasons.append('no exported reporter CDS')
    if not any(r['model_eligible'] for r in reporters):
        reasons.append('no model-eligible reporter CDS (length)')

    admit = not reasons
    return {
        'accession': accession, 'admit': admit, 'reasons': reasons,
        'length': len(record.seq), 'n_cds': n_cds,
        'n_exported': len(exported_rows), 'n_labelled': len(labelled),
        'n_reporters': len(reporters),
        'archetype': _infer_archetype([r for r in cds_rows if r['exported']], len(record.seq)),
        'reporter_families': sorted({r['family'] for r in reporters}),
        'cds': cds_rows,
    }


def audit(sources: Path, exclusions: set[str], out: Path):
    records = []
    for gb_path in sorted(sources.glob('*.gb')):
        try:
            records.append(audit_record(gb_path, exclusions))
        except Exception as error:  # noqa: BLE001 - record-level isolation, reported individually
            records.append({'accession': gb_path.stem, 'admit': False,
                            'reasons': [f'audit error: {type(error).__name__}: {error}'], 'cds': []})
    admitted = [r for r in records if r['admit']]
    # Transparent duplicate accounting by reporter-protein identity.
    protein_groups: dict[str, list[str]] = {}
    for r in admitted:
        for row in r['cds']:
            if row['role'] == 'reporter' and row['protein_sha256']:
                protein_groups.setdefault(row['protein_sha256'], []).append(r['accession'])
    family_counts: dict[str, int] = {}
    role_counts: dict[str, int] = {}
    for r in admitted:
        for row in r['cds']:
            family_counts[row['family']] = family_counts.get(row['family'], 0) + 1
            role_counts[row['role'] or 'unlabelled'] = role_counts.get(row['role'] or 'unlabelled', 0) + 1
    summary = {
        'records_audited': len(records),
        'records_admitted': len(admitted),
        'admitted_accessions': sorted(r['accession'] for r in admitted),
        'archetypes': {a: sum(1 for r in admitted if r.get('archetype') == a)
                       for a in sorted({r.get('archetype') for r in admitted})},
        'reporter_families_present': sorted({f for r in admitted for f in r.get('reporter_families', [])}),
        'all_family_component_counts': dict(sorted(family_counts.items())),
        'role_component_counts': dict(sorted(role_counts.items())),
        'distinct_reporter_proteins': len(protein_groups),
        'repeated_reporter_proteins': {sha: accs for sha, accs in sorted(protein_groups.items()) if len(accs) > 1},
        'rejection_reason_counts': _reason_counts(records),
    }
    write_json(out, {'tool': TOOL, 'summary': summary, 'records': records}, exist_ok=True)
    return summary


def _reason_counts(records):
    counts: dict[str, int] = {}
    for r in records:
        if r['admit']:
            continue
        for reason in r['reasons']:
            key = re.sub(r'\d+', 'N', reason)
            counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items(), key=lambda kv: -kv[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)

    d = sub.add_parser('discover')
    d.add_argument('--exclude', type=Path, action='append', default=[])
    d.add_argument('--out', type=Path, required=True)

    f = sub.add_parser('fetch')
    f.add_argument('--candidates', type=Path, required=True)
    f.add_argument('--sources', type=Path, required=True)
    f.add_argument('--limit', type=int, default=None)
    f.add_argument('--out', type=Path, required=True)

    a = sub.add_parser('audit')
    a.add_argument('--sources', type=Path, required=True)
    a.add_argument('--exclude', type=Path, action='append', default=[])
    a.add_argument('--out', type=Path, required=True)

    args = parser.parse_args()
    if args.stage == 'discover':
        payload = discover(load_exclusions(args.exclude), args.out)
        print(f'discovered {len(payload["candidates"])} candidate accessions')
    elif args.stage == 'fetch':
        candidates = read_json(args.candidates)['candidates']
        if args.limit:
            candidates = candidates[:args.limit]
        manifest = fetch(candidates, args.sources)
        write_json(args.out, {'tool': TOOL, 'fetched': manifest}, exist_ok=True)
        print(f'fetched {len(manifest)} records into {args.sources}')
    elif args.stage == 'audit':
        summary = audit(args.sources, load_exclusions(args.exclude), args.out)
        print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
