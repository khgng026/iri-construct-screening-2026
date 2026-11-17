"""Repeat the reviewed panel-v2 E1 helper audit without writing biological inputs."""
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from iri_constructs.annotation_quality import audit_quality, summarize_quality
from iri_constructs.annotations import json_safe, normalize_record
from iri_constructs.dataset_curation import parse_genbank

PANEL_HASH = 'f9111dd17ec028fdc170552aacc8af7275fccaf46d1e9712e3f7a681bdc0d70b'
REVIEW = ROOT / 'docs/coordination/PANEL_V2_AI_REVIEW_2026_11_13.md'


def main():
    panel = json.loads((ROOT / 'data/manifests/prospective-2026-11-12/panel.v2.json').read_text())
    proposed = panel.pop('proposal_sha256')
    digest = hashlib.sha256(json.dumps(panel, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if proposed != PANEL_HASH or digest != PANEL_HASH:
        raise ValueError('Panel checksum mismatch')
    review = REVIEW.read_bytes()
    if b'E1: CLEAR for the two exact sources' not in review:
        raise ValueError('Missing recorded E1 decision')
    results = []
    for item in panel['new_E1_only_sources']:
        acc = item['accession']
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode(
            {'db': 'nuccore', 'id': acc, 'rettype': 'gbwithparts', 'retmode': 'text'})
        with urlopen(Request(url, headers={'User-Agent': 'IRI-project-research/1.0'}), timeout=30) as response:
            data = response.read(1_500_001)
        if len(data) > 1_500_000 or hashlib.sha256(data).hexdigest() != item['original_sha256']:
            raise ValueError('Original checksum mismatch for ' + acc)
        record, _ = parse_genbank(data, acc)
        base, version = acc.rsplit('.', 1)
        row = {'record_id': 'ncbi-' + acc, 'source': 'NCBI', 'source_url': item['source_url'],
               'accession': base, 'accession_version': version, 'archetype': item['archetype'],
               'notes': 'Reviewed for in-memory annotation preservation only.'}
        doc, audit, extracted = normalize_record(record, row, item['original_sha256'])
        result = audit_quality(record, doc, json_safe(audit), extracted,
                               max_residues=panel['criteria']['max_residues'],
                               source_sha256=item['original_sha256'])
        if result['status'] == 'failed':
            raise ValueError('E1 audit failed for ' + acc + ': ' + str(result['failed_checks']))
        results.append(result)
    report = {
        'run_kind': 'panel_v2_in_memory_E1_annotation_quality_helper',
        'panel_canonical_sha256': PANEL_HASH,
        'review_file': str(REVIEW.relative_to(ROOT)), 'review_sha256': hashlib.sha256(review).hexdigest(),
        'timestamp_utc': datetime.now(timezone.utc).isoformat(), 'python': platform.python_version(),
        'biopython': __import__('Bio').__version__, 'max_residues': panel['criteria']['max_residues'],
        'method': 'In-memory exact NCBI source, normalize_record then audit_quality; no raw source, qualifier, protein or sidecar files saved. Not the production CLI.',
        'records': [{key: r[key] for key in ('construct_id', 'status', 'quality_stratum',
                    'source_checksum_binding', 'counts', 'metrics', 'source_preservation_audit_status',
                    'failed_checks')} | {'source_preservation_counts': r['source_preservation_audit'].get('counts')}
                    for r in results],
        'summary': summarize_quality(results),
        'limitations': ['Source-annotation agreement is not independent biological truth.',
                        'Not an independent retrieval, decomposition benefit, model inference or fusion-boundary experiment.',
                        'Exact source bytes were processed in memory and not persisted.',
                        'This helper-level run does not reproduce the full production CLI output or its file-based admission path.'],
    }
    # Require an explicit fresh destination; never overwrite the original evidence.
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python scripts/panel_v2_e1_in_memory.py <new-output.json>')
    output = Path(sys.argv[1])
    with output.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print('E1 source-free audit written:', output)


if __name__ == '__main__':
    main()
