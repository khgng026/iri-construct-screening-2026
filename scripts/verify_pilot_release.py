"""Verify the exact one-source E1 release and its HOLD boundaries; no downloads."""
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.annotations import annotated_translation, feature_location
from iri_constructs.validate_manifest import validate_manifest


def read(path):
    return json.loads((ROOT / path).read_text(encoding='utf-8'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    from Bio import SeqIO
    base = 'data/manifests/pilot_execution_2026_11_06/'
    release = read(base + 'release.json')
    selected = {r['record_id'] for r in read('data/manifests/pilot_shortlist_v1_1.json')['records']}
    require(len(selected) == 10 and {r['record_id'] for r in release['records']} == selected,
            'Release must account for exactly the approved ten identities')
    approved = []
    for row in release['records']:
        detail = read(base + row['manifest'])
        require(detail['record_id'] == row['record_id'] and detail['status'] == row['status'],
                'Per-record status or identity differs from release')
        if row['status'] == 'HOLD':
            require(detail['checksum_sha256'] is None and not detail['downloaded_locally']
                    and not detail['source_admitted_for_annotation_audit']
                    and not detail['cloud_storage_approved'], 'HOLD record acquired/admitted by implication')
        else:
            require(row['status'] == 'ADMITTED_E1_ONLY' and row['record_id'] == 'ncbi-J01749.1',
                    'Unexpected source admitted')
            approved.append(detail)
    require(len(approved) == release['admitted_e1_source_count'] == 1
            and release['held_source_count'] == 9 and not release['retrieval_panel_ready'],
            'Source counts or retrieval readiness changed')
    detail = approved[0]
    source = ROOT / detail['source_path']
    raw = source.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    require(digest == detail['exact_file']['file_sha256']
            == '7660b0fe95a91d492b102192dc72924508a69fda762988aa0923b3d5cc2073f6', 'Original-byte checksum mismatch')
    require(len(raw) == detail['exact_file']['file_bytes'] and raw.rstrip().endswith(b'//'),
            'Incomplete source file')
    record = SeqIO.read(source, 'genbank')
    require(record.id == 'J01749.1' and len(record) == 4361 and len(record.features) == 50,
            'Exact accession, length or feature count changed')
    cds = [f for f in record.features if f.type == 'CDS']
    require(len(cds) == len(detail['cds_evidence']) == 3, 'CDS count changed')
    for f, evidence in zip(cds, detail['cds_evidence']):
        seq, state, _ = annotated_translation(record, f, feature_location(f, len(record))[2])
        require(state == 'known' and seq is not None, 'Source translation no longer validates')
        require(hashlib.sha256(seq.encode('ascii')).hexdigest() == evidence['translation_sha256'],
                'Annotated protein fingerprint changed')
    manifest = ROOT / (base + 'admitted_e1_sources.csv')
    require(not validate_manifest(manifest), 'T02 manifest metadata invalid')
    with manifest.open(encoding='utf-8', newline='') as handle:
        rows = list(csv.DictReader(handle))
    require(len(rows) == 1 and rows[0]['checksum_sha256'] == digest, 'Unexpected admission row')
    inputs = read('config/pilot/inputs.verified_e1_2026_11_06.json')
    require(len(inputs['records']) == 1 and not inputs['families'] and not inputs['reference_ids'],
            'Unreviewed retrieval panel introduced')
    report = read('docs/runs/pilot-2026-11-06/E1/report.json')
    quality = report['quality_report']
    counts = quality['overall']['counts']
    require(report['status'] == 'ok' and not report['inference_executed']
            and counts['observable_features'] == counts['retained_observable_features'] == 48
            and counts['preserved_observable_qualifiers'] == 48
            and counts['source_aligned_exported_cds'] == 3
            and counts['erroneous_invention_count'] == 0, 'E1 evidence changed')
    for experiment, status in [('E2', 'blocked_insufficient_parent_panel'),
                               ('E3', 'blocked_insufficient_panel'), ('E4', 'blocked_no_measurements')]:
        result = read(f'docs/runs/pilot-2026-11-06/{experiment}/report.json')
        require(result['status'] == status and not result.get('inference_executed')
                and not result.get('scientific_scores_produced'), 'Blocked experiment promoted to science')
    index = read('docs/runs/pilot-2026-11-06/artifact_checksums.json')
    for path, expected in index['sha256'].items():
        require(hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected,
                'Published run evidence digest mismatch: ' + path)
    print('Verified: 1 exact E1 source, 9 HOLD records, preserved run evidence, zero retrieval scores.')


if __name__ == '__main__':
    main()
