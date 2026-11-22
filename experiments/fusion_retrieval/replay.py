"""Verify packaged historical evidence; no models, sequences or downloads."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def main():
    manifest = json.loads((ROOT / 'MANIFEST.json').read_text())
    for entry in manifest['files']:
        actual = hashlib.sha256((ROOT / entry['file']).read_bytes()).hexdigest()
        if actual != entry['sha256']:
            raise ValueError('Changed source artifact: ' + entry['file'])
    rows = json.loads((ROOT / 'results/fusion-example.json').read_text())
    count = 0
    print('Saved example: whole versus windows, five TOTAL returned hits')
    for row in rows:
        if row['budget'] != 5 or row['aggregation'] == 'per_segment_diagnostic':
            continue
        returned = row['returned_accessions']
        if len(returned) != len(set(returned)) or len(returned) > row['budget']:
            raise ValueError('Invalid output budget')
        if row['threshold'] != 1.0 or row['n_targets'] != 1:
            raise ValueError('Example requires exact single-target condition')
        rank = returned.index(row['payload']) + 1 if row['payload'] in returned else None
        if rank != row['target_rank_in_output'] or bool(rank) != bool(row['hit']):
            raise ValueError('Saved rank disagrees with actual returned identities')
        print(f"{row['arm']:12} {row['aggregation']:12} target rank {rank}: {', '.join(returned)}")
        count += 1
    if not count:
        raise ValueError('No example rows checked')
    print(f"Verified {len(manifest['files'])} source hashes and {count} actual output lists.")
    print('Replay of saved evidence, not fresh inference or raw cosine recovery.')

if __name__ == '__main__':
    main()
