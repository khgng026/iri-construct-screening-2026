"""Acquire the pre-declared seeded sample of public NCBI vector records.

The candidate pool is one fixed Entrez query; the sample is a seeded draw from
its accession-sorted ID list, made before any record content is inspected.
Original records go to an ignored private cache. A source-free acquisition
manifest (accessions, hashes, retrieval facts) is written for Git.

    python experiments/component_study_v1/fetch_ncbi_sample.py

This is the only step that uses the network. It is not part of the toolkit.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import time
import urllib.parse
import urllib.request

ROOT = Path(__file__).resolve().parents[2]
BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
QUERY = ('("cloning vector"[Title] OR "expression vector"[Title] OR "shuttle vector"[Title]) '
         'AND "complete sequence"[Title] AND "vectors"[Organism] AND 2000:15000[SLEN] NOT gbdiv_pat[PROP]')
SEED, SAMPLE_SIZE, BATCH = 20261117, 200, 20
TERMS = ('NCBI states that it places no restrictions on the use or distribution of GenBank data, while noting '
         'that submitters may hold rights and NCBI cannot transfer unrestricted rights '
         '(https://www.ncbi.nlm.nih.gov/home/about/policies/). Used here as a private research cache with '
         'attribution; no sequence or qualifier text is committed or published.')


def get(endpoint: str, **params) -> bytes:
    params.setdefault('tool', 'iri-constructs-research')
    url = BASE + endpoint + '?' + urllib.parse.urlencode(params)
    for attempt in range(5):          # bounded retries; stay well under 3 requests/second
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                data = response.read()
            time.sleep(0.5)
            return data
        except Exception:
            time.sleep(3 * (attempt + 1))
    raise SystemExit('NCBI request failed after bounded retries: ' + endpoint)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--cache', type=Path, default=ROOT / 'data' / 'raw' / 'component_study_v1' / 'ncbi')
    parser.add_argument('--manifest', type=Path, default=ROOT / 'docs' / 'results' / 'component-study-v1' / 'acquisition.json')
    args = parser.parse_args(argv)
    args.cache.mkdir(parents=True, exist_ok=True)
    ids, count, start = [], None, 0
    while count is None or start < count:
        page = json.loads(get('esearch.fcgi', db='nuccore', term=QUERY, retmax=5000, retstart=start,
                              retmode='json', sort='accession'))['esearchresult']
        count = int(page['count'])
        ids.extend(page['idlist'])
        start += 5000
    ids = sorted(set(ids), key=int)
    sample = sorted(random.Random(SEED).sample(ids, SAMPLE_SIZE), key=int)
    files = []
    for number, offset in enumerate(range(0, len(sample), BATCH), 1):
        path = args.cache / f'batch_{number:02d}.gb'
        if not path.exists():
            data = get('efetch.fcgi', db='nuccore', id=','.join(sample[offset:offset + BATCH]),
                       rettype='gbwithparts', retmode='text')
            if not data.lstrip().startswith(b'LOCUS'):
                raise SystemExit('Unexpected EFetch response; nothing cached for this batch')
            path.write_bytes(data)
        data = path.read_bytes()
        versions = [line.split()[1] for line in data.decode('utf-8', 'replace').splitlines() if line.startswith('VERSION')]
        files.append({'file': path.name, 'sha256': hashlib.sha256(data).hexdigest(), 'size_bytes': len(data),
                      'requested_gi': sample[offset:offset + BATCH], 'accession_versions': versions})
    retrieved = datetime.now(timezone.utc).date().isoformat()
    manifest = {
        'schema_version': 'component-study-acquisition-1', 'source': 'NCBI GenBank (nuccore) via E-utilities EFetch',
        'query': QUERY, 'query_result_count': len(ids), 'seed': SEED, 'sample_size': SAMPLE_SIZE,
        'sampling_rule': 'random.Random(seed).sample(numerically sorted unique GI list, sample_size)',
        'retrieved_on': retrieved, 'rettype': 'gbwithparts', 'terms_note': TERMS,
        'storage': 'Original records are in an ignored private cache; this manifest is source-free.',
        'files': files, 'records': sum(len(f['accession_versions']) for f in files)}
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    provenance = {f['file']: {'source': 'NCBI GenBank', 'source_url': BASE + 'efetch.fcgi?db=nuccore&rettype=gbwithparts',
                              'retrieved_on': retrieved, 'terms_note': TERMS} for f in files}
    (args.cache / 'provenance.json').write_text(json.dumps(provenance, indent=2, sort_keys=True) + '\n')
    print(f"pool {len(ids)}; sampled {len(sample)}; cached {manifest['records']} records in {len(files)} files")


if __name__ == '__main__':
    main()
