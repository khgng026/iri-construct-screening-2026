#!/usr/bin/env python3
"""Acquire only the prespecified gene candidates. No model inference or approval."""
import concurrent.futures
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.t05_data import sha256


def main():
    selection = json.loads((ROOT/'config/t05_expanded_selection.json').read_text())
    genes = [g for names in selection['groups'].values() for g in names]
    raw = ROOT/'data/raw/t05-expanded'
    raw.mkdir(parents=True, exist_ok=True)
    def fetch(chunk):
        index, names = chunk
        path = raw/f'candidate-page-{index:02d}.json'
        url = 'https://www.ebi.ac.uk/proteins/api/proteins?' + urlencode({
            'offset':0, 'size':100, 'taxid':selection['taxon_id'],
            'gene':','.join(names), 'reviewed':'true'})
        if path.exists() and json.loads(path.with_suffix('.provenance.json').read_text())['url'] != url:
            raise ValueError('Cached acquisition request differs; retain separately before retry')
        if not path.exists():
            with urlopen(Request(url, headers={'Accept':'application/json'}), timeout=40) as response:
                data = response.read()
                headers = dict(response.headers)
            json.loads(data)
            path.write_bytes(data)
            path.with_suffix('.provenance.json').write_text(json.dumps({
                'url':url, 'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),
                'headers':headers, 'sha256':sha256(data)},indent=2)+'\n')
        docs = json.loads(path.read_bytes())
        if len(docs) >= 100:
            raise ValueError('Page may be truncated; pagination review required')
        return path, docs
    # Proteins API limits the gene parameter to 50 characters, not 50 names.
    chunks = [(i,genes[j:j+6]) for i,j in enumerate(range(0,len(genes),6))]
    records, pages = {}, []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for path, docs in pool.map(fetch, chunks):
            pages.append(json.loads(path.with_suffix('.provenance.json').read_text()))
            for d in docs:
                matches = sorted(set(genes) & {g.get('name',{}).get('value') for g in d.get('gene',[])})
                if not matches:
                    continue  # API also matches synonyms; primary-gene rule is prespecified.
                accession = d['accession']
                record_path = raw/f'{accession}.json'
                # This is a JSON extraction, not original HTTP bytes. Retain page/hash too.
                data = (json.dumps(d,indent=2)+'\n').encode()
                if record_path.exists() and record_path.read_bytes() != data:
                    raise ValueError('Conflicting source versions; refusing overwrite')
                record_path.write_bytes(data)
                functions = [c.get('text',[]) for c in d.get('comments',[]) if c.get('type')=='FUNCTION']
                records[accession] = {'accession':accession,'genes':matches,
                    'name':d['protein'].get('recommendedName',{}).get('fullName',{}).get('value'),
                    'length':d['sequence']['length'],'type':d['info']['type'],
                    'taxon':d['organism']['taxonomy'],'functions':functions,
                    'raw_path':str(record_path.relative_to(ROOT)),
                    'source_page_path':str(path.relative_to(ROOT)),
                    'source_page_sha256':sha256(path.read_bytes())}
    audit = {'records':list(sorted(records.values(),key=lambda x:x['accession'])),
             'unresolved_genes':sorted(set(genes)-{g for r in records.values() for g in r['genes']}),
             'source_pages':pages,'review_status':'pending'}
    (raw/'candidate-review.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(f'Acquired {len(records)} candidate accessions; unresolved genes: {audit["unresolved_genes"]}. Review pending.')


if __name__ == '__main__':
    main()
