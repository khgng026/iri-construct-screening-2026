"""Acquire/verify a bounded reviewed UniProt reference proposal, never score it."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
from urllib.error import HTTPError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.reference_dataset import (build_source_queries, file_hash, freeze_proposal,
                                              summarize_reference)
from iri_constructs.runtime import InputError, local_output, write_json


def fetch(url):
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={'User-Agent': 'IRI-academic-source-evidence/1.0',
                                               'Accept': 'application/json'}), timeout=45) as response:
                data = response.read(5_000_001)
            if len(data) > 5_000_000:
                raise InputError('Reference source response exceeds size bound')
            return data
        except HTTPError as exc:
            if exc.code not in (429, 503) or attempt == 2:
                raise
            value = exc.headers.get('Retry-After', '10')
            # Never retry before an upstream server's requested delay.
            if not value.isdigit() or int(value) > 60:
                raise InputError('Reference provider requires a later retry; cached work is preserved') from exc
            time.sleep(max(1, int(value)))
    raise InputError('Reference retrieval did not complete')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue', type=Path, default=ROOT/'data/manifests/week2/references/catalogue.json')
    parser.add_argument('--queries', type=Path, help='Optional unscored query-proposals JSON with queries list')
    parser.add_argument('--query-specification', type=Path, help='Explicit source-protein family/evidence proposals')
    parser.add_argument('--unit-provenance', type=Path, help='Preparatory normalized unit provenance bound to source bytes')
    parser.add_argument('--source-manifest', type=Path, help='Exact source admission/evidence proposal bound to this panel')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--download', action='store_true', help='Permit anonymous UniProt source download; no inference')
    args = parser.parse_args(argv)
    if args.query_specification and (args.queries or not args.unit_provenance or not args.source_manifest):
        parser.error('--query-specification requires --unit-provenance and --source-manifest; do not also pass --queries')
    output = local_output(args.output)
    raw = output/'raw_references'
    raw.mkdir(parents=True, exist_ok=True)
    spec = json.loads(args.catalogue.read_text())
    refs, attempts = [], []
    for row in spec['references']:
        accession = row['accession']
        path = raw/(accession+'.json')
        record = {'accession': accession, 'status': 'blocked'}
        try:
            meta = path.with_suffix('.retrieval.json')
            if path.exists() != meta.exists():
                raise InputError('Incomplete cached reference pair; preserve the original and choose a new output directory')
            if path.exists() and meta.exists():
                retrieved = json.loads(meta.read_text())
                if file_hash(path) != retrieved['original_json_sha256']:
                    raise InputError('Cached reference bytes changed after retrieval')
            else:
                if not args.download:
                    raise InputError('Missing cached reference; explicit download flag required')
                data = fetch(f'https://rest.uniprot.org/uniprotkb/{accession}.json')
                # Validate response before storing metadata; raw source stays ignored.
                json.loads(data)
                path.write_bytes(data)
                retrieved = {'retrieved_at_utc': datetime.now(timezone.utc).isoformat(),
                             'original_json_sha256': hashlib.sha256(data).hexdigest()}
                write_json(meta, retrieved)
                time.sleep(0.25)
            summary = summarize_reference(json.loads(path.read_bytes()), row,
                raw_sha256=retrieved['original_json_sha256'],
                retrieved_at_utc=retrieved['retrieved_at_utc'])
            refs.append(summary)
            record.update(status='validated_reference_proposal', original_json_sha256=file_hash(path))
        except (OSError, ValueError, InputError) as exc:
            # No raw provider error bodies, sequence contents or credentials in reports.
            record['error_type'] = type(exc).__name__
            record['reason'] = str(exc) if isinstance(exc, InputError) else 'Source acquisition/validation failed'
        attempts.append(record)
        print(accession, record['status'], flush=True)
    query_bundle = None
    if args.query_specification:
        query_bundle = build_source_queries(json.loads(args.source_manifest.read_text()),
            json.loads(args.unit_provenance.read_text()), json.loads(args.query_specification.read_text()))
        write_json(output/'query_proposals.json', query_bundle)
        queries = query_bundle['queries']
    else:
        queries = json.loads(args.queries.read_text())['queries'] if args.queries else []
    proposal = freeze_proposal(refs, queries,
        source_manifest_sha256=file_hash(args.source_manifest) if args.source_manifest else None,
        catalogue_spec_sha256=file_hash(args.catalogue))
    proposal['source_rights_evidence'] = spec['source_rights_evidence']
    if query_bundle:
        proposal['parent_component_coverage'] = query_bundle['parent_component_coverage']
        proposal['query_specification_sha256'] = file_hash(args.query_specification)
        proposal['unit_provenance_sha256'] = file_hash(args.unit_provenance)
    # The proposal fingerprint includes the complete final metadata, not only roles.
    from iri_constructs.reference_dataset import canonical_hash
    proposal.pop('proposal_sha256')
    proposal['proposal_sha256'] = canonical_hash(proposal)
    write_json(output/'reference_proposal.json', proposal)
    write_json(output/'acquisition.json', {'attempts': attempts, 'scoring_performed': False})
    return 0 if len(refs) == len(spec['references']) else 1


if __name__ == '__main__':
    raise SystemExit(main())
