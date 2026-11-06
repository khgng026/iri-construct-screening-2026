"""Frozen benign input handling. Never include sequence contents in errors/logs."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
from urllib.request import Request

CANONICAL_AA = frozenset('ACDEFGHIKLMNPQRSTVWY')


def source_record(document: dict) -> dict:
    """Normalize only required fields from official UniProt and EBI Proteins JSON."""
    if 'results' in document:
        if len(document['results']) != 1:
            raise ValueError('Expected one source record')
        document = document['results'][0]
    if 'primaryAccession' in document:
        return document
    if document.get('info', {}).get('type') != 'Swiss-Prot':
        raise ValueError('Source is not reviewed Swiss-Prot')
    return {
        'primaryAccession': document['accession'],
        'entryType': 'UniProtKB reviewed (Swiss-Prot)',
        'entryAudit': {'entryVersion': document['info']['version'],
                       'sequenceVersion': document['sequence']['version'],
                       'lastAnnotationUpdateDate': document['info']['modified'],
                       'lastSequenceUpdateDate': document['sequence']['modified']},
        'organism': {'taxonId': document['organism']['taxonomy'],
                     'scientificName': next(x['value'] for x in document['organism']['names'] if x['type'] == 'scientific')},
        'proteinDescription': document['protein'],
        'sequence': {'value': document['sequence']['sequence'], 'length': document['sequence']['length']},
        'features': document.get('features', [])}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def validate_protein(sequence: str, max_residues: int = 1022) -> None:
    if not isinstance(sequence, str) or not sequence:
        raise ValueError('Protein must be a nonempty string')
    if len(sequence) > max_residues:
        raise ValueError(f'Protein exceeds operational limit of {max_residues} residues; not truncated')
    if not set(sequence) <= CANONICAL_AA:
        raise ValueError('Protein contains unsupported characters; not modified')


def local_data_path(root: Path, relative: str) -> Path:
    """Only ignored data paths may supply sequences; block escapes and symlinks out."""
    candidate = (root / relative).resolve()
    allowed = [(root / 'data' / part).resolve() for part in ('raw', 'processed')]
    if not any(candidate.is_relative_to(base) and candidate != base for base in allowed):
        raise ValueError('Sequence path must remain inside local raw/processed data')
    return candidate


def read_single_fasta(path: Path, expected_id: str) -> str:
    lines = path.read_text().splitlines()
    if not lines or not lines[0].startswith('>'):
        raise ValueError('Missing FASTA header')
    if lines[0][1:] != expected_id or sum(x.startswith('>') for x in lines) != 1:
        raise ValueError('Expected exactly one FASTA record with matching ID')
    sequence = ''.join(lines[1:])
    validate_protein(sequence)
    return sequence


def load_catalogue(root: Path, manifest: dict, *, download: bool = False) -> dict[str, str]:
    """Use exact frozen bytes or fail. Downloads never update the manifest."""
    sequences = {}
    seen_hashes = set()
    if manifest.get('manifest_version') != 't05-benign-v1':
        raise ValueError('Unsupported catalogue manifest version')
    for record in manifest['records']:
        accession = record['accession']
        if accession in sequences:
            raise ValueError('Duplicate catalogue accession')
        if any(record.get(gate) != 'approved' for gate in
               ('review_status', 'licensing_status', 'safety_review')):
            raise ValueError('Catalogue contains an unapproved record')
        if record.get('annotation_status') != 'known':
            raise ValueError('Catalogue annotations are not reviewed')
        license_path = local_data_path(root, record['license_evidence_path'])
        if not license_path.exists() and download:
            with urlopen('https://rest.uniprot.org/help/license', timeout=30) as response:
                license_bytes = response.read()
            if sha256(license_bytes) != record['license_evidence_sha256']:
                raise ValueError('License evidence changed; review required')
            license_path.parent.mkdir(parents=True, exist_ok=True)
            license_path.write_bytes(license_bytes)
        if sha256(license_path.read_bytes()) != record['license_evidence_sha256']:
            raise ValueError('License evidence checksum mismatch')
        raw_path = local_data_path(root, record['raw_path'])
        if not raw_path.exists():
            if not download:
                raise ValueError('Frozen source data missing; run prepare --download')
            request = Request(record['download_url'], headers={'Accept': 'application/json'})
            with urlopen(request, timeout=30) as response:
                raw = response.read()
            if sha256(raw) != record['checksum_sha256']:
                raise ValueError('Downloaded source differs from frozen checksum; review required')
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            raw_path.write_bytes(raw)
        raw = raw_path.read_bytes()
        if sha256(raw) != record['checksum_sha256']:
            raise ValueError('Source file checksum mismatch')
        document = source_record(json.loads(raw))
        if (document['primaryAccession'] != accession or
                document['entryType'] != 'UniProtKB reviewed (Swiss-Prot)' or
                document['entryAudit']['sequenceVersion'] != record['sequence_version'] or
                document['entryAudit']['entryVersion'] != record['entry_version'] or
                document['organism']['taxonId'] != record['taxon_id']):
            raise ValueError('Source identity/version differs from review')
        sequence = document['sequence']['value']
        validate_protein(sequence)
        digest = sha256(sequence.encode('ascii'))
        if digest != record['sequence_sha256'] or len(sequence) != record['sequence_length_aa']:
            raise ValueError('Sequence checksum or length mismatch')
        if digest in seen_hashes:
            raise ValueError('Duplicate catalogue sequence')
        seen_hashes.add(digest)
        ref = record['sequence_ref']
        fasta_path = local_data_path(root, ref['path'])
        expected = f'>{accession}\n{sequence}\n'.encode('ascii')
        if fasta_path.exists() and fasta_path.read_bytes() != expected:
            raise ValueError('Existing processed FASTA differs; refusing overwrite')
        fasta_path.parent.mkdir(parents=True, exist_ok=True)
        fasta_path.write_bytes(expected)
        if read_single_fasta(fasta_path, ref['record_id']) != sequence:
            raise ValueError('Processed sequence mismatch')
        if ref['sha256'] != digest or ref['length'] != len(sequence) or ref['format'] != 'fasta':
            raise ValueError('Sequence reference metadata mismatch')
        sequences[accession] = sequence
    if not sequences:
        raise ValueError('Empty reference catalogue')
    return sequences


def build_query(components: list[str], target: str, sequences: dict[str, str]):
    if components.count(target) != 1:
        raise ValueError('Each query must contain exactly one expected target')
    query = ''.join(sequences[a] for a in components)
    validate_protein(query)
    start = sum(len(sequences[a]) for a in components[:components.index(target)])
    end = start + len(sequences[target])
    oracle = query[start:end]
    if oracle != sequences[target]:
        raise ValueError('Oracle interval does not recover the target')
    return query, oracle, (start, end)
