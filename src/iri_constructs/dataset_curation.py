"""Bounded public-source acquisition and source-free Week 2 review proposals.

This is preparatory curation, not T02 scientific admission. Original files remain
in ignored storage. Successful acquisition never changes an approval status.
"""
from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
import hashlib
from io import StringIO
import json
import os
from pathlib import Path
import re
import time
from urllib.error import HTTPError
from urllib.parse import urlencode, urlparse
from urllib.request import Request, urlopen

from .annotations import annotated_translation, feature_location
from .runtime import InputError, ROOT, digest_file, local_output, write_json

MAX_BYTES = 1_500_000
LICENSE_UUID = 'd6c69ca7-8be4-4bc0-b4a8-d3ae1d428aa6'
NCBI_POLICY = 'https://www.ncbi.nlm.nih.gov/home/about/policies/'
IGEM_LICENSE = 'https://api.registry.igem.org/v1/licenses/' + LICENSE_UUID


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def retry_after_seconds(value: str, now=None) -> float | None:
    """Honor both RFC HTTP-date and nonnegative delay forms, never shorten them."""
    if value.isdigit():
        return float(value)
    try:
        date = parsedate_to_datetime(value)
        return max(0.0, (date - (now or datetime.now(timezone.utc))).total_seconds())
    except (TypeError, ValueError, OverflowError):
        return None


def safe_error(exc: Exception) -> str:
    if isinstance(exc, HTTPError):
        return f'HTTP_{exc.code}'
    if isinstance(exc, InputError):
        return str(exc)
    return type(exc).__name__


class SourceClient:
    """Cache original bytes with digests; cache corruption is never silently repaired.

    Long Retry-After responses produce a resumable HOLD instead of retrying early.
    This avoids blocking a worker for minutes while respecting provider policy.
    """
    def __init__(self, folder: Path, offline=False, opener=None, sleeper=None,
                 igem_interval=13.0, max_retries=2):
        self.folder = Path(folder)
        self.folder.mkdir(parents=True, exist_ok=True)
        self.offline = offline
        self.opener = opener or urlopen
        self.sleeper = sleeper or time.sleep
        self.igem_interval = max(0.0, igem_interval)
        self.max_retries = max_retries
        self._last = {}

    def get(self, url: str, filename: str) -> tuple[bytes, dict]:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+', filename):
            raise InputError('Invalid cache filename')
        if urlparse(url).scheme != 'https':
            raise InputError('Public source retrieval requires HTTPS')
        path = self.folder / filename
        index = path.with_name(path.name + '.retrieval.json')
        if path.exists() or index.exists():
            if not (path.is_file() and index.is_file()):
                raise InputError('Incomplete cache pair; preserve evidence and choose a new output directory')
            try:
                meta = json.loads(index.read_text(encoding='utf-8'))
            except (ValueError, OSError) as exc:
                raise InputError('Invalid cache metadata') from exc
            data = path.read_bytes()
            if meta.get('url') != url or meta.get('sha256') != sha256(data) or meta.get('bytes') != len(data):
                raise InputError('Cached source URL, byte count or digest mismatch')
            return data, {**meta, 'cache_hit': True}
        if self.offline:
            raise InputError('Exact original source not present in verified local cache')
        host = urlparse(url).hostname
        delay = self.igem_interval if host == 'api.registry.igem.org' else 0.4
        for attempt in range(self.max_retries + 1):
            elapsed = time.monotonic() - self._last.get(host, -float('inf'))
            if elapsed < delay:
                self.sleeper(delay - elapsed)
            self._last[host] = time.monotonic()
            try:
                request = Request(url, headers={
                    'User-Agent': 'IRI-Constructs/0.3 (public benign source review; no authentication)',
                    'Accept': 'application/json, text/plain, text/html;q=0.5',
                })
                with self.opener(request, timeout=35) as response:
                    data = response.read(MAX_BYTES + 1)
                    final_url = response.geturl()
                    content_type = response.headers.get('Content-Type')
                if len(data) > MAX_BYTES:
                    raise InputError('Source response exceeds bounded size')
                meta = {'url': url, 'final_url': final_url, 'sha256': sha256(data),
                        'bytes': len(data), 'retrieved_at_utc': utc_now(),
                        'content_type': content_type, 'cache_hit': False}
                path.write_bytes(data)
                write_json(index, meta)
                return data, meta
            except HTTPError as exc:
                if exc.code not in (429, 503) or attempt == self.max_retries:
                    raise
                retry = retry_after_seconds(exc.headers.get('Retry-After', ''))
                if retry is not None and retry > 60:
                    raise InputError(f'Provider Retry-After requires {int(retry + 0.999)} seconds; retry later using cache') from None
                self.sleeper(retry if retry is not None else min(30, (attempt + 1) * 10))
        raise InputError('No original source response obtained')


def parse_genbank(data: bytes, accession: str) -> tuple[object, dict]:
    """Read one complete annotated original, require exact version and DNA alphabet."""
    from .bench_inputs import validate_sequence
    try:
        text = data.decode('utf-8')
        if not text.lstrip().startswith('LOCUS') or not text.rstrip().endswith('//'):
            raise InputError('Incomplete or non-GenBank original source')
        try:
            from Bio import SeqIO
        except ImportError as exc:
            raise InputError('Biopython is required to inspect a complete GenBank original') from exc
        record = SeqIO.read(StringIO(text), 'genbank')
        if record.id != accession:
            raise InputError('Original GenBank accession.version mismatch')
        validate_sequence(str(record.seq).upper(), 'dna')
        if not record.features:
            raise InputError('GenBank feature table absent')
    except InputError:
        raise
    except Exception as exc:
        raise InputError('GenBank original source failed parsing or sequence integrity') from exc
    cds = []
    for ordinal, feature in enumerate(record.features, 1):
        if feature.type != 'CDS':
            continue
        _, segments, issue = feature_location(feature, len(record))
        seq, state, reason = annotated_translation(record, feature, issue)
        cds.append({'source_ordinal_1based': ordinal,
                    'protein_accession_version': feature.qualifiers.get('protein_id', [None])[0],
                    'product': feature.qualifiers.get('product', [None])[0],
                    'gene': feature.qualifiers.get('gene', [None])[0],
                    'location': str(feature.location), 'segments': segments,
                    'translation_state': state, 'translation_reason': reason,
                    'translation_length_aa': len(seq) if seq else None,
                    'translation_sha256': sha256(seq.encode()) if seq else None})
    notice_lines = [line.strip() for line in text.splitlines()
                    if re.search(r'copyright|all rights reserved|licen[sc]e|redistribut|patent|permission', line, re.I)]
    references = [{'authors': r.authors, 'title': r.title, 'journal': r.journal,
                   'pubmed_id': r.pubmed_id} for r in record.annotations.get('references', [])]
    return record, {'record_title': record.description, 'record_locus_date': record.annotations.get('date'),
                    'topology': record.annotations.get('topology', 'unknown'),
                    'sequence_length_nt': len(record), 'sequence_sha256': sha256(str(record.seq).upper().encode()),
                    'feature_count': len(record.features),
                    'feature_types': dict(sorted(Counter(f.type for f in record.features).items())),
                    'cds_count': len(cds), 'validated_source_translation_count': sum(c['translation_state'] == 'known' for c in cds),
                    'cds_evidence': cds, 'source_rights_notice_lines': notice_lines,
                    'attribution_references': references}


def parse_igem(data: bytes, annotations: bytes, expected_id: str, license_data: dict) -> tuple[dict, dict]:
    """Inspect original API data; never manufacture CDS translations from DNA."""
    from .bench_inputs import validate_sequence
    try:
        part, features = json.loads(data), json.loads(annotations)
        if part.get('name') != expected_id or part.get('status') != 'published':
            raise InputError('iGEM identity or published status mismatch')
        sequence = part.get('sequence')
        if not isinstance(sequence, str):
            raise InputError('iGEM deposited sequence absent')
        validate_sequence(sequence.upper(), 'dna')
        if not isinstance(features, dict) or not isinstance(features.get('data'), list):
            raise InputError('iGEM annotation response shape unsupported')
        rows = features['data']
        if (type(features.get('total')) is not int or features['total'] != len(rows)
                or type(features.get('page')) is not int or features['page'] != 1):
            raise InputError('iGEM annotation pagination incomplete')
        if part.get('licenseUUID') != license_data.get('uuid') or license_data.get('spdxID') != 'cc-by-4.0':
            raise InputError('iGEM exact-record license is not established by resolved license metadata')
        for feature in rows:
            locations = feature.get('locations')
            if not isinstance(locations, list) or not locations:
                raise InputError('iGEM source feature has no documented locations')
            for loc in locations:
                if (not isinstance(loc.get('start'), int) or isinstance(loc['start'], bool)
                    or not isinstance(loc.get('end'), int) or isinstance(loc['end'], bool)
                    or not 0 <= loc['start'] <= loc['end'] <= len(sequence)):
                    raise InputError('iGEM source feature coordinates invalid or unsupported')
    except InputError:
        raise
    except Exception as exc:
        raise InputError('iGEM original source failed parsing or integrity') from exc
    cds = [f for f in rows if f.get('role', {}).get('accession') == 'SO:0000316']
    return part, {'record_title': part.get('title'), 'description_evidence': part.get('description'),
                  'source_type': 'complete_plasmid' if part.get('topology', {}).get('label') == 'Circular' else 'cassette_or_part',
                  'topology': part.get('topology', {}).get('label'), 'sequence_length_nt': len(sequence),
                  'sequence_sha256': sha256(sequence.upper().encode()),
                  'feature_count': len(rows), 'feature_types': dict(sorted(Counter(f.get('role', {}).get('label', 'unknown') for f in rows).items())),
                  'source_cds_feature_count': len(cds), 'validated_source_translation_count': 0,
                  'cds_evidence': [{'source_feature_id': f.get('uuid'), 'label': f.get('label'),
                                    'locations_as_deposited': f.get('locations'), 'strand': f.get('strand'),
                                    'translation_state': 'unknown',
                                    'translation_reason': 'No source-provided translation qualifier in this API feature; no protein synthesized'} for f in cds],
                  'attribution_source_as_deposited': part.get('source'),
                  'license_uuid': part.get('licenseUUID'),
                  'source_created_at': part.get('audit', {}).get('created'),
                  'source_updated_at': part.get('audit', {}).get('updated'),
                  'annotation_coordinate_scope': 'Original raw intervals retained as deposited; coordinate indexing convention unverified, no normalized coordinate claim until documented',
                  'source_annotations_caveat': 'CDS role and label are deposited annotations, not validated translation, family/function evidence or reviewed fusion boundaries'}


def lineage_groups(source_id: str) -> list[str]:
    """Conservative prospective links; exact component copies are also computed."""
    mappings = {
        'pGL3-vector-lineage': ['ncbi-U47295.2', 'ncbi-U47296.2'],
        'pGL4-vector-lineage': ['ncbi-AY864931.1', 'ncbi-AY738222.1'],
        'pCellFree-vector-lineage': ['ncbi-KJ541668.1', 'ncbi-KJ541670.1'],
        'GfasPurple-construct-cassette': ['BBa_25NPQIRU', 'BBa_25UIP3DW'],
        'TyrBm-construct-cassette': ['BBa_25ULGTE2', 'BBa_25ECYPMT'],
        'YF1-FixJ-reporter-lineage': ['BBa_252QAUN3', 'BBa_25B5R4ZL'],
    }
    return [name for name, ids in mappings.items() if source_id in ids]


def source_scope_hold(accession: str, evidence: dict) -> tuple[str | None, str | None]:
    """Conservative source eligibility findings, never an inferred hazard classifier."""
    if any(re.search(r'\bccdb\b|\btoxin\b', (c.get('product') or '') + ' ' + (c.get('gene') or ''), re.I)
           for c in evidence.get('cds_evidence', [])):
        return ('T02_source_annotated_toxin_component_excluded',
                'Source-annotated toxin-based selection component conflicts with current T02 toxin-exclusion rule; do not include solely because the vector title is ordinary cloning.')
    context_holds = {
        'KJ170897.1': 'Original cited study concerns virulence in uropathogenic Proteus mirabilis; reporter-vector identity does not itself establish hazard, but this source association needs specialist/scope review under T02.',
        'PP098727.1': 'Original cited study concerns Type III effectors in bacterial infection; fluorescent reporter annotation does not resolve source-context eligibility under T02 and requires specialist/scope review.',
    }
    if accession in context_holds:
        return 'HOLD_context_review_source_publication_association', context_holds[accession]
    return None, None


def curate(output: Path, *, offline=False, source='all', only_ids=(), igem_interval=13.0) -> dict:
    """Acquire bounded preselected pool and write source-free evidence + local QA inputs."""
    output = local_output(output)
    output.mkdir(parents=True, exist_ok=True)
    raw = output / 'raw_sources'
    client = SourceClient(raw, offline=offline, igem_interval=igem_interval)
    proposals, inputs, rights = [], [], {}
    for name, url in [('ncbi_policy', NCBI_POLICY), ('igem_license', IGEM_LICENSE)]:
        try:
            data, meta = client.get(url, name + '.original')
            rights[name] = {**meta, 'path': str((raw / (name + '.original')).relative_to(ROOT))}
            if name == 'igem_license':
                rights[name]['metadata'] = json.loads(data)
        except Exception as exc:
            rights[name] = {'url': url, 'status': 'HOLD', 'reason': safe_error(exc)}
    catalog = json.loads((ROOT / 'data/manifests/open_archive_candidates_2026_11_06.json').read_text())['candidates']
    catalog = [{**r, 'source_id': 'ncbi-' + r['accession_version'], 'source': 'NCBI'} for r in catalog]
    catalog.append({'accession_version': 'J01749.1', 'source_id': 'ncbi-J01749.1', 'source': 'NCBI',
                    'proposed_role': 'backbone_heavy', 'source_url': 'https://www.ncbi.nlm.nih.gov/nuccore/J01749.1',
                    'caveat': 'Historical E1 source independently fetched; this does not create a new experimental observation'})
    igem = json.loads((ROOT / 'data/manifests/igem_additional_ten_2026_11_10.json').read_text())['candidates']
    catalog += [{**r, 'source_id': r['record_id'], 'source': 'iGEM'} for r in igem]
    for candidate in catalog:
        sid, provider = candidate['source_id'], candidate['source']
        if (source != 'all' and provider.lower() != source) or (only_ids and sid not in only_ids):
            continue
        row = {'source_id': sid, 'source': provider, 'source_url': candidate['source_url'],
               'candidate_selected': True, 'source_admitted_for_scoring': False,
               'cloud_copy_permitted_by_scientific_review': False, 'family_panel_approved': False,
               'run_authorized': False, 'reviewer_role': 'Codex evidence preparer; not delegated clearance',
               'related_lineage_proposals': lineage_groups(sid), 'selection_caveat': candidate.get('caveat')}
        try:
            if provider == 'NCBI':
                acc = candidate['accession_version']
                if not re.fullmatch(r'[A-Z]{1,4}\d{5,8}\.\d+', acc):
                    raise InputError('Invalid versioned NCBI accession')
                url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode(
                    {'db': 'nuccore', 'id': acc, 'rettype': 'gbwithparts', 'retmode': 'text'})
                name = 'NCBI_' + acc + '.gb'
                data, meta = client.get(url, name)
                record, evidence = parse_genbank(data, acc)
                row.update(evidence, accession_version=acc, source_type='complete_plasmid',
                           source_retrieval=meta, original_path=str((raw / name).relative_to(ROOT)))
                cds_count = evidence['cds_count']
                archetype = 'multi_cds' if cds_count > 1 else 'single_cds'
                if acc in ('J01749.1', 'JX560321.2'):
                    archetype = 'backbone_heavy'
                if acc in ('AF298787.1', 'KJ541668.1', 'KJ541670.1'):
                    row['architecture_caveat'] = 'Fusion-ready vector is not an independently validated completed fusion; no reviewed fusion boundary asserted'
                row['rights_evidence'] = {'policy_key': 'ncbi_policy', 'local_analysis_proposal': 'permitted_under_NCBI_molecular_data_usage_policy',
                    'private_copy_proposal': 'permitted_under_NCBI_molecular_data_usage_policy_subject_to_submitter_IP_and_exact_record_notices',
                    'limitations': 'NCBI imposes no distribution restrictions on molecular data and does not accept submitter-requested restrictions; submitter intellectual-property caveat remains. Not CC0, patent clearance or permission for public project release.',
                    'exact_source_conflicting_notice_found': bool(evidence['source_rights_notice_lines'])}
                if 'sha256' not in rights.get('ncbi_policy', {}):
                    raise InputError('NCBI policy evidence unavailable')
                hold_reason, hold_basis = source_scope_hold(acc, evidence)
                if hold_reason:
                    row['source_exclusion_basis'] = hold_basis
                    raise InputError(hold_reason)
                input_row = {'source_id': sid, 'path': 'raw_sources/' + name, 'format': 'genbank',
                    'sequence_id': acc, 'checksum_sha256': sha256(data), 'molecule_type': 'dna',
                    'source': 'NCBI', 'accession': acc.rsplit('.', 1)[0], 'accession_version': acc.rsplit('.', 1)[1]}
                row['benign_scope_proposal'] = 'Ordinary cloning/expression/reporter vector in preselected scope; source CDS products and original references retained for case-specific review. Antibiotic markers and viral-derived regulatory elements, where present, are not omitted.'
            else:
                sid = candidate['record_id']
                if not re.fullmatch(r'BBa_[A-Za-z0-9]+', sid):
                    raise InputError('Invalid iGEM record identifier')
                data, meta = client.get(candidate['source_api'], sid + '.json')
                annotations, annmeta = client.get(candidate['annotation_api'], sid + '.features.json')
                part, evidence = parse_igem(data, annotations, sid, rights.get('igem_license', {}).get('metadata', {}))
                row.update(evidence, source_retrieval=meta, annotation_retrieval=annmeta,
                           original_path=str((raw / (sid + '.json')).relative_to(ROOT)),
                           annotation_original_path=str((raw / (sid + '.features.json')).relative_to(ROOT)))
                count = evidence['source_cds_feature_count']
                archetype = 'multi_cds' if count > 1 else 'single_cds'
                if candidate['proposed_archetype'] == 'backbone_heavy' and evidence['source_type'] == 'complete_plasmid':
                    archetype = 'backbone_heavy'
                if candidate['proposed_archetype'] == 'tagged_fusion':
                    row['architecture_caveat'] = 'Source description documents fusion/tag intent; precise full-fusion CDS translation and amino-acid boundaries unreviewed. Not counted as validated tagged_fusion coverage.'
                row['rights_evidence'] = {'policy_key': 'igem_license', 'license': 'CC-BY-4.0',
                    'local_analysis_proposal': 'permitted_with_attribution', 'private_copy_proposal': 'permitted_with_attribution',
                    'attribution': {'provider': 'iGEM Registry', 'part_id': sid, 'url': candidate['source_url'],
                                    'title': part.get('title'), 'source_credit_as_deposited': part.get('source'),
                                    'license_url': rights['igem_license']['metadata']['url']},
                    'limitations': 'Preserve creator/source credits and modification notice on derived artifacts; API status is not experimental validation. No public release in this task.'}
                row['benign_scope_proposal'] = 'Source-deposited ordinary reporter, pigment or coating description retained verbatim for review; no experimental safety guarantee. Case-specific review is required before scoring.'
                if sid == 'BBa_25K8XDBD':
                    row['benign_scope_proposal'] += ' Engineered BslA-dCBM coating fusion is held from any scored study pending specialist/scope confirmation; no hazardous function inferred.'
                input_row = {'source_id': sid, 'path': 'raw_sources/' + sid + '.json', 'format': 'igem-json',
                    'sequence_id': sid, 'checksum_sha256': sha256(data), 'molecule_type': 'dna',
                    'annotations_path': 'raw_sources/' + sid + '.features.json',
                    'annotations_checksum_sha256': sha256(annotations), 'source': 'iGEM',
                    'accession': sid, 'accession_version': part.get('audit', {}).get('updated', 'unversioned')}
            row.update(status='READY_FOR_REVIEW', eligible_for_preparatory_QA=True, proposed_archetype=archetype,
                       scoring_status='HOLD_pending_actual_delegated_source_and_stage_clearance')
            input_row.update(source_url=row['source_url'], archetype=archetype, source_type=row['source_type'],
                             rights_evidence=row['rights_evidence'], benign_evidence=row['benign_scope_proposal'],
                             related_lineage_proposals=row['related_lineage_proposals'],
                             source_admitted_for_scoring=False)
            inputs.append(input_row)
        except Exception as exc:
            row.update(status='HOLD', eligible_for_preparatory_QA=False, failure_class=safe_error(exc))
        proposals.append(row)
        print(sid, row['status'], row.get('feature_count', 'unknown'), 'features', flush=True)
        # Persist progress in case provider cooldowns interrupt a larger run.
        write_json(output / 'source_review.partial.json', {'records': proposals, 'rights_evidence': rights})
    copies = {}
    for row in proposals:
        if row.get('sequence_sha256'):
            copies.setdefault(row['sequence_sha256'], []).append(row['source_id'])
    component_copies = {}
    for row in proposals:
        for cds in row.get('cds_evidence', []):
            if cds.get('translation_sha256'):
                component_copies.setdefault(cds['translation_sha256'], []).append({
                    'source_id': row['source_id'], 'protein_accession_version': cds.get('protein_accession_version'),
                    'source_ordinal_1based': cds.get('source_ordinal_1based')})
    report = {'schema_version': 'week2-source-review-proposal-1', 'generated_at_utc': utc_now(),
        'scope': 'Existing 11 GenBank + ten iGEM selection and independently fetched historical pBR322; no automatic scored-panel expansion',
        'evidence_kind': 'preparatory_source_integrity_and_annotation_review',
        'scientific_clearance_issued': False, 'original_source_bytes_committed': False,
        'record_count': len(proposals), 'preparatory_QA_eligible_count': len(inputs),
        'new_scoring_admission_count': 0, 'counts_by_proposed_archetype': dict(sorted(Counter(r['archetype'] for r in inputs).items())),
        'rights_evidence': rights, 'records': proposals,
        'identical_parent_sequence_groups': [ids for ids in copies.values() if len(ids) > 1],
        'identical_validated_component_groups': [ids for ids in component_copies.values() if len(ids) > 1],
        'limitations': ['Candidate count is not independent sample size; shared components and vector lineage must be excluded prospectively.',
            'Preparatory QA is allowed by current mission; scoring needs separate actual delegated scientific clearance.',
            'Source feature/translation agreement does not establish biological truth or protein function.',
            'iGEM annotations without provided translations remain unsupported for protein export; no ORFs or missing sequences are inferred.',
            '150-250 target deliberately deferred until rights, annotation quality, family evidence and independent reference coverage justify expansion.']}
    write_json(output / 'source_review.json', report)
    write_json(output / 'inputs.preparatory.json', {'schema_version': 'week2-source-inputs-1',
        'evidence_kind': 'preparatory_QA_only', 'scientific_scoring_authorized': False, 'records': inputs})
    return report


def publish_review(output: Path, destination: Path) -> None:
    """Copy only source-free review/provenance manifests; never original source bytes."""
    output = local_output(output)
    destination = Path(destination).resolve()
    if not destination.is_relative_to(ROOT / 'data/manifests/week2'):
        raise InputError('Published Week 2 review must stay under data/manifests/week2/')
    destination.mkdir(parents=True, exist_ok=True)
    review = json.loads((output / 'source_review.json').read_text())
    inputs = json.loads((output / 'inputs.preparatory.json').read_text())
    for row in inputs['records']:
        for key in ('path', 'annotations_path'):
            if key in row:
                row[key] = Path(os.path.relpath((output / row[key]).resolve(), destination)).as_posix()
    write_json(destination / 'source_review.json', review)
    write_json(destination / 'inputs.preparatory.json', inputs)
