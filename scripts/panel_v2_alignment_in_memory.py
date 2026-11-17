"""Run the specifically reviewed E3 panel-v2 alignment-only sensitivity in memory."""
import hashlib
import json
import resource
import sys
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.annotations import normalize_record
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.dataset_curation import parse_genbank
from iri_constructs.pilot_alignment import alignment_scores
from iri_constructs.pilot_scoring import _hash, validate_panel
from iri_constructs.reference_tiers import score_reference_tier

FREEZE_SHA = '00c487200404e58f4a6aefd029b6db8df9d70aa3adac7a20201f242d80910d4d'
PANEL_SHA = 'f9111dd17ec028fdc170552aacc8af7275fccaf46d1e9712e3f7a681bdc0d70b'
E3_SHA = 'f13e363599c4173ea7973301a696600d7dcf15c742e3d572e8c2b82ea87f1a53'


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def read_bytes(url):
    with urlopen(Request(url, headers={'User-Agent': 'IRI-project-research/1.0'}), timeout=30) as response:
        data = response.read(1_500_001)
    require(len(data) <= 1_500_000, 'Provider response exceeds size limit')
    return data


def main():
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python scripts/panel_v2_alignment_in_memory.py <new-output.json>')
    output = Path(sys.argv[1])
    require(not output.exists(), 'Refusing to overwrite result')
    freeze_path = ROOT / 'docs/results/panel-v2-2026-11-13/E3.alignment.freeze.json'
    freeze_bytes = freeze_path.read_bytes()
    require(hashlib.sha256(freeze_bytes).hexdigest() == FREEZE_SHA, 'Execution freeze changed')
    freeze = json.loads(freeze_bytes)
    review = (ROOT / 'docs/coordination/PANEL_V2_ALIGNMENT_REVIEW_2026_11_13.md').read_bytes()
    require(hashlib.sha256(review).hexdigest() == freeze['review_sha256'], 'Review changed')
    require(b'CLEAR_FOR_E3_PANEL_V2_ALIGNMENT_DESCRIPTIVE_ONLY' in review, 'E3 alignment is not cleared')
    proposal = json.loads((ROOT / 'data/manifests/prospective-2026-11-12/panel.v2.json').read_text())
    proposal_hash = proposal.pop('proposal_sha256')
    require(proposal_hash == PANEL_SHA == _hash(proposal), 'Proposal changed')
    old = (ROOT / 'docs/results/cleared-2026-11-11/E3.json').read_bytes()
    require(hashlib.sha256(old).hexdigest() == E3_SHA, 'Historical E3 evidence changed')
    panel = freeze['panel']
    validate_panel(panel)
    require(panel['eligible_reference_ids'] == proposal['eligible_reference_ids'], 'Changed eligibility')
    original = json.loads(old)['panel']
    require(set(panel['reference_ids']) == set(original['reference_ids']) | {'uniprot-O58462'}, 'Changed reference set')
    for qid in original['query_ids']:
        require(set(panel['eligible_reference_ids'][qid]) == set(original['eligible_reference_ids'][qid]) | {'uniprot-O58462'}, 'Changed original exclusions')
    # Input identities and exact provider bytes are checked before any scoring.
    proteins = {}
    source_queries = {}
    for query in proposal['queries']:
        source_queries.setdefault(query['source_id'], []).append(query)
    for source_id, queries in source_queries.items():
        acc = source_id.removeprefix('ncbi-')
        url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + urlencode(
            {'db': 'nuccore', 'id': acc, 'rettype': 'gbwithparts', 'retmode': 'text'})
        data = read_bytes(url)
        source_sha = hashlib.sha256(data).hexdigest()
        require(all(source_sha == q['source_sha256'] for q in queries), 'Changed query source')
        record, _ = parse_genbank(data, acc)
        base, version = acc.rsplit('.', 1)
        row = {'record_id': source_id, 'source': 'NCBI', 'source_url': queries[0]['source_url'],
               'accession': base, 'accession_version': version, 'archetype': 'multi_cds', 'notes': 'Scoped AI-assisted review.'}
        _, _, extracted = normalize_record(record, row, source_sha)
        for q in queries:
            qid = q['query_id']
            unit = qid.split('/', 1)[1]
            require(unit in extracted, 'Missing reviewed query unit')
            protein = ProteinInput(qid, source_id, extracted[unit], source_sha, q['feature_ids'][0])
            require(protein.sequence_sha256 == q['protein_sequence_sha256'], 'Changed query protein')
            proteins[qid] = protein
    refs = json.loads((ROOT / 'data/manifests/week2/references/reference_proposal.json').read_text())
    ref_proposals = {r['reference_id']: (r['accession'], r['original_json_sha256'], r['protein_sequence_sha256']) for r in refs['references']}
    new = proposal['new_reference']
    ref_proposals[new['reference_id']] = (new['accession'], new['original_sha256'], new['sequence_sha256'])
    require(set(ref_proposals) == set(panel['reference_ids']), 'Reference definitions changed')
    for rid in panel['reference_ids']:
        acc, original_sha, sequence_sha = ref_proposals[rid]
        data = read_bytes('https://rest.uniprot.org/uniprotkb/' + acc + '.json')
        require(hashlib.sha256(data).hexdigest() == original_sha, 'Changed reference bytes: ' + rid)
        reference = json.loads(data)
        require(reference['primaryAccession'] == acc and reference['entryType'] == 'UniProtKB reviewed (Swiss-Prot)', 'Invalid reviewed reference')
        protein = ProteinInput(rid, rid, reference['sequence']['value'], original_sha)
        require(protein.sequence_sha256 == sequence_sha, 'Changed reference sequence: ' + rid)
        proteins[rid] = protein
    for iid, fp in panel['input_fingerprints'].items():
        require(iid in proteins, 'Missing frozen protein')
        actual = proteins[iid].provenance()
        require(all(actual[k] == fp[k] for k in ('construct_id', 'source_sha256', 'sequence_sha256', 'feature_id')), 'Changed frozen protein')
    require(len(proteins) <= 64 and sum(len(p.sequence) for p in proteins.values()) <= 65536, 'Input budget exceeded')
    labels = {iid: {'input_id': iid, 'family': fp['family'], 'independence_group': fp['independence_group'],
                    'annotation_status': 'known', 'review_status': 'approved'} for iid, fp in panel['input_fingerprints'].items()}
    resource.setrlimit(resource.RLIMIT_AS, (6 * 1024**3, 6 * 1024**3))
    started = time.monotonic()
    scores = alignment_scores(proteins, panel, max_residues=1022)
    require(scores['status'] == 'ok' and time.monotonic() - started < 1800, 'Alignment did not complete within bounds')
    raw = score_reference_tier(scores['raw_scores'], panel, labels)
    normalized = score_reference_tier(scores['normalized_scores'], panel, labels)
    require(raw['scored_queries'] == len(panel['query_ids']), 'Incomplete raw scoring')
    require(normalized['scored_queries'] == len(panel['query_ids']), 'Incomplete normalized scoring')
    report = {'stage': 'E3_panel_v2_alignment_only', 'freeze_sha256': FREEZE_SHA, 'panel_sha256': panel['panel_sha256'],
              'raw_score_matrix_sha256': _hash(scores['raw_scores']), 'normalization_score_matrix_sha256': _hash(scores['normalized_scores']),
              'raw_primary': raw, 'normalized_sensitivity': normalized, 'configuration': scores['configuration'],
              'resources': scores['resources'], 'pairwise_score_calls': scores['pairwise_score_calls'],
              'reference_self_score_calls': scores['reference_self_score_calls'],
              'limitations': ['Alignment only; no matched ESM-2 execution.', 'Four queries, three parents, two dependence clusters; descriptive only.',
                              'Reference composition sensitivity, not a new independent-query benchmark.', 'No raw source, protein or score matrix is persisted.']}
    with output.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print('Source-free E3 alignment result written:', output)


if __name__ == '__main__':
    main()
