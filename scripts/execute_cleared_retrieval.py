"""Execute the externally cleared descriptive E3/E2 using existing scorers."""
import argparse, hashlib, json, os, sys, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.runtime import digest_file, write_json, run_metadata, Measure, InputError
from iri_constructs.bench_inputs import load_proteins, ProteinInput
from iri_constructs.extension_config import base_scoring_config
from iri_constructs.pilot_scoring import freeze_panel, _hash
from iri_constructs.reference_tiers import score_reference_tier
from iri_constructs.pilot_alignment import alignment_scores
from iri_constructs.pilot_experiments import _CachedAdapter, _vectors
from iri_constructs.extension_experiments import _vector_scores
from iri_constructs.model_adapters import LocalAdapter, prepare_device, ESM2_REVISION
from iri_constructs.multi_cds import load_multi_cds_inputs, freeze_multi_cds_panel, produce_multi_cds_scores, score_multi_cds, compare_multi_cds_conditions

EXPECTED_REVIEW_FREEZE_SHA256 = 'c90d063f7f17c9f929af2fddb4e3fe9bc719c9a4cae5c36f66b238f968e47c2e'

def verify_original_freeze():
    import subprocess
    freeze = json.loads((ROOT / 'docs/week2/review-freeze.json').read_text())
    require(_hash({k:v for k,v in freeze.items() if k != 'freeze_sha256'}) == freeze['freeze_sha256'], 'Changed original freeze')
    require(freeze['freeze_sha256'] == EXPECTED_REVIEW_FREEZE_SHA256, 'Review freeze differs from approved scope')
    for name, record in freeze['artifacts'].items():
        path = ROOT / name
        # README is a maintained entry point; its reviewed version remains in Git.
        if path.is_file() and name != 'README.md':
            require(digest_file(path) == record['sha256'], 'Changed frozen artifact: ' + name)
        else:
            require(name == 'README.md' or name.startswith('docs/runs/week2-2026-11-11/'), 'Missing scientific input: ' + name)
            original = subprocess.check_output(['git','show','3fb426231347873269c655c25cc67bdd4ce1ebca:' + name], cwd=ROOT)
            require(hashlib.sha256(original).hexdigest() == record['sha256'], 'Archived evidence mismatch')

def require(condition, message):
    if not condition:
        raise InputError(message)

def validate_authorization(auth, stage):
    expected = {'E2': 'CLEAR_FOR_E2_SUPPORTED_SUBSET_DESCRIPTIVE', 'E3': 'CLEAR_FOR_E3_TIER1_DESCRIPTIVE_ALIGNMENT_ESM2'}
    require(stage in expected, 'Unapproved stage')
    require(auth.get('review_comment_id') == 6800553713 and auth.get('review_comment_url') == URL, 'Wrong reviewer decision')
    require(auth.get('decisions', {}).get(stage) == expected[stage], 'Stage is not cleared')
    require(auth.get('review_freeze_sha256') == 'c90d063f7f17c9f929af2fddb4e3fe9bc719c9a4cae5c36f66b238f968e47c2e', 'Wrong approved evidence scope')

def immutable_e3():
    import subprocess
    raw = subprocess.check_output(['git', 'show', '3fb426231347873269c655c25cc67bdd4ce1ebca:docs/results/cleared-2026-11-11/E3.json'], cwd=ROOT)
    require(hashlib.sha256(raw).hexdigest() == 'f13e363599c4173ea7973301a696600d7dcf15c742e3d572e8c2b82ea87f1a53', 'Immutable E3 evidence mismatch')
    return json.loads(raw)

def verify_score_cache(path):
    hashes = json.loads((path / 'cache-checksums.json').read_text())
    for name in ['report.json', 'execution-freeze.json', 'esm2-scores.json']:
        if digest_file(path / name) != hashes[name]:
            raise InputError('Cached score evidence changed')

    anchor = immutable_e3()
    require(digest_file(path / 'execution-freeze.json') == anchor['execution_freeze_sha256'], 'Cached execution/source configuration differs from immutable evidence')
    require(digest_file(path / 'report.json') == anchor['local_report_sha256'], 'Cached report differs from immutable evidence')
    scores = json.loads((path / 'esm2-scores.json').read_text())
    canonical = hashlib.sha256(json.dumps(scores, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    require(canonical == anchor['reproducibility']['canonical_score_matrix_sha256']['esm2-scores.json'], 'Score matrix differs from immutable experimental evidence')
    frozen = json.loads((path / 'execution-freeze.json').read_text())
    require(frozen['panel'] == anchor['panel'], 'Cached panel/configuration differs from immutable evidence')
    report = json.loads((path / 'report.json').read_text())
    require(report['methods']['esm2']['model'] == anchor['methods']['esm2']['model'], 'Cached model differs from immutable evidence')
    for iid, fp in frozen['input_fingerprints'].items():
        expected = anchor['panel']['input_fingerprints'].get(iid, {})
        require(all(fp.get(k) == expected.get(k) for k in ['construct_id','source_sha256','sequence_sha256','feature_id']), 'Cached source fingerprint mismatch')

def freeze_model(metadata, path):
    require(metadata['resolved_revision'] == ESM2_REVISION, 'Cleared execution validation failed')
    require(metadata['weight_sha256'] == '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189', 'Cleared execution validation failed')
    write_json(path, metadata)
URL = 'https://github.com/khgng026/iri-construct-screening-2026/pull/11#issuecomment-6800553713'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--stage', choices=['E3', 'E2'], required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--reuse-e3', type=Path)
    a = p.parse_args()
    if a.reuse_e3 and a.stage != 'E2':
        p.error('--reuse-e3 is an E2 replay only')
    a.output.mkdir(parents=True, exist_ok=False)
    proposal = ROOT / 'data/manifests/week2/references/reference_proposal.json'
    require(digest_file(proposal) == '627fb59b5e8b6a0b307447d2e22355fc75449bb02c21c3063468004948928349', 'Cleared execution validation failed')
    verify_original_freeze()
    refs = json.loads(proposal.read_text())
    auth = json.loads((ROOT / 'docs/oversight/EXPERIMENT_RUN_AUTHORIZATION.json').read_text())
    require(auth['review_comment_id'] == 6800553713, 'Cleared execution validation failed')
    validate_authorization(auth, a.stage)
    cfg = json.loads((ROOT / f'config/week2/{a.stage}.proposed.json').read_text())
    cfg.update(models=['alignment', 'esm2'] if a.stage == 'E3' else ['esm2'], scientific_execution_clearance=auth['decisions'][a.stage], clearance_ref=URL)
    queries = [q for q in refs['query_proposals'] if q['proposal_role'] == 'PROPOSED_SUPPORTED_PENDING_REVIEW']
    if a.stage == 'E2':
        queries = [q for q in queries if q['parent_id'] == 'ncbi-AF298787.1']
    inp = json.loads((ROOT / 'data/manifests/execution_2026_11_11/E1.new.inputs.json').read_text())
    inp['records'] = [r for r in inp['records'] if r['record_id'] in {q['parent_id'] for q in queries}]
    inp['families'] = [dict(input_id=q['query_id'], family=q['proposed_family'], annotation_status='known', review_status='approved', evidence=URL, independence_group=q['lineage_group']) for q in queries]
    inp['reference_ids'] = []
    if a.stage == 'E2':
        inp['parent_ids'] = ['ncbi-AF298787.1']
    ip = ROOT / f'data/manifests/execution_2026_11_11/{a.stage}.inputs.json'
    write_json(ip, inp)
    if a.stage == 'E3':
        _, proteins, prov, _ = load_proteins(ip)
        labels = {v['input_id']: v for v in inp['families']}
        parents = None
    else:
        _, proteins, labels, parents, prov = load_multi_cds_inputs(ip)
    for q in queries:
        require(proteins[q['query_id']].sequence_sha256 == q['protein_sequence_sha256'], 'Cleared execution validation failed')
        require(proteins[q['query_id']].source_sha256 == q['source_sha256'], 'Cleared execution validation failed')
    archetypes = {r['record_id']: r['archetype'] for r in prov['sources']}
    for q in queries:
        labels[q['query_id']]['archetype'] = archetypes[q['parent_id']]
    reference_ids = []
    for r in refs['references']:
        raw = ROOT / f"results/week2/references/raw_references/{r['accession']}.json"
        require(digest_file(raw) == r['original_json_sha256'], 'Cleared execution validation failed')
        doc = json.loads(raw.read_text())
        seq = doc['sequence']['value']
        require(hashlib.sha256(seq.encode()).hexdigest() == r['protein_sequence_sha256'], 'Cleared execution validation failed')
        iid = r['reference_id']
        reference_ids.append(iid)
        proteins[iid] = ProteinInput(iid, iid, seq, r['original_json_sha256'])
        labels[iid] = dict(input_id=iid, family=r['proposed_family'], annotation_status='known', review_status='approved', evidence=URL, independence_group=r['lineage_group'], archetype='unknown')
    require(len(reference_ids) == 17, 'Cleared execution validation failed')
    if a.stage == 'E3':
        scoring = base_scoring_config(json.loads((ROOT / 'config/week2/E3.proposed.json').read_text()))
        scoring['models'] = ['alignment', 'esm2']
        panel = freeze_panel(scoring, proteins, labels, [q['query_id'] for q in queries], reference_ids)
        allowed = {r['query_id']: r['eligible_reference_ids_proposed'] for r in refs['eligibility_proposals']}
        for q in queries:
            qid = q['query_id']
            require(set(allowed[qid]) <= set(panel['eligible_reference_ids'][qid]), 'Cleared execution validation failed')
            panel['eligible_reference_ids'][qid] = allowed[qid]
        panel.pop('panel_sha256')
        panel['panel_sha256'] = _hash(panel)
    else:
        panel = freeze_multi_cds_panel(cfg, proteins, labels, parents, reference_ids)
        expected = next((r['eligible_reference_ids_proposed'] for r in refs['eligibility_proposals'] if r['query_id'] == queries[0]['query_id']))
        require(panel['parents']['ncbi-AF298787.1']['eligible_reference_ids'] == expected, 'Cleared execution validation failed')
    ids = set(reference_ids) | {q['query_id'] for q in queries}
    require(len(ids) <= 64 and sum((len(proteins[i].sequence) for i in ids)) <= 65536, 'Cleared execution validation failed')
    write_json(a.output / 'execution-freeze.json', {'protocol': cfg, 'panel': panel, 'clearance': URL, 'proposal_sha256': digest_file(proposal), 'input_sha256': digest_file(ip), 'source_provenance': prov, 'input_fingerprints': {i: proteins[i].provenance() for i in sorted(ids)}})
    report = run_metadata(a.stage + '_cleared_descriptive', cfg)
    report.update(clearance=URL, evidence_kind='admitted_research', panel_sha256=panel['panel_sha256'], methods={}, inference_executed=False)
    if not a.reuse_e3:
        import torch
        torch.set_num_threads(2)
        torch.set_num_interop_threads(3)
    if a.stage == 'E3':
        start = time.process_time()
        with Measure() as measure:
            align = alignment_scores(proteins, panel, max_residues=1022)
        report['methods']['alignment_raw'] = {'result': score_reference_tier(align['raw_scores'], panel, labels), 'costs': dict(measure.data, cpu_seconds=time.process_time() - start, alignment_calls=align['pairwise_score_calls'] + align['reference_self_score_calls'])}
        write_json(a.output / 'alignment-scores.json', align)
        start = time.process_time()
        with Measure() as load:
            adapter = LocalAdapter('esm2', prepare_device('cpu', cfg['seed']), ESM2_REVISION, offline=True)
        freeze_model(adapter.metadata, a.output / 'model-freeze.json')
        cached = _CachedAdapter(adapter)
        with Measure() as measure:
            vectors, details = _vectors(cached, proteins, ids)
            scores = _vector_scores(vectors, panel)
        result = score_reference_tier(scores, panel, labels)
        require(result['scored_queries'] == len(queries), 'Incomplete cleared query support')
        report['methods']['esm2'] = {'result': result, 'costs': dict(measure.data, cold_load_seconds=load.data['wall_seconds'], cpu_seconds=time.process_time() - start, model_calls=sum((e['status'] == 'ok' and (not e['cache_hit']) for e in cached.events)), processed_residues=sum((e['length_aa'] for e in cached.events if e['status'] == 'ok' and (not e['cache_hit']))), cache_hits=sum((e['cache_hit'] for e in cached.events))), 'model': adapter.metadata, 'inputs': details}
        write_json(a.output / 'esm2-scores.json', scores)
        adapter.close()
    elif a.reuse_e3:
        verify_score_cache(a.reuse_e3)
        original = json.loads((a.reuse_e3 / 'report.json').read_text())
        frozen = json.loads((a.reuse_e3 / 'execution-freeze.json').read_text())
        require(original['clearance'] == URL and original['methods']['esm2']['model']['resolved_revision'] == ESM2_REVISION, 'Cleared execution validation failed')
        require(original['methods']['esm2']['model']['weight_sha256'] == '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189', 'Cleared execution validation failed')
        for iid in ids:
            require(frozen['input_fingerprints'][iid] == proteins[iid].provenance(), 'Changed cached source identity')
        scores = json.loads((a.reuse_e3 / 'esm2-scores.json').read_text())
        start = time.perf_counter()
        result = score_multi_cds({q['query_id']: scores[q['query_id']] for q in queries}, panel, labels)
        report['methods']['esm2'] = {'result': result, 'costs': {'model_calls': 0, 'wall_seconds': time.perf_counter() - start, 'shared_execution': str(a.reuse_e3)}, 'cache_note': 'Reuses exact same-model E3 unit/reference scores; no new cold experiment or independent observation.'}
    else:
        start = time.process_time()
        execution = produce_multi_cds_scores('esm2', proteins, panel, cfg, device='cpu', allow_download=False, freeze_model=lambda m: freeze_model(m, a.output / 'model-freeze.json'))
        require(execution['method_metadata']['weight_sha256'] == '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189', 'Cleared execution validation failed')
        result = score_multi_cds(execution['unit_scores'], panel, labels)
        result['paired_condition_comparison'] = compare_multi_cds_conditions(result, panel, cfg)
        execution['costs']['cpu_seconds'] = time.process_time() - start
        report['methods']['esm2'] = {'result': result, **{k: v for k, v in execution.items() if k != 'unit_scores'}}
        write_json(a.output / 'unit-scores.json', execution['unit_scores'])
    if a.stage == 'E2':
        require(all((c['matched']['scored_parents'] == 1 for c in report['methods']['esm2']['result']['conditions'].values())), 'Cleared execution validation failed')
    report.update(status='completed_descriptive_partial_method_pilot' if a.stage == 'E3' else 'completed_supported_subset_descriptive', inference_executed=not bool(a.reuse_e3), scientific_scores_produced=True, confidence_interval=None, confidence_interval_reason='fewer_than_eight_independent_parents', limitations=['Molecular-family retrieval only; no functional validation or population/model-winner claim.', 'ESM-C, stricter tiers and Renilla remain HOLD; E4 is separately cleared for matched original E3 only.', 'E2 full documented-family endpoint unresolved; AmpR retained as unknown source component.'])
    write_json(a.output / 'report.json', report)
    if a.stage == 'E3':
        write_json(a.output / 'cache-checksums.json', {name: digest_file(a.output / name) for name in ['report.json', 'execution-freeze.json', 'esm2-scores.json']})
    print(report['status'])
if __name__ == '__main__':
    main()
