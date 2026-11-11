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

def verify_score_cache(path):
    hashes = json.loads((path / 'cache-checksums.json').read_text())
    for name in ['report.json', 'execution-freeze.json', 'esm2-scores.json']:
        if digest_file(path / name) != hashes[name]:
            raise InputError('Cached score evidence changed')

def freeze_model(metadata, path):
    assert metadata['resolved_revision'] == ESM2_REVISION
    assert metadata['weight_sha256'] == '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189'
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
    assert digest_file(proposal) == '627fb59b5e8b6a0b307447d2e22355fc75449bb02c21c3063468004948928349'
    from iri_constructs.review_freeze import verify_review_freeze
    verify_review_freeze(ROOT, json.loads((ROOT / 'docs/week2/review-freeze.json').read_text()))
    refs = json.loads(proposal.read_text())
    auth = json.loads((ROOT / 'docs/oversight/EXPERIMENT_RUN_AUTHORIZATION.json').read_text())
    assert auth['review_comment_id'] == 6800553713
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
        assert proteins[q['query_id']].sequence_sha256 == q['protein_sequence_sha256']
        assert proteins[q['query_id']].source_sha256 == q['source_sha256']
    reference_ids = []
    for r in refs['references']:
        raw = ROOT / f"results/week2/references/raw_references/{r['accession']}.json"
        assert digest_file(raw) == r['original_json_sha256']
        doc = json.loads(raw.read_text())
        seq = doc['sequence']['value']
        assert hashlib.sha256(seq.encode()).hexdigest() == r['protein_sequence_sha256']
        iid = r['reference_id']
        reference_ids.append(iid)
        proteins[iid] = ProteinInput(iid, iid, seq, r['original_json_sha256'])
        labels[iid] = dict(input_id=iid, family=r['proposed_family'], annotation_status='known', review_status='approved', evidence=URL, independence_group=r['lineage_group'], archetype='unknown')
    assert len(reference_ids) == 17
    if a.stage == 'E3':
        scoring = base_scoring_config(json.loads((ROOT / 'config/week2/E3.proposed.json').read_text()))
        scoring['models'] = ['alignment', 'esm2']
        panel = freeze_panel(scoring, proteins, labels, [q['query_id'] for q in queries], reference_ids)
        allowed = {r['query_id']: r['eligible_reference_ids_proposed'] for r in refs['eligibility_proposals']}
        for q in queries:
            qid = q['query_id']
            assert set(allowed[qid]) <= set(panel['eligible_reference_ids'][qid])
            panel['eligible_reference_ids'][qid] = allowed[qid]
        panel.pop('panel_sha256')
        panel['panel_sha256'] = _hash(panel)
    else:
        panel = freeze_multi_cds_panel(cfg, proteins, labels, parents, reference_ids)
        expected = next((r['eligible_reference_ids_proposed'] for r in refs['eligibility_proposals'] if r['query_id'] == queries[0]['query_id']))
        assert panel['parents']['ncbi-AF298787.1']['eligible_reference_ids'] == expected
    ids = set(reference_ids) | {q['query_id'] for q in queries}
    assert len(ids) <= 64 and sum((len(proteins[i].sequence) for i in ids)) <= 65536
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
        assert adapter.metadata['weight_sha256'] == '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189'
        cached = _CachedAdapter(adapter)
        with Measure() as measure:
            vectors, details = _vectors(cached, proteins, ids)
            scores = _vector_scores(vectors, panel)
        result = score_reference_tier(scores, panel, labels)
        assert result['scored_queries'] == len(queries), 'Incomplete cleared query support'
        report['methods']['esm2'] = {'result': result, 'costs': dict(measure.data, cold_load_seconds=load.data['wall_seconds'], cpu_seconds=time.process_time() - start, model_calls=sum((e['status'] == 'ok' and (not e['cache_hit']) for e in cached.events)), processed_residues=sum((e['length_aa'] for e in cached.events if e['status'] == 'ok' and (not e['cache_hit']))), cache_hits=sum((e['cache_hit'] for e in cached.events))), 'model': adapter.metadata, 'inputs': details}
        write_json(a.output / 'esm2-scores.json', scores)
        adapter.close()
    elif a.reuse_e3:
        verify_score_cache(a.reuse_e3)
        original = json.loads((a.reuse_e3 / 'report.json').read_text())
        frozen = json.loads((a.reuse_e3 / 'execution-freeze.json').read_text())
        assert original['clearance'] == URL and original['methods']['esm2']['model']['resolved_revision'] == ESM2_REVISION
        assert original['methods']['esm2']['model']['weight_sha256'] == '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189'
        for iid in ids:
            assert frozen['input_fingerprints'][iid] == proteins[iid].provenance(), 'Changed cached source identity'
        scores = json.loads((a.reuse_e3 / 'esm2-scores.json').read_text())
        start = time.perf_counter()
        result = score_multi_cds({q['query_id']: scores[q['query_id']] for q in queries}, panel, labels)
        report['methods']['esm2'] = {'result': result, 'costs': {'model_calls': 0, 'wall_seconds': time.perf_counter() - start, 'shared_execution': str(a.reuse_e3)}, 'cache_note': 'Reuses exact same-model E3 unit/reference scores; no new cold experiment or independent observation.'}
    else:
        start = time.process_time()
        execution = produce_multi_cds_scores('esm2', proteins, panel, cfg, device='cpu', allow_download=False, freeze_model=lambda m: freeze_model(m, a.output / 'model-freeze.json'))
        assert execution['method_metadata']['weight_sha256'] == '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189'
        result = score_multi_cds(execution['unit_scores'], panel, labels)
        result['paired_condition_comparison'] = compare_multi_cds_conditions(result, panel, cfg)
        execution['costs']['cpu_seconds'] = time.process_time() - start
        report['methods']['esm2'] = {'result': result, **{k: v for k, v in execution.items() if k != 'unit_scores'}}
        write_json(a.output / 'unit-scores.json', execution['unit_scores'])
    if a.stage == 'E2':
        assert all((c['matched']['scored_parents'] == 1 for c in report['methods']['esm2']['result']['conditions'].values()))
    report.update(status='completed_descriptive_partial_method_pilot' if a.stage == 'E3' else 'completed_supported_subset_descriptive', inference_executed=not bool(a.reuse_e3), scientific_scores_produced=True, confidence_interval=None, confidence_interval_reason='fewer_than_eight_independent_parents', limitations=['Molecular-family retrieval only; no functional validation or population/model-winner claim.', 'ESM-C, strictertiers,Renilla and formalE4 remainHOLD.', 'E2 full documented-family endpoint unresolved; AmpR retained as unknown source component.'])
    write_json(a.output / 'report.json', report)
    if a.stage == 'E3':
        write_json(a.output / 'cache-checksums.json', {name: digest_file(a.output / name) for name in ['report.json', 'execution-freeze.json', 'esm2-scores.json']})
    print(report['status'])
if __name__ == '__main__':
    main()
