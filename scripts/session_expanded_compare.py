#!/usr/bin/env python3
"""Compare completed nested-catalogue experiments on an unchanged query panel."""
import argparse
import csv
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash
from iri_constructs.session_benchmark import load_suite
from session_experiment import check_deadline
from session_segment_experiment import validate_source_suite

MODEL_ROSTER = {'facebook/esm2_t6_8M_UR50D', 'facebook/esm2_t30_150M_UR50D', 'facebook/esm2_t33_650M_UR50D',
                'biohub/ESMC-300M', 'biohub/ESMC-600M'}
MODEL_CODE_PATHS = {'scripts/session_expanded_model.py', 'src/iri_constructs/session_expansion.py',
                    'src/iri_constructs/session_models.py', 'src/iri_constructs/session_benchmark.py',
                    'src/iri_constructs/session_store.py', 'scripts/session_experiment.py',
                    'scripts/session_catalogue_analysis.py', 'scripts/session_segment_experiment.py', 'scripts/session_independent_audit.py'}
ALIGNMENT_CODE_PATHS = {'scripts/session_expanded_alignment.py', 'src/iri_constructs/session_expansion.py',
                        'src/iri_constructs/session_alignment.py', 'src/iri_constructs/session_store.py',
                        'src/iri_constructs/session_benchmark.py', 'scripts/session_alignment_analysis.py',
                        'scripts/session_catalogue_analysis.py', 'scripts/session_independent_audit.py'}


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def require_code_binding(identity, expected_paths):
    require(set(identity['code']) == expected_paths, 'Source scientific dependency roster differs')
    for path, expected in identity['code'].items():
        require(file_hash(ROOT / path) == expected, 'Source scientific implementation changed: ' + path)


def completed(path):
    identity = read(path / 'identity.json')
    store = UnitStore(path, identity)
    descriptors = sorted(path.glob('completion-*.json'), key=lambda p: p.stat().st_mtime)
    require(descriptors, 'Source has no completed output')
    descriptor = descriptors[-1]
    completion = read(descriptor)
    require(completion['identity_hash'] == store.identity_hash and
            len(completion['units']) == completion['expected_count'], 'Invalid completion identity or count')
    for key, unit in completion['units'].items():
        folder = (path / unit['path']).resolve()
        require(folder.is_relative_to(path.resolve()), 'Completion unit escapes source')
        verified = store.verify(folder)
        require(verified['key'] == key and 'record.json' in verified['files'], 'Completion key or scientific record binding differs')
        require(file_hash(folder / 'complete.json') == unit['descriptor_sha256'], 'Completion unit descriptor changed')
    require('summary' in completion['units'], 'Missing source summary')
    record = read(path / completion['units']['summary']['path'] / 'record.json')
    return identity, store.identity_hash, file_hash(descriptor), record, set(completion['units'])


def require_common_design(identities, query_ids):
    require(identities, 'No sources')
    expected = set(query_ids)
    require(len(expected) == len(query_ids) == 1146, 'Expected frozen full query panel')
    shared = ('expansion_manifest_sha256', 'retrieval_design_sha256', 'subsets', 'target_aliases')
    for identity in identities:
        ids = identity['query_ids']
        require(len(ids) == len(set(ids)) and set(ids) == expected, 'Source query coverage differs; no silent intersection')
        for key in shared:
            require(identity[key] == identities[0][key], 'Source catalogue or design differs: ' + key)
        require(set(identity['subsets']) == {'101', '250', '500', '1000'}, 'Expected all four nested catalogue sizes')
        prior = set()
        for size in ('101', '250', '500', '1000'):
            entries = identity['subsets'][size]
            require(len(entries) == len(set(entries)) == int(size) and prior <= set(entries), 'Catalogue is not a unique nested prefix')
            prior = set(entries)


def headline_rows(summary, label, job, mode):
    rows = []
    expected_arms = ({'whole', 'isolated_mixture', 'contextual_component', 'isolated_component'}
                     if mode == 'expanded_model_retrieval_v1' else {'local_alignment_raw', 'local_alignment_reference_selfscore'})
    require(set(summary['summaries']) == {'101', '250', '500', '1000'}, 'Summary omits a required catalogue')
    for size, arms in summary['summaries'].items():
        require(set(arms) == expected_arms, 'Summary omits or substitutes a required retrieval arm')
        for arm, scopes in arms.items():
            require(set(scopes) == {'clean', 'combined'}, 'Clean and combined populations must remain separate')
            for scope, aggregate in scopes.items():
                nq, nc = (16, 16) if scope == 'clean' else (1130, 3520)
                require(aggregate['query_count'] == nq and aggregate['component_count'] == nc, 'Summary population changed')
                exact = aggregate['kinds']['exact']
                expected_counts = {'1': (16, 16)} if scope == 'clean' else {'2': (230, 460), '3': (540, 1620), '4': (360, 1440)}
                require(set(exact['by_component_count']) == set(expected_counts) and exact['original_components'] == nc and
                        exact['parent_queries'] == nq and exact['query_metrics']['original_queries'] == nq,
                        'Overall or component-count strata differ from the frozen panel')
                for count, (stratum_queries, stratum_components) in expected_counts.items():
                    values = exact['by_component_count'][count]
                    require(values['parent_queries'] == stratum_queries and values['original_components'] == stratum_components and
                            values['query_metrics']['original_queries'] == stratum_queries, 'Component-count support differs')
                groups = [('all', exact)] + [(count, values) for count, values in exact['by_component_count'].items()]
                for count, values in groups:
                    require(values['scorable_components'] == values['original_components'] and values['unscorable_components'] == 0,
                            'Exact catalogue diagnostic unexpectedly has absent labels')
                    query = values['query_metrics']
                    require(query['excluded_unscorable_queries'] == 0 and query['all_components_scorable_queries'] == query['original_queries'],
                            'Exact query denominator unexpectedly excludes examples')
                    rows.append({'label': label, 'job': job, 'mode': mode, 'catalogue_size': int(size), 'arm': arm,
                                 'scope': scope, 'component_count': count,
                                 'budget': 'k per component (known boundaries)' if arm in ('contextual_component', 'isolated_component') else 'k total per query',
                                 **{key: values[key] for key in ('original_components', 'scorable_components', 'hits_at_1', 'hits_at_5', 'hits_at_10',
                                                               'recall_at_1', 'recall_at_5', 'recall_at_10', 'per_target_macro_recall_at_5')},
                                 **{key: query[key] for key in ('original_queries', 'all_components_scorable_queries',
                                                               'all_components_hit_at_1', 'all_components_hit_at_5', 'all_components_hit_at_10',
                                                               'any_component_hit_at_5')}})
    return rows


def plot(rows, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sizes = [101, 250, 500, 1000]
    selected = [r for r in rows if r['scope'] == 'combined' and r['component_count'] == 'all']
    labels = sorted({r['label'] for r in selected if r['mode'] == 'expanded_model_retrieval_v1'})
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.0), layout='constrained')
    for label in labels:
        for ax, arm in zip(axes, ('whole', 'contextual_component')):
            by_size = {r['catalogue_size']: r['recall_at_5'] for r in selected if r['label'] == label and r['arm'] == arm}
            require(set(by_size) == set(sizes), 'Plot is missing a catalogue endpoint')
            ax.plot(sizes, [by_size[s] for s in sizes], marker='o', linewidth=2, label=label)
    for arm, label, style in [('local_alignment_raw', 'SW raw', '--'),
                              ('local_alignment_reference_selfscore', 'SW self-normalized (exploratory)', ':')]:
        by_size = {r['catalogue_size']: r['recall_at_5'] for r in selected if r['arm'] == arm}
        require(set(by_size) == set(sizes), 'Alignment plot is missing a catalogue endpoint')
        axes[0].plot(sizes, [by_size[s] for s in sizes], marker='s', linestyle=style, color='black', label=label)
    for ax, title in zip(axes, ('Whole-query retrieval: 5 total returns', 'Context oracle: 5 returns per known component')):
        ax.set_title(title, fontsize=11)
        ax.set_xticks(sizes, ['101', '250', '500', '1,000'])
        ax.set_xlabel('Unique reference sequences')
        ax.set_ylabel('Exact-component recall')
        ax.set_ylim(0, 1.03)
        ax.grid(axis='y', alpha=.2)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8, loc='best')
    fig.suptitle('Catalogue competition on the same 1,130 joined queries', fontsize=15)
    fig.text(.5, -.035, '3,520 component observations. Nested reviewed catalogue; size and composition change together. Oracle budget is larger.',
             ha='center', fontsize=9)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def run(args):
    from iri_constructs.session_expansion import load_expansion
    session = args.session.resolve()
    check_deadline(session)
    started = time.monotonic()
    suite = load_suite(ROOT, session / 'inputs')
    expansion = load_expansion(ROOT, session, session / 'expansion/milestone-1000/admitted/expanded_catalogue_v1.json')
    retrieval_design_sha256 = file_hash(session / 'inputs/expanded_retrieval_design_v1.json')
    require(len(args.model_job) == len(args.parent_job) == args.expected_models == len(MODEL_ROSTER),
            'Explicit model/parent roster is incomplete')
    sources, identities, rows = [], [], []
    model_ids = []
    for job_arg, parent_arg in zip(args.model_job, args.parent_job):
        job, parent = job_arg.resolve(), parent_arg.resolve()
        identity, identity_hash, completion_hash, summary, units = completed(job)
        require(identity['mode'] == 'expanded_model_retrieval_v1' and units == {'query-' + q for q in identity['query_ids']} | {'summary'},
                'Wrong or incomplete expanded model source')
        require_code_binding(identity, MODEL_CODE_PATHS)
        parent_identity = read(parent / 'identity.json')
        require(digest(canonical(parent_identity)) == identity['parent_identity'] and parent_identity['mode'] == 'model', 'Wrong model parent')
        validate_source_suite(parent_identity, suite, suite['verification']['challenge_manifest_sha256'])
        require(parent_identity['query_ids'] == identity['query_ids'], 'Expanded source omitted or changed parent query coverage')
        parent_completions = [p for p in parent.glob('completion-*.json') if file_hash(p) == identity['parent_completion_sha256']]
        require(len(parent_completions) == 1, 'Pinned model parent completion is missing or ambiguous')
        method = parent_identity['method']
        require(method['device'] == 'cpu' and method['dtype'] == 'float32' and method['layer'] == method['num_layers'],
                'This report compares only full-depth CPU FP32 sources')
        model_id = method['model_id']
        require(model_id in MODEL_ROSTER and parent_identity['config'] ==
                read(session / 'inputs/configs' / (model_id.split('/')[-1] + '-cpu.json')),
                'Source model configuration differs from its frozen primary configuration')
        require(model_id not in model_ids, 'Repeated model in the explicit roster')
        model_ids.append(model_id)
        require(summary['independent_acceptance']['all_rankings_evaluations_summaries_checked'] is True, 'Missing independent acceptance')
        rows.extend(headline_rows(summary, model_id.split('/')[-1], job.name, identity['mode']))
        identities.append(identity)
        sources.append({'job': str(job), 'identity_hash': identity_hash, 'completion_sha256': completion_hash,
                        'parent_job': str(parent), 'parent_identity_hash': identity['parent_identity'], 'method': method})
    require(set(model_ids) == MODEL_ROSTER, 'A required model is missing or substituted in the fixed five-model roster')
    job = args.alignment_job.resolve()
    identity, identity_hash, completion_hash, summary, units = completed(job)
    require(identity['mode'] == 'expanded_alignment_retrieval_v1' and
            summary['all_query_metrics_and_summary_denominators_independently_checked'] is True, 'Wrong or unaudited expanded alignment source')
    require_code_binding(identity, ALIGNMENT_CODE_PATHS)
    require(identity['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'] and
            identity['query_sequence_sha256'] == {q['query_id']: q['joined_sequence_sha256'] for q in suite['queries']},
            'Expanded alignment challenge or query sequences differ')
    expected_hashes = {a: digest(sequence.encode()) for a, sequence in expansion['sequences'].items()}
    shards = identity['new_reference_shards']
    flat = [a for shard in shards for a in shard]
    expected_new = [a for a in expansion['manifest']['new_representatives_hash_order']
                    if a in expected_hashes and a not in expansion['subsets']['101']]
    require(identity['reference_sequence_sha256'] == expected_hashes and flat == expected_new and
            len(flat) == len(set(flat)) == 899 and all(1 <= len(shard) <= 32 for shard in shards),
            'Expanded alignment shards or canonical reference identities are incomplete')
    expected_units = {'summary', 'reference-selfscores'} | {'query-' + q for q in identity['query_ids']} | {
        f'raw-{q}-{i:03d}' for q in identity['query_ids'] for i in range(len(identity['new_reference_shards']))}
    require(units == expected_units, 'Expanded alignment completion omits or adds query/shard evidence')
    require(identity['normalization_design_sha256'] == file_hash(session / 'audit/alignment_normalization_design_v1.json'),
            'Exploratory alignment normalization changed')
    rows.extend(headline_rows(summary, 'Smith-Waterman', job.name, identity['mode']))
    identities.append(identity)
    sources.append({'job': str(job), 'identity_hash': identity_hash, 'completion_sha256': completion_hash, 'method': identity['method']})
    require_common_design(identities, [q['query_id'] for q in suite['queries']])
    require(identities[0]['expansion_manifest_sha256'] == expansion['manifest_sha256'] and
            identities[0]['retrieval_design_sha256'] == retrieval_design_sha256 and
            identities[0]['subsets'] == expansion['subsets'] and identities[0]['target_aliases'] == expansion['target_aliases'],
            'Reported catalogue or labels differ from independently reloaded reviewed expansion')
    identity = {'mode': 'expanded_common_panel_comparison_v1', 'sources': sources,
                'expected_model_count': args.expected_models, 'challenge_sha256': suite['verification']['challenge_manifest_sha256'],
                'shared_design': {k: identities[0][k] for k in ('expansion_manifest_sha256', 'retrieval_design_sha256', 'subsets', 'target_aliases')},
                'code': {name: file_hash(ROOT / name) for name in ('scripts/session_expanded_compare.py', 'src/iri_constructs/session_store.py',
                                                                  'src/iri_constructs/session_benchmark.py', 'src/iri_constructs/session_expansion.py',
                                                                  'scripts/session_segment_experiment.py')}}
    store = UnitStore(session / 'analysis' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    report = {'model_ids': model_ids, 'rows': rows,
              'interpretation': 'Descriptive exact-reference retrieval on identical frozen queries and nested reviewed catalogues. Size and competitor composition are confounded; not a training holdout or function validation.',
              'budgets': 'Whole/mixture/alignment use k total; contextual/isolated oracle use k per supplied component. Do not compare these as equal total budgets.',
              'verification': 'Every source completion unit checksum verified; original experiment independently checked rankings, evaluations and summaries. No inference or scoring rerun here.'}
    old = store.read('report')
    if old is None:
        store.publish('report', report)
    else:
        require(old == report, 'Regenerated comparison differs')
    exports = store.read('exports')
    if exports is None:
        table = store.root / 'catalogue_comparison.csv'
        with table.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        paths = [table]
        if args.plots:
            picture = store.root / 'catalogue_competition.png'
            plot(rows, picture)
            paths.append(picture)
        store.publish('exports', {p.name: {'sha256': file_hash(p), 'bytes': p.stat().st_size} for p in paths})
    else:
        for name, descriptor in exports.items():
            require(file_hash(store.root / name) == descriptor['sha256'], 'Saved export changed')
        require(not args.plots or 'catalogue_competition.png' in exports, 'Existing completion did not request plots')
    store.complete(['report', 'exports'], {'elapsed_seconds': time.monotonic() - started})
    print(json.dumps({'status': 'complete', 'job': str(store.root)}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--model-job', type=Path, action='append', required=True)
    parser.add_argument('--parent-job', type=Path, action='append', required=True)
    parser.add_argument('--expected-models', type=int, required=True)
    parser.add_argument('--alignment-job', type=Path, required=True)
    parser.add_argument('--name', required=True)
    parser.add_argument('--plots', action='store_true')
    run(parser.parse_args())
