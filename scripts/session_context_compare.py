#!/usr/bin/env python3
"""Compare saved join-distance profiles on the identical frozen repeat panel."""
import argparse
from collections import Counter
import json
from pathlib import Path
import statistics
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, canonical, digest, file_hash
from iri_constructs.session_benchmark import load_suite
from session_experiment import check_deadline, public_query
from session_context_profile_audit import completed_job, independent_summaries, check_close


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def require_common_panel(identities, expected):
    wanted = set(expected)
    require(len(expected) == len(wanted) == 32, 'Frozen repeat panel must contain 32 unique queries')
    for identity in identities:
        ids = identity['query_ids']
        require(len(ids) == len(set(ids)) and wanted <= set(ids), 'Source omits or duplicates required common queries')
    return sorted(wanted)


def require_shared_policy(identities, design, design_sha256, sequence_hashes):
    require(identities, 'Missing source profiles')
    for identity in identities:
        require(identity['design'] == design and identity['design_sha256'] == design_sha256,
                'Profile design or join geometry differs from the frozen design')
        require(identity['quantiles'] == [0., .05, .25, .5, .75, .95, 1.] and identity['tolerance'] == 1e-6,
                'Profile quantiles or tolerance differ')
        require(identity['numerical_policy'] == 'per-residue comparisons and reductions float64; stored means float32' and
                identity['storage_policy'] == 'standalone residue states only; joined matrices and per-residue deltas never persisted',
                'Profile numerical or storage policy differs')
        require(identity['sequence_sha256'] == sequence_hashes, 'Profile canonical target sequences differ')
        require(identity['code_hashes'] == identities[0]['code_hashes'] and
                identity['code_hashes']['scripts/session_context_profile.py'] ==
                '6c7cf52848981adf1b707cc25948c9ea5850e55fc9fa205c2b73fadd8666e18d',
                'Source profile implementations differ from the reviewed comparison')


def contrast(rows):
    return {'group_target_pairs': len(rows), 'targets': len({r['accession'] for r in rows}),
            'matched_query_component_observations': sum(r['matched_query_count'] for r in rows),
            'cosine_far_minus_near_mean': statistics.fmean(r['cosine_far_minus_near'] for r in rows) if rows else None,
            'relative_l2_far_minus_near_mean': statistics.fmean(r['relative_l2_far_minus_near'] for r in rows) if rows else None,
            'positive_cosine_difference_pairs': sum(r['cosine_far_minus_near'] > 0 for r in rows),
            'negative_l2_difference_pairs': sum(r['relative_l2_far_minus_near'] < 0 for r in rows)}


def plot(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    bins = ['0-7', '8-31', '32-63', '64-127', '128+']
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6), layout='constrained')
    for source in report['sources']:
        rows = {r['bin']: r for r in source['common_summary']['by_scope_suite_bin'] if r['scope'] == 'combined' and r['suite'] == 'all'}
        for ax, metric in zip(axes, ('cosine', 'relative_l2')):
            ax.plot(range(5), [rows[b][metric]['group_target_macro_mean'] for b in bins], marker='o', linewidth=2,
                    label=source['label'])
            ax.set_xticks(range(5), bins)
            ax.set_xlabel('Residues from nearest artificial join')
            ax.grid(axis='y', alpha=.2)
            ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('Cosine: joined versus isolated residue state')
    axes[1].set_ylabel('Relative L2 change from isolated state')
    axes[0].legend(frameon=False)
    support = report['combined_bin_support']
    fig.suptitle('Context changes are largest near artificial joins', fontsize=15)
    fig.text(.5, -.035, 'Same 16 combined queries; clean 16 excluded. Group-target means; bin support: ' +
             ', '.join(str(support[b]['group_target_pairs']) for b in bins) + '. Descriptive, correlated observations.',
             ha='center', fontsize=9)
    fig.savefig(path, dpi=180, bbox_inches='tight')
    plt.close(fig)


def run(args):
    session = args.session.resolve()
    check_deadline(session)
    started = time.monotonic()
    suite = load_suite(ROOT, session / 'inputs')
    query_map = {q['query_id']: q for q in suite['queries']}
    wanted = suite['manifest']['repeat_query_ids']
    paths = [p.resolve() for p in args.profile_job]
    require(len(paths) == len(set(paths)) == 2, 'This contrast requires exactly two distinct source profiles')
    identities = [read(path / 'identity.json') for path in paths]
    design_path = session / 'inputs/context_profile_design_v1.json'
    design_sha256 = file_hash(design_path)
    require(design_sha256 == '5b475b60b21447f91698bf1c4ed1d33fd10b891c0afb2de2931181d46e9b0ec1',
            'Frozen profile design changed')
    require_shared_policy(identities, read(design_path), design_sha256,
                          {a: digest(suite['sequences'][a].encode()) for a in suite['manifest']['targets']})
    ids = require_common_panel(identities, wanted)
    sources = []
    for path, identity in zip(paths, identities):
        require(identity['mode'] == 'context_profile_v1' and identity['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'],
                'Profile source mode or challenge differs')
        require(identity['subset'] in ('full', 'repeat'), 'Unknown profile coverage')
        expected_ids = set(query_map) if identity['subset'] == 'full' else set(wanted)
        require(set(identity['query_ids']) == expected_ids, 'Profile coverage does not match its declared panel')
        require(identity['targets'] == sorted(suite['manifest']['targets']), 'Standalone target panel differs')
        expected = {'query-' + q for q in identity['query_ids']} | {'standalone-' + a for a in identity['targets']} | {'summary'}
        source = completed_job(path, expected)
        common_records = []
        for qid in ids:
            record = read(source['units']['query-' + qid] / 'record.json')
            require(record['query'] == public_query(query_map[qid]), 'Common query metadata differs')
            require(record['scope'] == ('clean' if query_map[qid]['component_count'] == 1 else 'combined'),
                    'Clean/combined profile scope differs')
            common_records.append(record)
        bins, matched = independent_summaries(common_records)
        if identity['subset'] == 'repeat':
            native = read(source['units']['summary'] / 'record.json')
            require(native['query_count'] == 32 and native['standalone_target_count'] == 16, 'Repeat summary support differs')
            check_close(native['summary']['by_scope_suite_bin'], bins, 'Independent repeat bin summary')
            check_close(native['summary']['matched_far128plus_vs_near0to7_by_group_target'], matched, 'Independent repeat paired contrast')
        sources.append({'job': str(path), 'identity_sha256': source['identity_hash'],
                        'completion_sha256': source['completion_sha256'], 'method': identity['method'],
                        'label': identity['config']['model_id'].split('/')[-1], 'original_query_count': len(identity['query_ids']),
                        'common_summary': {'by_scope_suite_bin': bins, 'matched_group_target_contrasts': matched},
                        'matched_near_far': contrast(matched)})
    common_support = []
    for source in sources:
        common_support.append({r['bin']: {k: r[k] for k in ('queries', 'groups', 'target_count', 'group_target_pairs',
                                 'component_bin_observations', 'residue_observations')} for r in source['common_summary']['by_scope_suite_bin']
                               if r['scope'] == 'combined' and r['suite'] == 'all'})
    require(common_support[0] == common_support[1], 'Compared distance bins have different support')
    matched_keys = [{(r['group_id'], r['accession'], tuple(r['query_ids'])) for r in source['common_summary']['matched_group_target_contrasts']}
                    for source in sources]
    require(matched_keys[0] == matched_keys[1], 'Matched near/far populations differ between models')
    identity = {'mode': 'context_common_panel_comparison_v1', 'source_evidence': [{k: s[k] for k in ('job', 'identity_sha256', 'completion_sha256')} for s in sources],
                'challenge_sha256': suite['verification']['challenge_manifest_sha256'], 'common_query_ids': ids,
                'code_hashes': {name: file_hash(ROOT / name) for name in ('scripts/session_context_compare.py',
                               'scripts/session_context_profile_audit.py', 'src/iri_constructs/session_store.py')}}
    store = UnitStore(session / 'analysis' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    report = {'sources': sources, 'common_query_ids': ids, 'query_count_by_suite': dict(Counter(query_map[q]['suite'] for q in ids)),
              'combined_bin_support': common_support[0], 'summary_reconstruction': 'Independent math.fsum grouping and sorted quantiles from checksum-verified saved scalar observations',
              'limitation': 'No fresh residue inference or independent per-residue scalar reconstruction here. Correlated group/target observations; distance-bin availability varies with component length. Cross-model cosine scales are not calibrated causal measures.',
              'interpretation': 'Same-panel descriptive context perturbation; no causal decomposition of retrieval loss.'}
    old = store.read('report')
    if old is None:
        store.publish('report', report)
    else:
        require(old == report, 'Regenerated common-panel report differs')
    expected = ['report']
    if args.plots:
        path = store.root / 'common_repeat_context_profile.png'
        export = store.read('figure')
        if export is None:
            plot(report, path)
            store.publish('figure', {'path': path.name, 'sha256': file_hash(path), 'bytes': path.stat().st_size})
        else:
            require(file_hash(path) == export['sha256'], 'Saved figure changed')
        expected.append('figure')
    store.complete(expected, {'elapsed_seconds': time.monotonic() - started})
    print(json.dumps({'status': 'complete', 'job': str(store.root)}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--profile-job', type=Path, action='append', required=True)
    parser.add_argument('--name', required=True)
    parser.add_argument('--plots', action='store_true')
    run(parser.parse_args())
