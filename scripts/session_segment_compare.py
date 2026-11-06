#!/usr/bin/env python3
"""Descriptive held-out segmentation comparison from immutable completed jobs.

No inference, alignment, selection, or recalibration occurs here. The original
950 combined test queries and 16 clean controls remain separate populations.
Family transfer is a separate annotation-labelled, reference-removed analysis.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.session_store import UnitStore, atomic_bytes, canonical, digest, file_hash
from iri_constructs.session_benchmark import evaluate_query, evaluate_oracle_query, filter_ranking
from iri_constructs.session_segmentation import per_segment_union
from session_catalogue_analysis import (SummaryAccumulator, latest_completion, source_arms,
    subset_denominators, denominator_check, validate_family_freeze)
from session_compare_models import csv_bytes, display_name, primary_arm, comparison_role
from session_experiment import check_deadline, public_query
from session_independent_audit import verify_unit, validate_ranking
from session_segment_experiment import (setup, expected_keys, prediction_config, validate_proposal,
    validate_retrieval_record, validate_source_suite)
from session_family_segmentation import (expected_transfer_keys, mask_for_query, unit_key,
    validate_transfer_record, validate_parent_development_method)

KS = (1, 5, 10)
CORE_MODELS = {'facebook/esm2_t6_8M_UR50D', 'facebook/esm2_t30_150M_UR50D',
               'facebook/esm2_t33_650M_UR50D', 'biohub/ESMC-300M'}
ORACLES = {'contextual_component', 'isolated_component'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def expected_source_keys(identity):
    mode = identity['mode']
    if mode.startswith('family_automatic_'):
        return expected_transfer_keys(identity)
    if mode == 'alignment_normalization':
        ids = identity['query_ids']
        require(ids and len(ids) == len(set(ids)), 'Empty or duplicate normalized query IDs')
        return {'query-' + q for q in ids} | {'summary', 'reference-selfscores'}
    return expected_keys(identity)


class Snapshot:
    """Pin a completion descriptor once; never look up latest loose attempts."""
    def __init__(self, path):
        path = Path(path).resolve()
        self.root = path.parent if path.is_file() else path
        self.completion_path = path if path.is_file() else latest_completion(self.root)
        self.identity = read(self.root / 'identity.json')
        self.identity_hash = digest(canonical(self.identity))
        self.completion_sha256 = file_hash(self.completion_path)
        descriptor = read(self.completion_path)
        expected = expected_source_keys(self.identity)
        require(descriptor['identity_hash'] == self.identity_hash
                and descriptor['expected_count'] == len(expected)
                and set(descriptor['units']) == expected,
                'Completed snapshot lacks exact expected unit coverage')
        self.paths = {key: verify_unit(self.root, key, entry, self.identity_hash)
                      for key, entry in descriptor['units'].items()}
        self.index = descriptor['units']

    def record(self, key):
        require(key in self.paths, 'Requested unit is outside pinned completion')
        return read(self.paths[key] / 'record.json')

    def reference(self, key, pointer=''):
        return {'job_identity_hash': self.identity_hash, 'completion_sha256': self.completion_sha256,
                'unit_key': key, 'unit_path': str(self.paths[key]),
                'unit_descriptor_sha256': self.index[key]['descriptor_sha256'], 'json_pointer': pointer}

    def provenance(self):
        return {'job_id': self.root.name, 'path': str(self.root), 'identity_hash': self.identity_hash,
                'identity_file_sha256': file_hash(self.root / 'identity.json'),
                'completion_path': str(self.completion_path), 'completion_sha256': self.completion_sha256,
                'mode': self.identity['mode'], 'source_query_count': len(self.identity['query_ids'])}


def validate_panel(source_ids, panel_ids, *, exact=False):
    require(source_ids and len(source_ids) == len(set(source_ids)), 'Duplicate or empty source query IDs')
    require(len(panel_ids) == len(set(panel_ids)), 'Duplicate frozen panel query IDs')
    require(set(panel_ids) <= set(source_ids), 'Source is missing held-out panel queries')
    if exact:
        require(set(source_ids) == set(panel_ids), 'Segmentation source contains queries outside frozen test panel')
    return {'source_queries': len(source_ids), 'compared_queries': len(panel_ids),
            'source_queries_outside_comparison': len(set(source_ids) - set(panel_ids)),
            'missing_panel_queries': 0}


def validate_parent(child, parent, *, development_method=None):
    require(parent.identity['mode'] == 'model', 'Segmentation parent is not a core model')
    require(child.identity['parent_model_identity'] == parent.identity_hash
            and child.identity['parent_completion_sha256'] == parent.completion_sha256,
            'Segmentation parent snapshot differs')
    require(child.identity['model'] == parent.identity['method'], 'Segmentation and whole-query representations differ')
    if development_method is not None:
        require(child.identity['model'] == development_method, 'Family transfer differs from exact frozen development method')


def scopes(query):
    return ['clean' if len(query['components']) == 1 else 'combined', 'suite:' + query['suite']]


def component_strata(component, count):
    return [('overall', 'all'), ('component_count', str(count)), ('fraction', component['fraction_bin']),
            ('position', component['position']), ('target', component['accession'])]


def return_volume(rankings, query, k, *, acceptable=None, budget='k total per original query'):
    """Count review volume without treating unlabelled candidates as functional negatives."""
    lists = list(rankings.values()) if isinstance(rankings, dict) else rankings
    if budget == 'k total per original query':
        lists = [rankings]
    union = per_segment_union(lists, k)
    present = set(union['accessions'])
    labels = set(query['components']) if acceptable is None else set().union(
        *(set(acceptable.get(a, [])) for a in query['components']))
    return {'k': k, 'retrieval_budget': budget, 'retrieval_lists': len(lists),
            'returned_hit_occurrences': union['returned_hit_occurrences'],
            'unique_returned_hits': union['unique_returned_hits'],
            'duplicate_returned_hits': union['duplicate_returned_hits'],
            'outside_exact_labels': len(present - set(query['components'])),
            'outside_acceptable_labels': len(present - labels),
            'outside_label_interpretation': 'Not a verified functional false-positive count.'}


def union_observation(evaluation, rankings, query, *, acceptable=None):
    """Membership recovery at k per predicted segment; no fictitious union rank."""
    kind = 'labelled' if acceptable is not None else 'exact'
    observations = {}
    for k in KS:
        union = per_segment_union(rankings, k)
        present = set(union['accessions'])
        components = []
        for component in evaluation['components']:
            a = component['accession']
            surviving = set(component[kind]['accepted_surviving'])
            components.append({**{field: component[field] for field in (
                'accession', 'component_index', 'position', 'fraction', 'fraction_bin')},
                'scorable': bool(surviving), 'hit': int(bool(surviving & present)) if surviving else None})
        n = sum(c['scorable'] for c in components)
        hits = sum(c['hit'] for c in components if c['scorable'])
        full = n == len(components)
        observations[str(k)] = {'components': components, 'scorable_components': n,
            'component_hits': hits, 'all_components_scorable': full,
            'any_component_hit': int(hits > 0) if full else None,
            'all_components_hit': int(hits == len(components)) if full else None,
            **return_volume(rankings, query, k, acceptable=acceptable, budget='k per predicted segment')}
    return {'query_id': query['query_id'], 'group_id': query['group_id'], 'suite': query['suite'],
            'component_count': len(query['components']), 'label_kind': kind, 'by_k': observations}


def union_tables(observations, metadata):
    require(len(observations) == len({o['query_id'] for o in observations}), 'Duplicate union query observations')
    components, queries = [], []
    for scope in sorted({s for o in observations for s in scopes({'components': [None]*o['component_count'], 'suite': o['suite']})}):
        rows = [o for o in observations if scope in scopes({'components': [None]*o['component_count'], 'suite': o['suite']})]
        for k in KS:
            component_groups = defaultdict(list)
            query_groups = defaultdict(list)
            for observation in rows:
                value = observation['by_k'][str(k)]
                for dimension, key in [('overall', 'all'), ('component_count', str(observation['component_count']))]:
                    query_groups[(dimension, key)].append((observation, value))
                for c in value['components']:
                    for dimension, key in component_strata(c, observation['component_count']):
                        component_groups[(dimension, key)].append((observation, c))
            base = {**metadata, 'scope': scope, 'k': k, 'label_kind': rows[0]['label_kind'],
                    'retrieval_budget': 'k per predicted segment'}
            for (dimension, key), values in sorted(component_groups.items()):
                scored = [(q, c) for q, c in values if c['scorable']]
                target = defaultdict(list)
                for _, c in scored:
                    target[c['accession']].append(c['hit'])
                hits = sum(c['hit'] for _, c in scored)
                components.append({**base, 'stratum_dimension': dimension, 'stratum': key,
                    'original_components': len(values), 'scorable_components': len(scored),
                    'unscorable_components': len(values)-len(scored), 'component_hits': hits,
                    'component_recall': hits/len(scored) if scored else None,
                    'per_target_macro_component_recall': math.fsum(sum(v)/len(v) for v in target.values())/len(target) if target else None,
                    'parent_queries': len({q['query_id'] for q, _ in values}),
                    'distinct_groups': len({q['group_id'] for q, _ in values})})
            for (dimension, key), values in sorted(query_groups.items()):
                fully = [v for _, v in values if v['all_components_scorable']]
                any_hits = sum(v['any_component_hit'] for v in fully)
                all_hits = sum(v['all_components_hit'] for v in fully)
                queries.append({**base, 'stratum_dimension': dimension, 'stratum': key,
                    'original_queries': len(values), 'all_components_scorable_queries': len(fully),
                    'excluded_unscorable_queries': len(values)-len(fully),
                    'any_component_hit_count': any_hits, 'all_components_hit_count': all_hits,
                    'any_component_hit': any_hits/len(fully) if fully else None,
                    'all_components_hit': all_hits/len(fully) if fully else None,
                    **{field: sum(v[field] for _, v in values) for field in (
                        'returned_hit_occurrences', 'unique_returned_hits', 'duplicate_returned_hits',
                        'outside_exact_labels', 'outside_acceptable_labels')}})
    return components, queries


def boundary_rows(query, proposal, metadata):
    evaluation = proposal['boundary_evaluation']
    # Geometry was independently regenerated by validate_proposal before this call.
    sample = evaluate_query([], query, {})
    result = []
    for boundary, component in zip(evaluation['components'], sample['components']):
        require(boundary['accession'] == component['accession'], 'Boundary/component ordering differs')
        result.append({**metadata, 'query_id': query['query_id'], 'group_id': query['group_id'],
            'suite': query['suite'], 'component_count': len(query['components']),
            'scope': 'clean' if len(query['components']) == 1 else 'combined',
            **{k: component[k] for k in ('accession', 'position', 'fraction', 'fraction_bin')}, **boundary})
    return result


def resource_row(query, proposal, record, metadata):
    scan, encoded = proposal['logical_candidate_scan_residues'], proposal['logical_encoded_segment_residues']
    require(proposal['logical_total_queried_residues'] == scan + encoded, 'Work must include scanning AND encoding')
    return {**metadata, 'query_id': query['query_id'], 'group_id': query['group_id'], 'suite': query['suite'],
            'component_count': len(query['components']), 'scope': 'clean' if len(query['components']) == 1 else 'combined',
            'predicted_segments': proposal['prediction']['segment_count'],
            'cap_bound': int(proposal['prediction']['cap_bound']),
            'whole_query_fallback': int(proposal['prediction'].get('whole_query_fallback_used', False)),
            'logical_candidate_scan_residues': scan, 'logical_encoded_segment_residues': encoded,
            'logical_total_queried_residues': scan + encoded,
            'prediction_unit_seconds': proposal['seconds'], 'retrieval_unit_seconds': record['seconds'],
            'cache_counts': record.get('cache_counts', {'cache_hits': record.get('cache_hits', 0),
                                                      'new_encodes': record.get('new_encodes', 0)}),
            'ranking_reused_within_query': record.get('ranking_reused_within_query', False),
            'prediction_cumulative_cache_counts': proposal.get('mask_cumulative_cache_counts', {
                'source_cache_hits': proposal.get('query_cumulative_cache_hits', 0),
                'fresh_alignment_calls': proposal.get('query_cumulative_fresh_alignment_calls', 0)}),
            'runtime_interpretation': 'Observed saved-unit wall time; cache/restart dependent, excludes model loading. Shared prediction costs are not additive across models.'}


def summarize_diagnostics(boundaries, resources):
    boundary_summary, resource_summary = [], []
    bgroups, rgroups = defaultdict(list), defaultdict(list)
    identity_fields = ('job_id', 'configuration_id', 'policy', 'scope')
    for row in boundaries:
        component = {**row, 'accession': row['accession']}
        for dimension, key in component_strata(component, row['component_count']):
            bgroups[tuple(row[x] for x in identity_fields) + (dimension, key)].append(row)
    for key, rows in sorted(bgroups.items()):
        boundary_summary.append({**dict(zip(identity_fields, key[:4])), 'stratum_dimension': key[4], 'stratum': key[5],
            'component_observations': len(rows), 'queries': len({r['query_id'] for r in rows}),
            'mean_best_interval_iou': math.fsum(r['best_interval_iou'] for r in rows)/len(rows),
            'mean_covered_fraction_by_union': math.fsum(r['covered_fraction_by_union'] for r in rows)/len(rows),
            'exact_boundary_match_count': sum(r['exact_boundary_match'] for r in rows),
            'exact_boundary_match_fraction': sum(r['exact_boundary_match'] for r in rows)/len(rows)})
    for row in resources:
        for dimension, key in [('overall', 'all'), ('component_count', str(row['component_count']))]:
            rgroups[tuple(row[x] for x in identity_fields) + (dimension, key)].append(row)
    for key, rows in sorted(rgroups.items()):
        cache = defaultdict(int)
        for row in rows:
            for name, value in row['cache_counts'].items():
                cache[name] += value
        resource_summary.append({**dict(zip(identity_fields, key[:4])), 'stratum_dimension': key[4], 'stratum': key[5],
            'queries': len(rows), **{field: sum(r[field] for r in rows) for field in (
                'predicted_segments', 'cap_bound', 'whole_query_fallback', 'logical_candidate_scan_residues',
                'logical_encoded_segment_residues', 'logical_total_queried_residues',
                'prediction_unit_seconds', 'retrieval_unit_seconds')}, 'cache_counts': dict(cache),
            'prediction_cost_scope': 'One prediction source/configuration; reused costs across models must not be added.',
            'runtime_interpretation': rows[0]['runtime_interpretation']})
    return boundary_summary, resource_summary


def summarize_volumes(rows):
    groups = defaultdict(list)
    fields = ('job_id', 'arm', 'policy', 'scope', 'retrieval_budget', 'k')
    for row in rows:
        for dimension, key in [('overall', 'all'), ('component_count', str(row['component_count']))]:
            groups[tuple(row[f] for f in fields) + (dimension, key)].append(row)
    output = []
    for key, values in sorted(groups.items()):
        counts = {field: sum(r[field] for r in values) for field in (
            'returned_hit_occurrences', 'unique_returned_hits', 'duplicate_returned_hits',
            'outside_exact_labels', 'outside_acceptable_labels')}
        output.append({**dict(zip(fields, key[:6])), 'stratum_dimension': key[6], 'stratum': key[7],
            'queries': len(values), **counts, **{'mean_' + k: v/len(values) for k, v in counts.items()},
            'outside_label_interpretation': 'Not a verified functional false-positive count.'})
    return output


def summarize_baseline(rows):
    groups = defaultdict(list)
    for row in rows:
        for dimension, value in [('overall', 'all'), ('component_count', str(row['component_count']))]:
            groups[(row['job_id'], row['arm'], row['policy'], row['scope'], dimension, value)].append(row)
    result = []
    for key, values in sorted(groups.items()):
        n = len(values); components = sum(r['component_count'] for r in values)
        scorable = sum(r['scorable_components'] for r in values)
        fully = [r for r in values if r['all_components_scorable']]
        result.append({'job_id': key[0], 'arm': key[1], 'policy': key[2], 'scope': key[3],
            'stratum_dimension': key[4], 'stratum': key[5], 'queries': n, 'components': components,
            'label_kind': values[0]['label_kind'], 'scorable_components': scorable,
            'unscorable_components': components-scorable, 'all_components_scorable_queries': len(fully),
            'excluded_unscorable_queries': n-len(fully),
            'distinct_groups': len({r['group_id'] for r in values}),
            'retrieval_budget': 'k total per original query',
            'contrast': 'selected_segmenter_minus_same_method_whole',
            **{f'component_hit_difference_at_{k}': sum(r[f'component_hit_difference_at_{k}'] for r in values) for k in KS},
            **{f'component_recall_difference_at_{k}': sum(r[f'component_hit_difference_at_{k}'] for r in values)/scorable if scorable else None for k in KS},
            **{f'all_hit_rate_difference_at_{k}': sum(r[f'all_hit_difference_at_{k}'] for r in fully)/len(fully) if fully else None for k in KS},
            **{f'any_hit_rate_difference_at_{k}': sum(r[f'any_hit_difference_at_{k}'] for r in fully)/len(fully) if fully else None for k in KS}})
    return result


class ReportBook:
    def __init__(self):
        self.summaries = {}
        self.metadata = {}
        self.unions = defaultdict(list)
        self.tables = {name: [] for name in ('source_coverage', 'per_query', 'component_observations',
            'returned_volume', 'boundary_observations', 'resource_observations', 'matched_whole_baseline', 'core_runtime')}

    def add(self, metadata, query, evaluation, rankings, reference, *, acceptable=None):
        key = (metadata['job_id'], metadata['arm'], metadata['policy'])
        self.metadata[key] = metadata
        bundle = self.summaries.setdefault(key, {})
        for scope in scopes(query):
            bundle.setdefault(scope, SummaryAccumulator()).add(evaluation)
        self.tables['per_query'].append({**metadata, 'query_id': query['query_id'], 'group_id': query['group_id'],
            'suite': query['suite'], 'component_count': len(query['components']),
            'retrieval_budget': evaluation['retrieval_budget'], 'query_metrics': evaluation['query_metrics'],
            'source': reference})
        for c in evaluation['components']:
            for kind in evaluation['query_metrics']:
                metric = c[kind]
                self.tables['component_observations'].append({**metadata, 'query_id': query['query_id'],
                    'group_id': query['group_id'], 'suite': query['suite'], 'component_count': len(query['components']),
                    'retrieval_budget': evaluation['retrieval_budget'], 'label_kind': kind,
                    **{name: c[name] for name in ('accession', 'component_index', 'fraction', 'fraction_bin', 'position')},
                    **{name: metric[name] for name in ('scorable', 'rank', 'hit_at_1', 'hit_at_5', 'hit_at_10')}})
        for k in KS:
            self.tables['returned_volume'].append({**metadata, 'query_id': query['query_id'], 'group_id': query['group_id'],
                'suite': query['suite'], 'component_count': len(query['components']),
                'scope': 'clean' if len(query['components']) == 1 else 'combined',
                **return_volume(rankings, query, k, acceptable=acceptable, budget=evaluation['retrieval_budget'])})

    def add_segments(self, metadata, query, record, proposal, segment_rankings, reference, *, acceptable=None):
        self.add(metadata, query, record['evaluation'], record['ranking'], reference, acceptable=acceptable)
        key = (metadata['job_id'], metadata['arm'], metadata['policy'])
        self.unions[key].append(union_observation(record['evaluation'], segment_rankings, query, acceptable=acceptable))
        self.tables['boundary_observations'] += boundary_rows(query, proposal, metadata)
        self.tables['resource_observations'].append(resource_row(query, proposal, record, metadata))

    def result(self):
        tables = {key: list(values) for key, values in self.tables.items()}
        tables.update(component_metrics=[], query_metrics=[], headlines=[], union_component_metrics=[], union_query_metrics=[], union_per_query=[])
        summaries = []
        for key, bundle in sorted(self.summaries.items()):
            metadata = self.metadata[key]
            for scope, accumulator in sorted(bundle.items()):
                summary = accumulator.result()
                summaries.append({**metadata, 'scope': scope, 'summary': summary})
                for kind, metrics in summary['kinds'].items():
                    base = {**metadata, 'scope': scope, 'label_kind': kind, 'retrieval_budget': summary['retrieval_budget']}
                    entries = [('overall', 'all', metrics)]
                    for field, dimension in [('by_component_count', 'component_count'), ('by_fraction', 'fraction'),
                                             ('by_position', 'position'), ('by_target', 'target')]:
                        entries += [(dimension, k, v) for k, v in metrics[field].items()]
                    for dimension, stratum, value in entries:
                        row = {**base, 'stratum_dimension': dimension, 'stratum': stratum,
                               **{k: v for k, v in value.items() if not k.startswith('by_') and k != 'query_metrics'}}
                        tables['component_metrics'].append(row)
                        if 'query_metrics' in value:
                            tables['query_metrics'].append({**base, 'stratum_dimension': dimension, 'stratum': stratum,
                                                           **value['query_metrics']})
                        if (scope == 'combined' and kind == ('labelled' if metadata['policy'] != 'original' else 'exact')
                            and dimension == 'overall' and metadata['role'] != 'oracle_control'):
                            tables['headlines'].append({**row, **value['query_metrics']})
            if key in self.unions:
                crows, qrows = union_tables(self.unions[key], metadata)
                tables['union_component_metrics'] += crows
                tables['union_query_metrics'] += qrows
                tables['union_per_query'] += [{**metadata, **o, 'retrieval_budget': 'k per predicted segment'} for o in self.unions[key]]
        tables['boundary_metrics'], tables['resource_metrics'] = summarize_diagnostics(
            tables['boundary_observations'], tables['resource_observations'])
        tables['returned_volume_metrics'] = summarize_volumes(tables['returned_volume'])
        tables['matched_whole_baseline_metrics'] = summarize_baseline(tables['matched_whole_baseline'])
        # Keep family outputs physically separate from exact-reference headline tables.
        family, primary = {}, {}
        for name, rows in tables.items():
            family[name] = [r for r in rows if r.get('policy', 'original') != 'original']
            primary[name] = [r for r in rows if r.get('policy', 'original') == 'original']
        return {'summaries': [r for r in summaries if r['policy'] == 'original'],
                'family_summaries': [r for r in summaries if r['policy'] != 'original'],
                'tables': primary, 'family_tables': family,
                'interpretation': 'Descriptive held-out retrieval; correlated orders and proteins are not independent replicates. No p-values, calibration, biological-function validation, or training holdout is claimed.'}


def metadata(source, arm, *, config=None, policy='original'):
    mode = source.identity['mode']
    role = ('oracle_control' if arm in ORACLES or arm == 'isolated_mixture' else
            'exploratory_family_transfer' if policy != 'original' else
            'selected_segmentation' if config is not None else comparison_role(mode))
    return {'job_id': source.root.name, 'mode': mode, 'arm': arm, 'policy': policy, 'role': role,
            'configuration_id': config['configuration_id'] if config else '',
            'segmenter': config['segmenter'] if config else '',
            'label': (source.identity['model']['model_id'] + ' / ' + config['configuration_id']) if config else display_name(source.identity)}


def cached_family_rankings(record, query, identity, session, verified):
    expected = {'mode': 'family_transfer_segment_cache_v1', 'model': identity['model'],
                'pooling': identity['pooling'], 'catalogue_sequence_sha256': identity['catalogue_sequence_sha256'],
                'code_hashes': identity['code_hashes']}
    cache_hash = digest(canonical(expected))
    root = (session / 'caches' / ('family-segments-' + cache_hash[:16])).resolve()
    require(read(root / 'identity.json') == expected, 'Family segment cache representation changed')
    references, spans = record['segment_ranking_references'], record['prediction']['segments']
    require(len(references) == len(spans), 'Family segment cache coverage differs')
    result = []
    for evidence, span in zip(references, spans):
        sequence = query['sequence'][span['start0']:span['end0']]
        sha = digest(sequence.encode())
        require(evidence['sequence_sha256'] == sha and evidence['length'] == len(sequence)
                and evidence['cache_identity_hash'] == cache_hash
                and evidence['ranking_pointer'] == '/full_catalogue_ranking', 'Family segment reference differs')
        path = Path(evidence['unit_path']).resolve()
        require(path.is_relative_to(root.resolve()), 'Family segment reference escaped cache')
        key = (str(path), evidence['descriptor_sha256'])
        if key not in verified:
            unit = verify_unit(root, sha, {'path': str(path.relative_to(root)),
                'descriptor_sha256': evidence['descriptor_sha256']}, cache_hash)
            value = read(unit / 'record.json')
            require(value['sequence_sha256'] == sha and value['length'] == len(sequence)
                    and value['representation_identity_hash'] == digest(canonical(identity['model'])),
                    'Family cache record identity differs')
            validate_ranking(value['full_catalogue_ranking'], identity['catalogue_accessions'])
            verified[key] = value['full_catalogue_ranking']
        result.append(verified[key])
    return result


def matched_baseline(metadata_row, query, segment_eval, whole_eval, *, kind='exact'):
    a, b = segment_eval['query_metrics'][kind], whole_eval['query_metrics'][kind]
    require(segment_eval['query_id'] == whole_eval['query_id'] == query['query_id'], 'Paired baseline query differs')
    require(segment_eval['retrieval_budget'] == whole_eval['retrieval_budget'] == 'k total per original query',
            'Paired baseline retrieval budgets differ')
    require(a['scorable_components'] == b['scorable_components']
            and a['all_components_scorable'] == b['all_components_scorable']
            and [c[kind]['accepted_surviving'] for c in segment_eval['components']]
                == [c[kind]['accepted_surviving'] for c in whole_eval['components']],
            'Paired baseline surviving labels/denominators differ')
    full = a['all_components_scorable']
    return {**metadata_row, 'query_id': query['query_id'], 'group_id': query['group_id'], 'suite': query['suite'],
            'component_count': len(query['components']), 'scope': 'clean' if len(query['components']) == 1 else 'combined',
            'label_kind': kind, 'scorable_components': a['scorable_components'],
            'unscorable_components': a['unscorable_components'], 'all_components_scorable': full,
            'retrieval_budget': 'k total per original query', 'contrast': 'selected_segmenter_minus_same_method_whole',
            **{f'component_hit_difference_at_{k}': a[f'component_hits_at_{k}'] - b[f'component_hits_at_{k}'] for k in KS},
            **{f'all_hit_difference_at_{k}': a[f'all_components_hit_at_{k}'] - b[f'all_components_hit_at_{k}'] if full else None for k in KS},
            **{f'any_hit_difference_at_{k}': a[f'any_component_hit_at_{k}'] - b[f'any_component_hit_at_{k}'] if full else None for k in KS}}


def derive_family_baseline(identity, record, query, suite, labels, policy):
    """Mask the native primary arm with the same frozen all-component policy."""
    require(record['query'] == public_query(query), 'Family baseline query differs from frozen panel')
    arm = primary_arm(identity['mode'])
    original, pointer = source_arms(identity, record)[arm]
    require(isinstance(original, list), 'Family baseline must use the native single-list primary arm')
    validate_ranking(original, sorted(suite['sequences']))
    mask = mask_for_query(query, suite, labels, policy)
    ranking = filter_ranking(original, mask['included_accessions'])
    validate_ranking(ranking, mask['included_accessions'])
    evaluation = evaluate_query(ranking, query, suite['equivalence'], labels['acceptable_matches'])
    frozen = next(row for row in labels['removal_policies'][policy]['query_scoring_status']
                  if row['query_id'] == query['query_id'])
    scored = evaluation['query_metrics']['labelled']
    require(scored['scorable_components'] == frozen['scorable_components']
            and scored['all_components_scorable'] == frozen['all_components_scorable'],
            'Family baseline denominators differ from frozen panel')
    require(evaluation['query_metrics']['exact']['scorable_components'] == 0
            and evaluation['query_metrics']['duplicate_equivalent']['scorable_components'] == 0,
            'Family baseline mask left a query exact or duplicate reference')
    return {'arm': arm, 'ranking': ranking, 'evaluation': evaluation, 'common_mask': mask,
            'source_pointer': pointer,
            'regeneration': 'Read complete immutable native primary ranking; retain common-mask accessions, preserve score/accession order, renumber and recompute ties; evaluate frozen annotation labels.'}


def add_family_core_baselines(book, sources, queries, suite, labels, *, check=lambda: None):
    """Add all supplied core methods on the exact requested panel; return paired whole evaluations."""
    ids = [q['query_id'] for q in queries]
    paired = {}
    for source in sources:
        validate_panel(source.identity['query_ids'], ids)
        for policy in ('exact_removed', 'near90_removed'):
            arm = primary_arm(source.identity['mode'])
            meta = {**metadata(source, arm, policy=policy), 'role': 'family_retrieval_baseline',
                    'baseline_method_role': comparison_role(source.identity['mode'])}
            book.tables['source_coverage'].append({**meta, **validate_panel(source.identity['query_ids'], ids),
                'combined_queries': sum(len(q['components']) > 1 for q in queries),
                'clean_queries': sum(len(q['components']) == 1 for q in queries)})
            for query in queries:
                check()
                key = 'query-' + query['query_id']
                derived = derive_family_baseline(source.identity, source.record(key), query, suite, labels, policy)
                reference = {**source.reference(key, derived['source_pointer']),
                    'common_mask': derived['common_mask'], 'regeneration': derived['regeneration']}
                book.add(meta, query, derived['evaluation'], derived['ranking'], reference,
                         acceptable=labels['acceptable_matches'])
                if arm == 'whole':
                    paired[(source.identity_hash, query['query_id'], policy)] = derived['evaluation']
            expected = subset_denominators(labels, policy, ids)
            bundle = book.summaries[(source.root.name, arm, policy)]
            if 'combined' in bundle: denominator_check(bundle['combined'].result(), expected['combined'])
            for suite_name, denominator in expected['by_suite'].items():
                denominator_check(bundle['suite:' + suite_name].result(), denominator)
    return paired


def run(args):
    session = args.session.resolve()
    check_deadline(session)
    design_path = args.design or session / 'inputs/segmentation_design_v1.json'
    selection_path = args.selection or session / 'inputs/segmentation_selection_v1.json'
    session, suite, design, queries, configs = setup(SimpleNamespace(
        session=session, design=design_path, selection=selection_path, partition='test'))
    require(len(design['test_query_ids']) == 950 and len(design['clean_sanity_query_ids']) == 16,
            'Expected frozen 950 combined / 16 clean test panel')
    require(len(configs) == 2 and len({c['segmenter'] for c in configs}) == 2, 'Require both frozen segmenter winners')
    ids = [q['query_id'] for q in queries]
    selection = read(selection_path)
    dev_method = read(Path(selection['development_job']) / 'identity.json')['model']
    cores = [Snapshot(p) for p in args.core_job]
    segments = [Snapshot(p) for p in args.segment_job]
    spans = [Snapshot(p) for p in args.spans_job]
    families = [Snapshot(p) for p in args.family_job]
    family_spans = [Snapshot(p) for p in args.family_spans_job]
    sources = cores + segments + spans + families + family_spans
    require(len({s.identity_hash for s in sources}) == len(sources), 'Duplicate source job identity')
    models = [s for s in cores if s.identity['mode'] == 'model']
    require(len(models) == 4 and {s.identity['method']['model_id'] for s in models} == CORE_MODELS,
            'Require the four original core model representations')
    require(len(cores) == 6 and sum(s.identity['mode'] == 'alignment' for s in cores) == 1
            and sum(s.identity['mode'] == 'alignment_normalization' for s in cores) == 1,
            'Require four core models, raw SW, and separately labelled normalized SW')
    raw = next(s for s in cores if s.identity['mode'] == 'alignment')
    norm = next(s for s in cores if s.identity['mode'] == 'alignment_normalization')
    require(norm.identity['source_identity_hash'] == raw.identity_hash
            and norm.identity['source_completion_sha256'] == raw.completion_sha256,
            'Normalized alignment must derive from the supplied raw completion')
    book = ReportBook()
    for source in sources:
        validate_source_suite(source.identity, suite, suite['verification']['challenge_manifest_sha256'])
    for source in cores:
        coverage = validate_panel(source.identity['query_ids'], ids)
        book.tables['source_coverage'].append({**metadata(source, primary_arm(source.identity['mode'])), **coverage})
    parent_by_hash = {s.identity_hash: s for s in models}
    spans_by_hash = {s.identity_hash: s for s in spans}
    family_spans_by_hash = {s.identity_hash: s for s in family_spans}
    require(any(s.identity['model'] == dev_method for s in segments), 'Missing exact development-method segmentation comparison')
    baseline = {}
    for source in cores:
        for query in queries:
            check_deadline(session)
            key = 'query-' + query['query_id']
            record = source.record(key)
            require(record['query'] == public_query(query), 'Core query metadata differs')
            if source.identity['mode'] != 'alignment_normalization':
                book.tables['core_runtime'].append({**metadata(source, primary_arm(source.identity['mode'])),
                    'query_id': query['query_id'], 'suite': query['suite'], 'component_count': len(query['components']),
                    'scope': 'clean' if len(query['components']) == 1 else 'combined', 'query_unit_seconds': record['seconds'],
                    'runtime_interpretation': 'Saved native query unit time includes its controls; excludes catalogue/model loading. Observed cache-dependent execution, not a fresh speed benchmark.'})
            for arm, (ranking, pointer) in source_arms(source.identity, record).items():
                if isinstance(ranking, dict):
                    for rows in ranking.values(): validate_ranking(rows, sorted(suite['sequences']))
                    evaluation = evaluate_oracle_query(ranking, query, suite['equivalence'])
                else:
                    validate_ranking(ranking, sorted(suite['sequences']))
                    evaluation = evaluate_query(ranking, query, suite['equivalence'])
                expected_key = 'accession105' if source.identity['mode'] == 'alignment_normalization' else arm
                require(record['evaluations'][expected_key] == evaluation, 'Core saved/native reevaluation differs')
                book.add(metadata(source, arm), query, evaluation, ranking, source.reference(key, pointer))
                if arm == 'whole': baseline[(source.identity_hash, query['query_id'])] = evaluation
    used_spans = set()
    for source in segments:
        require(source.identity['mode'] == 'automatic_segment_retrieval_v1' and source.identity['partition'] == 'test',
                'Require completed held-out segment retrieval')
        require(source.identity['configs'] == configs and source.identity['selection_sha256'] == file_hash(selection_path)
                and source.identity['design_sha256'] == file_hash(design_path), 'Segmentation winners or freeze differ')
        parent = parent_by_hash.get(source.identity['parent_model_identity'])
        require(parent is not None, 'Segmentation parent omitted from supplied core jobs')
        validate_parent(source, parent)
        prediction = spans_by_hash.get(source.identity['predictions_identity'])
        require(prediction is not None and prediction.identity['mode'] == 'automatic_spans_v1'
                and prediction.completion_sha256 == source.identity['predictions_completion_sha256'],
                'Segmentation prediction snapshot omitted or changed')
        require(prediction.identity['partition'] == 'test' and prediction.identity['design_sha256'] == file_hash(design_path)
                and prediction.identity['selection_sha256'] == file_hash(selection_path), 'Prediction freeze differs')
        require(prediction.identity['alignment_identity'] == raw.identity_hash
                and prediction.identity['alignment_completion_sha256'] == raw.completion_sha256,
                'Primary prediction raw-alignment source differs')
        expected_predictions = {digest(canonical(prediction_config(c)))[:12]: prediction_config(c) for c in configs}
        require(prediction.identity['configs'] == expected_predictions, 'Prediction configurations differ from frozen winners')
        used_spans.add(prediction.identity_hash)
        coverage = validate_panel(source.identity['query_ids'], ids, exact=True)
        validate_panel(prediction.identity['query_ids'], ids, exact=True)
        book.tables['source_coverage'].append({**metadata(source, 'selected_segmentation', config=configs[0]), **coverage})
        for query in queries:
            check_deadline(session)
            for cfg in configs:
                key = query['query_id'] + '--' + cfg['configuration_id']
                prediction_key = query['query_id'] + '--' + digest(canonical(prediction_config(cfg)))[:12]
                proposal, record = prediction.record(prediction_key), source.record(key)
                validate_proposal(proposal, query, prediction_config(cfg))
                validate_retrieval_record(record, query, cfg, proposal, sorted(suite['sequences']), suite['equivalence'])
                require(record['prediction_unit'] == prediction_key, 'Retrieval prediction unit differs')
                meta = metadata(source, cfg['configuration_id'], config=cfg)
                ref = {**source.reference(key, '/ranking'), 'prediction_source': prediction.reference(prediction_key),
                       'parent_source': parent.reference('query-' + query['query_id'], '/rankings/whole')}
                book.add_segments(meta, query, record, proposal, record['segment_rankings'], ref)
                book.tables['matched_whole_baseline'].append({**matched_baseline(
                    meta, query, record['evaluation'], baseline[(parent.identity_hash, query['query_id'])]),
                    'baseline_job_id': parent.root.name, 'baseline_identity_hash': parent.identity_hash})
    require(used_spans == set(spans_by_hash), 'Unused or missing primary prediction snapshot')
    family_evidence = None
    if families:
        family_design_path = args.family_design or session / 'inputs/family_segmentation_design_v1.json'
        family_design = read(family_design_path)
        labels_path = session / 'inputs' / family_design['family_labels_filename']
        labels = read(labels_path)
        audit = validate_family_freeze(labels, suite, labels_path, session / 'audit/independent_family_label_review.json')
        require(family_design['status'] == 'prospective_frozen_before_transfer_scoring'
                and family_design['retuning_allowed'] is False
                and family_design['primary_design_sha256'] == file_hash(design_path)
                and family_design['family_labels_sha256'] == file_hash(labels_path), 'Family transfer freeze differs')
        family_evidence = {'design_sha256': file_hash(family_design_path), 'labels_sha256': file_hash(labels_path), 'label_audit': audit}
        # Whole-model/alignment baselines use identical masks and test queries.
        # Only the transferred segmenter is restricted to the exact dev 8M method.
        family_baselines = add_family_core_baselines(book, cores, queries, suite, labels,
                                                    check=lambda: check_deadline(session))
        verified_cache, used_family_spans = {}, set()
        for source in families:
            require(source.identity['mode'] == 'family_automatic_retrieval_v1'
                    and source.identity['configs'] == configs, 'Family configurations changed')
            parent = parent_by_hash.get(source.identity['parent_model_identity'])
            require(parent is not None, 'Family parent omitted')
            validate_parent(source, parent, development_method=dev_method)
            validate_parent_development_method(parent.identity['method'], selection_path)
            prediction = family_spans_by_hash.get(source.identity['prediction_identity'])
            require(prediction is not None and prediction.identity['mode'] == 'family_automatic_spans_v1'
                    and prediction.completion_sha256 == source.identity['prediction_completion_sha256'], 'Family prediction snapshot differs')
            require(prediction.identity['raw_alignment_identity'] == raw.identity_hash
                    and prediction.identity['raw_alignment_completion_sha256'] == raw.completion_sha256,
                    'Family prediction raw-alignment source differs')
            family_ids = list(design['test_query_ids'])
            if source.identity['include_clean_sanity']: family_ids += design['clean_sanity_query_ids']
            family_queries = [q for q in queries if q['query_id'] in set(family_ids)]
            for obj in (source, prediction):
                require(obj.identity['selection_sha256'] == file_hash(selection_path)
                        and obj.identity['transfer_design_sha256'] == file_hash(family_design_path)
                        and obj.identity['primary_design_sha256'] == file_hash(design_path)
                        and obj.identity['family_labels_sha256'] == file_hash(labels_path), 'Family source freeze differs')
                require(obj.identity['include_clean_sanity'] == source.identity['include_clean_sanity'],
                        'Family prediction/retrieval clean coverage differs')
                validate_panel(obj.identity['query_ids'], family_ids, exact=True)
            used_family_spans.add(prediction.identity_hash)
            book.tables['source_coverage'].append({**metadata(source, 'family_segmentation', config=configs[0], policy='family_panel'),
                **validate_panel(source.identity['query_ids'], family_ids, exact=True),
                'combined_queries': 950, 'clean_queries': len(family_ids)-950})
            for query in family_queries:
                check_deadline(session)
                for policy in source.identity['policies']:
                    mask = mask_for_query(query, suite, labels, policy)
                    for cfg in configs:
                        key = unit_key(query['query_id'], policy, cfg['configuration_id'])
                        proposal, record = prediction.record(key), source.record(key)
                        require(record['prediction_unit_descriptor_sha256'] == prediction.index[key]['descriptor_sha256'],
                                'Family prediction unit descriptor differs')
                        rankings = cached_family_rankings(record, query, source.identity, session, verified_cache)
                        validate_transfer_record(record, query, cfg, proposal, mask, labels, suite['equivalence'], rankings)
                        masked = [filter_ranking(r, mask['included_accessions']) for r in rankings]
                        ref = {**source.reference(key, '/ranking'), 'prediction_source': prediction.reference(key),
                               'common_mask': mask, 'segment_rankings': record['segment_ranking_references'],
                               'regeneration': 'Filter/rerank each complete cached segment list by common mask BEFORE frozen fusion.'}
                        meta = metadata(source, cfg['configuration_id'], config=cfg, policy=policy)
                        book.add_segments(meta,
                                          query, record, proposal, masked, ref, acceptable=labels['acceptable_matches'])
                        book.tables['matched_whole_baseline'].append({**matched_baseline(meta, query,
                            record['evaluation'], family_baselines[(parent.identity_hash, query['query_id'], policy)],
                            kind='labelled'), 'baseline_job_id': parent.root.name,
                            'baseline_identity_hash': parent.identity_hash})
            for policy in source.identity['policies']:
                expected = subset_denominators(labels, policy, family_ids)
                for cfg in configs:
                    bundle = book.summaries[(source.root.name, cfg['configuration_id'], policy)]
                    denominator_check(bundle['combined'].result(), expected['combined'])
                    for suite_name, denominator in expected['by_suite'].items():
                        denominator_check(bundle['suite:' + suite_name].result(), denominator)
        require(used_family_spans == set(family_spans_by_hash), 'Unused or missing family prediction snapshot')
    else:
        require(not family_spans, 'Family prediction source supplied without retrieval output')
    report = book.result()
    if families:
        baseline_headlines = [r for r in report['family_tables']['headlines'] if r['role'] == 'family_retrieval_baseline']
        require(len(baseline_headlines) == 2*len(cores)
                and {(r['job_id'], r['policy']) for r in baseline_headlines}
                    == {(s.root.name, p) for s in cores for p in ('exact_removed', 'near90_removed')},
                'Family report lacks a core baseline/policy headline')
    report.update(frozen_panel={'combined_query_ids': design['test_query_ids'], 'clean_query_ids': design['clean_sanity_query_ids']},
                  selection_evidence={'selection_sha256': file_hash(selection_path), 'selection': selection,
                                      'exact_development_method': dev_method, 'selected_configurations': configs},
                  family_evidence=family_evidence, sources=[s.provenance() for s in sources],
                  validation_scope='Saved-result regeneration and report consistency; not an independent implementation audit.',
                  family_interpretation='Exploratory transfer of unchanged original-dev winners. Acceptable source-annotation family membership is distinct from whole-protein function; unavailable surviving labels are unscorable, and any/all query recovery requires every component to be scorable.',
                  no_recalibration=True, cross_model_scores_compared=False)
    identity = {'mode': 'heldout_segment_comparison_v1', 'sources': report['sources'],
                'selection_sha256': file_hash(selection_path), 'design_sha256': file_hash(design_path),
                'family_evidence': family_evidence, 'query_ids': ids,
                'code_hashes': {str(p.relative_to(ROOT)): file_hash(p) for p in [Path(__file__),
                    ROOT/'scripts/session_segment_experiment.py', ROOT/'scripts/session_family_segmentation.py',
                    ROOT/'scripts/session_catalogue_analysis.py', ROOT/'scripts/session_compare_models.py',
                    ROOT/'scripts/session_independent_audit.py', ROOT/'src/iri_constructs/session_store.py',
                    ROOT/'src/iri_constructs/session_benchmark.py', ROOT/'src/iri_constructs/session_segmentation.py']}}
    store = UnitStore(session / 'analysis' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    existing = store.read('report')
    if existing is None: store.publish('report', report)
    else: require(existing == report, 'Regenerated report differs from immutable result')
    artifacts = {}
    for namespace in ('tables', 'family_tables'):
        for name, rows in report[namespace].items():
            relative = namespace + '/' + name + '.csv'
            path = store.root / relative
            data = csv_bytes(rows)
            if path.exists(): require(path.read_bytes() == data, 'Existing table bytes differ')
            else: atomic_bytes(path, data)
            artifacts[relative] = {'sha256': file_hash(path), 'bytes': path.stat().st_size}
    exports = {'artifacts': artifacts, 'report_unit_descriptor_sha256': file_hash(store.find('report') / 'complete.json')}
    existing = store.read('exports')
    if existing is None: store.publish('exports', exports)
    else: require(existing == exports, 'Existing export evidence differs')
    store.complete(['report', 'exports'], {'query_count': len(ids), 'combined_queries': 950, 'clean_queries': 16})
    print(json.dumps({'status': 'complete', 'job': str(store.root), 'source_count': len(sources)}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--core-job', type=Path, action='append', required=True)
    parser.add_argument('--segment-job', type=Path, action='append', required=True)
    parser.add_argument('--spans-job', type=Path, action='append', required=True)
    parser.add_argument('--family-job', type=Path, action='append', default=[])
    parser.add_argument('--family-spans-job', type=Path, action='append', default=[])
    parser.add_argument('--design', type=Path)
    parser.add_argument('--selection', type=Path)
    parser.add_argument('--family-design', type=Path)
    parser.add_argument('--name', default='heldout-segment-comparison-v1')
    run(parser.parse_args())


if __name__ == '__main__':
    main()
