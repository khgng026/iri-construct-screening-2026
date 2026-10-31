#!/usr/bin/env python3
"""Freeze source-family labels, then filter complete expanded rankings only."""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'src'), str(ROOT / 'scripts')]
from iri_constructs.session_store import UnitStore, atomic_json, canonical, digest, file_hash
from iri_constructs.session_expansion import load_expansion
from iri_constructs.session_benchmark import load_suite, filter_ranking
from session_context_profile_audit import completed_job
from session_independent_audit import validate_ranking
from session_experiment import public_query, check_deadline

SIZES = ('101', '250', '500', '1000')
KS = (1, 5, 10)
MODELS = {'facebook/esm2_t6_8M_UR50D', 'facebook/esm2_t30_150M_UR50D',
          'facebook/esm2_t33_650M_UR50D', 'biohub/ESMC-300M', 'biohub/ESMC-600M'}
ARMS = ('local_alignment_raw', 'local_alignment_reference_selfscore')
MANIFEST_SHA = 'f9cc5513f130f2c1f6b5907e0104f4967a3549f9007eba3424ac8a05b90366f2'
LEGACY_LABEL_SHA = '92a2bb852426948fb627f4029e46af7c97dcdbafbf7ad579f57fa5e7a5155d94'
RETRIEVAL_SHA = '8cf5f2a2184d24465324a589afb5016cddd1421e78055b63c9759f3a08fd1cfc'
ORIGINAL_FAMILY_DESIGN_SHA = '4ae9b36df27d5616b5750fc138a754129df802c70e368dbcf20e7212deabdde7'
FAMILY_LABEL_SHA = '5868ebcdd2e42570679d06b8e5cb06a8d89906d91d6840952a36d4f1bf409883'


def require(value, message):
    if not value: raise ValueError(message)


def read(path): return json.loads(Path(path).read_text())


def source_code():
    names = ('scripts/session_expanded_family.py', 'scripts/session_context_profile_audit.py',
             'scripts/session_independent_audit.py', 'scripts/session_experiment.py',
             'src/iri_constructs/session_expansion.py', 'src/iri_constructs/session_store.py',
             'src/iri_constructs/session_benchmark.py', 'src/iri_constructs/t05_data.py',
             'src/iri_constructs/t05_composition.py')
    return {name: file_hash(ROOT / name) for name in names}


def first_family(text):
    first = text.split('.', 1)[0].strip()
    return first if first.startswith('Belongs to the ') and first.endswith(' family') else None


def tsv_family(text):
    require(not text or text.startswith('SIMILARITY:'), 'Unexpected SIMILARITY transport prefix')
    statements = {value for comment in text.split('SIMILARITY:')[1:]
                  if (value := first_family(comment)) is not None}
    require(len(statements) <= 1, 'Ambiguous multiple first-family annotations')
    return next(iter(statements), None)


def pointer(document, path):
    for part in path.strip('/').split('/'):
        part = part.replace('~1', '/').replace('~0', '~')
        document = document[int(part)] if isinstance(document, list) else document[part]
    return document


def annotations(session, expansion):
    legacy_path = session / 'inputs/family_labels_v1.json'
    require(file_hash(legacy_path) == LEGACY_LABEL_SHA, 'Original family label snapshot changed')
    legacy = read(legacy_path); output = {}
    for accession in expansion['subsets']['1000']:
        record = expansion['manifest']['catalogue'][accession]
        if record['origin'] == 'original_reviewed_105':
            alias_statements, evidence = [], []
            for alias in record['aliases']:
                original = legacy['source_evidence'][alias]
                raw_path = ROOT / original['raw_path']
                require(file_hash(raw_path) == original['raw_sha256'], 'Original annotation source changed')
                raw = read(raw_path); statements = set()
                for entry in original['annotation_evidence']:
                    node = pointer(raw, entry['json_pointer'])
                    require(node['value'] == entry['full_statement'], 'Original annotation pointer differs')
                    statement = first_family(node['value'])
                    require(statement == entry['family_statement'] and statement is not None, 'Original first-family rule differs')
                    statements.add(statement)
                require(len(statements) <= 1, 'Original representative has ambiguous families')
                alias_statements.append(statements); evidence.append(original)
            require(all(x == alias_statements[0] for x in alias_statements), 'Exact aliases have conflicting family labels')
            statement = next(iter(alias_statements[0]), None)
            source = {'kind': 'original_json', 'aliases': evidence}
        else:
            metadata = record['metadata']; statement = tsv_family(metadata['cc_similarity'])
            source = {'kind': 'reviewed_expansion_tsv', 'full_similarity_text': metadata['cc_similarity'],
                      'manifest_json_pointer': '/catalogue/' + accession + '/metadata/cc_similarity',
                      'metadata_row_sha256': record['metadata_row_sha256'],
                      'source_page': record['review']['source_page'], 'source_page_sha256': record['review']['source_page_sha256'],
                      'source_url': record['review']['source_url'], 'source_response': record['source_response'],
                      'source_response_sha256': record['source_response_sha256'],
                      'entry_version': metadata['version'], 'sequence_version': metadata['sequence_version']}
        output[accession] = {'family_statement': statement, 'label_status': 'source_annotated' if statement else 'unlabelled',
                             'sequence_sha256': record['sequence_sha256'], 'aliases': record['aliases'], 'source': source}
    return output


def query_mask(query, subset, target_aliases, sequence_hashes):
    hashes = {sequence_hashes[target_aliases[a]] for a in query['components']}
    removed = sorted(a for a in subset if sequence_hashes[a] in hashes)
    included = sorted(set(subset) - set(removed))
    require(included and all(target_aliases[a] in removed for a in query['components']), 'Exact sequence removal failed')
    return {'removed_representatives': removed, 'remaining_count': len(included),
            'included_representatives_sha256': digest(canonical(included))}


def denominator(rows, stratum):
    total = sum(len(row['components']) for row in rows)
    flags = [[c['native_scorable'] if stratum == 'native' else c['common_scorable'] for c in row['components']] for row in rows]
    n = sum(sum(values) for values in flags); fully = sum(all(values) for values in flags)
    return {'original_queries': len(rows), 'original_components': total, 'scorable_components': n,
            'unscorable_components': total - n, 'all_components_scorable_queries': fully,
            'excluded_unscorable_queries': len(rows) - fully}


def build_labels(session, expansion, suite):
    annotation = annotations(session, expansion)
    aliases = expansion['target_aliases']; hashes = {a: row['sequence_sha256'] for a, row in annotation.items()}
    targets = sorted(suite['manifest']['targets'])
    require(all(annotation[aliases[t]]['family_statement'] for t in targets), 'Target lacks its frozen explicit family')
    accepted = {size: {t: sorted(a for a in expansion['subsets'][size]
                         if annotation[a]['family_statement'] == annotation[aliases[t]]['family_statement']) for t in targets}
                for size in SIZES}
    rows = {}
    for query in suite['queries']:
        by_prefix = {}
        for size in SIZES:
            mask = query_mask(query, expansion['subsets'][size], aliases, hashes)
            removed = set(mask['removed_representatives'])
            components = [{'accession': a, 'accepted_survivors': sorted(set(accepted[size][a]) - removed)} for a in query['components']]
            by_prefix[size] = {'mask': mask, 'components': components}
        common = {a for a in query['components'] if all(next(c for c in by_prefix[size]['components'] if c['accession'] == a)['accepted_survivors'] for size in SIZES)}
        for row in by_prefix.values():
            for c in row['components']:
                c.update(native_scorable=bool(c['accepted_survivors']), common_scorable=c['accession'] in common)
        rows[query['query_id']] = {'query': public_query(query), 'prefixes': by_prefix}
    primary = read(session / 'inputs/segmentation_design_v1.json')
    cohorts = {'clean16': sorted(q['query_id'] for q in suite['queries'] if len(q['components']) == 1),
               'combined1130': sorted(q['query_id'] for q in suite['queries'] if len(q['components']) > 1),
               'heldout_combined950': sorted(primary['test_query_ids'])}
    require({k: len(v) for k, v in cohorts.items()} == {'clean16': 16, 'combined1130': 1130, 'heldout_combined950': 950}, 'Frozen cohort counts differ')
    require(set(cohorts['heldout_combined950']) <= set(cohorts['combined1130']), 'Heldout cohort is not a combined subset')
    counts = {size: {stratum: {cohort: denominator([rows[q]['prefixes'][size] for q in ids], stratum)
                             for cohort, ids in cohorts.items()} for stratum in ('native', 'common11')} for size in SIZES}
    common_targets = sorted({c['accession'] for row in rows.values() for c in row['prefixes']['101']['components'] if c['common_scorable']})
    require(len(common_targets) == 11, 'Expected the independently reviewed common 11-target stratum')
    return {'version': 'expanded-source-family-labels-v1', 'status': 'frozen_before_expanded_family_outcome_evaluation',
            'expansion_manifest_sha256': expansion['manifest_sha256'], 'legacy_family_labels_sha256': LEGACY_LABEL_SHA,
            'challenge_sha256': suite['verification']['challenge_manifest_sha256'],
            'heldout_design_sha256': file_hash(session / 'inputs/segmentation_design_v1.json'),
            'rule': 'Exact first literal-period family sentence of original SIMILARITY comment; preserve case and spelling; TSV marker stripped; subfamily text is evidence only',
            'missing_labels': 'Unlabelled is not a verified biological negative; only source-family retrieval is evaluated',
            'annotations': annotation, 'accepted_representatives_by_prefix': accepted,
            'subsets': expansion['subsets'], 'target_aliases': aliases, 'common_scorable_targets': common_targets,
            'cohorts': cohorts, 'query_masks_and_survivors': rows, 'denominators': counts,
            'strata': {'native': 'All original components; scorable if at least one family-labelled representative survives',
                       'common11': 'All original counts retained; only components scorable at every prefix are eligible; any/all query recovery requires every original component eligible'},
            'budget': 'k total unique representatives per original query', 'budgets': list(KS), 'near90': 'not included'}


def parent_metadata(path):
    path = Path(path).resolve(); identity = read(path / 'identity.json'); identity_hash = digest(canonical(identity))
    candidates = [(read(p)['completed_utc'], p) for p in path.glob('completion-*.json')]
    require(candidates, 'Core parent is incomplete')
    completion_path = max(candidates)[1]; completion = read(completion_path)
    ids = identity['query_ids']; require(len(ids) == len(set(ids)) == 1146, 'Parent query panel is incomplete')
    expected = {'query-' + q for q in ids} | {'summary'}
    if identity['mode'] == 'model': expected |= {'catalogue-' + a for a in identity['catalogue_accessions']}
    else: require(identity['mode'] == 'alignment', 'Unsupported core parent')
    require(completion['identity_hash'] == identity_hash and completion['expected_count'] == len(expected) and
            set(completion['units']) == expected, 'Parent completion identity/count differs')
    return {'path': str(path), 'identity_sha256': identity_hash, 'completion_sha256': file_hash(completion_path),
            'identity': identity}


def freeze(args):
    session = args.session.resolve(); suite = load_suite(ROOT, session / 'inputs')
    manifest_path = session / 'expansion/milestone-1000/admitted/expanded_catalogue_v1.json'
    require(file_hash(manifest_path) == MANIFEST_SHA, 'Wrong frozen 1000 catalogue')
    expansion = load_expansion(ROOT, session, manifest_path); labels = build_labels(session, expansion, suite)
    parents = [parent_metadata(p) for p in args.parent_job]
    require(len(parents) == 5 and {p['identity']['method']['model_id'] for p in parents} == MODELS, 'Require the five fixed core model parents')
    roster = {}
    for parent in parents:
        identity = parent['identity']; method = identity['method']
        require(identity['mode'] == 'model' and method['device'] == 'cpu' and method['dtype'] == 'float32' and method['layer'] == method['num_layers'], 'Wrong parent representation')
        require(identity['challenge_sha256'] == suite['verification']['challenge_manifest_sha256'] and
                set(identity['query_ids']) == {q['query_id'] for q in suite['queries']} and
                identity['catalogue_sequence_sha256'] == {a: digest(seq.encode()) for a, seq in suite['sequences'].items()}, 'Parent canonical suite differs')
        roster[method['model_id']] = {k: v for k, v in parent.items() if k != 'identity'} | {'method': method}
    raw = parent_metadata(session / 'jobs/alignment-full-cc35ef35e5cc')
    require(raw['identity']['mode'] == 'alignment', 'Wrong original alignment source')
    labels_path = session / 'inputs/expanded_family_labels_v1.json'
    require(not labels_path.exists(), 'Preserve the already frozen expanded labels')
    atomic_json(labels_path, labels)
    design = {'version': 'expanded-exact-removal-family-v1', 'status': 'frozen_before_expanded_family_outcome_evaluation',
              'basis': 'Exploratory extension after original core/family and expanded-exact observations; metadata-only label inspection precedes expanded-family outcomes',
              'manifest_path': str(manifest_path.relative_to(session)), 'manifest_sha256': MANIFEST_SHA,
              'labels_path': str(labels_path.relative_to(session)), 'labels_sha256': file_hash(labels_path),
              'retrieval_design_sha256': RETRIEVAL_SHA, 'normalization_design_sha256': file_hash(session / 'audit/alignment_normalization_design_v1.json'),
              'challenge_sha256': suite['verification']['challenge_manifest_sha256'], 'prefixes': list(SIZES),
              'model_roster': roster, 'model_arms': ['whole'], 'alignment_arms': list(ARMS),
              'raw_alignment_source': {k: v for k, v in raw.items() if k != 'identity'} | {'method': raw['identity']['method']},
              'ranking_implementation_sha256': {name: file_hash(ROOT / name) for name in ('scripts/session_expanded_model.py', 'scripts/session_expanded_alignment.py')},
              'cohorts': labels['cohorts'], 'strata': labels['strata'], 'budgets': list(KS), 'headline_k': 5,
              'removal_policy': 'One common union removes canonical exact-sequence equivalents of every original component before ranking evaluation',
              'near90': 'not included', 'code_hashes': source_code(),
              'interpretation': 'Source-annotation family recovery only; catalogue size/composition and label coverage change; no functional or pretrained-holdout claims'}
    design_path = session / 'inputs/expanded_family_design_v1.json'
    require(not design_path.exists(), 'Preserve the already frozen expanded family design')
    atomic_json(design_path, design)
    print(json.dumps({'status': 'labels_and_design_frozen_without_outcomes', 'design_sha256': file_hash(design_path),
                      'labels_sha256': file_hash(labels_path), 'runner_sha256': file_hash(__file__)}), flush=True)


def freeze_execution(args):
    """Pin stricter source-admission code without revising scientific artifacts."""
    session = args.session.resolve(); design_path = session / 'inputs/expanded_family_design_v1.json'
    require(file_hash(design_path) == ORIGINAL_FAMILY_DESIGN_SHA, 'Preserve the original scientific design')
    design = read(design_path)
    require(file_hash(session / design['labels_path']) == design['labels_sha256'] == FAMILY_LABEL_SHA, 'Preserve original frozen labels')
    manifest = read(session / design['manifest_path'])
    require(file_hash(session / design['manifest_path']) == MANIFEST_SHA, 'Expansion manifest changed')
    labels = read(session / design['labels_path'])
    new_order = [a for a in manifest['new_representatives_hash_order'] if a in set(labels['subsets']['1000']) - set(labels['subsets']['101'])]
    require(len(new_order) == len(set(new_order)) == 899 and set(new_order) == set(labels['subsets']['1000']) - set(labels['subsets']['101']), 'Frozen new-reference order differs')
    model_source = read(session / 'jobs/expanded1000-8m-rank-v1-f61a400723f8/identity.json')
    model_code = model_source['code']
    require(model_source['expansion_manifest_sha256'] == MANIFEST_SHA and
            model_code['scripts/session_expanded_model.py'] == design['ranking_implementation_sha256']['scripts/session_expanded_model.py'],
            'Existing accepted model producer differs')
    require(all(file_hash(ROOT / name) == sha for name, sha in model_code.items()), 'A model scientific dependency changed since accepted expanded8M scoring')
    alignment_names = ('scripts/session_expanded_alignment.py', 'src/iri_constructs/session_expansion.py',
                       'src/iri_constructs/session_alignment.py', 'src/iri_constructs/session_store.py',
                       'src/iri_constructs/session_benchmark.py', 'scripts/session_alignment_analysis.py',
                       'scripts/session_catalogue_analysis.py', 'scripts/session_independent_audit.py')
    alignment_code = {name: file_hash(ROOT / name) for name in alignment_names}
    require(alignment_code['scripts/session_expanded_alignment.py'] == design['ranking_implementation_sha256']['scripts/session_expanded_alignment.py'], 'Original alignment wrapper changed')
    raw_identity = read(Path(design['raw_alignment_source']['path']) / 'identity.json')
    require(digest(canonical(raw_identity)) == design['raw_alignment_source']['identity_sha256'] and
            alignment_code['src/iri_constructs/session_alignment.py'] == raw_identity['code_hashes']['src/iri_constructs/session_alignment.py'],
            'Original independently audited alignment scoring implementation changed')
    consumer = source_code(); old_consumer = design['code_hashes']
    require(consumer == {**old_consumer, 'scripts/session_expanded_family.py': file_hash(__file__)}, 'Unrelated consumer dependency changed')
    record = {'version': 'expanded-family-source-admission-v2', 'status': 'frozen_before_expanded_family_outcome_evaluation',
              'scientific_design_sha256': ORIGINAL_FAMILY_DESIGN_SHA, 'labels_sha256': FAMILY_LABEL_SHA,
              'original_runner_sha256': old_consumer['scripts/session_expanded_family.py'], 'runner_sha256': file_hash(__file__),
              'consumer_code_hashes': consumer, 'expanded_model_code_hashes': model_code,
              'expanded_alignment_code_hashes': alignment_code, 'new_reference_order': new_order,
              'manifest_sha256': MANIFEST_SHA, 'source_model_code_basis_identity_sha256': digest(canonical(model_source)),
              'change': 'Source admission only: exact899-reference partition and canonical hashes; every scientific dependency pinned; unchanged labels/cohorts/masks/budgets/methods'}
    output = session / 'inputs/expanded_family_execution_v2.json'
    require(not output.exists(), 'Preserve existing execution freeze')
    atomic_json(output, record)
    print(json.dumps({'status': 'source_admission_freeze_without_outcomes', 'execution_sha256': file_hash(output),
                      'runner_sha256': file_hash(__file__), 'design_sha256': ORIGINAL_FAMILY_DESIGN_SHA, 'labels_sha256': FAMILY_LABEL_SHA}), flush=True)


def evaluate_filtered(query, ranking, label_row, stratum):
    positions = {row['accession']: row for row in ranking}; components = []
    for label in label_row['components']:
        surviving = label['accepted_survivors']
        require(set(surviving) <= set(positions), 'Frozen surviving label missing from complete ranking')
        eligible = label['native_scorable'] if stratum == 'native' else label['common_scorable']
        best = min((positions[a] for a in surviving), key=lambda r: r['rank']) if eligible else None
        components.append({'accession': label['accession'], 'scorable': bool(eligible), 'accepted_survivors': surviving,
                           'best_accession': best['accession'] if best else None, 'rank': best['rank'] if best else None,
                           'hit': {str(k): int(best['rank'] <= k) if best else None for k in KS}})
    require([c['accession'] for c in components] == query['components'], 'Frozen component order differs')
    fully = all(c['scorable'] for c in components)
    return {'components': components, 'all_components_scorable': fully,
            'query_hits': {str(k): {'any': int(any(c['hit'][str(k)] for c in components)) if fully else None,
                                   'all': int(all(c['hit'][str(k)] for c in components)) if fully else None} for k in KS}}


class Accumulator:
    def __init__(self):
        self.queries = self.components = self.scorable = self.fully = 0
        self.hits = {k: 0 for k in KS}; self.any = {k: 0 for k in KS}; self.all = {k: 0 for k in KS}
        self.targets = defaultdict(lambda: {'n': 0, 'hits': {k: 0 for k in KS}})
        self.groups = defaultdict(lambda: {'n': 0, 'hits': {k: 0 for k in KS}})

    def add(self, query, evaluation):
        self.queries += 1; self.components += len(evaluation['components']); self.fully += evaluation['all_components_scorable']
        for c in evaluation['components']:
            if not c['scorable']: continue
            self.scorable += 1; target = self.targets[c['accession']]; group = self.groups[query['group_id']]
            target['n'] += 1; group['n'] += 1
            for k in KS:
                hit = c['hit'][str(k)]; self.hits[k] += hit; target['hits'][k] += hit; group['hits'][k] += hit
        if evaluation['all_components_scorable']:
            for k in KS:
                self.any[k] += evaluation['query_hits'][str(k)]['any']; self.all[k] += evaluation['query_hits'][str(k)]['all']

    def result(self, expected):
        counts = {'original_queries': self.queries, 'original_components': self.components, 'scorable_components': self.scorable,
                  'unscorable_components': self.components - self.scorable, 'all_components_scorable_queries': self.fully,
                  'excluded_unscorable_queries': self.queries - self.fully}
        require(counts == expected, 'Aggregated family denominators differ from frozen metadata')
        return {**counts, 'scorable_target_count': len(self.targets), 'scorable_group_count': len(self.groups),
                'by_k': {str(k): {'component_hits': self.hits[k], 'component_recall': self.hits[k] / self.scorable if self.scorable else None,
                                 'per_target_macro_recall': math.fsum(v['hits'][k]/v['n'] for v in self.targets.values()) / len(self.targets) if self.targets else None,
                                 'group_macro_component_recall': math.fsum(v['hits'][k]/v['n'] for v in self.groups.values()) / len(self.groups) if self.groups else None,
                                 'any_component_hit_count': self.any[k], 'all_components_hit_count': self.all[k],
                                 'any_component_hit': self.any[k] / self.fully if self.fully else None,
                                 'all_components_hit': self.all[k] / self.fully if self.fully else None} for k in KS}}


def validate_alignment_partition(identity, reference_order, reference_hashes):
    require(reference_order and len(reference_order) == len(set(reference_order)), 'Frozen new-reference order is empty or duplicated')
    shards = identity['new_reference_shards']
    require(isinstance(shards, list) and shards and all(isinstance(shard, list) and 1 <= len(shard) <= 32 for shard in shards),
            'Alignment shard sizes must be nonempty and bounded1..32')
    flattened = [accession for shard in shards for accession in shard]
    require(len(flattened) == len(set(flattened)) and flattened == reference_order,
            'Alignment shards do not reproduce the exact frozen new-reference order/coverage')
    require(identity['reference_sequence_sha256'] == reference_hashes, 'Canonical alignment reference sequence hashes differ')


def validate_source_code(identity, execution):
    mode = identity['mode']
    require(mode in ('expanded_model_retrieval_v1', 'expanded_alignment_retrieval_v1'), 'Unknown expanded source mode')
    key = 'expanded_model_code_hashes' if mode == 'expanded_model_retrieval_v1' else 'expanded_alignment_code_hashes'
    require(identity['code'] == execution[key], 'Expanded source scientific dependency dictionary differs')


def expected_units(identity, reference_order=None, reference_hashes=None):
    ids = identity['query_ids']; require(len(ids) == len(set(ids)) == 1146, 'Expanded source must cover all1146 queries')
    expected = {'query-' + q for q in ids} | {'summary'}
    if identity['mode'] == 'expanded_alignment_retrieval_v1':
        require(reference_order is not None and reference_hashes is not None, 'Alignment evidence needs independently frozen reference expectations')
        validate_alignment_partition(identity, reference_order, reference_hashes)
        expected |= {'reference-selfscores'} | {f'raw-{q}-{i:03d}' for q in ids for i in range(len(identity['new_reference_shards']))}
    else: require(identity['mode'] == 'expanded_model_retrieval_v1', 'Unknown expanded source mode')
    return expected


def open_sources(args, session, design, labels, suite, execution):
    require(len(args.model_job) == 5 and args.alignment_job is not None, 'Require all five model sources and expanded alignment')
    sources = []; seen = set(); ids = {q['query_id'] for q in suite['queries']}
    for path in [*args.model_job, args.alignment_job]:
        path = path.resolve(); require(path.is_relative_to(session / 'jobs'), 'Expanded source outside owned jobs')
        identity = read(path / 'identity.json')
        require(set(identity['query_ids']) == ids and identity['expansion_manifest_sha256'] == design['manifest_sha256'] and
                identity['retrieval_design_sha256'] == design['retrieval_design_sha256'] and identity['subsets'] == labels['subsets'] and
                identity['target_aliases'] == labels['target_aliases'], 'Expanded source panel/catalogue/design differs')
        validate_source_code(identity, execution)
        expected = expected_units(identity, execution['new_reference_order'],
                                  {a: labels['annotations'][a]['sequence_sha256'] for a in labels['subsets']['1000']})
        source = completed_job(path, expected)
        require(source['identity_hash'] not in seen, 'Duplicate expanded source'); seen.add(source['identity_hash'])
        if identity['mode'] == 'expanded_model_retrieval_v1':
            require(identity['code']['scripts/session_expanded_model.py'] == design['ranking_implementation_sha256']['scripts/session_expanded_model.py'], 'Expanded model ranking implementation differs')
            matches = [(model, parent) for model, parent in design['model_roster'].items() if identity['parent_identity'] == parent['identity_sha256']]
            require(len(matches) == 1, 'Expanded model does not match frozen roster')
            model, parent = matches[0]
            require(identity['parent_completion_sha256'] == parent['completion_sha256'], 'Expanded parent completion changed')
            require(digest(canonical(read(Path(parent['path']) / 'identity.json'))) == parent['identity_sha256'], 'Parent identity bytes changed')
            require(any(file_hash(p) == parent['completion_sha256'] for p in Path(parent['path']).glob('completion-*.json')), 'Pinned parent completion missing')
            source.update(label=model, arms=['whole'])
        else:
            raw = design['raw_alignment_source']
            require(identity['code']['scripts/session_expanded_alignment.py'] == design['ranking_implementation_sha256']['scripts/session_expanded_alignment.py'], 'Expanded alignment implementation differs')
            require(identity['source_raw_identity'] == raw['identity_sha256'] and identity['source_raw_completion_sha256'] == raw['completion_sha256'] and
                    identity['normalization_design_sha256'] == design['normalization_design_sha256'] and identity['method'] == raw['method'] and
                    identity['challenge_sha256'] == design['challenge_sha256'] and
                    identity['query_sequence_sha256'] == {q['query_id']: digest(q['sequence'].encode()) for q in suite['queries']}, 'Expanded alignment method/source changed')
            source.update(label='Smith-Waterman', arms=list(ARMS))
        sources.append(source)
    require({s['label'] for s in sources if s['identity']['mode'] == 'expanded_model_retrieval_v1'} == MODELS and
            sum(s['identity']['mode'] == 'expanded_alignment_retrieval_v1' for s in sources) == 1, 'Source roster incomplete or repeated')
    return sorted(sources, key=lambda source: source['label'])


def run(args):
    session = args.session.resolve(); check_deadline(session); started = time.monotonic()
    design_path = session / 'inputs/expanded_family_design_v1.json'
    require(file_hash(__file__) == args.expected_code_sha256 and file_hash(design_path) == args.expected_design_sha256, 'Explicit frozen invocation changed')
    design = read(design_path)
    require(file_hash(design_path) == ORIGINAL_FAMILY_DESIGN_SHA and design['near90'] == 'not included' and design['budgets'] == list(KS), 'Original scientific design/settings changed')
    require(args.execution_freeze is not None and args.expected_execution_sha256 and
            file_hash(args.execution_freeze) == args.expected_execution_sha256, 'Explicit execution/code freeze is required')
    execution = read(args.execution_freeze)
    require(execution['version'] == 'expanded-family-source-admission-v2' and
            execution['status'] == 'frozen_before_expanded_family_outcome_evaluation' and
            execution['scientific_design_sha256'] == ORIGINAL_FAMILY_DESIGN_SHA and execution['labels_sha256'] == FAMILY_LABEL_SHA and
            execution['runner_sha256'] == file_hash(__file__) and execution['manifest_sha256'] == MANIFEST_SHA and
            execution['original_runner_sha256'] == design['code_hashes']['scripts/session_expanded_family.py'] and
            execution['consumer_code_hashes'] == source_code() == {**design['code_hashes'], 'scripts/session_expanded_family.py': file_hash(__file__)},
            'Execution amendment changes unrelated scientific inputs or dependencies')
    require(file_hash(session / 'inputs/expanded_retrieval_design_v1.json') == design['retrieval_design_sha256'] == RETRIEVAL_SHA and
            file_hash(session / 'audit/alignment_normalization_design_v1.json') == design['normalization_design_sha256'], 'Original retrieval or normalization design changed')
    labels_path = session / design['labels_path']; require(file_hash(labels_path) == design['labels_sha256'], 'Frozen labels changed')
    labels = read(labels_path); suite = load_suite(ROOT, session / 'inputs')
    require(suite['verification']['challenge_manifest_sha256'] == design['challenge_sha256'], 'Challenge changed')
    expansion = load_expansion(ROOT, session, session / design['manifest_path'])
    require(expansion['manifest_sha256'] == design['manifest_sha256'] and build_labels(session, expansion, suite) == labels,
            'Independent metadata/label/mask reconstruction differs from freeze')
    new_order = [a for a in expansion['manifest']['new_representatives_hash_order'] if a in set(labels['subsets']['1000']) - set(labels['subsets']['101'])]
    require(len(new_order) == len(set(new_order)) == 899 and execution['new_reference_order'] == new_order,
            'Execution freeze does not bind the exact899 new representatives')
    sources = open_sources(args, session, design, labels, suite, execution)
    evidence = [{'path': str(s['path']), 'identity_sha256': s['identity_hash'], 'completion_sha256': s['completion_sha256'], 'label': s['label'], 'arms': s['arms']} for s in sources]
    identity = {'mode': 'expanded_exact_family_filtering_v1', 'design_sha256': file_hash(design_path), 'labels_sha256': file_hash(labels_path),
                'execution_freeze_sha256': file_hash(args.execution_freeze),
                'sources': evidence, 'query_ids': sorted(labels['query_masks_and_survivors']), 'code_hashes': source_code()}
    require(Path(args.name).name == args.name and args.name not in ('', '.', '..'), 'Unsafe job name')
    store = UnitStore(session / 'analysis' / (args.name + '-' + digest(canonical(identity))[:12]), identity, session / 'job_ledger.jsonl')
    cohorts = {name: set(ids) for name, ids in labels['cohorts'].items()}; accumulators = defaultdict(Accumulator); expected = []
    for source_index, source in enumerate(sources):
        for qid, frozen in sorted(labels['query_masks_and_survivors'].items()):
            check_deadline(session)
            query = frozen['query']; original = read(source['units']['query-' + qid] / 'record.json')
            require(original['query'] == query, 'Expanded query metadata differs from freeze')
            expected_query = {**query, 'components': [labels['target_aliases'][a] for a in query['components']],
                              'spans': [{**span, 'accession': labels['target_aliases'][span['accession']]} for span in query['spans']]}
            require(original['representative_label_query'] == expected_query and set(original['rankings']) == set(SIZES), 'Representative query or prefix coverage differs')
            record = {'query_id': qid, 'source': evidence[source_index], 'source_unit_descriptor_sha256': source['unit_descriptor_sha256']['query-' + qid], 'prefixes': {}}
            for size in SIZES:
                row = frozen['prefixes'][size]; mask = row['mask']; included = sorted(set(labels['subsets'][size]) - set(mask['removed_representatives']))
                require(len(included) == mask['remaining_count'] and digest(canonical(included)) == mask['included_representatives_sha256'], 'Frozen query mask changed')
                arms = {}; included_set = set(included)
                for arm in source['arms']:
                    native = original['rankings'][size][arm]; validate_ranking(native, labels['subsets'][size])
                    ranking = filter_ranking(native, included); validate_ranking(ranking, included)
                    require([r['accession'] for r in ranking] == [r['accession'] for r in native if r['accession'] in included_set], 'Filtering changed score order')
                    evaluations = {stratum: evaluate_filtered(query, ranking, row, stratum) for stratum in ('native', 'common11')}
                    arms[arm] = {'evaluations': evaluations, 'top_returns': [{k: r[k] for k in ('accession', 'rank', 'score')} for r in ranking[:10]],
                                 'returned_unique_counts': {str(k): len(ranking[:k]) for k in KS}}
                    for stratum, evaluation in evaluations.items():
                        for cohort, ids in cohorts.items():
                            if qid in ids: accumulators[(source['label'], arm, size, stratum, cohort)].add(query, evaluation)
                record['prefixes'][size] = {'mask': mask, 'arms': arms}
            key = f'source-{source_index}-query-{qid}'; expected.append(key); old = store.read(key)
            if old is None: store.publish(key, record)
            else: require(old == record, 'Recomputed derived family query differs')
    summaries = []
    for (method, arm, size, stratum, cohort), accumulator in sorted(accumulators.items()):
        summaries.append({'method': method, 'arm': arm, 'prefix': size, 'stratum': stratum, 'cohort': cohort,
                          'retrieval_budget': 'k total unique representatives per original query',
                          **accumulator.result(labels['denominators'][size][stratum][cohort])})
    require(len(summaries) == 7 * 4 * 2 * 3, 'Missing method/prefix/stratum/cohort summary')
    report = {'summaries': summaries, 'sources': evidence, 'design_sha256': file_hash(design_path), 'labels_sha256': file_hash(labels_path),
              'execution_freeze_sha256': file_hash(args.execution_freeze),
              'interpretation': design['interpretation'], 'near90': 'not included',
              'scope': 'Complete saved-ranking filtering only; no model or fresh alignment scoring; common11 keeps all original component/query counts and excludes ineligible targets from recovery denominators'}
    old = store.read('summary')
    if old is None: store.publish('summary', report)
    else: require(old == report, 'Recomputed family summary differs')
    expected.append('summary'); store.complete(expected, {'elapsed_seconds': time.monotonic()-started})
    print(json.dumps({'status': 'complete', 'job': str(store.root), 'summary_rows': len(summaries)}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('freeze', 'freeze-execution', 'run')); parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--parent-job', action='append', type=Path, default=[])
    parser.add_argument('--model-job', action='append', type=Path, default=[]); parser.add_argument('--alignment-job', type=Path)
    parser.add_argument('--expected-code-sha256'); parser.add_argument('--expected-design-sha256'); parser.add_argument('--name', default='expanded-family-exact-v1')
    parser.add_argument('--execution-freeze', type=Path); parser.add_argument('--expected-execution-sha256')
    args = parser.parse_args(); {'freeze': freeze, 'freeze-execution': freeze_execution, 'run': run}[args.phase](args)
