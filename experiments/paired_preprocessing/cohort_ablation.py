"""Leave-close-references-out ablation: does the alignment/model gap close as identity falls?

For each documented component we take its validated protein as the query and score it against
the reference panel with each method. We then remove references that are too similar to the
query (identity >= threshold T) and recompute the top-family call on the thinned panel. As T
falls, the panel keeps only ever-more-distant references, manufacturing the low-identity
regime the biosecurity literature identifies as where sequence alignment fails
(Wittmann et al., Science 2025; Front. Bioeng. Biotechnol. 2026). The headline is the
crossover identity: the retained-reference identity below which a foundation model's recall
exceeds conventional alignment, if any.

Benign and defensive: queries are the cohort's public benign components; the regime is reached
by withholding close references from naturally occurring homologs, never by designing divergent
sequences.

    score     encode references and distinct component proteins once per model; store, for
              every (component, reference) pair, the BLOSUM identity and the normalized BLOSUM
              score, the mean-pooled cosine and the residue-level local score
    analyze    per identity threshold, recompute each method's top-family recall over the
              components that still retain a same-family reference, with group bootstrap

Methods: normalized local BLOSUM62 (the locked conventional comparator), mean-pooled ESM-2
cosine, and residue-level ESM-2 local alignment (the related-family study's scorer).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

from .cohort_scorers import Encoder, MODELS, _mean_vector, _normalized_rows, _residue_scorer, sha
from .cohort_study import ROOT, identity_coverage, read, write

THRESHOLDS = [1.01, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.35, 0.30, 0.25, 0.20]
BOOTSTRAP = {'replicates': 4000, 'seed': 20261123, 'interval': 0.95}


def _components(run: Path):
    """One query per corroborated component: its validated protein, family and lineage group."""
    construct_map = {c['id']: c for c in read(run / 'prepared' / 'holdout' / 'evaluator' / 'construct-map.private.json')}
    units = read(run / 'prepared' / 'holdout' / 'evaluator' / 'annotation-units.json')
    key = {c['id']: c for c in read(run / 'prepared' / 'holdout' / 'evaluator' / 'answer-key.corroborated.json')['constructs']}
    rows = []
    for cid, construct in key.items():
        meta = {c['id']: c for c in construct_map[cid]['components']}
        unit_by_feature = {u['id'].rsplit('-', 1)[-1]: u for u in units[cid]}
        for component in construct['components']:
            feature = component['id'].split('-', 1)[0]
            unit = unit_by_feature.get(feature)
            if unit is None:
                continue
            info = meta[component['id']]
            rows.append({'construct': cid, 'component': component['id'], 'family': component['family'],
                         'lineage_group': construct_map[cid]['lineage_group'], 'protein': unit['protein'],
                         'protein_sha256': info['protein_sha256'], 'fusion_class': info['fusion_class']})
    return rows


def score(run: Path, model: str, out: Path, threads: int = 4):
    from iri_constructs.metrics import cosine

    rows = _components(run)
    references = read(run / 'prepared' / 'references.private.json')
    encoder = Encoder(model, threads)
    residue_score = _residue_scorer()

    started = time.perf_counter()
    ref_states = {r['accession']: encoder.states(r['sequence']) for r in references}
    ref_mean = {a: _mean_vector(s) for a, s in ref_states.items()}
    ref_rows = {a: _normalized_rows(s) for a, s in ref_states.items()}
    ref_family = {r['accession']: r['family'] for r in references}
    ref_seq = {r['accession']: r['sequence'] for r in references}

    # Score each distinct component protein once; identity needs the BLOSUM aligner.
    from Bio.Align import PairwiseAligner, substitution_matrices
    aligner = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                              open_gap_score=-10.0, extend_gap_score=-0.5, epsilon=1e-6)
    self_scores = {a: aligner.score(s, s) for a, s in ref_seq.items()}

    cache = {}
    distinct = {r['protein_sha256']: r['protein'] for r in rows}
    for i, (digest, protein) in enumerate(distinct.items()):
        states = encoder.states(protein)
        q_mean, q_rows = _mean_vector(states), _normalized_rows(states)
        per_ref = {}
        for acc in ref_family:
            ic = identity_coverage(aligner, protein, ref_seq[acc])
            per_ref[acc] = {'identity': round(ic['identity'], 4),
                            'blosum': round(aligner.score(protein, ref_seq[acc]) / self_scores[acc], 6),
                            'mean': round(cosine(q_mean.tolist(), ref_mean[acc].tolist()), 6),
                            'residue': round(residue_score(q_rows, ref_rows[acc]), 6)}
        cache[digest] = per_ref
        if (i + 1) % 50 == 0:
            print(f'  scored {i + 1}/{len(distinct)} distinct proteins in {time.perf_counter() - started:.0f}s')
    write(out, {'model': encoder.metadata, 'reference_family': ref_family,
                'components': [{k: r[k] for k in ('construct', 'component', 'family', 'lineage_group',
                                                  'protein_sha256', 'fusion_class')} for r in rows],
                'scores_by_protein': cache, 'encode_seconds': time.perf_counter() - started})
    print(f'{model}: scored {len(distinct)} distinct proteins against {len(references)} references')


def _top_family(per_ref, ref_family, method, keep):
    best = {}
    for acc in keep:
        score = per_ref[acc][method]
        family = ref_family[acc]
        if family not in best or score > best[family]:
            best[family] = score
    if not best:
        return None
    return max(best, key=lambda f: (best[f], f))


def _interval(values):
    alpha = (1 - BOOTSTRAP['interval']) / 2
    return [float(np.quantile(values, alpha)), float(np.quantile(values, 1 - alpha))]


def _group_recall(rows, flags):
    groups = defaultdict(list)
    for row, hit in zip(rows, flags, strict=True):
        groups[row['lineage_group']].append(hit)
    per_group = {g: float(np.mean(v)) for g, v in groups.items()}
    return (sum(flags), len(flags), float(np.mean(list(per_group.values()))), per_group)


def analyze(score_paths, out: Path):
    """Recall vs identity threshold for each method, on the components scorable at each threshold."""
    loaded = [read(Path(p)) for p in score_paths]
    runs = {doc['model']['model']: doc for doc in loaded}
    base = loaded[0]
    ref_family = base['reference_family']
    components = base['components']
    # BLOSUM identity/score are model-independent, so they come from any file; the ESM
    # columns are per model. One method per (model, representation).
    model_methods = {f'{model}_{name}': (model, name) for model in runs for name in ('mean', 'residue')}
    report = {'thresholds': THRESHOLDS, 'bootstrap': BOOTSTRAP, 'components': len(components),
              'models': {m: runs[m]['model'] for m in runs}, 'by_threshold': {}}
    rng = np.random.default_rng(BOOTSTRAP['seed'])
    for threshold in THRESHOLDS:
        scorable, flags = [], defaultdict(list)
        for row in components:
            per_ref = base['scores_by_protein'][row['protein_sha256']]
            keep = [a for a in ref_family if per_ref[a]['identity'] < threshold]
            if not any(ref_family[a] == row['family'] for a in keep):
                continue  # no same-family reference survives; outside the denominator at this T
            scorable.append(row)
            flags['blosum62'].append(_top_family(per_ref, ref_family, 'blosum', keep) == row['family'])
            for label, (model, col) in model_methods.items():
                pr = runs[model]['scores_by_protein'][row['protein_sha256']]
                flags[label].append(_top_family(pr, ref_family, col, keep) == row['family'])
        entry = {'scorable_components': len(scorable),
                 'lineage_groups': len({r['lineage_group'] for r in scorable}),
                 'max_retained_identity_below': threshold, 'recall': {}}
        groups = sorted({r['lineage_group'] for r in scorable})
        gindex = {g: i for i, g in enumerate(groups)}
        draws = rng.integers(0, len(groups), size=(BOOTSTRAP['replicates'], len(groups))) if groups else None
        for method, hits in flags.items():
            total, n, eq, per_group = _group_recall(scorable, hits)
            rec = {'recovered': total, 'eligible': n, 'pooled_recall': total / n if n else None,
                   'equal_group_recall': eq}
            if groups and len(groups) > 1:
                gvec = np.array([per_group[g] for g in groups])
                rec['equal_group_interval'] = _interval(gvec[draws].mean(1))
            entry['recall'][method] = rec
        report['by_threshold'][f'{threshold:.2f}'] = entry
    write(out, report)
    for t in THRESHOLDS:
        e = report['by_threshold'][f'{t:.2f}']
        line = ' '.join(f"{m}={e['recall'][m]['equal_group_recall']*100:.0f}%" for m in sorted(e['recall']))
        print(f"T<{t:.2f}: n={e['scorable_components']:4d} G={e['lineage_groups']:3d}  {line}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    s = sub.add_parser('score')
    s.add_argument('--run', type=Path, default=ROOT / 'results' / 'paired-cohort-v2')
    s.add_argument('--model', choices=sorted(MODELS), required=True)
    s.add_argument('--threads', type=int, default=4)
    s.add_argument('--out', type=Path, required=True)
    a = sub.add_parser('analyze')
    a.add_argument('--scores', type=Path, nargs='+', required=True)
    a.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.stage == 'score':
        score(args.run, args.model, args.out, args.threads)
    else:
        analyze(args.scores, args.out)


if __name__ == '__main__':
    main()
