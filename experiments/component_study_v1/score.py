"""Score every distinct eligible unit against the frozen reference panel.

One stage per method, same queries and references:

    python experiments/component_study_v1/score.py --method esm2
    python experiments/component_study_v1/score.py --method alignment

Reuses the project's existing ESM-2 adapter, pooling/cosine helpers and
BLOSUM62 alignment scorer. Score matrices stay private; a source-free record of
model identity, work and resources is written for Git. No condition, merge rule
or label is evaluated here.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import resource
import time
from types import SimpleNamespace

from iri_constructs.runtime import run_metadata  # noqa: F401  (kept importable; not used for paths)

ROOT = Path(__file__).resolve().parents[2]
SEED, THREADS = 20261117, 2


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def canonical(value) -> str:
    return sha(json.dumps(value, sort_keys=True, separators=(',', ':')))


def load(args):
    freeze = json.loads(args.freeze.read_text())
    private = json.loads(args.private.read_text())
    if private['panel_sha256'] != freeze['panel_sha256']:
        raise SystemExit('Private panel does not belong to the committed freeze')
    references = {r: private['references'][r]['sequence'] for r in freeze['reference_ids']}
    frozen = {r['reference_id']: r['sequence_sha256'] for f in freeze['families'].values() for r in f['references']}
    if any(sha(s) != frozen[r] for r, s in references.items()):
        raise SystemExit('Reference sequence changed after the freeze')
    sequences, current = {}, None
    for line in (args.bundle / 'proteins.fasta').read_text().splitlines():
        if line.startswith('>'):
            current = line[1:].split()[0]
            sequences[current] = ''
        else:
            sequences[current] += line
    queries = {}
    for unit in freeze['units']:
        if unit['eligible']:
            if sha(sequences[unit['unit_id']]) != unit['sequence_sha256']:
                raise SystemExit('Unit sequence changed after the freeze')
            queries[unit['sequence_sha256']] = sequences[unit['unit_id']]
    return freeze, references, dict(sorted(queries.items()))


def esm2(references, queries):
    import torch
    from iri_constructs.metrics import cosine, mean_pool
    from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device
    torch.set_num_threads(THREADS)
    device = prepare_device('cpu', SEED)
    started = time.perf_counter()
    adapter = LocalAdapter('esm2', device, ESM2_REVISION, 1022, offline=True)
    load_seconds = time.perf_counter() - started
    vectors, stages = {}, {}
    for stage, items in (('references', references), ('queries', queries)):
        started = time.perf_counter()
        for key, sequence in items.items():
            vectors[(stage, key)] = mean_pool(adapter.residues(sequence))
        stages[stage] = {'model_calls': len(items), 'residues': sum(map(len, items.values())),
                         'seconds': round(time.perf_counter() - started, 3)}
    metadata = dict(adapter.metadata)
    adapter.close()
    started = time.perf_counter()
    scores = {q: {r: cosine(vectors[('queries', q)], vectors[('references', r)]) for r in references} for q in queries}
    return scores, {'model': metadata, 'load_seconds': round(load_seconds, 3), 'stages': stages,
                    'cosine_seconds': round(time.perf_counter() - started, 3), 'seed': SEED, 'torch_threads': THREADS}


def alignment(references, queries):
    from iri_constructs.pilot_alignment import alignment_scores
    proteins = {k: SimpleNamespace(sequence=s) for k, s in {**references, **queries}.items()}
    result = alignment_scores(proteins, {'query_ids': list(queries), 'eligible_reference_ids': {q: list(references) for q in queries}})
    if result['failures']:
        raise SystemExit('Alignment failed for some queries')
    return result['raw_scores'], {'configuration': result['configuration'], 'pairwise_score_calls': result['pairwise_score_calls'],
                                  'seconds': round(result['resources']['wall_seconds'], 3),
                                  'score_used': 'raw local alignment score (the historical primary); normalized values are not used'}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--method', choices=('esm2', 'alignment'), required=True)
    parser.add_argument('--bundle', type=Path, default=ROOT / 'results' / 'component_study_v1' / 'bundle')
    parser.add_argument('--private', type=Path, default=ROOT / 'results' / 'component_study_v1' / 'panel.private.json')
    parser.add_argument('--freeze', type=Path, default=ROOT / 'docs' / 'results' / 'component-study-v1' / 'panel.freeze.json')
    parser.add_argument('--out', type=Path, default=ROOT / 'results' / 'component_study_v1')
    parser.add_argument('--record', type=Path, default=ROOT / 'docs' / 'results' / 'component-study-v1')
    args = parser.parse_args(argv)
    started = time.perf_counter()
    freeze, references, queries = load(args)
    scores, work = (esm2 if args.method == 'esm2' else alignment)(references, queries)
    (args.out / f'scores.{args.method}.json').write_text(json.dumps({'panel_sha256': freeze['panel_sha256'], 'scores': scores}, sort_keys=True))
    record = {'schema_version': 'component-study-scoring-1', 'method': args.method, 'panel_sha256': freeze['panel_sha256'],
              'distinct_query_sequences': len(queries), 'references': len(references), 'score_pairs': len(queries) * len(references),
              'score_matrix_sha256': canonical(scores), 'work': work,
              'stage_wall_seconds': round(time.perf_counter() - started, 3),
              'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'bounds': {'wall_seconds': 1800, 'rss_gib': 6}, 'device': 'cpu',
              'note': 'Wall time and RSS describe this run on this machine only; OS caches were not controlled.'}
    (args.record / f'scoring.{args.method}.json').write_text(json.dumps(record, indent=1, sort_keys=True) + '\n')
    print(json.dumps({k: record[k] for k in ('method', 'distinct_query_sequences', 'references', 'stage_wall_seconds', 'peak_rss_kib', 'score_matrix_sha256')}))


if __name__ == '__main__':
    main()
