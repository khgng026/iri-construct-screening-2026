"""Exploratory scorers on the locked cohort v2 inputs: residue-level ESM-2 and ESM-2 150M.

The preprocessing arms, selected inputs, references, caps and evaluator key are those of
the locked holdout. Only the characterization step changes. Stages:

    collect   re-derive every arm's selected inputs, check them against the locked
              traces' fingerprints, and list the distinct proteins to score
    score     encode references and distinct inputs once per model and rank families
              by mean-pooled cosine and by residue-level local alignment
    trace     write verifier traces for each arm from the stored rankings

The residue scorer and its settings are those fixed in the related-family study
(``related_study._residue_local_score``; gaps -0.15/-0.02; centred cosine; per
shorter-protein residue). The best reference score per family ranks families.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from contextlib import contextmanager
from pathlib import Path

from .cohort_study import (AMINO_ACIDS, MAX_RESIDUES, MIN_RESIDUES, PYRODIGAL_VERSION, ROOT, candidate_units,
                           file_digest, read, run, scorer_protocol, write)

MODELS = {
    'esm2_8m': {'repository': 'facebook/esm2_t6_8M_UR50D',
                'revision': 'c731040fcd8d73dceaa04b0a8e6329b345b0f5df',
                'weight_sha256': '24c5fa474c48f3b754b86efe752d5f189d2bcd88190fa2270fc92b2ef3034189'},
    'esm2_150m': {'repository': 'facebook/esm2_t30_150M_UR50D',
                  'revision': 'a695f6045e2e32885fa60af20c13cb35398ce30c',
                  'weight_sha256': 'c3f1da8aea53bddd32c246c86168c23b9fd72341fb9db9a94436f855f5053566'},
}
# Scorer name -> (model, representation). esm2_8m_mean is an equivalence check against the locked traces.
SCORERS = {'esm2_8m_mean': ('esm2_8m', 'mean'), 'esm2_8m_residue': ('esm2_8m', 'residue'),
           'esm2_150m_mean': ('esm2_150m', 'mean'), 'esm2_150m_residue': ('esm2_150m', 'residue')}
PARTITIONS = {
    'holdout': {'inputs_key': 'holdout', 'fasta': 'prepared/holdout/method_inputs/constructs.fasta',
                'units': {'annotation_reference': 'prepared/holdout/evaluator/annotation-units.json'},
                'arms': ('minimal', 'start_ranked', 'gene_caller', 'annotation_reference'), 'locked': 'holdout'},
    'holdout_fusion': {'inputs_key': 'holdout_fusion', 'fasta': 'prepared/holdout/fusion/constructs.fasta',
                       'units': {'annotation_reference': 'prepared/holdout/fusion/annotation-units.json',
                                 'documented_regions': 'prepared/holdout/fusion/documented-regions.json'},
                       'arms': ('gene_caller', 'gene_caller_windows', 'annotation_reference', 'documented_regions'),
                       'locked': 'holdout-fusion'},
}


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


# ------------------------------------------------------------------ collect

def _records(run_dir: Path, spec, arm):
    from Bio import SeqIO

    if arm in spec['units']:
        return read(run_dir / spec['units'][arm])
    return {r.id: str(r.seq) for r in SeqIO.parse(run_dir / spec['fasta'], 'fasta')}


def selected_inputs(run_dir: Path, partition: str, arm: str):
    """Selected inputs per cap, rebuilt with the locked runner's rules and checked against its fingerprint."""
    spec = PARTITIONS[partition]
    protocol = scorer_protocol(spec['inputs_key'], 'esm2')
    finder = None
    if arm.startswith('gene_caller'):
        import pyrodigal
        if pyrodigal.__version__ != PYRODIGAL_VERSION:
            raise ValueError(f'requires pyrodigal=={PYRODIGAL_VERSION}')
        finder = pyrodigal.GeneFinder(meta=True, closed=False, mask=False)
    records = _records(run_dir, spec, arm)
    pools = {key: [u for u in candidate_units(arm, original, finder)
                   if MIN_RESIDUES <= len(u['protein']) <= MAX_RESIDUES and not set(u['protein']) - AMINO_ACIDS]
             for key, original in sorted(records.items())}
    fingerprint = hashlib.sha256()
    for cap in protocol['query_caps']:
        for key in sorted(records):
            encoded = json.dumps((cap, key, pools[key][:cap]), sort_keys=True, separators=(',', ':')).encode()
            fingerprint.update(encoded + b'\n')
    locked = read(run_dir / spec['locked'] / 'esm2' / f'{arm}.json')['selected_inputs_sha256']
    if fingerprint.hexdigest() != locked:
        raise ValueError(f'{partition}/{arm}: selected inputs differ from the locked trace')
    cap = max(protocol['query_caps'])
    return {u['protein'] for pool in pools.values() for u in pool[:cap]}


def collect(run_dir: Path, out: Path):
    proteins = set()
    for partition, spec in PARTITIONS.items():
        for arm in spec['arms']:
            found = selected_inputs(run_dir, partition, arm)
            print(f'{partition}/{arm}: fingerprint matches locked trace; {len(found)} distinct inputs')
            proteins |= found
    ordered = sorted(proteins, key=lambda p: (-len(p), p))
    write(out, {'proteins': ordered})
    print(f'{len(ordered)} distinct inputs, {sum(map(len, ordered))} residues')


# -------------------------------------------------------------------- score

class Encoder:
    """Final-layer residue states, special tokens excluded: the toolkit adapter's ESM-2 procedure."""

    def __init__(self, model: str, threads: int):
        import torch
        from huggingface_hub import snapshot_download
        from transformers import AutoModel, AutoTokenizer

        torch.set_num_threads(threads)
        spec = MODELS[model]
        snapshot = Path(snapshot_download(spec['repository'], revision=spec['revision'], local_files_only=True,
                                          allow_patterns=['config.json', 'tokenizer_config.json', 'vocab.txt',
                                                          'special_tokens_map.json', 'model.safetensors']))
        if snapshot.name != spec['revision'] or file_digest(snapshot / 'model.safetensors') != spec['weight_sha256']:
            raise ValueError(f'{model}: snapshot or weights differ from the pinned checkpoint')
        torch.manual_seed(20261122)
        self.tokenizer = AutoTokenizer.from_pretrained(snapshot, local_files_only=True)
        self.model, loading = AutoModel.from_pretrained(snapshot, local_files_only=True, add_pooling_layer=False,
                                                        output_loading_info=True)
        if loading.get('missing_keys') or loading.get('mismatched_keys') or any(
                not k.startswith('lm_head.') for k in loading.get('unexpected_keys', [])):
            raise ValueError(f'{model}: checkpoint parameters do not match the encoder')
        self.model = self.model.float().eval()
        self.metadata = {**spec, 'model': model, 'torch_threads': threads,
                         'representation': 'final layer; BOS/EOS excluded; float32'}

    def states(self, sequence: str):
        import numpy as np
        import torch

        with torch.inference_mode():
            tokens = self.tokenizer(sequence, return_tensors='pt', truncation=False, return_special_tokens_mask=True)
            special = tokens.pop('special_tokens_mask')[0].bool()
            if self.tokenizer.unk_token_id in tokens['input_ids'][0].tolist():
                raise ValueError('unknown residue token')
            keep = ~special & tokens['attention_mask'][0].bool()
            states = self.model(**tokens).last_hidden_state[0][keep]
        if states.shape[0] != len(sequence) or not torch.isfinite(states).all():
            raise ValueError('invalid residue representation')
        return states.float().numpy().astype(np.float32)


def _mean_vector(states):
    """Mean-pooled vector matching the locked runner (exact summation via the toolkit helper)."""
    import numpy as np
    from iri_constructs.metrics import mean_pool

    return np.asarray(mean_pool(states.tolist()), dtype=np.float64)


def _normalized_rows(states):
    import numpy as np

    return states / np.linalg.norm(states, axis=1, keepdims=True)


def _residue_scorer():
    from Bio.Align import PairwiseAligner

    from experiments.benign_protein_pilot.related_study import LOCAL_EXTEND, LOCAL_OPEN, _residue_local_score

    aligner = PairwiseAligner(mode='local', open_gap_score=LOCAL_OPEN, extend_gap_score=LOCAL_EXTEND)

    def score(query_rows, reference_rows):
        return _residue_local_score(query_rows, reference_rows, aligner)

    return score


def _rank(reference_scores):
    from .reporter_pilot import _rank_families

    return _rank_families(reference_scores)


def score(proteins_path: Path, references_path: Path, model: str, out: Path, threads: int = 2):
    """Encode the references and distinct inputs once, then rank families by mean cosine and residue alignment."""
    import numpy as np
    from iri_constructs.metrics import cosine

    proteins = read(proteins_path)['proteins']
    references = read(references_path)
    encoder = Encoder(model, threads)
    residue_score = _residue_scorer()

    started = time.perf_counter()
    ref_states = {r['accession']: encoder.states(r['sequence']) for r in references}
    ref_mean = {acc: _mean_vector(st) for acc, st in ref_states.items()}
    ref_rows = {acc: _normalized_rows(st) for acc, st in ref_states.items()}
    ref_family = {r['accession']: r['family'] for r in references}
    encode_seconds = time.perf_counter() - started

    mean_rank, residue_rank = {}, {}
    mean_seconds = residue_seconds = 0.0
    for protein in proteins:
        states = encoder.states(protein)
        t = time.perf_counter()
        vector = _mean_vector(states)
        mean_rank[sha(protein)] = _rank((ref_family[a], cosine(vector.tolist(), ref_mean[a].tolist()))
                                        for a in ref_family)
        mean_seconds += time.perf_counter() - t
        rows = _normalized_rows(states)
        t = time.perf_counter()
        residue_rank[sha(protein)] = _rank((ref_family[a], residue_score(rows, ref_rows[a])) for a in ref_family)
        residue_seconds += time.perf_counter() - t
    write(out, {'model': encoder.metadata, 'reference_count': len(references),
                'reference_residues': sum(len(r['sequence']) for r in references),
                'proteins': len(proteins), 'input_residues': sum(map(len, proteins)),
                'encode_seconds': encode_seconds, 'mean_scoring_seconds': mean_seconds,
                'residue_scoring_seconds': residue_seconds,
                'mean_rankings': mean_rank, 'residue_rankings': residue_rank})
    print(f'{model}: encoded {len(references)} refs + {len(proteins)} inputs in {encode_seconds:.0f}s; '
          f'mean {mean_seconds:.0f}s, residue {residue_seconds:.0f}s')


# -------------------------------------------------------------------- trace

@contextmanager
def _precomputed(rankings, representation):
    """Reporter-interface characterizer backed by stored family rankings for one representation."""
    table = rankings[f'{representation}_rankings']

    def rank_families(protein):
        return table[sha(protein)]

    yield rank_families, {'model': rankings['model'], 'representation': representation}


def trace(run_dir: Path, rankings_path: Path, scorer: str, out_dir: Path):
    """Write verifier traces for every arm of both partitions under one exploratory scorer."""
    rankings = read(rankings_path)
    model, representation = SCORERS[scorer]
    if rankings['model']['model'] != model:
        raise ValueError(f'{scorer}: rankings are for {rankings["model"]["model"]}, not {model}')
    references = run_dir / 'prepared' / 'references.private.json'
    for partition, spec in PARTITIONS.items():
        out = out_dir / partition / scorer
        out.mkdir(parents=True, exist_ok=True)
        protocol = scorer_protocol(spec['inputs_key'], 'esm2')
        protocol['characterizer'] = scorer
        protocol['characterization_id'] = f'{scorer}-cohort-v2'
        write(out / 'protocol.json', protocol)
        for arm in spec['arms']:
            inputs = run_dir / (spec['units'][arm] if arm in spec['units'] else spec['fasta'])
            run(arm, inputs, references, read(out / 'protocol.json'), out / f'{arm}.json',
                characterizer_factory=lambda name, refs: _precomputed(rankings, representation))
        print(f'{partition}/{scorer}: wrote {len(spec["arms"])} arms')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    c = sub.add_parser('collect')
    c.add_argument('--run', type=Path, default=ROOT / 'results' / 'paired-cohort-v2')
    c.add_argument('--out', type=Path, required=True)
    s = sub.add_parser('score')
    s.add_argument('--proteins', type=Path, required=True)
    s.add_argument('--references', type=Path, required=True)
    s.add_argument('--model', choices=sorted(MODELS), required=True)
    s.add_argument('--threads', type=int, default=2)
    s.add_argument('--out', type=Path, required=True)
    t = sub.add_parser('trace')
    t.add_argument('--run', type=Path, default=ROOT / 'results' / 'paired-cohort-v2')
    t.add_argument('--rankings', type=Path, required=True)
    t.add_argument('--scorer', choices=sorted(SCORERS), required=True)
    t.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.stage == 'collect':
        collect(args.run, args.out)
    elif args.stage == 'score':
        score(args.proteins, args.references, args.model, args.out, args.threads)
    else:
        trace(args.run, args.rankings, args.scorer, args.out_dir)


if __name__ == '__main__':
    main()