"""Fixed composition design, duplicate-aware evaluation and validated cache I/O."""
from collections import defaultdict
import itertools
import json
from pathlib import Path

from .t05_data import load_catalogue, local_data_path, sha256, build_query


def canonical_hash(value):
    return sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode())


def verify_freeze(root, freeze):
    for relative, digest in freeze['file_sha256'].items():
        path = (root/relative).resolve()
        if not path.is_relative_to(root.resolve()) or sha256(path.read_bytes()) != digest:
            raise ValueError('Frozen input hash mismatch')


def load_reviewed(root, manifest):
    if manifest['manifest_version'] != 't05-composition-v1':
        raise ValueError('Unsupported composition manifest')
    sequences = {}
    for record in manifest['records']:
        accession = record['accession']
        if accession in sequences:
            raise ValueError('Duplicate accession')
        page = local_data_path(root,record['source_page_path'])
        if sha256(page.read_bytes()) != record['source_page_sha256']:
            raise ValueError('Original source page checksum mismatch')
        # Original validator enforces source versions, permissions and FASTA refs.
        # Single-record validation intentionally permits duplicates ACROSS accessions.
        sequences.update(load_catalogue(root,{'manifest_version':'t05-benign-v1','records':[record]}))
    if not sequences:
        raise ValueError('Empty catalogue')
    if duplicate_groups(sequences) != manifest['duplicate_sequence_groups']:
        raise ValueError('Duplicate-sequence audit differs from frozen manifest')
    return sequences


def duplicate_groups(sequences):
    groups = defaultdict(list)
    for accession,sequence in sorted(sequences.items()):
        groups[sha256(sequence.encode('ascii'))].append(accession)
    return sorted([g for g in groups.values() if len(g)>1])


def make_plan(targets, sequences, limit=1022):
    if len(targets) != len(set(targets)) or any(t not in sequences for t in targets):
        raise ValueError('Targets must be unique catalogue IDs')
    pairs, excluded = [], []
    for a,b in itertools.combinations(sorted(targets),2):
        n = len(sequences[a])+len(sequences[b])
        item = {'pair_id':a+'__'+b,'components':[a,b],'length_aa':n,
                'smaller_fraction':min(len(sequences[a]),len(sequences[b]))/n}
        if n > limit:
            excluded.append({**item,'reason':'combined length exceeds 1022; no truncation'})
        else:
            pairs.append(item)
    queries = [{'query_id':'clean_'+t,'components':[t],'pair_id':None,'artificial':False} for t in targets]
    for p in pairs:
        for order in (p['components'],p['components'][::-1]):
            queries.append({'query_id':'join_'+'__'.join(order),'pair_id':p['pair_id'],
                            'components':order,'artificial':True})
    subset_pairs = set()
    if pairs:
        for key in (lambda p:(p['length_aa'],p['pair_id']),
                    lambda p:(-p['length_aa'],p['pair_id']),
                    lambda p:(p['smaller_fraction'],p['pair_id']),
                    lambda p:(-p['smaller_fraction'],p['pair_id'])):
            subset_pairs.add(min(pairs,key=key)['pair_id'])
    return {'version':'t05-composition-queries-v1','targets':targets,'pairs':pairs,
            'excluded_pairs':excluded,'queries':queries,
            'repeat_query_ids':[q['query_id'] for q in queries if not q['artificial'] or q['pair_id'] in subset_pairs],
            'repeat_pair_ids':sorted(subset_pairs),'fraction_bins':[0,.25,.40,.60,.75,1],
            'near_tie_tolerance':1e-6}


def fraction_bin(fraction):
    if not 0 < fraction <= 1:
        raise ValueError('Invalid fraction')
    for lo,hi in ((0,.25),(.25,.4),(.4,.6),(.6,.75),(.75,1)):
        if lo <= fraction < hi or (hi==1 and fraction==1):
            return f'{lo:.2f}–{hi:.2f}'


def query_recall(ranking, expected, equivalence):
    if not expected or len(expected) != len(set(expected)):
        raise ValueError('Expected components must be distinct')
    if not set(expected) <= {r['accession'] for r in ranking}:
        raise ValueError('Expected component absent from rankings')
    result = {}
    for k in (1,5,10):
        found = {r['accession'] for r in ranking[:k]}
        hits = len(set(expected)&found)
        result[f'recall_at_{k}'] = hits/len(expected)
        result[f'any_component_at_{k}'] = int(hits>0)
        result[f'all_components_at_{k}'] = int(hits==len(expected))
        result[f'sequence_equivalent_recall_at_{k}'] = sum(bool(set(equivalence[a])&found) for a in expected)/len(expected)
    return result


def component_metrics(ranking, accession, equivalence, clean=None, tolerance=1e-6):
    item = next(r for r in ranking if r['accession']==accession)
    score, rank = item['cosine_similarity'], item['rank']
    tied = [r['rank'] for r in ranking if r['cosine_similarity']==score]
    near = [r['rank'] for r in ranking if abs(r['cosine_similarity']-score)<=tolerance]
    eqrank = min(r['rank'] for r in ranking if r['accession'] in equivalence[accession])
    result = {**item,'reciprocal_rank':1/rank,'tie_best_rank':min(tied),'tie_worst_rank':max(tied),
              'near_tie_count_1e6':len(near),'near_tie_best_rank':min(near),'near_tie_worst_rank':max(near),
              'sequence_equivalent_best_rank':eqrank,'duplicate_sequence_count':len(equivalence[accession])}
    for k in (1,5,10):
        result[f'hit_at_{k}'] = int(rank<=k)
        result[f'sequence_equivalent_hit_at_{k}'] = int(eqrank<=k)
    base = clean or result
    result.update(clean_rank=base['rank'],rank_change=rank-base['rank'],
                  clean_cosine=base['cosine_similarity'],cosine_change=score-base['cosine_similarity'],
                  reciprocal_rank_change=1/rank-base['reciprocal_rank'])
    return result


def equivalence_map(sequences):
    groups = defaultdict(list)
    for a,s in sequences.items():
        groups[sha256(s.encode())].append(a)
    return {a:sorted(groups[sha256(s.encode())]) for a,s in sequences.items()}


def cache_write(directory, key_data, vectors):
    import numpy as np
    directory.mkdir(parents=True,exist_ok=True)
    key = canonical_hash(key_data)
    array_path, meta_path = directory/(key+'.npz'), directory/(key+'.json')
    if array_path.exists() or meta_path.exists():
        raise ValueError('Refusing to overwrite existing cache')
    np.savez_compressed(array_path,vectors=vectors)
    meta_path.write_text(json.dumps({'key':key,'key_data':key_data,
                                    'array_sha256':sha256(array_path.read_bytes())},indent=2)+'\n')
    return key


def cache_read(directory, key_data, expected_shape):
    import numpy as np
    key = canonical_hash(key_data)
    array_path, meta_path = directory/(key+'.npz'), directory/(key+'.json')
    if not array_path.exists() and not meta_path.exists():
        return None
    metadata = json.loads(meta_path.read_text())
    if (metadata['key'] != key or metadata['key_data'] != key_data or
            metadata['array_sha256'] != sha256(array_path.read_bytes())):
        raise ValueError('Cache provenance/checksum mismatch')
    with np.load(array_path,allow_pickle=False) as data:
        result = data['vectors'].copy()
    if result.shape != expected_shape or result.dtype != np.float32 or not np.isfinite(result).all():
        raise ValueError('Invalid cached embedding array')
    return result
