"""Checkpointed fresh geometry/alignment and original-checkpoint ESM inference.

Run ONLY through the session resource guard. This imports Nishit's pinned code
without editing it. New score caches are separate from every historical store.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'sources/nishit-working/src'))
from pilot.decompose import best_hit_segments, oracle_span, tile_segments
from pilot.embed import Esm2Embedder
from budget_compare import atomic_json, digest


def configuration():
    d = json.loads((ROOT/'configs/frozen.json').read_text())
    if digest({k:v for k,v in d.items() if k != 'freeze_sha256'}) != d['freeze_sha256']:
        raise ValueError('Frozen design changed')
    return d


def save_array(path, array):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp.npy')
    np.save(tmp, array)
    tmp.replace(path)


def seq_key(seq):
    return hashlib.sha256(seq.encode()).hexdigest()


def file_hash(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def implementation_signature():
    paths=[Path(__file__),ROOT/'sources/nishit-working/src/pilot/embed.py',
           ROOT/'sources/nishit-working/src/pilot/decompose.py',ROOT/'sources/nishit-working/src/pilot/align.py']
    return digest({str(p.relative_to(ROOT)):file_hash(p) for p in paths})


def valid_array(path,shape,cfg,signature):
    meta_path=path.with_suffix('.json')
    if not path.exists() or not meta_path.exists():return None
    meta=json.loads(meta_path.read_text())
    if (meta.get('sequence_sha256')!=path.stem or meta.get('freeze_sha256')!=cfg['freeze_sha256'] or
        meta.get('implementation_signature')!=signature or meta.get('array_sha256')!=file_hash(path)):
        raise ValueError('Cache provenance mismatch: '+str(path))
    value=np.load(path)
    if value.shape!=shape or not np.isfinite(value).all():raise ValueError('Bad cache shape/value')
    return value


class BioAlignmentScorer:
    """Same SW scoring objective as upstream Parasail, different exact backend.

    Biopython gap length1 pays open only; set -11/-1 explicitly, because its
    scoring='blastp' preset instead sets -12/-1. Upstream span locator remains
    unmodified and is therefore -12/-1, exactly as in Nishit's own workflow.
    """
    def __init__(self, cfg):
        self.cfg,self.signature=cfg,implementation_signature()
        self.refs = [r['sequence'] for r in cfg['references']]
        self.path = ROOT/'results/alignment_scores'/cfg['freeze_sha256']
        self.al = Align.PairwiseAligner()
        self.al.mode = 'local'
        self.al.substitution_matrix = substitution_matrices.load('BLOSUM62')
        self.al.open_gap_score = -11
        self.al.extend_gap_score = -1
        self.fresh, self.reused, self.seconds = 0, 0, 0.0

    def score_all(self, query):
        path = self.path/(seq_key(query)+'.npy')
        x=valid_array(path,(len(self.refs),),self.cfg,self.signature)
        if x is not None:
            self.reused += 1
            return x
        t = time.perf_counter()
        x = np.array([self.al.score(query,s) for s in self.refs], dtype=np.float32)
        seconds = time.perf_counter()-t
        save_array(path,x)
        atomic_json(path.with_suffix('.json'), {'sequence_sha256':seq_key(query), 'seconds':seconds,
                    'fresh':True, 'backend':'Bio.PairwiseAligner', 'score_parameters':'local BLOSUM62 -11/-1',
                    'freeze_sha256':self.cfg['freeze_sha256'],'implementation_signature':self.signature,
                    'array_sha256':file_hash(path)})
        self.fresh += 1
        self.seconds += seconds
        return x


def span_record(s):
    return {'start':s.start-1,'end':s.end,'sequence_sha256':seq_key(s.seq),'sequence':s.seq,
            'removed':list(getattr(s,'removed',()))}


def geometry(cfg, subset):
    scorer = BioAlignmentScorer(cfg)
    refs = [x['accession'] for x in cfg['references']]
    refseqs = [x['sequence'] for x in cfg['references']]
    byacc = dict(zip(refs,refseqs))
    idx = {a:i for i,a in enumerate(refs)}
    run_start = time.perf_counter()
    for n,item in enumerate(cfg['items']):
        if subset != 'all' and item['split'] != subset:
            continue
        path = ROOT/'results/geometry'/(digest(item['construct_id'])[:24]+'.json')
        if path.exists():
            v=json.loads(path.read_text())
            if v['freeze_sha256'] != cfg['freeze_sha256'] or v['implementation_signature']!=scorer.signature:
                raise ValueError('Stale geometry')
            print('reuse geometry',item['construct_id'],flush=True)
            continue
        seq=item['sequence']; payload=item['payload_accession']; t=time.perf_counter()
        oracle=oracle_span(seq,byacc[payload])
        if not oracle or oracle.seq != byacc[payload] or oracle.length != int(item['oracle_len']):
            raise ValueError('Synthetic oracle does not match exact documented target: '+item['construct_id'])
        spans={'whole':[{'start':0,'end':len(seq),'sequence':seq,'sequence_sha256':seq_key(seq),'removed':[]}],
               'windows':[span_record(s) for s in tile_segments(seq,cfg['window'],cfg['step'],cfg['min_segment'])],
               'oracle':[span_record(oracle)]}
        conditions=[]
        for thr in cfg['thresholds']:
            hs=[h for h in cfg['homologs'] if h['payload_accession']==payload]
            mask=sorted({idx[h['homolog_accession']] for h in hs if float(h['identity'])>thr}) if thr<1 else []
            targets=sorted({idx[h['homolog_accession']] for h in hs if float(h['identity'])<=thr}) if thr<1 else [idx[payload]]
            bil=max((float(h['identity']) for h in hs if float(h['identity'])<=thr),default=None) if thr<1 else 1.0
            bhs_start=time.perf_counter()
            result=best_hit_segments(seq,refseqs,scorer,sum(map(len,refseqs)),refs,
                                    max_evalue=cfg['max_evalue'],min_len=cfg['min_segment'],exclude_idx=np.array(mask,dtype=int))
            bhs_seconds=time.perf_counter()-bhs_start
            bhs=[span_record(s) for s in result.segments]
            fallback=not bhs
            if fallback:
                bhs=spans['whole']
            bhs=sorted(bhs,key=lambda s:(s['start'],s['end']))
            condition={'threshold':thr,'mask':mask,'targets':targets,'best_identity_left':bil,
                       'bhs':bhs,'bhs_fallback':fallback,'bhs_stop':result.stop_reason,'bhs_floor':result.floor,
                       'bhs_trace':[asdict(s) for s in result.trace], 'bhs_seconds':bhs_seconds}
            conditions.append(condition)
        # Score all emitted unmodified query subsequences once, checkpoint by hash.
        for s in [s for ss in spans.values() for s in ss]+[s for c in conditions for s in c['bhs']]:
            if seq[s['start']:s['end']] != s['sequence']:
                raise ValueError('Coordinates do not reconstruct source sequence')
            scorer.score_all(s['sequence'])
        atomic_json(path,{'freeze_sha256':cfg['freeze_sha256'],'item_id':item['construct_id'],
                         'implementation_signature':scorer.signature,
                         'spans':spans,'conditions':conditions,'seconds':time.perf_counter()-t,
                         'input_sha256':item['sequence_sha256']})
        print(json.dumps({'completed':n+1,'total':len(cfg['items']),'item':item['construct_id'],
              'seconds':round(time.perf_counter()-t,2),'fresh_align_rows':scorer.fresh}),flush=True)
    atomic_json(ROOT/f'results/geometry_run_{subset}.json',{'freeze_sha256':cfg['freeze_sha256'],
                'wall_seconds':time.perf_counter()-run_start,'fresh_score_rows':scorer.fresh,
                'reused_score_rows':scorer.reused,'fresh_alignment_seconds':scorer.seconds,
                'process_maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})


def embeddings(cfg,subset,limit):
    import torch
    torch.set_num_threads(cfg['threads'])
    torch.set_num_interop_threads(1)
    os.environ['HF_HUB_OFFLINE']='1'
    os.environ['TRANSFORMERS_OFFLINE']='1'
    modelpath=Path.home()/'.cache/huggingface/hub/models--facebook--esm2_t33_650M_UR50D/snapshots'/cfg['model_revision']
    signature=implementation_signature()
    receipt_path=ROOT/'evidence'/('inference-source-receipt-'+signature[:16]+'.json')
    if not receipt_path.exists():
        sources=[modelpath/'model.safetensors',modelpath/'config.json',modelpath/'tokenizer_config.json',
                 ROOT.parents[1]/'iri-reviewed-core-inputs-2026-10-31.zip',
                 ROOT.parents[1]/'iri-construct-screening-2026/data/manifests/t05_benign_catalogue.json']
        atomic_json(receipt_path,{'implementation_signature':signature,'freeze_sha256':cfg['freeze_sha256'],
                    'files':{str(p):file_hash(p) for p in sources},'torch':torch.__version__})
    embedder=Esm2Embedder(str(modelpath),device='cpu',batch_size=1,max_length=cfg['max_length'])
    seqs={seq_key(r['sequence']):r['sequence'] for r in cfg['references']}
    for item in cfg['items']:
        if subset!='all' and item['split']!=subset:
            continue
        path=ROOT/'results/geometry'/(digest(item['construct_id'])[:24]+'.json')
        g=json.loads(path.read_text())
        if g['freeze_sha256'] != cfg['freeze_sha256'] or g['implementation_signature']!=signature:
            raise ValueError('Stale geometry')
        for s in [s for ss in g['spans'].values() for s in ss]+[s for c in g['conditions'] for s in c['bhs']]:
            seqs[s['sequence_sha256']]=s['sequence']
    folder=ROOT/'results/embeddings'/cfg['freeze_sha256']
    fresh=reused=0; start=time.perf_counter()
    # Length ordering bounds peak early; all keys are sequence hashes.
    for k,seq in sorted(seqs.items(),key=lambda x:(len(x[1]),x[0])):
        path=folder/(k+'.npy')
        v=valid_array(path,(1280,),cfg,signature)
        if v is not None:
            if abs(np.linalg.norm(v)-1)>1e-5:
                raise ValueError('Bad cached representation')
            reused+=1
            continue
        if len(seq)>cfg['max_query_residues']:
            raise ValueError('No silent truncation')
        t=time.perf_counter(); vector=embedder.embed([seq])[0]
        if vector.shape!=(1280,) or not np.isfinite(vector).all():
            raise ValueError('Invalid fresh representation')
        save_array(path,vector)
        meta={'sequence_sha256':k,'residues':len(seq),'fresh':True,'seconds':time.perf_counter()-t,
              'checkpoint':cfg['model_revision'],'freeze_sha256':cfg['freeze_sha256'],
              'implementation_signature':signature,'array_sha256':file_hash(path),
              'dtype':str(vector.dtype),'process_maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
        atomic_json(path.with_suffix('.json'),meta)
        fresh+=1
        print(json.dumps({'fresh':fresh,'reused':reused,'total':len(seqs),**meta}),flush=True)
        if limit and fresh>=limit:
            break
    atomic_json(ROOT/f'results/embedding_run_{subset}_{limit or "full"}.json',
                {'freeze_sha256':cfg['freeze_sha256'],'fresh':fresh,'reused':reused,'required':len(seqs),
                 'wall_seconds':time.perf_counter()-start,'process_maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 'complete':fresh+reused==len(seqs),'torch':torch.__version__,'device':'cpu'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['geometry','embed']);p.add_argument('--split',default='all',choices=['all','development','evaluation']);p.add_argument('--limit',type=int,default=0)
    a=p.parse_args();c=configuration()
    if a.mode=='geometry':geometry(c,a.split)
    else:embeddings(c,a.split,a.limit)
