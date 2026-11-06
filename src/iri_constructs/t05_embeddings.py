"""Small, fixed ESM-2 encoder and cosine ranking. No training or segmentation."""
from __future__ import annotations

from time import perf_counter

from .t05_data import validate_protein


def mean_pool(hidden, attention_mask, special_tokens_mask):
    import torch
    if hidden.ndim != 3 or hidden.shape[:2] != attention_mask.shape or special_tokens_mask.shape != attention_mask.shape:
        raise ValueError('Incompatible pooling shapes')
    mask = attention_mask.bool() & ~special_tokens_mask.bool()
    counts = mask.sum(dim=1)
    if (counts == 0).any():
        raise ValueError('No actual residues to pool')
    if not torch.isfinite(hidden).all():
        raise ValueError('Nonfinite hidden state')
    pooled = hidden.masked_fill(~mask.unsqueeze(-1), 0).sum(dim=1) / counts.unsqueeze(-1)
    return pooled, mask


def rank_cosine(query, references, accessions):
    import torch
    if query.ndim != 1 or references.ndim != 2 or references.shape != (len(accessions), query.shape[0]):
        raise ValueError('Incompatible ranking shapes')
    if not accessions or len(set(accessions)) != len(accessions):
        raise ValueError('Catalogue IDs must be unique and nonempty')
    if not torch.isfinite(query).all() or not torch.isfinite(references).all():
        raise ValueError('Nonfinite ranking vectors')
    norms = references.norm(dim=1)
    qnorm = query.norm()
    if qnorm <= 0 or (norms <= 0).any():
        raise ValueError('Cannot rank a zero vector')
    scores = (references @ query / (norms * qnorm)).tolist()
    order = sorted(range(len(scores)), key=lambda i: (-scores[i], accessions[i]))
    return [{'rank': rank, 'accession': accessions[i], 'cosine_similarity': scores[i],
             'exact_score_tie_count': sum(score == scores[i] for score in scores)}
            for rank, i in enumerate(order, 1)]


class Encoder:
    def __init__(self, config, *, offline=True):
        import torch
        from transformers import EsmModel, EsmTokenizer
        if config['device'] != 'cpu' or config['dtype'] != 'float32':
            raise ValueError('This frozen pilot requires CPU float32')
        if config['max_residues'] != 1022 or config['tokenizer_truncation'] is not False:
            raise ValueError('Unsupported input policy')
        torch.set_num_threads(config['threads'])
        torch.manual_seed(config['seed'])
        torch.use_deterministic_algorithms(True)
        self.config = config
        options = dict(revision=config['revision'], local_files_only=offline)
        started = perf_counter()
        self.tokenizer = EsmTokenizer.from_pretrained(config['model_id'], **options)
        self.model = EsmModel.from_pretrained(config['model_id'], use_safetensors=True,
                                             add_pooling_layer=False, **options).float().eval()
        self.load_seconds = perf_counter() - started

    def encode(self, sequences):
        import torch
        if not sequences:
            raise ValueError('Empty input batch')
        for sequence in sequences:
            validate_protein(sequence, self.config['max_residues'])
        started = perf_counter()
        inputs = self.tokenizer(sequences, return_tensors='pt', padding=True,
                                truncation=False, return_special_tokens_mask=True)
        special = inputs.pop('special_tokens_mask')
        mask = inputs['attention_mask'].bool() & ~special.bool()
        counts = mask.sum(dim=1).tolist()
        if counts != [len(s) for s in sequences] or (inputs['input_ids'] == self.tokenizer.unk_token_id).any():
            raise ValueError('Tokenizer did not preserve one token per residue')
        token_seconds = perf_counter() - started
        with torch.inference_mode():
            started = perf_counter()
            hidden = self.model(**inputs).last_hidden_state
            forward_seconds = perf_counter() - started
            started = perf_counter()
            pooled, mask = mean_pool(hidden, inputs['attention_mask'], special)
            pool_seconds = perf_counter() - started
        if pooled.requires_grad or hidden.requires_grad:
            raise ValueError('Inference unexpectedly tracks gradients')
        return pooled, [hidden[i, mask[i]].clone() for i in range(len(sequences))], {
            'hidden_shape': list(hidden.shape), 'pooled_shape': list(pooled.shape),
            'residue_counts': counts, 'token_counts': inputs['attention_mask'].sum(dim=1).tolist(),
            'tokenize_seconds': token_seconds, 'forward_seconds': forward_seconds,
            'pool_seconds': pool_seconds, 'inference_mode': True, 'model_training': self.model.training,
            'truncated': False}
