"""Local inference adapters. No training, generation, or remote inference APIs."""
from __future__ import annotations

import importlib.metadata
import random
import re
from pathlib import Path
from typing import Protocol

from .bench_inputs import validate_sequence
from .runtime import InputError, digest_file

MODEL_REPOS = {
    'esm2': 'facebook/esm2_t6_8M_UR50D',
    'esmc': 'EvolutionaryScale/esmc-300m-2024-12',
    'esm3': 'EvolutionaryScale/esm3-sm-open-v1',
}
ESM2_REVISION = 'c731040fcd8d73dceaa04b0a8e6329b345b0f5df'


def _validate_snapshot_revision(snapshot: Path, requested: str) -> None:
    if re.fullmatch('[0-9a-f]{40}', snapshot.name) is None:
        raise InputError('Model cache did not resolve an immutable snapshot revision')
    if re.fullmatch('[0-9a-f]{40}', requested) and snapshot.name != requested:
        raise InputError('Resolved model snapshot differs from the pinned requested revision')


class Adapter(Protocol):
    metadata: dict
    device: str
    def residues(self, sequence: str) -> list[list[float]]: ...
    def close(self) -> None: ...


def prepare_device(device: str, seed: int) -> str:
    import torch
    if device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    if device == 'cuda' and not torch.cuda.is_available():
        raise InputError('CUDA requested but unavailable; use --device cpu')
    random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.benchmark = False
    torch.use_deterministic_algorithms(True)
    return device


class LocalAdapter:
    def __init__(self, name: str, device: str, revision: str, max_residues: int = 1022,
                 offline: bool = True):
        if name not in MODEL_REPOS:
            raise InputError('Unknown model; choose esm2, esmc or esm3')
        if not 1 <= max_residues <= 1022:
            raise InputError('max_residues must be in 1..1022 (shared experimental cap)')
        self.name, self.device, self.max_residues = name, device, max_residues
        self.model = None
        self.metadata = {'model': name, 'repository': MODEL_REPOS[name], 'requested_revision': revision,
                         'max_residues': max_residues, 'device': device, 'offline': offline,
                         'representation': 'Final adapter output; residue mean; BOS/EOS/PAD excluded',
                         'input_tracks': 'sequence only', 'precision_policy': 'float32 weights; no autocast',
                         'truncation': False}
        if name == 'esm2':
            from transformers import AutoModel, AutoTokenizer
            from huggingface_hub import snapshot_download
            snapshot = Path(snapshot_download(MODEL_REPOS[name], revision=revision,
                local_files_only=offline, allow_patterns=['config.json', 'tokenizer_config.json',
                    'tokenizer.json', 'vocab.txt', 'special_tokens_map.json', 'model.safetensors']))
            _validate_snapshot_revision(snapshot, revision)
            self.tokenizer = AutoTokenizer.from_pretrained(snapshot, local_files_only=True)
            self.model, loading = AutoModel.from_pretrained(snapshot, local_files_only=True,
                add_pooling_layer=False, output_loading_info=True)
            if loading.get('missing_keys') or loading.get('mismatched_keys'):
                raise InputError('Checkpoint has missing/mismatched parameters; pretrained encoder load refused')
            if any(not key.startswith('lm_head.') for key in loading.get('unexpected_keys', [])):
                raise InputError('Checkpoint contains unexpected encoder parameters')
            self.model = self.model.float().to(device).eval()
            self.metadata.update(resolved_revision=snapshot.name,
                weight_sha256=digest_file(snapshot / 'model.safetensors'),
                unused_checkpoint_heads='Masked language model head excluded; no additional pooler initialized')
        else:
            if importlib.metadata.version('esm') != '3.2.1':
                raise InputError('ESM-C/ESM-3 adapters require esm==3.2.1 in a separate model environment')
            from huggingface_hub import snapshot_download
            from unittest.mock import patch
            import esm.pretrained as pretrained
            weight = 'esmc_300m_2024_12_v0.pth' if name == 'esmc' else 'esm3_sm_open_v1.pth'
            snapshot = Path(snapshot_download(MODEL_REPOS[name], revision=revision,
                                             local_files_only=offline,
                                             allow_patterns=['data/weights/' + weight]))
            _validate_snapshot_revision(snapshot, revision)
            self.metadata['resolved_revision'] = snapshot.name
            self.metadata['weight_sha256'] = digest_file(snapshot / 'data' / 'weights' / weight)
            # SDK 3.2.1 builders otherwise resolve mutable main internally. Bind just the
            # synchronous load to the selected snapshot; sequence inference needs no decoders.
            with patch.object(pretrained, 'data_root', lambda _: snapshot):
                if name == 'esmc':
                    self.model = pretrained.ESMC_300M_202412(device=device, use_flash_attn=False)
                    self.tokenizer = self.model.tokenizer
                else:
                    self.model = pretrained.ESM3_sm_open_v0(device=device)
                    self.tokenizer = self.model.tokenizers.sequence
            self.model = self.model.float().to(device).eval()

        import torch
        self.metadata.update(torch_num_threads=torch.get_num_threads(),
                             torch_num_interop_threads=torch.get_num_interop_threads())

    def residues(self, sequence: str) -> list[list[float]]:
        import torch
        validate_sequence(sequence, 'protein')
        if len(sequence) > self.max_residues:
            raise InputError('Protein exceeds shared residue cap; no truncation was performed')
        with torch.inference_mode():
            if self.name == 'esm2':
                tokens = self.tokenizer(sequence, return_tensors='pt', truncation=False,
                                        return_special_tokens_mask=True)
                special = tokens.pop('special_tokens_mask')[0].bool()
                if self.tokenizer.unk_token_id in tokens['input_ids'][0].tolist():
                    raise InputError('Tokenizer produced an unknown residue token')
                keep = (~special & tokens['attention_mask'][0].bool()).to(self.device)
                tokens = {k: v.to(self.device) for k, v in tokens.items()}
                full = self.model(**tokens).last_hidden_state
                if full.ndim != 3 or full.shape[0] != 1 or full.shape[1] != len(keep):
                    raise InputError('Model/token shape mismatch; no pooling performed')
                states = full[0][keep]
            else:
                from esm.sdk.api import ESMProtein
                encoded = self.model.encode(ESMProtein(sequence=sequence))
                ids = encoded.sequence
                if ids is None or ids.ndim != 1 or len(ids) != len(sequence) + 2:
                    raise InputError('SDK token count differs from residue count plus BOS/EOS')
                if ids[0].item() != self.tokenizer.bos_token_id or ids[-1].item() != self.tokenizer.eos_token_id:
                    raise InputError('SDK BOS/EOS token convention changed')
                if self.tokenizer.unk_token_id in ids[1:-1].tolist():
                    raise InputError('SDK tokenizer produced an unknown residue token')
                full = self.model(sequence_tokens=ids.unsqueeze(0)).embeddings
                if full is None or full.ndim != 3 or full.shape[0] != 1 or full.shape[1] != len(ids):
                    raise InputError('SDK model/token shape mismatch; no pooling performed')
                states = full[0, 1:-1]
            if states.ndim != 2 or states.shape[0] != len(sequence) or not torch.isfinite(states).all():
                raise InputError('Model returned invalid residue representations; no pooling performed')
            return states.float().cpu().tolist()

    def close(self):
        self.model = None
        if self.device.startswith('cuda'):
            import torch
            torch.cuda.empty_cache()
