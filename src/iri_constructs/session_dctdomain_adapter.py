"""Offline HF ESM2 bridge for pinned DCTdomain; not an admitted reproduction.

The caller must establish numerical equivalence and resource admission first.
This shim supplies DCTdomain's original Embedding class with layers 15/21 and
learned contacts. Original window merging, RecCut and DCT code stay upstream.
No downloads, model generation, training, environment installation, or C++ build.
"""
from __future__ import annotations

import inspect
import re

from .session_models import load_encoder, sha256_file, validate_sequence, validate_token_mapping


MODEL_ID = "facebook/esm2_t30_150M_UR50D"
REVISION = "a695f6045e2e32885fa60af20c13cb35398ce30c"
WEIGHT_SHA256 = "c3f1da8aea53bddd32c246c86168c23b9fd72341fb9db9a94436f855f5053566"
DCT_REVISION = "17c1a423a337e808e2734071154495cfa56af295"
EMBEDDING_SHA256 = "dd8ace63fd2e7722e3241c9d71d5ff9d6717cf0320240024e9847127604eef75"
LAYERS = (15, 21)


def validate_dct_config(config):
    """Freeze the published method settings; no silent alternative model/device."""
    c = dict(config)
    required = {"model_id": MODEL_ID, "revision": REVISION, "device": "cpu",
                "dtype": "float32", "attn_implementation": "eager",
                "local_files_only": True, "layer": 15}
    for name, value in required.items():
        if name in c and c[name] != value:
            raise ValueError(f"DCTdomain adapter requires {name}={value}")
        c[name] = value
    expected = dict(c.get("expected_file_sha256", {}))
    if expected.get("model.safetensors", WEIGHT_SHA256) != WEIGHT_SHA256:
        raise ValueError("DCTdomain adapter requires the admitted 150M checkpoint hash")
    expected["model.safetensors"] = WEIGHT_SHA256
    c["expected_file_sha256"] = expected
    return c


class DCTdomainHFModel:
    """Minimal duck-typed replacement for upstream embedding.Model, batch one.

    Use the unchanged pinned ``Embedding(...).embed_seq(adapter, 'cpu', [15,21],
    500)`` or ``encode_with_upstream``. This captures only two hidden layers but
    temporarily needs all attention heads for the official learned contact head.
    Do not substitute the primary last-layer mean-pooling cache.
    """

    def __init__(self, config):
        self.encoder = load_encoder(validate_dct_config(config))
        self.torch = self.encoder.torch
        self.model = self.encoder.model
        self.tokenizer = self.encoder.tokenizer
        c = self.model.config
        if (c.num_hidden_layers, c.hidden_size, c.num_attention_heads) != (30, 640, 20):
            self.close()
            raise ValueError("Unexpected ESM2 150M architecture")
        from safetensors import safe_open
        weight_path = self.encoder._cached_file("model.safetensors")
        with safe_open(str(weight_path), framework="pt", device="cpu") as source:
            keys = set(source.keys())
            if not {"esm.contact_head.regression.weight", "esm.contact_head.regression.bias"} <= keys:
                self.close()
                raise ValueError("Checkpoint lacks the learned contact head")
            for suffix in ("weight", "bias"):
                expected = source.get_tensor(f"esm.contact_head.regression.{suffix}").float()
                actual = getattr(self.model.contact_head.regression, suffix).detach().cpu()
                if not self.torch.equal(actual, expected):
                    self.close()
                    raise ValueError("Loaded contact head differs from checkpoint")
        if tuple(self.model.contact_head.regression.weight.shape) != (1, 600):
            self.close()
            raise ValueError("Unexpected contact-head feature layout")
        self.metadata = dict(self.encoder.metadata, comparator="DCTdomain_HF_adapter",
            adapter_source_sha256=sha256_file(__file__), upstream_revision=DCT_REVISION,
            layers=list(LAYERS), representation="after blocks 15 and 21; no final norm",
            contacts="official EsmContactPredictionHead; all 30 layers x 20 heads",
            max_window_residues=500, overlap_residues=200,
            reccut_contact_multiplier=2.6, quantization_dimensions=[3, 80, 3, 80],
            numerical_equivalence="pending independent fair-esm comparison",
            implementation_status="adapted comparator; admission pending")

    def esm_tokenizer(self, batch):
        if len(batch) != 1:
            raise ValueError("DCTdomain adapter supports batch one only")
        pid, sequence = batch[0]
        validate_sequence(sequence, 500)
        tokens = self.tokenizer(sequence, return_tensors="pt", padding=False,
                                truncation=False, return_special_tokens_mask=True)
        validate_token_mapping(sequence, tokens["input_ids"][0].numpy(),
            tokens["attention_mask"][0].numpy(),
            tokens["special_tokens_mask"][0].numpy(), self.tokenizer)
        return [pid], [sequence], tokens["input_ids"]

    def esm_encoder(self, tokens, repr_layers, return_contacts):
        """Single eager forward; identical contact postprocessing to HF API."""
        torch = self.torch
        if tuple(repr_layers) != LAYERS or return_contacts is not True:
            raise ValueError("Only the pinned DCTdomain layer/contact contract is supported")
        if self.encoder.closed or tokens.device.type != "cpu":
            raise ValueError("Closed encoder or unsupported device")
        if tokens.ndim != 2 or tokens.shape[0] != 1 or not 3 <= tokens.shape[1] <= 502:
            raise ValueError("Invalid single-window token shape")
        if (int(tokens[0, 0]) != self.tokenizer.cls_token_id or
                int(tokens[0, -1]) != self.tokenizer.eos_token_id):
            raise ValueError("Missing boundary special tokens")
        if (tokens == self.tokenizer.pad_token_id).any():
            raise ValueError("Padding is not supported in the single-sequence method")
        attention_mask = tokens.ne(self.tokenizer.pad_token_id).long()
        captured = {}
        handles = []
        def capture(layer):
            def hook(module, args, output):
                if layer in captured:
                    raise ValueError("Repeated layer capture")
                captured[layer] = output[0] if isinstance(output, tuple) else output
            return hook
        try:
            for layer in LAYERS:
                handles.append(self.model.encoder.layer[layer - 1].register_forward_hook(capture(layer)))
            # Upstream stitching mutates returned slices outside this scope.
            # no_grad matches fair-esm here; inference tensors forbid that edit.
            with torch.no_grad():
                output = self.model(input_ids=tokens, attention_mask=attention_mask,
                    output_hidden_states=False, output_attentions=True, return_dict=True)
                if set(captured) != set(LAYERS) or len(output.attentions) != 30:
                    raise ValueError("Missing intermediate representations or attention heads")
                for hidden in captured.values():
                    if tuple(hidden.shape) != (1, tokens.shape[1], 640) or not torch.isfinite(hidden).all():
                        raise ValueError("Invalid intermediate representation")
                # Exactly the two padding masks in EsmModel.predict_contacts.
                attention = torch.stack(output.attentions, dim=1)
                if tuple(attention.shape) != (1, 30, 20, tokens.shape[1], tokens.shape[1]):
                    raise ValueError("Unexpected layer/head attention layout")
                attention *= attention_mask.unsqueeze(1).unsqueeze(2).unsqueeze(3)
                attention *= attention_mask.unsqueeze(1).unsqueeze(2).unsqueeze(4)
                contacts = self.model.contact_head(tokens, attention)
                if (tuple(contacts.shape) != (1, tokens.shape[1] - 2, tokens.shape[1] - 2)
                        or not torch.isfinite(contacts).all() or
                        (contacts < 0).any() or (contacts > 1).any()):
                    raise ValueError("Invalid actual-residue contact map")
                return {"representations": dict(captured), "contacts": contacts}
        finally:
            for handle in handles:
                handle.remove()
            captured.clear()

    def encode_with_upstream(self, sequence, pid, embedding_class):
        """Reuse exact upstream500/200 split and stitching; no algorithm rewrite."""
        import numpy as np
        validate_sequence(sequence, self.encoder.config["max_residues"])
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,48}", pid):
            raise ValueError("Use a short safe identifier for upstream temporary filenames")
        source = inspect.getsourcefile(embedding_class)
        if not source or sha256_file(source) != EMBEDDING_SHA256:
            raise ValueError("The pinned unchanged DCTdomain Embedding class is required")
        embedding = embedding_class(pid=pid, seq=sequence)
        embedding.embed_seq(self, "cpu", list(LAYERS), 500)
        if tuple(embedding.embed) != LAYERS:
            raise ValueError("Layer insertion order changed, affecting DCT concatenation")
        for matrix in embedding.embed.values():
            if matrix.shape != (len(sequence), 640) or not np.isfinite(matrix).all():
                raise ValueError("Invalid stitched representations")
        if embedding.contacts.shape != (len(sequence), len(sequence)) or not np.isfinite(embedding.contacts).all():
            raise ValueError("Invalid stitched contact map")
        return embedding

    def close(self):
        self.model = None
        self.encoder.close()


def hf_esm2_to_fair_state(state, num_layers=30):
    """Inverse of HF's official ESM2 converter, using existing tensor objects.

    Source: transformers commit fe861e578f50dc9c06de33cd361d2f625017e624,
    src/transformers/models/esm/convert_esm.py (Apache-2.0). No transposition,
    quantization, regenerated rotary frequency, or checkpoint serialization.
    Input is the raw safetensors mapping, not an EsmModel state missing LM keys.
    Caller constructs pinned official fair-esm ESM2 and load_state_dict(strict=True).
    This mapping alone is not a numerical equivalence result.
    """
    if type(num_layers) is not int or num_layers < 1:
        raise ValueError("Invalid layer count")
    mapped, consumed = {}, set()
    def assign(target, source):
        if source not in state:
            raise ValueError(f"Missing source tensor: {source}")
        mapped[target] = state[source]
        consumed.add(source)
    assign("embed_tokens.weight", "esm.embeddings.word_embeddings.weight")
    for suffix in ("weight", "bias"):
        assign(f"emb_layer_norm_after.{suffix}", f"esm.encoder.emb_layer_norm_after.{suffix}")
        assign(f"contact_head.regression.{suffix}", f"esm.contact_head.regression.{suffix}")
        assign(f"lm_head.dense.{suffix}", f"lm_head.dense.{suffix}")
        assign(f"lm_head.layer_norm.{suffix}", f"lm_head.layer_norm.{suffix}")
    assign("lm_head.bias", "lm_head.bias")
    mapped["lm_head.weight"] = mapped["embed_tokens.weight"]
    # HF safetensors usually removes these tied aliases. If present, verify ties.
    for alias, target in (("lm_head.decoder.weight", "embed_tokens.weight"),
                          ("lm_head.decoder.bias", "lm_head.bias")):
        if alias in state:
            if not state[alias].equal(mapped[target]):
                raise ValueError(f"Unexpected untied tensor: {alias}")
            consumed.add(alias)
    pairs = {
        "self_attn.q_proj": "attention.self.query",
        "self_attn.k_proj": "attention.self.key",
        "self_attn.v_proj": "attention.self.value",
        "self_attn.out_proj": "attention.output.dense",
        "self_attn_layer_norm": "attention.LayerNorm",
        "final_layer_norm": "LayerNorm", "fc1": "intermediate.dense", "fc2": "output.dense",
    }
    for layer in range(num_layers):
        for target, source in pairs.items():
            for suffix in ("weight", "bias"):
                assign(f"layers.{layer}.{target}.{suffix}", f"esm.encoder.layer.{layer}.{source}.{suffix}")
        assign(f"layers.{layer}.self_attn.rot_emb.inv_freq",
               f"esm.encoder.layer.{layer}.attention.self.rotary_embeddings.inv_freq")
    unused = set(state) - consumed
    if unused:
        raise ValueError(f"Unexpected unmapped source tensors: {sorted(unused)}")
    return mapped
