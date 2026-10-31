"""Isolated, offline protein representation adapters for the 2026-10-31 session.

No generation, training, downloading, or changes to the frozen T05 encoder. The
caller owns resource admission and supplies a pinned, already acquired checkpoint.
``encode`` returns one float32 row per actual residue, without pooling.
"""
from __future__ import annotations

import gc
import hashlib
import importlib.metadata
import inspect
from pathlib import Path
import re
from time import perf_counter


MODEL_FAMILIES = {
    "facebook/esm2_t6_8M_UR50D": "esm2",
    "facebook/esm2_t30_150M_UR50D": "esm2",
    "facebook/esm2_t33_650M_UR50D": "esm2",
    "biohub/ESMC-300M": "esmc",
    "biohub/ESMC-600M": "esmc",
    "biohub/esm3-sm-open-v1": "esm3",
}
CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")
ESM3_WEIGHT_FILE = "data/weights/esm3_sm_open_v1.pth"


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_config(config):
    """Resolve explicit defaults and reject unsupported or ambiguous configurations."""
    c = dict(config)
    model_id = c.get("model_id")
    if model_id not in MODEL_FAMILIES:
        raise ValueError("Model is outside the authorized checkpoint list")
    family = MODEL_FAMILIES[model_id]
    if c.get("family", family) != family:
        raise ValueError("Model family does not match model identity")
    c["family"] = family
    if not re.fullmatch(r"[0-9a-f]{40}", c.get("revision", "")):
        raise ValueError("An immutable 40-character checkpoint revision is required")
    if c.get("local_files_only", True) is not True:
        raise ValueError("Adapters are offline; acquisition belongs to the lead")
    c["local_files_only"] = True
    for key, default in (("device", "cpu"), ("dtype", "float32"),
                         ("threads", 4), ("seed", 0), ("max_residues", 1022),
                         ("layer", -1), ("tokenizer_truncation", False)):
        c.setdefault(key, default)
    if c["device"] not in ("cpu", "mps") or c["dtype"] != "float32":
        raise ValueError("Only separately admitted CPU/MPS float32 configurations are supported")
    if type(c["threads"]) is not int or not 1 <= c["threads"] <= 4:
        raise ValueError("CPU threads must be an integer between one and four")
    if type(c["max_residues"]) is not int or not 1 <= c["max_residues"] <= 1022:
        raise ValueError("Shared session maximum is 1022 residues")
    if type(c["layer"]) is not int or c["layer"] < -1:
        raise ValueError("Layer must be -1 (final) or a nonnegative integer")
    if c["tokenizer_truncation"] is not False:
        raise ValueError("Truncation is forbidden")
    if family in ("esmc", "esm3") and c.get("attn_implementation", "sdpa") != "sdpa":
        raise ValueError("This adapter uses the installed official SDPA route")
    if family == "esm3" and c["device"] != "cpu":
        raise ValueError("ESM3 MPS is not admitted by this adapter")
    if family == "esm3":
        c.setdefault("representation_variant", "official_embeddings_prenorm")
        if c["representation_variant"] not in ("official_embeddings_prenorm", "post_final_norm"):
            raise ValueError("Unknown ESM3 representation variant")
    return c


def validate_sequence(sequence, max_residues):
    # Reject, never normalize or trim. Error messages deliberately omit sequences.
    if not isinstance(sequence, str) or not sequence or len(sequence) > max_residues:
        raise ValueError("Invalid sequence length or type")
    if set(sequence) - CANONICAL:
        raise ValueError("Only reviewed canonical amino-acid strings are accepted")


def validate_token_mapping(sequence, input_ids, attention_mask, special_mask, tokenizer):
    """Check both cardinality and ordered residue identities, not just token count."""
    import numpy as np
    ids = np.asarray(input_ids)
    attention = np.asarray(attention_mask, dtype=bool)
    special = np.asarray(special_mask, dtype=bool)
    if ids.ndim != 1 or ids.shape != attention.shape or ids.shape != special.shape:
        raise ValueError("Incompatible token/mask shapes")
    mask = attention & ~special
    positions = np.flatnonzero(mask)
    if len(ids) != len(sequence) + 2 or positions.tolist() != list(range(1, len(sequence) + 1)):
        raise ValueError("Tokenizer did not preserve contiguous one-token-per-residue mapping")
    expected = tokenizer.convert_tokens_to_ids(list(sequence))
    if tokenizer.unk_token_id in ids or ids[mask].tolist() != expected:
        raise ValueError("Tokenizer changed residue identities or produced an unknown token")
    if ids[0] != tokenizer.cls_token_id or ids[-1] != tokenizer.eos_token_id:
        raise ValueError("Unexpected start/end special tokens")
    return mask


def _package_versions():
    result = {}
    for name in ("torch", "transformers", "numpy", "esm", "safetensors", "tokenizers", "huggingface-hub"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    return result


class ResidueEncoder:
    """One model copy, batch one, with exact residue mapping and bounded state capture."""

    def __init__(self, config):
        self.config = validate_config(config)
        self.closed = False
        self.last_encode_metadata = None
        self.model = None
        self._source_files = {}
        self._weight_files = {}
        started = perf_counter()
        import torch
        self.torch = torch
        torch.set_num_threads(self.config["threads"])
        torch.manual_seed(self.config["seed"])
        torch.use_deterministic_algorithms(True)
        if self.config["device"] == "mps" and not torch.backends.mps.is_available():
            raise RuntimeError("MPS requested but unavailable; no silent device fallback")
        if self.config["family"] == "esm2":
            self._load_esm2()
        elif self.config["family"] == "esmc":
            self._load_esmc()
        else:
            self._load_esm3()
        self.model.eval()
        self.num_layers = (48 if self.config["family"] == "esm3"
                           else int(self.model.config.num_hidden_layers))
        layer = self.config["layer"]
        if layer > self.num_layers:
            self.close()
            raise ValueError("Requested representation layer exceeds model depth")
        self.layer = self.num_layers if layer == -1 else layer
        self.load_seconds = perf_counter() - started
        self.metadata = {
            "adapter_version": 1, "config": self.config,
            "model_id": self.config["model_id"], "revision": self.config["revision"],
            "device": self.config["device"], "dtype": self.config["dtype"],
            "representation": (("official_output_embeddings_pre_final_norm"
                                if self.config["family"] == "esm3" and self.config["representation_variant"] == "official_embeddings_prenorm"
                                else "last_layer_post_final_norm") if self.layer == self.num_layers
                               else f"hidden_state_after_{self.layer}_blocks"),
            "layer": self.layer, "num_layers": self.num_layers,
            "pooling": "none; actual-residue states returned",
            "batch_size": 1, "max_residues": self.config["max_residues"],
            "tokenizer": type(self.tokenizer).__module__ + "." + type(self.tokenizer).__name__,
            "tokenizer_vocab_sha256": hashlib.sha256(repr(sorted(self.tokenizer.get_vocab().items())).encode()).hexdigest(),
            "adapter_source_sha256": sha256_file(__file__),
            "library_source_sha256": self._source_files,
            "checkpoint_files": self._weight_files,
            "packages": _package_versions(), "load_seconds": self.load_seconds,
            "extra_tracks": "absent; official missing-track defaults" if self.config["family"] == "esm3" else None,
            "attention_implementation": self._attention_implementation,
            "deterministic_algorithms": True, "truncated": False,
        }

    def _record_source(self, obj):
        path = inspect.getsourcefile(obj)
        if path:
            self._source_files[str(Path(path).resolve())] = sha256_file(path)

    def _cached_file(self, filename):
        from huggingface_hub import hf_hub_download
        path = Path(hf_hub_download(self.config["model_id"], filename,
                    revision=self.config["revision"], local_files_only=True))
        digest = sha256_file(path)
        expected = self.config.get("expected_file_sha256", {}).get(filename)
        if expected and digest != expected:
            raise ValueError(f"Checkpoint checksum mismatch for {filename}")
        self._weight_files[filename] = {"sha256": digest, "bytes": path.stat().st_size}
        return path

    def _load_esm2(self):
        from transformers import EsmModel, EsmTokenizer
        self._record_source(EsmModel)
        self._record_source(EsmTokenizer)
        self._cached_file("model.safetensors")
        self._cached_file("config.json")
        self._cached_file("vocab.txt")
        options = {"revision": self.config["revision"], "local_files_only": True}
        self.tokenizer = EsmTokenizer.from_pretrained(self.config["model_id"], **options)
        extra = {}
        if "attn_implementation" in self.config:
            extra["attn_implementation"] = self.config["attn_implementation"]
        self.model = EsmModel.from_pretrained(self.config["model_id"],
                         use_safetensors=True, add_pooling_layer=False, **options, **extra)
        self.model = self.model.to(device=self.config["device"], dtype=self.torch.float32).eval()
        self._attention_implementation = getattr(self.model.config, "_attn_implementation", "library_default")

    def _load_esmc(self):
        from esm.models.esmc.model import EsmcModel
        from esm.models.esmc.tokenizer import EsmcTokenizer
        self._record_source(EsmcModel)
        self._record_source(EsmcTokenizer)
        self._cached_file("model.safetensors")
        self._cached_file("config.json")
        self.tokenizer = EsmcTokenizer()
        self.model = EsmcModel.from_pretrained(self.config["model_id"],
                         revision=self.config["revision"], local_files_only=True,
                         device=self.config["device"], dtype=self.torch.float32,
                         attn_implementation="sdpa")
        self._attention_implementation = "sdpa"

    def _load_esm3(self):
        # Match the official factory architecture, but resolve one immutable local
        # file instead of its unpinned full-snapshot/optional-decoder download path.
        from accelerate import init_empty_weights
        from esm.models.esm3 import ESM3
        from esm.tokenization import TokenizerCollection
        from esm.tokenization.sequence_tokenizer import EsmSequenceTokenizer
        from esm.tokenization.structure_tokenizer import StructureTokenizer
        from esm.tokenization.ss_tokenizer import SecondaryStructureTokenizer
        from esm.tokenization.sasa_tokenizer import SASADiscretizingTokenizer
        from esm.tokenization.function_tokenizer import InterProQuantizedTokenizer
        from esm.tokenization.residue_tokenizer import ResidueAnnotationsTokenizer
        from esm.utils.constants import esm3 as constants
        self._record_source(ESM3)
        self._record_source(EsmSequenceTokenizer)
        path = self._cached_file(ESM3_WEIGHT_FILE)
        snapshot = path.parent.parent.parent
        self.tokenizer = EsmSequenceTokenizer()
        tokenizers = TokenizerCollection(
            sequence=self.tokenizer, structure=StructureTokenizer(),
            secondary_structure=SecondaryStructureTokenizer(kind="ss8"),
            sasa=SASADiscretizingTokenizer(), function=InterProQuantizedTokenizer(),
            # Explicit path prevents data_root() from requesting any snapshot.
            # This deferred table is never read for absent annotation tracks.
            residue_annotations=ResidueAnnotationsTokenizer(csv_path=str(snapshot / constants.RESID_CSV)))
        def forbidden_decoder(*args, **kwargs):
            raise RuntimeError("Decoder/structure/function services are outside this sequence-only adapter")
        with init_empty_weights():
            self.model = ESM3(d_model=1536, n_heads=24, v_heads=256, n_layers=48,
                structure_encoder_fn=forbidden_decoder, structure_decoder_fn=forbidden_decoder,
                function_decoder_fn=forbidden_decoder, tokenizers=tokenizers).eval()
        state = self.torch.load(path, map_location="cpu", weights_only=True, mmap=True)
        self.model.load_state_dict(state, strict=True, assign=True)
        del state
        self.model = self.model.to(device="cpu", dtype=self.torch.float32).eval()
        self._attention_implementation = "sdpa"

    def _synchronize(self):
        if self.config["device"] == "mps":
            self.torch.mps.synchronize()

    def _capture_module(self):
        """A single hook avoids retaining every layer during sensitivity runs."""
        if self.config["family"] == "esm2":
            return self.model.embeddings if self.layer == 0 else self.model.encoder.layer[self.layer - 1]
        if self.config["family"] == "esmc":
            return self.model.embed if self.layer == 0 else self.model.transformer.blocks[self.layer - 1]
        if self.layer == self.num_layers:
            return self.model.transformer.norm
        return self.model.encoder if self.layer == 0 else self.model.transformer.blocks[self.layer - 1]

    def encode(self, sequence):
        if self.closed:
            raise RuntimeError("Encoder is closed")
        import numpy as np
        torch = self.torch
        validate_sequence(sequence, self.config["max_residues"])
        started = perf_counter()
        inputs = self.tokenizer(sequence, return_tensors="pt", padding=False,
                                truncation=False, return_special_tokens_mask=True)
        special = inputs.pop("special_tokens_mask")
        ids = inputs["input_ids"][0].numpy()
        attention = inputs["attention_mask"][0].numpy()
        mask = validate_token_mapping(sequence, ids, attention, special[0].numpy(), self.tokenizer)
        inputs = {key: value.to(self.config["device"]) for key, value in inputs.items()
                  if key in ("input_ids", "attention_mask")}
        token_seconds = perf_counter() - started
        captured = []
        handle = None
        requires_hook = self.layer != self.num_layers or (
            self.config["family"] == "esm3" and self.config["representation_variant"] == "post_final_norm")
        if requires_hook:
            def capture(module, args, output):
                captured.append(output[0] if isinstance(output, tuple) else output)
            handle = self._capture_module().register_forward_hook(capture)
        try:
            self._synchronize()
            started = perf_counter()
            with torch.inference_mode():
                if self.config["family"] == "esm3":
                    # Match ESM3.logits()'s sequence-only conditioning exactly.
                    output = self.model(sequence_tokens=inputs["input_ids"],
                        average_plddt=torch.tensor(1.0, device=self.config["device"]),
                        output_attentions=False)
                elif self.config["family"] == "esmc":
                    output = self.model(**inputs, output_hidden_states=False, output_attentions=False,
                                        return_dict=True, compute_sae=False)
                else:
                    output = self.model(**inputs, output_hidden_states=False, output_attentions=False,
                                        return_dict=True)
                self._synchronize()
                forward_seconds = perf_counter() - started
                if requires_hook and len(captured) != 1:
                    raise ValueError("Expected exactly one representation capture")
                hidden = (captured[0] if requires_hook else output.embeddings
                          if self.config["family"] == "esm3" else output.last_hidden_state)
                if hidden.ndim != 3 or hidden.shape[0] != 1 or hidden.shape[1] != len(sequence) + 2:
                    raise ValueError("Model hidden-state shape does not match exact token mapping")
                if hidden.requires_grad or not torch.isfinite(hidden).all():
                    raise ValueError("Nonfinite hidden states or unexpected gradient tracking")
                residue = hidden[0, torch.as_tensor(mask, device=hidden.device)]
                result = residue.detach().to(device="cpu", dtype=torch.float32).numpy().copy()
                if result.shape[0] != len(sequence) or not np.isfinite(result).all():
                    raise ValueError("Invalid actual-residue representation")
                self.last_encode_metadata = {
                    "residue_count": len(sequence), "token_count": len(ids),
                    "hidden_shape": list(hidden.shape), "residue_shape": list(result.shape),
                    "tokenize_seconds": token_seconds, "forward_seconds": forward_seconds,
                    "device": self.config["device"], "dtype": "float32", "layer": self.layer,
                    "truncated": False, "model_training": self.model.training,
                    "inference_mode": True,
                }
                return result
        finally:
            if handle is not None:
                handle.remove()
            captured.clear()

    def verify_official_equivalence(self, sequence, tolerance=1e-6):
        """ESM3 admission check against its public encode/logits API, no generation.

        This runs two short forwards on the same loaded model. The caller must
        include it inside the resource-guarded probe, not alongside another job.
        """
        import numpy as np
        if (self.config["family"] != "esm3" or self.layer != self.num_layers or
                self.config["representation_variant"] != "official_embeddings_prenorm"):
            raise ValueError("Equivalence check is for primary final-layer ESM3 embeddings")
        from esm.sdk.api import ESMProtein, LogitsConfig
        actual = self.encode(sequence)
        with self.torch.inference_mode():
            tokens = self.model.encode(ESMProtein(sequence=sequence))
            if any(getattr(tokens, name) is not None for name in
                   ("structure", "secondary_structure", "sasa", "function",
                    "residue_annotations", "coordinates")):
                raise ValueError("Official sequence-only encoding unexpectedly supplied extra tracks")
            official = self.model.logits(tokens, LogitsConfig(return_embeddings=True))
            expected = official.embeddings[0, 1:-1].detach().float().cpu().numpy()
        if expected.shape != actual.shape or not np.isfinite(expected).all():
            raise ValueError("Official ESM3 representation shape/finiteness mismatch")
        error = float(np.max(np.abs(actual.astype(np.float64) - expected.astype(np.float64))))
        result = {"interface": "ESM3.encode(ESMProtein(sequence)) then logits(return_embeddings=True)",
                  "max_abs_error": error, "tolerance": tolerance, "shape": list(actual.shape),
                  "extra_tracks_absent": True, "passed": error <= tolerance}
        if not result["passed"]:
            raise ValueError(f"ESM3 public-interface numerical disagreement: max_abs_error={error}")
        return result

    def close(self):
        self.model = None
        self.closed = True
        gc.collect()
        if self.config["device"] == "mps":
            self.torch.mps.empty_cache()


def load_encoder(config):
    return ResidueEncoder(config)
