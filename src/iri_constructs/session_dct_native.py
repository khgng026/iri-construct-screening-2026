"""Pinned native DCTdomain integration; all model imports occur on invocation.

Upstream ESM2, Embedding, Fingerprint and RecCut remain unchanged. The local
facades validate inputs/outputs and observe temporary files, never fit parameters.
"""
from __future__ import annotations
from contextlib import contextmanager
import gc
import importlib.util
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys
import tempfile

from .session_store import canonical, digest, file_hash
from .session_models import validate_sequence
from .session_dctdomain_adapter import MODEL_ID, REVISION, WEIGHT_SHA256, LAYERS, hf_esm2_to_fair_state


def require(value, message):
    if not value:
        raise ValueError(message)


def verify_sources(repo, session):
    manifest = json.loads((session / "audit/dct_source_acquisition.json").read_text())
    roots, hashes = {}, {}
    for family in ("fair_esm", "dctdomain"):
        root = (repo / manifest[family]["root"]).resolve()
        require(root.is_relative_to((session / "vendor").resolve()), "Vendor source outside owned session")
        roots[family] = root
        for row in manifest[family]["files"]:
            path = (root / row["path"]).resolve()
            require(path.is_relative_to(root) and path.stat().st_size == row["bytes"] and file_hash(path) == row["sha256"], "Pinned vendor source changed")
            hashes[str(path.relative_to(repo))] = row["sha256"]
    header = manifest["compatibility_header"]
    require(file_hash(repo / header["path"]) == header["sha256"], "Compatibility header changed")
    hashes[header["path"]] = header["sha256"]
    return roots, hashes


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class NativeDCTModel:
    """Official fair-esm forward using inverse-mapped local HF safetensors."""
    def __init__(self, fair_root, threads=4):
        require(type(threads) is int and 1 <= threads <= 4, "CPU thread bound exceeded")
        require("esm" not in sys.modules, "Start a fresh process before selecting native fair-esm")
        sys.path.insert(0, str(fair_root))
        import esm
        from esm.model.esm2 import ESM2
        import torch
        from huggingface_hub import hf_hub_download
        from safetensors.torch import load_file
        require(Path(esm.__file__).resolve() == (fair_root / "esm/__init__.py").resolve(), "Wrong esm namespace loaded")
        self.torch = torch
        torch.set_num_threads(threads)
        torch.manual_seed(0)
        torch.use_deterministic_algorithms(True)
        files = {}
        paths = {}
        for name in ("model.safetensors", "config.json", "vocab.txt"):
            path = Path(hf_hub_download(MODEL_ID, name, revision=REVISION, local_files_only=True))
            paths[name] = path
            files[name] = {"sha256": file_hash(path), "bytes": path.stat().st_size}
        require(files["model.safetensors"]["sha256"] == WEIGHT_SHA256, "Wrong local checkpoint")
        config = json.loads(paths["config.json"].read_text())
        expected = {"num_hidden_layers": 30, "hidden_size": 640, "num_attention_heads": 20,
                    "intermediate_size": 2560, "layer_norm_eps": 1e-5, "position_embedding_type": "rotary",
                    "token_dropout": True, "emb_layer_norm_before": False,
                    "attention_probs_dropout_prob": 0.0, "hidden_dropout_prob": 0.0}
        require(all(config.get(k) == v for k, v in expected.items()), "Unexpected checkpoint architecture/configuration")
        # No initialized second set of model weights and no serialized conversion.
        with torch.device("meta"):
            self.model = ESM2(num_layers=30, embed_dim=640, attention_heads=20,
                              alphabet="ESM-1b", token_dropout=True)
        state = load_file(str(paths["model.safetensors"]), device="cpu")
        require(all(not value.is_floating_point() or value.dtype == torch.float32 for value in state.values()), "Raw checkpoint is not FP32")
        mapped = hf_esm2_to_fair_state(state, num_layers=30)
        result = self.model.load_state_dict(mapped, strict=True, assign=True)
        require(not result.missing_keys and not result.unexpected_keys, "Strict native mapping incomplete")
        require(all(not tensor.is_meta and tensor.device.type == "cpu" for tensor in
                    [*self.model.parameters(), *self.model.buffers()]), "Unloaded/meta native tensor")
        norms = sorted({module.eps for module in self.model.modules() if isinstance(module, torch.nn.LayerNorm)})
        require(norms == [1e-5], "Native normalization epsilon differs")
        self.model.eval()
        self.alphabet = self.model.alphabet
        self.converter = self.alphabet.get_batch_converter()
        require(self.alphabet.all_toks == paths["vocab.txt"].read_text().splitlines(), "Native/HF token vocabulary differs")
        require(config["pad_token_id"] == self.alphabet.padding_idx and config["mask_token_id"] == self.alphabet.mask_idx,
                "Checkpoint special token identity differs")
        state_count = len(state)
        del mapped, state
        gc.collect()
        self.metadata = {"backend": "native_fair_esm", "model_id": MODEL_ID, "revision": REVISION,
            "checkpoint_files": files, "checkpoint_config": config, "strict_mapped_tensor_count": state_count,
            "token_vocabulary": list(self.alphabet.all_toks), "esm_source_version": esm.__version__,
            "esm_source_path": str(Path(esm.__file__).resolve()), "device": "cpu", "dtype": "float32",
            "threads": threads, "seed": 0, "layers": list(LAYERS), "batch_size": 1,
            "native_layer_norm_eps": norms, "native_token_dropout": self.model.token_dropout,
            "representation": "native blocks15/21 before final norm; original learned all-head contacts",
            "label": "Native fair-esm with inverse-mapped pinned HF checkpoint; no separately downloaded original checkpoint",
            "packages": {p: importlib.metadata.version(p) for p in ("torch", "numpy", "scipy", "safetensors", "huggingface-hub")}}

    def esm_tokenizer(self, batch):
        require(len(batch) == 1, "Native DCT requires batch one")
        sequence = batch[0][1]
        validate_sequence(sequence, 500)
        labels, strings, tokens = self.converter(batch)
        expected = [self.alphabet.cls_idx] + [self.alphabet.get_idx(a) for a in sequence] + [self.alphabet.eos_idx]
        require(tokens.shape == (1, len(sequence) + 2) and tokens[0].tolist() == expected, "Native residue token mapping changed")
        return labels, strings, tokens

    def esm_encoder(self, tokens, repr_layers, return_contacts):
        require(tuple(repr_layers) == LAYERS and return_contacts is True, "DCT representation contract changed")
        with self.torch.no_grad():
            result = self.model(tokens, repr_layers=list(LAYERS), return_contacts=True)
        # Drop LM logits and all-head attentions after official contact prediction.
        return {"representations": result["representations"], "contacts": result["contacts"]}

    def close(self):
        self.model = None
        gc.collect()


class ObservedModel:
    """Observe raw precision/windows on either backend without altering tensors."""
    def __init__(self, backend, capture=False):
        self.backend, self.capture = backend, capture
        self.torch = backend.torch
        self.windows, self.arrays = [], {}
        counts = {}
        for name, parameter in backend.model.named_parameters():
            if parameter.is_floating_point():
                require(parameter.dtype == self.torch.float32 and parameter.device.type == "cpu", "Wrong native floating parameter dtype/device")
                key = str(parameter.dtype) + "@" + str(parameter.device)
                counts[key] = counts.get(key, 0) + parameter.numel()
        require(counts and not backend.model.training, "Missing parameters or model in training mode")
        self.parameter_evidence = counts

    def esm_tokenizer(self, batch):
        return self.backend.esm_tokenizer(batch)

    def esm_encoder(self, tokens, repr_layers, return_contacts):
        torch = self.torch
        require(not torch.is_autocast_enabled("cpu"), "DCT autocast is forbidden")
        output = self.backend.esm_encoder(tokens, repr_layers, return_contacts)
        index = len(self.windows)
        evidence = {"window_index": index, "token_ids_sha256": digest(canonical(tokens[0].tolist())), "residue_count": tokens.shape[1] - 2,
                    "raw_dtype": "torch.float32", "device": "cpu", "autocast": False}
        if self.capture:
            evidence["token_ids"] = tokens[0].tolist()
        require(set(output["representations"]) == set(LAYERS), "Wrong returned layers")
        values = {f"layer{layer}": output["representations"][layer] for layer in LAYERS}
        values["contacts"] = output["contacts"]
        for name, tensor in values.items():
            expected = (1, tokens.shape[1] - 2, tokens.shape[1] - 2) if name == "contacts" else (1, tokens.shape[1], 640)
            require(tuple(tensor.shape) == expected and tensor.dtype == torch.float32 and tensor.device.type == "cpu"
                    and not tensor.requires_grad and torch.isfinite(tensor).all(), "Invalid raw DCT forward tensor")
            if self.capture:
                self.arrays[f"window{index}_{name}"] = tensor.detach().numpy().copy()
        self.windows.append(evidence)
        return output


def window_lengths(length):
    require(type(length) is int and length > 0, "Invalid input length")
    return [length] if length <= 500 else [min(500, length - start) for start in range(0, length, 300) if length - start > 200]


def parse_domains(domains, length):
    """Preserve interval order and discontiguous sets; exclude only appended global."""
    require(bool(domains) and len(set(domains)) == len(domains), "Empty/repeated domains")
    appended = len(domains) > 1 and domains[-1] == f"1-{length}"
    predicted = domains[:-1] if appended else domains
    require(len(predicted) > 0, "Missing predicted domains")
    parsed, partition = [], set()
    for number, text in enumerate(domains):
        spans, residues = [], set()
        for item in text.split(","):
            require(re.fullmatch(r"[1-9][0-9]*-[1-9][0-9]*", item) is not None, "Malformed domain interval")
            first, last = map(int, item.split("-"))
            require(1 <= first <= last <= length, "Domain endpoint outside input")
            current = set(range(first - 1, last))
            require(not current & residues, "Repeated residues within discontinuous domain")
            residues |= current
            spans.append({"start0": first - 1, "end0": last})
        whole_extra = appended and number == len(domains) - 1
        if not whole_extra:
            require(not partition & residues, "Predicted domains overlap")
            partition |= residues
        parsed.append({"domain": text, "segments": spans, "residue_count": len(residues),
                       "appended_whole_fingerprint": whole_extra})
    require(partition == set(range(length)), "Predicted domains do not partition input")
    return parsed


def ce_evidence(path, length, save_graph=False):
    """Read only the CON line; raw sequence stays inside owned temporary storage."""
    import numpy as np
    pairs, values = [], []
    with Path(path).open() as source:
        for line in source:
            if line.startswith("CON"):
                for item in line[3:].strip().split(","):
                    first, last, value = item.split()
                    first, last = int(first), int(last)
                    require(0 <= first < last < length and last - first >= 5, "Invalid selected contact endpoints")
                    require(re.fullmatch(r"[0-9]+\.[0-9]{6}", value) is not None, "Contact precision changed")
                    pairs.append((first, last))
                    values.append(float(value))
    pair_array = np.asarray(pairs, dtype=np.int32).reshape(-1, 2)
    rounded = np.asarray(values, dtype=np.float64)
    require(np.isfinite(rounded).all() and ((0 <= rounded) & (rounded <= 1)).all(), "Invalid contact probabilities")
    weights = (rounded * 100 + 0.5).astype(np.int32)
    graph = np.zeros((length, length), dtype=np.int32)
    for (i, j), value in zip(pair_array, weights):
        graph[i, j] = graph[j, i] = value
    for i in range(length):
        for gap in (1, 2, 3):
            if i + gap < length:
                graph[i, i + gap] = graph[i + gap, i] = 100
    arrays = {"ce_pairs": pair_array, "ce_rounded_values": rounded, "ce_graph_weights": weights}
    if save_graph:
        arrays["integer_graph"] = graph
    return {"selected_contact_count": len(pairs), "integer_graph_sha256": digest(graph.astype("<i4").tobytes()),
            "graph_shape": [length, length], "mandatory_adjacent_edges": "separation1,2,3;weight100"}, arrays


@contextmanager
def isolated_working_directory(base):
    base.mkdir(parents=True, exist_ok=True)
    old = Path.cwd()
    with tempfile.TemporaryDirectory(prefix="dct-", dir=base) as directory:
        try:
            os.chdir(directory)
            yield Path(directory)
        finally:
            os.chdir(old)


def fingerprint_embedding(embedding, fingerprint_class, temporary_root, capture=False):
    import numpy as np
    fp = fingerprint_class(pid=embedding.pid, seq=embedding.seq, embed=embedding.embed, contacts=embedding.contacts)
    original_writece = fp.writece
    evidence, contact_arrays = {}, {}
    def observe_writece(path, threshold):
        original_writece(path, threshold)
        record, arrays = ce_evidence(path, len(embedding.seq), save_graph=capture)
        evidence.update(record)
        if capture:
            contact_arrays.update(arrays)
    fp.writece = observe_writece
    try:
        with isolated_working_directory(temporary_root):
            fp.reccut(2.6)
    except Exception as exc:
        exc.iri_dct_stage = "original_reccut_or_contact_observation"
        raise
    finally:
        # Remove the observing closure's bound-method cycle before the next input.
        del fp.writece
    domains = list(fp.domains)
    try:
        parsed = parse_domains(domains, len(embedding.seq))
    except Exception as exc:
        exc.iri_dct_stage = "domain_partition_validation"
        raise
    try:
        with np.errstate(invalid="raise", divide="raise", over="raise"):
            fp.quantize([3, 80, 3, 80])
    except Exception as exc:
        exc.iri_dct_stage = "original_quantization_with_precast_error_check"
        raise
    require(fp.domains == domains and list(fp.quants) == domains, "Upstream quantization silently altered domains")
    fingerprints = np.stack([fp.quants[domain] for domain in domains])
    require(fingerprints.shape == (len(domains), 480) and np.issubdtype(fingerprints.dtype, np.integer)
            and ((0 <= fingerprints) & (fingerprints <= 127)).all(), "Invalid discrete fingerprint")
    return {"domains": parsed, "contact_evidence": evidence, "fingerprint_count": len(domains)}, {
        "fingerprints": fingerprints.astype(np.int16), **contact_arrays}


def boundary_metrics(domain_records, spans):
    predicted = [{x for s in d["segments"] for x in range(s["start0"], s["end0"])}
                 for d in domain_records if not d["appended_whole_fingerprint"]]
    rows = []
    for span in spans:
        truth = set(range(span["start0"], span["end0"]))
        comparisons = [{"predicted_domain_index": i, "intersection": len(truth & p),
                        "iou": len(truth & p) / len(truth | p), "component_coverage": len(truth & p) / len(truth),
                        "predicted_domain_precision": len(truth & p) / len(p)} for i, p in enumerate(predicted)]
        rows.append({"accession": span["accession"], "best_iou": max(x["iou"] for x in comparisons),
                     "best_single_domain_component_coverage": max(x["component_coverage"] for x in comparisons),
                     "domain_comparisons": comparisons})
    return {"predicted_domain_count": len(predicted), "components": rows,
            "semantics": "Residue-set unions, not interval hulls; appended whole fingerprint excluded; structural domains need not match constructed components."}


def exhaustive_distances(query, references):
    """Complete widened raw L1, with a deterministic native domain-pair order."""
    import numpy as np
    original = np.asarray(query)
    require(np.issubdtype(original.dtype, np.integer) and ((0 <= original) & (original <= 127)).all(), "Invalid discrete query values")
    q = original.astype(np.int64)
    require(q.ndim == 2 and q.shape[1] == 480 and q.shape[0] > 0, "Invalid query fingerprints")
    rows, minima = [], {}
    accessions = sorted(references)
    for accession_index, accession in enumerate(accessions):
        original = np.asarray(references[accession])
        require(np.issubdtype(original.dtype, np.integer) and ((0 <= original) & (original <= 127)).all(), "Invalid discrete reference values")
        reference = original.astype(np.int64)
        require(reference.ndim == 2 and reference.shape[1] == 480 and len(reference) > 0, "Invalid reference fingerprints")
        distance = np.abs(q[:, None, :] - reference[None, :, :]).sum(axis=2, dtype=np.int64)
        minima[accession] = int(distance.min())
        for qi in range(len(q)):
            for ri in range(len(reference)):
                rows.append((accession_index, qi, ri, int(distance[qi, ri])))
    rows.sort(key=lambda row: (row[3], accessions[row[0]], row[1], row[2]))
    return minima, np.asarray(rows, dtype=np.int64).reshape(-1, 4)
