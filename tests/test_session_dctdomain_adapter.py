"""Lightweight contract and exact key-map checks; no model/network required."""
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from iri_constructs.session_dctdomain_adapter import (
    MODEL_ID, REVISION, WEIGHT_SHA256, validate_dct_config, hf_esm2_to_fair_state,
    DCTdomainHFModel,
)


def source_state():
    # Distinct sentinels let the tests establish identity-preserving assignment.
    state = {}
    def put(key):
        state[key] = object()
    put("esm.embeddings.word_embeddings.weight")
    for stem in ("esm.encoder.emb_layer_norm_after", "esm.contact_head.regression",
                 "lm_head.dense", "lm_head.layer_norm"):
        for suffix in ("weight", "bias"):
            put(f"{stem}.{suffix}")
    put("lm_head.bias")
    for stem in ("attention.self.query", "attention.self.key", "attention.self.value",
                 "attention.output.dense", "attention.LayerNorm", "LayerNorm",
                 "intermediate.dense", "output.dense"):
        for suffix in ("weight", "bias"):
            put(f"esm.encoder.layer.0.{stem}.{suffix}")
    put("esm.encoder.layer.0.attention.self.rotary_embeddings.inv_freq")
    return state


class DCTConfigurationTests(unittest.TestCase):
    def test_frozen_identity_and_offline_eager_cpu_defaults(self):
        result = validate_dct_config({})
        self.assertEqual(result["model_id"], MODEL_ID)
        self.assertEqual(result["revision"], REVISION)
        self.assertEqual(result["attn_implementation"], "eager")
        self.assertEqual(result["expected_file_sha256"]["model.safetensors"], WEIGHT_SHA256)

    def test_different_method_configuration_rejected(self):
        for change in ({"model_id": "facebook/esm2_t6_8M_UR50D"}, {"revision": "main"},
                       {"device": "mps"}, {"dtype": "float16"}, {"local_files_only": False},
                       {"layer": -1}, {"attn_implementation": "sdpa"},
                       {"expected_file_sha256": {"model.safetensors": "0" * 64}}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate_dct_config(change)


class InverseConversionTests(unittest.TestCase):
    def test_preserves_rotary_contact_and_intermediate_tensor_objects(self):
        source = source_state()
        result = hf_esm2_to_fair_state(source, num_layers=1)
        expected = {
            "layers.0.self_attn.rot_emb.inv_freq": "esm.encoder.layer.0.attention.self.rotary_embeddings.inv_freq",
            "contact_head.regression.weight": "esm.contact_head.regression.weight",
            "layers.0.self_attn.q_proj.bias": "esm.encoder.layer.0.attention.self.query.bias",
            "layers.0.self_attn_layer_norm.weight": "esm.encoder.layer.0.attention.LayerNorm.weight",
            "layers.0.final_layer_norm.weight": "esm.encoder.layer.0.LayerNorm.weight",
            "layers.0.fc2.weight": "esm.encoder.layer.0.output.dense.weight",
        }
        for target, key in expected.items():
            self.assertIs(result[target], source[key])
        self.assertIs(result["lm_head.weight"], result["embed_tokens.weight"])
        self.assertEqual(len(result), len(source) + 1)

    def test_missing_contact_or_frequency_is_fatal(self):
        for key in ("esm.contact_head.regression.weight",
                    "esm.encoder.layer.0.attention.self.rotary_embeddings.inv_freq"):
            state = source_state()
            del state[key]
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, "Missing source tensor"):
                hf_esm2_to_fair_state(state, num_layers=1)

    def test_unknown_keys_are_not_silently_ignored(self):
        state = source_state()
        state["unexpected.weight"] = object()
        with self.assertRaisesRegex(ValueError, "Unexpected unmapped"):
            hf_esm2_to_fair_state(state, num_layers=1)

    def test_invalid_depth_and_incomplete_deeper_model_rejected(self):
        for depth in (0, -1, True):
            with self.subTest(depth=depth), self.assertRaises(ValueError):
                hf_esm2_to_fair_state(source_state(), num_layers=depth)
        with self.assertRaises(ValueError):
            hf_esm2_to_fair_state(source_state(), num_layers=2)


class ForwardContractTests(unittest.TestCase):
    def setUp(self):
        import torch
        self.torch = torch
        class TinyInterfaceModel(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.encoder = SimpleNamespace(layer=torch.nn.ModuleList(
                    [torch.nn.Identity() for _ in range(30)]))
                self.fail = False
            def forward(self, input_ids, **kwargs):
                hidden = input_ids.float().unsqueeze(-1).expand(-1, -1, 640).clone()
                for layer in self.encoder.layer:
                    hidden = layer(hidden)
                if self.fail:
                    raise RuntimeError("fixture failure after hooks")
                length = input_ids.shape[1]
                attention = torch.ones(1, 20, length, length) / length
                return SimpleNamespace(attentions=(attention,) * 30)
            def contact_head(self, tokens, attentions):
                return attentions.mean((1, 2))[:, 1:-1, 1:-1]
        self.adapter = DCTdomainHFModel.__new__(DCTdomainHFModel)
        self.adapter.encoder = SimpleNamespace(closed=False)
        self.adapter.torch = torch
        self.adapter.model = TinyInterfaceModel()
        self.adapter.tokenizer = SimpleNamespace(cls_token_id=0, eos_token_id=2, pad_token_id=1)
        self.tokens = torch.tensor([[0, 5, 23, 13, 2]])

    def test_returns_exact_requested_layers_and_mutable_no_grad_states(self):
        result = self.adapter.esm_encoder(self.tokens, [15, 21], True)
        self.assertEqual(tuple(result["representations"]), (15, 21))
        self.assertEqual(tuple(result["contacts"].shape), (1, 3, 3))
        self.assertFalse(result["representations"][15].requires_grad)
        result["representations"][15][:, 1, :] = 0  # upstream overlap operation must be legal
        self.assertTrue(all(not layer._forward_hooks for layer in self.adapter.model.encoder.layer))

    def test_failed_forward_removes_hooks_and_wrong_contract_rejected(self):
        self.adapter.model.fail = True
        with self.assertRaisesRegex(RuntimeError, "fixture failure"):
            self.adapter.esm_encoder(self.tokens, [15, 21], True)
        self.assertTrue(all(not layer._forward_hooks for layer in self.adapter.model.encoder.layer))
        with self.assertRaises(ValueError):
            self.adapter.esm_encoder(self.tokens, [21, 15], True)
        with self.assertRaises(ValueError):
            self.adapter.esm_encoder(self.tokens, [15, 21], False)


if __name__ == "__main__":
    unittest.main()
