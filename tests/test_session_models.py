"""Adapter checks with toy tokens and tiny tensors; no weights or network calls."""
import copy
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from iri_constructs.session_models import (
    ResidueEncoder, validate_config, validate_sequence, validate_token_mapping,
)


CONFIG = {"model_id": "facebook/esm2_t6_8M_UR50D",
          "revision": "c731040fcd8d73dceaa04b0a8e6329b345b0f5df"}


class ToyTokenizer:
    # Deliberately fictional test tokens, not source records.
    cls_token_id = 0
    eos_token_id = 2
    unk_token_id = 3
    mapping = {"A": 5, "C": 23, "D": 13}

    def convert_tokens_to_ids(self, tokens):
        return [self.mapping.get(token, self.unk_token_id) for token in tokens]

    def __call__(self, sequence, **kwargs):
        import torch
        assert kwargs["truncation"] is False and kwargs["padding"] is False
        tokens = [0] + self.convert_tokens_to_ids(sequence) + [2]
        return {"input_ids": torch.tensor([tokens]),
                "attention_mask": torch.ones((1, len(tokens)), dtype=torch.long),
                "special_tokens_mask": torch.tensor([[1] + [0] * len(sequence) + [1]])}


class ConfigurationTests(unittest.TestCase):
    def test_requires_immutable_revision_offline_and_correct_family(self):
        for change in ({"revision": "main"}, {"local_files_only": False},
                       {"family": "esmc"}, {"model_id": "custom/unapproved"}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate_config(dict(CONFIG, **change))

    def test_resource_and_representation_limits_are_not_silently_changed(self):
        for change in ({"threads": 5}, {"threads": True}, {"device": "cuda"},
                       {"dtype": "float16"}, {"max_residues": 2048},
                       {"tokenizer_truncation": True}, {"layer": -2}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                validate_config(dict(CONFIG, **change))
        self.assertEqual(validate_config(CONFIG)["layer"], -1)

    def test_esm3_representation_is_explicit_and_mps_rejected(self):
        config = dict(CONFIG, model_id="biohub/esm3-sm-open-v1")
        self.assertEqual(validate_config(config)["representation_variant"], "official_embeddings_prenorm")
        with self.assertRaises(ValueError):
            validate_config(dict(config, device="mps"))
        with self.assertRaises(ValueError):
            validate_config(dict(config, representation_variant="arbitrary"))

    def test_invalid_sequence_is_rejected_without_normalization(self):
        validate_sequence("ACD", 3)
        for sequence in ("", "ac", "AC ", "AX", "A*", "ACDA", ["A"]):
            with self.subTest(sequence=sequence), self.assertRaises(ValueError):
                validate_sequence(sequence, 3)


class MappingTests(unittest.TestCase):
    def test_one_row_per_residue_excludes_both_specials(self):
        mask = validate_token_mapping("ACD", [0, 5, 23, 13, 2],
                        [1, 1, 1, 1, 1], [1, 0, 0, 0, 1], ToyTokenizer())
        self.assertEqual(mask.tolist(), [False, True, True, True, False])

    def test_same_count_but_reordered_residues_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "residue identities"):
            validate_token_mapping("ACD", [0, 23, 5, 13, 2],
                        [1] * 5, [1, 0, 0, 0, 1], ToyTokenizer())

    def test_wrong_special_unknown_and_noncontiguous_mappings_are_rejected(self):
        cases = [([2, 5, 23, 13, 0], [1, 0, 0, 0, 1]),
                 ([0, 5, 3, 13, 2], [1, 0, 0, 0, 1]),
                 ([0, 5, 23, 13, 2], [0, 1, 0, 0, 1])]
        for ids, special in cases:
            with self.subTest(ids=ids, special=special), self.assertRaises(ValueError):
                validate_token_mapping("ACD", ids, [1] * 5, special, ToyTokenizer())


class TensorAdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import torch
        cls.torch = torch
        torch.set_num_threads(1)

    def make_encoder(self, *, family="esmc", layer=2, bad=False, esm3_variant="official_embeddings_prenorm"):
        torch = self.torch

        class ToyEmbedding(torch.nn.Module):
            def forward(self, ids):
                return torch.stack((ids.float(), ids.float() + 1), dim=-1)

        class ToyBlock(torch.nn.Module):
            def forward(self, x):
                return x + 1, None

        class ToyNorm(torch.nn.Module):
            def forward(self, x):
                return x * 2

        class ToyModel(torch.nn.Module):
            def __init__(self):
                super().__init__()
                self.embed = ToyEmbedding()
                self.transformer = torch.nn.Module()
                self.transformer.blocks = torch.nn.ModuleList([ToyBlock(), ToyBlock()])
                self.transformer.norm = ToyNorm()

            def forward(self, input_ids=None, sequence_tokens=None, **kwargs):
                assert kwargs.get("output_hidden_states", False) is False
                assert kwargs["output_attentions"] is False
                x = self.embed(input_ids if input_ids is not None else sequence_tokens)
                for block in self.transformer.blocks:
                    x, _ = block(x)
                final = self.transformer.norm(x)
                if bad:
                    final[0, 1, 0] = float("nan")
                return SimpleNamespace(last_hidden_state=final, embeddings=x)

        encoder = ResidueEncoder.__new__(ResidueEncoder)
        encoder.config = {"family": family, "device": "cpu", "max_residues": 1022,
                          "representation_variant": esm3_variant}
        encoder.torch = torch
        encoder.tokenizer = ToyTokenizer()
        encoder.model = ToyModel().eval()
        encoder.closed = False
        encoder.num_layers = 2
        encoder.layer = layer
        return encoder

    def test_final_states_shape_dtype_special_removal_and_metadata(self):
        import numpy as np
        encoder = self.make_encoder()
        result = encoder.encode("ACD")
        self.assertEqual(result.shape, (3, 2))
        self.assertEqual(result.dtype, np.float32)
        np.testing.assert_array_equal(result[:, 0], [14, 50, 30])
        self.assertEqual(encoder.last_encode_metadata["residue_count"], 3)
        self.assertFalse(encoder.last_encode_metadata["truncated"])
        encoder.close()
        with self.assertRaises(RuntimeError):
            encoder.encode("ACD")

    def test_intermediate_capture_uses_one_layer_and_removes_hook(self):
        import numpy as np
        encoder = self.make_encoder(layer=1)
        result = encoder.encode("ACD")
        np.testing.assert_array_equal(result[:, 0], [6, 24, 14])
        self.assertEqual(len(encoder.model.transformer.blocks[0]._forward_hooks), 0)

    def test_esm3_primary_is_official_embeddings_with_postnorm_separate(self):
        import numpy as np
        primary = self.make_encoder(family="esm3").encode("ACD")
        alternate = self.make_encoder(family="esm3", esm3_variant="post_final_norm").encode("ACD")
        np.testing.assert_array_equal(primary[:, 0], [7, 25, 15])
        np.testing.assert_array_equal(alternate, primary * 2)

    def test_nonfinite_states_fail_instead_of_propagating(self):
        with self.assertRaisesRegex(ValueError, "Nonfinite"):
            self.make_encoder(bad=True).encode("ACD")


if __name__ == "__main__":
    unittest.main()
