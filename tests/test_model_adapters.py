"""Optional tensor-level adapter contracts using numerical tensors, no pretrained weights."""
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
from types import ModuleType
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.model_adapters import LocalAdapter
from iri_constructs.runtime import InputError


@unittest.skipUnless(importlib.util.find_spec('torch'), 'Optional PyTorch adapter contract tests')
class AdapterContractTests(unittest.TestCase):
    def setUp(self):
        import torch
        self.torch = torch
        self.adapter = LocalAdapter.__new__(LocalAdapter)
        self.adapter.name, self.adapter.device, self.adapter.max_residues = 'esm2', 'cpu', 10
        class Tokenizer:
            unk_token_id = 99
            def __call__(_, sequence, **kwargs):
                # BOS, residues, EOS, PAD. Large special representations must not pool.
                n = len(sequence)
                return {'input_ids': torch.tensor([[0] + [4]*n + [2,1]]),
                        'attention_mask': torch.tensor([[1]*(n+2) + [0]]),
                        'special_tokens_mask': torch.tensor([[1] + [0]*n + [1,1]])}
        self.adapter.tokenizer = Tokenizer()
        self.adapter.model = lambda **inputs: SimpleNamespace(last_hidden_state=
            torch.arange(inputs['input_ids'].numel(), dtype=torch.float32).reshape(1,-1,1))

    def test_bos_eos_and_padding_excluded(self):
        self.assertEqual(self.adapter.residues('A' * 2), [[1.], [2.]])

    def test_length_and_alphabet_rejected_without_truncation(self):
        for sequence in ('A' * 11, 'X', ''):
            with self.assertRaises(InputError):
                self.adapter.residues(sequence)

    def test_token_count_mismatch_and_nonfinite_output_fail(self):
        for values in ([[[1.]]], [[[float('nan')]]*4]):
            self.adapter.model = lambda **kwargs: SimpleNamespace(last_hidden_state=self.torch.tensor(values))
            with self.assertRaises(InputError):
                self.adapter.residues('A' * 2)

    def test_sdk_bos_eos_contract_for_both_models(self):
        torch = self.torch
        api = ModuleType('esm.sdk.api')
        api.ESMProtein = lambda sequence: SimpleNamespace(sequence=sequence)
        class NumericalSdk:
            def encode(self, protein):
                return SimpleNamespace(sequence=torch.tensor([0] + [4]*len(protein.sequence) + [2]))
            def __call__(self, sequence_tokens):
                return SimpleNamespace(embeddings=torch.arange(sequence_tokens.numel(),dtype=torch.float32).reshape(1,-1,1))
        self.adapter.tokenizer = SimpleNamespace(bos_token_id=0,eos_token_id=2,unk_token_id=99)
        self.adapter.model = NumericalSdk()
        with patch.dict(sys.modules, {'esm': ModuleType('esm'), 'esm.sdk': ModuleType('esm.sdk'), 'esm.sdk.api': api}):
            for name in ('esmc', 'esm3'):
                self.adapter.name = name
                self.assertEqual(self.adapter.residues('A'*2), [[1.],[2.]])
            self.adapter.tokenizer.bos_token_id = 7
            with self.assertRaisesRegex(InputError, 'BOS/EOS'):
                self.adapter.residues('A'*2)
