"""v1.1 engineering guards do not mutate frozen v1.0 scientific contracts."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.extension_config import base_scoring_config,load_extension_config,validate_extension_config
from iri_constructs.pilot_scoring import validate_scoring_config
from iri_constructs.runtime import InputError


class ExtensionConfigurationTests(unittest.TestCase):
    def test_four_separate_versioned_designs(self):
        configs=[load_extension_config(ROOT/'config/extensions'/f'E{i}.json') for i in range(1,5)]
        self.assertEqual({c['extension'] for c in configs},{'E1','E2','E3','E4'})
        self.assertTrue(all(c['protocol_version']=='1.1' and c['retrieval_budget']==5 for c in configs))

    def test_invalid_budget_version_and_outcome_tuning_rejected(self):
        config=load_extension_config(ROOT/'config/extensions/E2.json')
        for change in ({'top_k':True},{'top_k':10},{'retrieval_budget':4},
                       {'protocol_version':'1.0'},{'design_status':'tuned_after_results'},
                       {'models':['alignment']},{'models':['esm2','esm2']},{'max_residues':2000}):
            with self.subTest(change=change),self.assertRaises(InputError):
                validate_extension_config({**config,**change})

    def test_extension_discriminator_is_checked(self):
        with self.assertRaises(InputError):
            load_extension_config(ROOT/'config/extensions/E2.json',expected='E3')

    def test_e3_adapter_reuses_existing_scorer_contract(self):
        config=load_extension_config(ROOT/'config/extensions/E3.json')
        original=copy.deepcopy(config)
        derived=base_scoring_config(config)
        validate_scoring_config(derived)
        self.assertEqual(config,original)
        self.assertEqual(derived['protocol_version'],'1.1')
        self.assertEqual(derived['extension'],'E3')

    def test_v1_0_frozen_hashes_still_match(self):
        freeze=json.loads((ROOT/'config/pilot/design_freeze.json').read_text(encoding='utf-8'))
        for path,expected in freeze['files'].items():
            with self.subTest(path=path):
                actual=hashlib.sha256((ROOT/path).read_text(encoding='utf-8').encode()).hexdigest()
                self.assertEqual(actual,expected)

    def test_v1_1_has_a_separate_design_freeze(self):
        freeze=json.loads((ROOT/'config/extensions/design_freeze.json').read_text(encoding='utf-8'))
        self.assertEqual(freeze['protocol_version'],'1.1')
        self.assertEqual(len(freeze['files']),5)
        for path,expected in freeze['files'].items():
            with self.subTest(path=path):
                self.assertEqual(hashlib.sha256((ROOT/path).read_text(encoding='utf-8').encode()).hexdigest(),expected)
