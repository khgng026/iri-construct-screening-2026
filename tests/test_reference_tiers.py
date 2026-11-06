"""Fictional metadata/pair policies only; no actual source approval or model findings."""
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.bench_inputs import ProteinInput
from iri_constructs.reference_tiers import (assert_same_tier_panel, freeze_reference_tiers,
                                          score_reference_tier)
from iri_constructs.runtime import InputError


class ReferenceTierTests(unittest.TestCase):
    def setUp(self):
        self.config=json.loads((ROOT/'config'/'extensions'/'E3.json').read_text())
        self.proteins={key:ProteinInput(key,'parent-'+key,aa*3,'fictional-software-fixture')
            for key,aa in [('q','A'),('copy','A'),('parent','C'),('linked','D'),
                           ('positive1','E'),('positive2','F'),('negative1','G'),('negative2','H')]}
        self.proteins['parent']=ProteinInput('parent','parent-q','C'*3,'fictional-software-fixture')
        self.labels={key:dict(input_id=key,family='other' if key.startswith('negative') else 'target',
            annotation_status='known',review_status='approved',evidence='Fictional unit-test metadata, not a real approval',
            independence_group='group-'+key) for key in self.proteins}
        self.labels['linked']['independence_group']='group-q'
        self.references=[key for key in self.proteins if key!='q']

    def policy(self,scope,pairs):
        return dict(schema_version='reviewed-relatedness-policy-1',policy_id='fictional-'+scope,
            policy_kind='reviewed_pair_exclusions',policy_scope=scope,review_status='approved',
            reviewer='Fictional software reviewer',reviewed_on='2026-11-05',declared_before_results=True,
            evidence='Fictional pair-policy assertions for code tests only',
            excluded_pairs=[dict(query_id=q,reference_id=r,reason='Fictional linked pair',
                evidence='Fictional test evidence, not biological research',
                query_sequence_sha256=self.proteins[q].sequence_sha256,
                reference_sequence_sha256=self.proteins[r].sequence_sha256) for q,r in pairs])

    def tiers(self,**kwargs):
        return freeze_reference_tiers(self.config,self.proteins,self.labels,['q'],self.references,**kwargs)

    def test_tier_partition_and_all_t07_guards_are_invariant(self):
        result=self.tiers()
        for name in ('tier1_t07_invariants','tier2_additional_linkage'):
            tier=result['tiers'][name]
            allowed=tier['panel']['eligible_reference_ids']['q']
            self.assertEqual(allowed,['negative1','negative2','positive1','positive2'])
            self.assertEqual(set(allowed)|set(tier['removals']['q']),set(self.references))
            self.assertFalse(set(allowed)&set(tier['removals']['q']))
            self.assertEqual(tier['reference_counts_by_query']['q']['positive_references'],2)
            self.assertEqual(tier['reference_counts_by_query']['q']['competing_references'],2)
        self.assertIn('redundant',result['tiers']['tier2_additional_linkage']['additional_policy_status'])
        self.assertIsNone(result['tiers']['tier3_reviewed_relatedness']['panel'])
        self.assertEqual(result,self.tiers())

    def test_declared_reviewed_pair_policies_produce_monotonic_attrition(self):
        link=self.policy('vector_linkage',[('q','positive1')])
        strict=self.policy('sequence_relatedness',[('q','positive2')])
        result=self.tiers(linkage_policy=link,relatedness_policy=strict)
        second=result['tiers']['tier2_additional_linkage']
        third=result['tiers']['tier3_reviewed_relatedness']
        self.assertNotIn('positive1',second['panel']['eligible_reference_ids']['q'])
        self.assertEqual(second['eligible_queries'],1)
        self.assertEqual(third['eligible_queries'],0)
        self.assertEqual(third['status'],'blocked_insufficient_panel')
        self.assertEqual(third['query_attrition_from_prior_tier'],['q'])
        self.assertEqual(third['reference_counts_by_query']['q']['positive_references'],0)

    def test_unsupported_or_unreviewed_strict_policy_is_blocked(self):
        base=self.policy('sequence_relatedness',[('q','positive1')])
        for change in ({'policy_kind':'automatic_identity_threshold'},{'review_status':'pending'},
                       {'reviewer':''},{'declared_before_results':False},{'reviewed_on':'not-a-date'}):
            with self.subTest(change=change):
                result=self.tiers(relatedness_policy={**base,**change})
                self.assertEqual(result['tiers']['tier1_t07_invariants']['status'],'frozen')
                self.assertEqual(result['tiers']['tier3_reviewed_relatedness']['status'],'blocked_relatedness_policy')

    def test_policy_fingerprint_and_role_id_mismatch_are_blocked(self):
        policy=self.policy('sequence_relatedness',[('q','positive1')])
        policy['excluded_pairs'][0]['query_sequence_sha256']='wrong'
        self.assertIsNone(self.tiers(relatedness_policy=policy)['tiers']['tier3_reviewed_relatedness']['panel'])
        policy=self.policy('vector_linkage',[('q','positive1')])
        policy['excluded_pairs'][0]['reference_id']='not-in-panel'
        result=self.tiers(linkage_policy=policy)
        self.assertIsNone(result['tiers']['tier2_additional_linkage']['panel'])
        self.assertIsNone(result['tiers']['tier3_reviewed_relatedness']['panel'])

    def test_exact_method_panel_parity_is_required(self):
        panel=self.tiers()['tiers']['tier1_t07_invariants']['panel']
        self.assertEqual(assert_same_tier_panel({'embedding':panel,'alignment':copy.deepcopy(panel)}),panel['panel_sha256'])
        other=self.tiers(linkage_policy=self.policy('vector_linkage',[('q','positive1')]))['tiers']['tier2_additional_linkage']['panel']
        with self.assertRaises(InputError):
            assert_same_tier_panel({'embedding':panel,'alignment':other})

    def test_wrong_family_hits_have_explicit_retrieval_denominators(self):
        tier=self.tiers()['tiers']['tier1_t07_invariants']
        scores={'q':{'positive1':1.0,'positive2':.9,'negative1':.2,'negative2':.1}}
        result=score_reference_tier(scores,tier,self.labels)
        self.assertEqual(result['false_family_hit_count'],2)
        self.assertEqual(result['returned_hit_denominator'],4)
        self.assertEqual(result['false_family_hit_fraction'],.5)
        self.assertEqual(result['queries_with_false_family_match_rate'],1.0)
        self.assertTrue(result['rows'][0]['budget_saturated'])

    def test_missing_scores_and_empty_panels_remain_abstentions(self):
        tier=self.tiers()['tiers']['tier1_t07_invariants']
        result=score_reference_tier({'q':{'positive1':1.0}},tier,self.labels)
        self.assertEqual(result['scored_queries'],0)
        self.assertIsNone(result['false_family_hit_fraction'])
        empty=freeze_reference_tiers(self.config,{}, {},[],[])
        self.assertEqual(empty['tiers']['tier1_t07_invariants']['status'],'blocked_insufficient_panel')
        blocked=score_reference_tier({},empty['tiers']['tier3_reviewed_relatedness'],{})
        self.assertIsNone(blocked['false_family_hit_fraction'])

    def test_unapproved_family_label_does_not_enter_any_tier(self):
        self.labels['positive1']['review_status']='pending'
        with self.assertRaises(InputError):
            self.tiers()
