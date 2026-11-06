import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from iri_constructs.t05_data import build_query, sha256
from iri_constructs.t05_composition import (make_plan, duplicate_groups, equivalence_map,
    query_recall, component_metrics, canonical_hash, verify_freeze, fraction_bin)
try:
    import numpy as np
except ImportError:
    np=None


class CompositionTests(unittest.TestCase):
    def test_all_pairs_cap_and_order_preserve_components(self):
        seq={'a':'A'*500,'b':'C'*522,'c':'D'*523}
        p=make_plan(['a','b','c'],seq)
        self.assertEqual([x['components'] for x in p['pairs']],[['a','b']])
        self.assertEqual(len(p['excluded_pairs']),2)
        self.assertEqual([q['components'] for q in p['queries'] if q['artificial']],[['a','b'],['b','a']])
        for order in [['a','b'],['b','a']]:
            for target in order:
                whole,oracle,(start,end)=build_query(order,target,seq)
                self.assertEqual(oracle,seq[target]); self.assertEqual(whole[start:end],seq[target])
        with self.assertRaises(ValueError): build_query(['a','c'],'a',seq)

    def test_two_relevant_components_recall_not_precision(self):
        ranks=[{'accession':a,'rank':i} for i,a in enumerate(['partner','offtarget','target'],1)]
        eq={a:[a] for a in ['partner','offtarget','target']}
        r=query_recall(ranks,['target','partner'],eq)
        self.assertEqual(r['recall_at_1'],.5); self.assertEqual(r['any_component_at_1'],1)
        self.assertEqual(r['all_components_at_1'],0); self.assertEqual(r['recall_at_5'],1)
        with self.assertRaises(ValueError):query_recall(ranks,['missing'],eq)

    def test_duplicate_aware_retrieval_does_not_erase_accession_failure(self):
        seq={'A':'ACD','B':'ACD','C':'EFG'}; eq=equivalence_map(seq)
        self.assertEqual(duplicate_groups(seq),[['A','B']])
        ranking=[{'accession':a,'rank':i,'cosine_similarity':s,'exact_score_tie_count':n}
                 for i,(a,s,n) in enumerate([('A',1.,2),('B',1.,2),('C',.5,1)],1)]
        m=component_metrics(ranking,'B',eq)
        self.assertEqual(m['rank'],2);self.assertEqual(m['hit_at_1'],0)
        self.assertEqual(m['sequence_equivalent_hit_at_1'],1)
        self.assertEqual((m['tie_best_rank'],m['tie_worst_rank']),(1,2))
        self.assertEqual(query_recall(ranking,['B'],eq)['recall_at_1'],0)

    def test_delta_sign_and_near_tie(self):
        rank=[{'accession':'A','rank':1,'cosine_similarity':.9,'exact_score_tie_count':1},
              {'accession':'B','rank':2,'cosine_similarity':.8999999,'exact_score_tie_count':1}]
        base={'rank':1,'cosine_similarity':1.,'reciprocal_rank':1.}
        m=component_metrics(rank,'B',{'B':['B']},base)
        self.assertEqual(m['rank_change'],1);self.assertLess(m['cosine_change'],0)
        self.assertEqual(m['near_tie_count_1e6'],2);self.assertEqual(m['tie_best_rank'],2)

    def test_frozen_input_tampering(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'config.json';p.write_text('{}')
            freeze={'file_sha256':{'config.json':sha256(p.read_bytes())}}
            verify_freeze(root,freeze);p.write_text('{"changed":1}')
            with self.assertRaises(ValueError):verify_freeze(root,freeze)

    def test_repeat_subset_selected_without_scores(self):
        seq={a:a*length for a,length in [('A',50),('C',70),('D',300),('E',700)]}
        p=make_plan(list(seq),seq)
        selected=[q for q in p['queries'] if q['query_id'] in p['repeat_query_ids']]
        self.assertEqual(sum(not q['artificial'] for q in selected),4)
        for pair in p['repeat_pair_ids']:
            self.assertEqual(sum(q['pair_id']==pair for q in selected),2)
        self.assertEqual(fraction_bin(.25),'0.25–0.40')
        self.assertEqual(fraction_bin(.4),'0.40–0.60')
        self.assertEqual(fraction_bin(.75),'0.75–1.00')

    def test_cache_key_sensitive_to_model_sequence_and_settings(self):
        original={'revision':'abc','sequence_hash':'123','pooling':'residue','threads':4}
        for k in original:
            altered={**original,k:str(original[k])+'changed'}
            self.assertNotEqual(canonical_hash(original),canonical_hash(altered))
        self.assertEqual(canonical_hash(original),canonical_hash(dict(reversed(list(original.items())))))

    @unittest.skipIf(np is None,'NumPy needed for cache check')
    def test_cache_roundtrip_corruption_and_no_overwrite(self):
        from iri_constructs.t05_composition import cache_write,cache_read
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);key={'revision':'test','sequences':{'A':'hash'}}
            arr=np.array([[1.,2.]],dtype=np.float32)
            self.assertIsNone(cache_read(root,key,(1,2)))
            digest=cache_write(root,key,arr)
            self.assertTrue(np.array_equal(arr,cache_read(root,key,(1,2))))
            with self.assertRaises(ValueError):cache_write(root,key,arr)
            self.assertIsNone(cache_read(root,{'revision':'other'},(1,2)))
            (root/(digest+'.npz')).write_bytes(b'corrupted')
            with self.assertRaises(ValueError):cache_read(root,key,(1,2))

    @unittest.skipIf(np is None,'NumPy needed for descriptive analysis checks')
    def test_spearman_ties_and_constant_outcome(self):
        sys.path.insert(0,str(ROOT/'scripts'))
        from t05_composition_analyze import averaged_ranks,spearman
        self.assertTrue(np.array_equal(averaged_ranks([2,1,2,4]),[2.5,1,2.5,4]))
        self.assertAlmostEqual(spearman([1,2,3],[3,2,1]),-1.)
        self.assertIsNone(spearman([1,2,3],[1,1,1]))

    @unittest.skipIf(np is None,'NumPy needed for descriptive analysis checks')
    def test_grouped_counts_do_not_call_reversals_independent_pairs(self):
        sys.path.insert(0,str(ROOT/'scripts'))
        from t05_composition_analyze import summarize
        rows=[]
        for qid in ['A_B','B_A']:
            for a in ['A','B']:
                rank=1 if a=='A' else 8
                rows.append({'accession':a,'pair_id':'A__B','query_id':qid,'fraction_in_parent':.5,
                    'rank':rank,'rank_change':rank-1,'reciprocal_rank':1/rank,
                    'cosine_similarity':.9,'cosine_change':-.1,'exact_score_tie_count':1,'near_tie_count_1e6':1,
                    **{f'hit_at_{k}':int(rank<=k) for k in (1,5,10)},
                    **{f'sequence_equivalent_hit_at_{k}':int(rank<=k) for k in (1,5,10)}})
        s=summarize(rows)
        self.assertEqual((s['component_instances'],s['distinct_pairs'],s['distinct_proteins'],s['ordered_queries']),(4,1,2,2))
        self.assertEqual(s['recall_at_1'],.5);self.assertEqual(s['recall_at_5'],.5)
        self.assertEqual(s['recall_at_10'],1.)


if __name__=='__main__':unittest.main()
