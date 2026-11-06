"""Geometry, deterministic fusion and dev-only selection; no model execution."""
import inspect
import unittest

from iri_constructs.session_benchmark import rank_scores
from iri_constructs.session_segmentation import (
    SegmentationConfig, alignment_guided_segments, fixed_windows, fuse_rankings,
    grouped_partition, interval_metrics, per_segment_union, select_configuration,
)


class SegmentationTest(unittest.TestCase):
    def test_windows_cover_termini_without_labels(self):
        for length in (32, 128, 129, 256, 1011):
            for width in (128, 256):
                result = fixed_windows(length, width)
                spans = result['segments']
                self.assertEqual(spans[0]['start0'], 0)
                self.assertEqual(spans[-1]['end0'], length)
                self.assertEqual(len({(s['start0'],s['end0']) for s in spans}), len(spans))
                self.assertTrue(all(a['end0'] >= b['start0'] for a,b in zip(spans,spans[1:])))
                self.assertTrue(all(s['end0']-s['start0'] <= width for s in spans))
        self.assertEqual(len(fixed_windows(1011,128)['segments']), 15)

    def test_iterative_geometry_and_compute_counts_use_no_target_boundaries(self):
        seen = []
        def ranker(fragment, catalogue):
            seen.append(len(fragment))
            values = {100: (200,[20,60]), 20: (80,[0,20]), 40: (100,[0,30])}
            score, span = values[len(fragment)]
            return {'ranking': [{'accession':'reference', 'rank':1, 'score':score,
                                 'query_span':span, 'reference_coverage':1., 'identity_fraction':1.}]}
        # Numeric dummy text tests geometry only; no biological sequences created.
        result = alignment_guided_segments('0123456789'*10, {'reference':'dummy'},
                                            SegmentationConfig(minimum_span=10,maximum_segments=2), ranker=ranker)
        self.assertEqual([(s['start0'],s['end0']) for s in result['segments']], [(20,60),(60,90)])
        self.assertEqual(seen, [100,20,40])
        self.assertEqual(result['total_queried_residues'], 160)
        self.assertEqual(result['alignment_comparisons'], 3)
        self.assertTrue(result['cap_bound'])
        self.assertEqual(result['unprocessed_intervals'], [[0,20],[90,100]])

    def test_filters_search_best_eligible_hit_and_explicit_fallback(self):
        def ranker(fragment, catalogue):
            return {'ranking': [
                {'accession':'A','rank':1,'score':200,'query_span':[0,len(fragment)],
                 'reference_coverage':.1,'identity_fraction':1.},
                {'accession':'B','rank':2,'score':100,'query_span':[0,len(fragment)],
                 'reference_coverage':1.,'identity_fraction':1.}]}
        args = ('0123456789', {'A':'dummy','B':'dummy'})
        result = alignment_guided_segments(*args, SegmentationConfig(minimum_span=5), ranker=ranker)
        self.assertEqual(result['segments'][0]['seed_hit']['accession'], 'B')
        self.assertFalse(result['whole_query_fallback_used'])
        empty = alignment_guided_segments(*args, SegmentationConfig(minimum_span=5,minimum_score=300), ranker=ranker)
        self.assertTrue(empty['whole_query_fallback_used'])
        self.assertEqual(empty['accepted_alignment_segment_count'], 0)
        self.assertEqual((empty['segments'][0]['start0'],empty['segments'][0]['end0']), (0,10))

    def test_prediction_signatures_have_no_label_or_boundary_inputs(self):
        params = set(inspect.signature(alignment_guided_segments).parameters)
        self.assertEqual(params, {'sequence','catalogue','config','ranker'})
        self.assertEqual(set(inspect.signature(fixed_windows).parameters), {'length','width','stride'})

    def test_fusion_is_label_free_complete_and_order_invariant(self):
        one = rank_scores([.9,.8,.1], ['A','B','C'], 'cosine_similarity')
        two = rank_scores([.1,.8,.9], ['A','B','C'], 'cosine_similarity')
        forward = fuse_rankings([one,two])
        backward = fuse_rankings([two,one])
        self.assertEqual(forward, backward)
        # Reciprocal-rank sums slightly favor [rank1,rank3] over [rank2,rank2].
        self.assertEqual([r['accession'] for r in forward['ranking']], ['A','C','B'])
        maximum = fuse_rankings([one,two], 'max_score')
        self.assertEqual([r['score'] for r in maximum['ranking']], [.9,.9,.8])
        larger = per_segment_union([one,two], 2)
        self.assertEqual(larger['accessions'], ['A','B','C'])
        self.assertEqual(larger['returned_hit_occurrences'], 4)
        self.assertEqual(larger['duplicate_returned_hits'], 1)
        with self.assertRaises(ValueError):
            fuse_rankings([one,two[:2]])
        with self.assertRaises(ValueError):
            fuse_rankings([rank_scores([1],['A'])], 'max_score')

    def test_interval_evaluation_union_does_not_double_count_overlap(self):
        result = interval_metrics([{'start0':0,'end0':60},{'start0':40,'end0':100}],
                                  [{'accession':'target','start0':20,'end0':80}])
        row = result['components'][0]
        self.assertEqual(row['covered_fraction_by_union'], 1)
        self.assertEqual(row['best_interval_iou'], .5)
        self.assertFalse(row['exact_boundary_match'])

    def test_all_permutations_stay_in_one_split(self):
        queries = [{'query_id':'ab','group_id':'AB','components':['A','B']},
                   {'query_id':'ba','group_id':'AB','components':['B','A']}]
        result = grouped_partition(queries)
        self.assertEqual(result['queries']['ab'], result['queries']['ba'])
        self.assertEqual(result, grouped_partition(queries[::-1]))
        queries[1]['group_id'] = 'BA'
        with self.assertRaises(ValueError):
            grouped_partition(queries)

    def test_dev_selection_uses_group_macro_and_rejects_test_or_missing_rows(self):
        design = {'selection_query_ids':['q1','q2','q3'],
                  'group_queries':{'g1':['q1','q2'],'g2':['q3'],'g3':['test']},
                  'partition':{'queries':{'q1':'dev','q2':'dev','q3':'dev','test':'test'}},
                  'candidate_configurations':[{'configuration_id':'a','segmenter':'windows'},
                                              {'configuration_id':'b','segmenter':'alignment'}]}
        def rows(values, work=100):
            return [{'query_id':q,'component_recall_at_5':v,'all_components_hit_at_5':int(v==1),
                     'total_queried_residues':work} for q,v in zip(['q1','q2','q3'],values)]
        candidates = [{'configuration_id':'a','rows':rows([1,1,0])},
                      {'configuration_id':'b','rows':rows([.5,.5,1])}]
        selected = select_configuration(candidates,design)
        self.assertEqual(selected['selected_configuration_id'], 'b')
        self.assertEqual(selected['development_scores'][0]['group_macro_component_recall_at_5'], .75)
        self.assertEqual(selected['selected_by_segmenter'], {'alignment':'b','windows':'a'})
        candidates[0]['rows'][-1]['query_id'] = 'test'
        with self.assertRaises(ValueError):
            select_configuration(candidates,design)


if __name__ == '__main__':
    unittest.main()
