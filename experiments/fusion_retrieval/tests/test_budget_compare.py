import sys
from pathlib import Path
import numpy as np
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from budget_compare import fixed_ranking, rank_lists, coverage


def test_family_crowding_without_any_target_labels():
    # Two component score lists favor different references. A total budget must
    # preserve each list's first hit before spending it on the first family.
    s=np.array([[1,.99,.98,.1,.05],[.1,.09,.08,.8,.79]])
    ids=['a','b','c','d','e']
    assert fixed_ranking(s,ids,list(range(5)),2,'max_score')==[0,1]
    assert fixed_ranking(s,ids,list(range(5)),2,'round_robin')==[0,3]


def test_exclusion_duplicates_and_one_total_budget():
    s=np.array([[9,8,7,6],[9,8,6,7]])
    ids=['a','b','c','d']
    for policy in ('max_score','round_robin'):
        r=fixed_ranking(s,ids,[1,2,3],3,policy)
        assert len(r)==len(set(r))==3 and 0 not in r


def test_score_row_permutation_and_ties():
    s=np.array([[1.,1.,1.]])
    assert fixed_ranking(s,['c','a','b'],[0,1,2],2)==[1,2]
    for policy in ('round_robin','max_score'):
        assert fixed_ranking(s,['c','a','b'],[0,1,2],1,policy)==[1]


def test_nonfinite_and_bad_dimensions_fail_closed():
    with unittest.TestCase().assertRaises(ValueError):rank_lists([[1,np.nan]],['a','b'],[0,1])
    with unittest.TestCase().assertRaises(ValueError):rank_lists([[1]],['a','b'],[0,1])
    with unittest.TestCase().assertRaises(ValueError):rank_lists([[1,2]],['a','a'],[0,1])


def test_zero_based_half_open_coverage():
    assert coverage([(0,5),(3,8)],10)==.8
    with unittest.TestCase().assertRaises(ValueError):coverage([(1,11)],10)


def load_tests(loader, tests, pattern):
    return unittest.TestSuite(unittest.FunctionTestCase(v) for k,v in globals().copy().items() if k.startswith('test_') and callable(v))


if __name__ == '__main__':
    unittest.main()
