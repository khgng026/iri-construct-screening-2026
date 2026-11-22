"""Apply existing answer-blind mergers to a complete chunk/reference score matrix."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent / 'scripts'))
from budget_compare import fixed_ranking

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scores', type=Path)
    parser.add_argument('--budget', type=int, default=5)
    parser.add_argument('--policy', choices=['max_score', 'round_robin'], default='round_robin')
    args = parser.parse_args()
    data = json.loads(args.scores.read_text())
    if set(data) != {'reference_ids', 'allowed', 'scores'}:
        raise ValueError('Input must contain only reference_ids, allowed, scores; no answer labels')
    ranking = fixed_ranking(np.asarray(data['scores'], dtype=float), data['reference_ids'],
                            data['allowed'], args.budget, args.policy)
    print(json.dumps({'policy':args.policy, 'budget':args.budget,
                      'returned_ids':[data['reference_ids'][i] for i in ranking]}, indent=2))

if __name__ == '__main__':
    main()
