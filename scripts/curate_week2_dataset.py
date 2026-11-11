"""Acquire bounded public source evidence for preparatory QA, never scoring approval."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.dataset_curation import curate, publish_review
from iri_constructs.runtime import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/week2/curation')
    parser.add_argument('--offline', action='store_true', help='Require verified cached original bytes; make no network requests')
    parser.add_argument('--source', choices=('all', 'ncbi', 'igem'), default='all')
    parser.add_argument('--publish-manifests', type=Path, help='Snapshot source-free review and rebased QA inputs under data/manifests/week2/')
    parser.add_argument('--only-id', action='append', default=[], help='Bound to named source IDs; repeat as needed')
    args = parser.parse_args()
    report = curate(args.output, offline=args.offline, source=args.source, only_ids=args.only_id)
    if args.publish_manifests:
        publish_review(args.output, args.publish_manifests)
    print('Preparatory QA eligible:', report['preparatory_QA_eligible_count'], '/', report['record_count'],
          '; newly scoring admitted:', report['new_scoring_admission_count'])
    return 0 if report['preparatory_QA_eligible_count'] == report['record_count'] and report['record_count'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
