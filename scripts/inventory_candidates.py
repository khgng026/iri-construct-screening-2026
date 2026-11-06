"""Validate sequence-free T06 metadata and produce an explicit manual-review queue."""
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from iri_constructs.candidate_inventory import (build_pilot_review_queue, load_candidate_inventory,
                                               review_queue, validate_pilot_review_queue)
from iri_constructs.config import configure_logging
from iri_constructs.runtime import InputError, digest_file, failure, local_output, run_metadata, write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, default=ROOT / 'data/manifests/t06_candidate_inventory.json')
    parser.add_argument('--output', type=Path, required=True, help='New ignored results/ or data/processed/ subdirectory')
    parser.add_argument('--check-pilot-queue', type=Path, help='Check a saved metadata-only pilot queue matches this inventory')
    parser.add_argument('--shortlist', type=Path, help='Optional versioned researcher candidate-selection decision; never source approval')
    args = parser.parse_args(argv)
    log = configure_logging('INFO')
    try:
        output = local_output(args.output)
        if output.exists() and any(output.iterdir()):
            raise InputError('Inventory output must be a new or empty ignored directory')
    except (InputError, OSError) as exc:
        parser.error(str(exc))
    report = run_metadata('t06_candidate_inventory', {'downloads': False, 'sequence_processing': False})
    report['inference_executed'] = False
    try:
        rows = load_candidate_inventory(args.inventory)
        queue = review_queue(rows)
        if args.check_pilot_queue is not None:
            validate_pilot_review_queue(args.check_pilot_queue, args.inventory, args.shortlist)
        write_json(output / 'review_queue.json', queue)
        write_json(output / 'pilot_review_queue.json', build_pilot_review_queue(args.inventory, args.shortlist))
        report.update(status='metadata_validated_pending_review', candidate_count=len(rows),
                      approved_for_scoring=0, proposed_pilot_count=queue['proposed_pilot_count'],
                      inventory_sha256=digest_file(args.inventory))
        code = 0
    except Exception as exc:
        report.update(status='blocked', **failure(exc))
        code = 1
    write_json(output / 'report.json', report)
    log.info('event=candidate_inventory_validation status=%s', report['status'])
    return code


if __name__ == '__main__':
    raise SystemExit(main())
