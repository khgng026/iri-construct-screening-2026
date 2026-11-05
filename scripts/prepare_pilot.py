"""Prepare a T02-reviewed local T08 pilot and audit source annotation preservation."""
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.config import configure_logging
from iri_constructs.pilot import prepare_pilot
from iri_constructs.runtime import InputError


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True, help='Reviewed t05-inputs-1 JSON; source paths relative to this file')
    parser.add_argument('--output', type=Path, required=True, help='New ignored results/ or data/processed/ subdirectory')
    args = parser.parse_args(argv)
    log = configure_logging('INFO')
    try:
        report = prepare_pilot(args.inputs, args.output)
    except InputError as exc:
        parser.error(str(exc))
    log.info('event=pilot_preparation_end status=%s parsed_records=%d scored_records=%d',
             report['status'], report.get('parsed_records', 0), report['scored_records'])
    print('Pilot report:', args.output / 'report.json')
    return 0 if report['status'] == 'ok' else 1


if __name__ == '__main__':
    raise SystemExit(main())
