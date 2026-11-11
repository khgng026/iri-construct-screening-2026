"""Create/verify an exact-file evidence freeze for a separate scientific review."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.review_freeze import build_review_freeze, verify_review_freeze
from iri_constructs.runtime import InputError, read_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact', type=Path, action='append', default=[],
                        help='Repository-relative file to bind; repeat for every relevant artifact')
    parser.add_argument('--output', type=Path, help='New evidence-freeze JSON; never overwritten')
    parser.add_argument('--verify', type=Path, help='Verify an existing freeze without approving execution')
    args = parser.parse_args(argv)
    try:
        if args.verify:
            if args.output or args.artifact:
                raise InputError('Verification takes only --verify')
            print(json.dumps(verify_review_freeze(ROOT, read_json(args.verify)), sort_keys=True))
            return 0
        if not args.output:
            raise InputError('--output is required when creating an evidence freeze')
        output = args.output.resolve()
        if not output.is_relative_to(ROOT):
            raise InputError('Review freeze output must remain inside the repository')
        freeze = build_review_freeze(ROOT, args.artifact)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open('x', encoding='utf-8') as handle:
            handle.write(json.dumps(freeze, indent=2, sort_keys=True, allow_nan=False) + '\n')
        print('Evidence frozen for separate review; scientific execution remains unauthorized:',
              freeze['freeze_sha256'])
        return 0
    except FileExistsError:
        print('blocked: output already exists; retain the original freeze', file=sys.stderr)
    except InputError as exc:
        print('blocked:', exc, file=sys.stderr)
    return 1


if __name__ == '__main__':
    raise SystemExit(main())
