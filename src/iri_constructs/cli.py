"""``iri-constructs`` command line. Exit status: 0 complete, 2 partial batch, 1 unusable."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from . import __version__


class _Parser(argparse.ArgumentParser):
    def error(self, message):  # argparse's default status 2 means "partial batch" here
        self.print_usage(sys.stderr)
        self.exit(1, f'{self.prog}: error: {message}\n')


def build_parser() -> argparse.ArgumentParser:
    parser = _Parser(prog='iri-constructs', description='Turn annotated GenBank constructs into traceable features, '
                     'validated protein units and diagnostics. No network access, model or approval manifest is used.')
    parser.add_argument('--version', action='version', version=f'iri-constructs {__version__}')
    commands = parser.add_subparsers(dest='command', required=True, parser_class=_Parser)
    pre = commands.add_parser('preprocess', help='process GenBank files into an output bundle',
        description='Write normalized.json, annotations.json, proteins.fasta, features.tsv, units.tsv, report.json '
                    'and run.json. Proteins are source /translation values that validate against the documented '
                    'coordinates, frame and table; everything else is reported with a reason code.',
        epilog='Exit status: 0 every record processed; 2 some records failed but usable output was written; '
               '1 invalid invocation, unsafe output location or no usable record.')
    pre.add_argument('inputs', nargs='+', metavar='INPUT', type=Path,
                     help='GenBank file (one or many records), nucleotide FASTA, or a directory of such files')
    pre.add_argument('--out', required=True, metavar='DIRECTORY', type=Path, help='new or empty output directory')
    pre.add_argument('--overwrite', action='store_true', help='replace an existing bundle in --out (never other files)')
    pre.add_argument('--max-model-residues', type=int, metavar='N',
                     help='also report which exported proteins exceed N residues (e.g. 1022 for ESM-2); nothing is truncated or dropped')
    pre.add_argument('--provenance', type=Path, metavar='JSON',
                     help='optional JSON object mapping an input file name or SHA-256 to source facts to carry into the outputs')
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    from .preprocess import PreprocessError, PreprocessOptions, preprocess_files
    try:
        provenance = None
        if args.provenance:
            try:
                provenance = json.loads(args.provenance.read_text(encoding='utf-8-sig'))
            except (OSError, ValueError) as exc:
                raise PreprocessError(f'Cannot read --provenance JSON ({type(exc).__name__})') from exc
            if not isinstance(provenance, dict) or any(not isinstance(v, dict) for v in provenance.values()):
                raise PreprocessError('--provenance must be a JSON object of objects')
        result = preprocess_files(args.inputs, args.out, options=PreprocessOptions(
            overwrite=args.overwrite, model_max_residues=args.max_model_residues, provenance=provenance))
    except PreprocessError as exc:
        print(f'iri-constructs: error: {exc}', file=sys.stderr)
        return 1
    s = result.summary()
    print(f"records: {s['records_ok']} ok, {s['records_failed']} failed, {s['records_skipped_duplicate']} duplicate "
          f"(from {s['inputs']} input file(s))")
    print(f"features: {s['features']} ({s['cds_features']} CDS, {s['noncoding_metadata_features']} kept as metadata)")
    withheld = ', '.join(f'{k} {v}' for k, v in s['cds_withheld_by_reason'].items())
    print(f"proteins exported: {s['proteins_exported']}; CDS withheld: {s['cds_withheld']}" + (f' ({withheld})' if withheld else ''))
    if s['model_compatibility']['evaluated']:
        m = s['model_compatibility']
        print(f"model limit {m['max_residues']} residues: {m['compatible']} compatible, {m['exceeds_max_residues']} too long (still exported)")
    for record in result.failed_records:
        where = record.input_name + (f'#{record.record_index}' if record.record_index else '')
        print(f'failed: {where}: {record.failure.code}: {record.failure.message}', file=sys.stderr)
    print(f'output: {result.output_dir}')
    if result.exit_code == 2:
        print('partial batch: successful records were written; see report.json for failures', file=sys.stderr)
    elif result.exit_code == 1:
        print('no usable record: report.json lists the failures', file=sys.stderr)
    return result.exit_code


if __name__ == '__main__':
    raise SystemExit(main())
