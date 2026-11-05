"""Demonstrate reviewed local FASTA/GenBank parsing without model dependencies."""
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.annotations import json_safe, normalize_record
from iri_constructs.bench_inputs import read_source, reviewed_sources
from iri_constructs.config import configure_logging
from iri_constructs.runtime import InputError, Measure, digest_file, failure, local_output, run_metadata, write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='New subdirectory of results/ or data/processed/')
    args = parser.parse_args(argv)
    log = configure_logging('INFO')
    try:
        output = local_output(args.output)
        if (output / 'report.json').exists():
            raise InputError('Output already has a report; choose a new directory')
    except InputError as exc:
        parser.error(str(exc))
    report = run_metadata('annotation_parsing', {})
    report['inference_executed'] = False
    with Measure() as measured:
        try:
            config, sources, provenance = reviewed_sources(args.inputs)
            documents, audits, proteins, seen = [], {}, {}, set()
            for source in sources:
                if source.molecule != 'dna':
                    raise InputError('Annotation demonstration requires DNA inputs; protein FASTA belongs in the probe')
                cid = source.review['record_id']
                if cid in seen:
                    raise InputError('Parse each construct once per annotation run')
                seen.add(cid)
                doc, audit, extracted = normalize_record(read_source(source), source.review, digest_file(source.path))
                documents.append(doc)
                audits[cid] = json_safe(audit)
                proteins.update({source.input_id + '/' + uid: sequence for uid, sequence in extracted.items()})
            write_json(output / 'normalized.json', documents)
            write_json(output / 'annotation-audit.json', audits)
            # Local, gitignored artifact only; never emit sequences to logs/reports.
            with (output / 'units.fasta').open('w', encoding='utf-8') as handle:
                for iid, sequence in sorted(proteins.items()):
                    handle.write('>' + iid + '\n' + sequence + '\n')
            report.update(status='ok', provenance=provenance, constructs=len(documents),
                          extracted_annotated_proteins=len(proteins), inference_executed=False)
            code = 0
        except Exception as exc:
            report.update(status='blocked', **failure(exc))
            code = 1
    report['resources_total'] = measured.data
    write_json(output / 'report.json', report)
    log.info('event=annotation_parse_end status=%s', report['status'])
    return code


if __name__ == '__main__':
    raise SystemExit(main())
