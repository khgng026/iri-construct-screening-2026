"""Prepare deterministic Week 2 source QA without granting scientific clearance."""
from pathlib import Path
import argparse
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.config import configure_logging
from iri_constructs.runtime import InputError, Measure, failure, local_output, run_metadata, write_json
from iri_constructs.week2_pipeline import run_preparatory_qa


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    log = configure_logging('INFO')
    try:
        output = local_output(args.output)
        if output.exists():
            raise InputError('Output already exists; choose a new Week 2 QA directory')
        execution = run_metadata('week2_preparatory_source_QA', {'input_manifest': str(args.inputs)})
        with Measure() as measured:
            report = run_preparatory_qa(args.inputs, output)
        execution.update(status=report['status'], evidence_kind='preparatory_source_QA',
                         inference_executed=False, scientific_scores_produced=False,
                         scientific_admission=False, resources_total=measured.data,
                         deterministic_content_sha256=report['content_sha256'])
        write_json(output / 'execution.json', execution)
        for row in report['records']:
            log.info('event=week2_source_qa source_id=%s status=%s failure_class=%s',
                     row['source_id'], row['status'], row.get('failure_class', 'none'))
        log.info('event=week2_preparatory_qa status=%s prepared=%s blocked=%s proteins=%s scientific_admission=false',
                 report['status'], report['prepared_records'], report['blocked_records'], report['exported_annotated_proteins'])
        return 0 if not report['blocked_records'] else 1
    except Exception as exc:
        controlled = failure(exc)
        log.error('event=week2_preparatory_qa status=blocked error_type=%s reason=%s',
                  controlled['error_type'], controlled['reason'])
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
