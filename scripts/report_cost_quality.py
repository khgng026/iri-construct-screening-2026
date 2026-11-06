"""Aggregate actual E2/E3/A-D measurements without launching new model inference."""
from pathlib import Path
import argparse
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from iri_constructs.cost_quality import records_from_report, write_cost_quality_report
from iri_constructs.extension_config import load_extension_config
from iri_constructs.runtime import InputError, read_json


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records',type=Path,nargs='*',default=[],help='Normalized record batches or conservative legacy reports')
    parser.add_argument('--protocol',type=Path,default=Path(__file__).resolve().parents[1]/'config'/'extensions'/'E4.json')
    parser.add_argument('--output',type=Path,required=True,help='New ignored local output directory')
    args=parser.parse_args(argv)
    try:
        config=load_extension_config(args.protocol,expected='E4')
        records=[]
        for path in args.records:
            records.extend(records_from_report(read_json(path)))
        result=write_cost_quality_report(records,args.output,protocol=config)
    except InputError as exc:
        parser.error(str(exc))
    print('Cost/quality report:',args.output/'report.json')
    return 0 if result['status']=='ok' else 1


if __name__=='__main__':
    raise SystemExit(main())
