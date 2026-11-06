"""Run the prospective E2 parent-level documented-CDS comparison."""
from pathlib import Path
import argparse
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from iri_constructs.config import configure_logging
from iri_constructs.multi_cds import run_multi_cds
from iri_constructs.runtime import InputError


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--protocol',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--device',choices=('cpu','cuda','auto'),default='cpu')
    parser.add_argument('--allow-download',action='store_true',help='Model weights only; no source-record downloads')
    args=parser.parse_args(argv)
    log=configure_logging('INFO')
    try:
        report=run_multi_cds(args.inputs,args.protocol,args.output,device=args.device,allow_download=args.allow_download)
    except InputError as exc:
        parser.error(str(exc))
    log.info('event=multi_cds_end status=%s scientific_scores=%s',report['status'],report['scientific_scores_produced'])
    print('E2 report:',args.output/'report.json')
    return 0 if report['status']=='ok' else 1


if __name__=='__main__':
    raise SystemExit(main())
