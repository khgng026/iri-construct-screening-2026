"""Record actual dependency/device/cache availability without reading sequences or downloading."""
from pathlib import Path
import argparse
import importlib.util
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.model_adapters import ESM2_REVISION, MODEL_REPOS, LocalAdapter, prepare_device
from iri_constructs.runtime import InputError, Measure, ROOT, failure, local_output, run_metadata, write_json


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results' / 't05_environment')
    parser.add_argument('--load-esm2', action='store_true', help='Also test offline cached ESM-2 loading; no sequence inference')
    args = parser.parse_args(argv)
    try:
        output = local_output(args.output)
    except InputError as exc:
        parser.error(str(exc))
    report = run_metadata('environment_preflight', {'downloads': False, 'load_cached_esm2': args.load_esm2})
    report['inference_executed'] = False
    report['availability'] = {m: importlib.util.find_spec(m) is not None for m in ('Bio','torch','transformers','esm','psutil')}
    try:
        import torch
        report['cuda_available'] = torch.cuda.is_available()
        report['cuda_device_count'] = torch.cuda.device_count()
        report['torch_build'] = torch.__version__
    except Exception as exc:
        report['torch_import_failure'] = failure(exc)
    try:
        from huggingface_hub import try_to_load_from_cache
        cached = try_to_load_from_cache(MODEL_REPOS['esm2'], 'model.safetensors', revision=ESM2_REVISION)
        report['pinned_esm2_safetensors_cached'] = isinstance(cached, str) and Path(cached).is_file()
    except Exception as exc:
        report['cache_check_failure'] = failure(exc)
    if args.load_esm2:
        adapter = None
        with Measure() as measured:
            try:
                device = prepare_device('cpu', 0)
                adapter = LocalAdapter('esm2', device, ESM2_REVISION, offline=True)
                report['esm2_load_check'] = {'status': 'loaded', 'model': adapter.metadata}
            except Exception as exc:
                report['esm2_load_check'] = {'status': 'blocked', **failure(exc)}
            finally:
                if adapter is not None:
                    adapter.close()
        report['esm2_load_check']['resources'] = measured.data
    report.update(status='environment_checked',
                  reviewed_input_status='No input data reviewed or approved by this command')
    write_json(output / 'report.json', report)
    print('Environment report:', output / 'report.json')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
