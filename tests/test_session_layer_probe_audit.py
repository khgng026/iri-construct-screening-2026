import copy
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_layer_probe_audit import validate_probe_record


class SourceProbeEvidenceTests(unittest.TestCase):
    def test_all_four_source_inputs_need_complete_raw_precision_evidence(self):
        query = {'query_id': 'canonical-fixture-id', 'length_aa': 67}
        config = {'layer': 3}
        record = {'query_id': query['query_id'], 'length': 67, 'times_seconds': [.1, .09, .08],
                  'actual_forward_outputs': [{'dtype': 'torch.float32', 'device': 'cpu'} for _ in range(3)],
                  'within_config_max_abs': 0., 'actual_parameter_dtypes': ['torch.float32'],
                  'encode': {'layer': 3, 'device': 'cpu', 'dtype': 'float32', 'truncated': False, 'model_training': False,
                             'inference_mode': True, 'residue_count': 67, 'token_count': 69}}
        validate_probe_record(record, query, config)
        mutations = [lambda x: x['actual_forward_outputs'].pop(),
                     lambda x: x['actual_forward_outputs'][0].update(dtype='torch.float16'),
                     lambda x: x['actual_forward_outputs'][2].update(device='mps'),
                     lambda x: x['times_seconds'].__setitem__(1, float('nan')),
                     lambda x: x.update(query_id='different-id'),
                     lambda x: x['encode'].update(layer=2),
                     lambda x: x['encode'].update(token_count=67),
                     lambda x: x.update(within_config_max_abs=-1.)]
        for mutation in mutations:
            malformed = copy.deepcopy(record)
            mutation(malformed)
            with self.assertRaises(ValueError):
                validate_probe_record(malformed, query, config)


if __name__ == '__main__':
    unittest.main()
