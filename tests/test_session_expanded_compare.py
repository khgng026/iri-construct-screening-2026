import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from session_expanded_compare import require_common_design, completed, headline_rows, require_code_binding
from iri_constructs.session_store import UnitStore, file_hash


class CommonExpandedPanelTests(unittest.TestCase):
    def test_omitted_or_changed_scientific_dependency_is_not_accepted(self):
        with mock.patch('session_expanded_compare.file_hash', return_value='current'):
            require_code_binding({'code': {'required.py': 'current'}}, {'required.py'})
            with self.assertRaises(ValueError):
                require_code_binding({'code': {}}, {'required.py'})
            with self.assertRaises(ValueError):
                require_code_binding({'code': {'required.py': 'changed'}}, {'required.py'})

    def test_valid_hashes_cannot_hide_a_swapped_unit_key(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            store = UnitStore(root, {'fixture': 'completion-key-binding'})
            store.publish('summary', {'fixture': 'summary'})
            store.publish('query-a', {'fixture': 'query'})
            store.complete(['summary', 'query-a'])
            path = next(root.glob('completion-*.json'))
            value = json.loads(path.read_text())
            value['units']['summary'], value['units']['query-a'] = value['units']['query-a'], value['units']['summary']
            path.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError, 'key or scientific record binding'):
                completed(root)

    def test_loose_scientific_record_is_not_authenticated_by_another_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            store = UnitStore(root, {'fixture': 'record-binding'})
            store.publish('summary', {'fixture': 'summary'})
            store.complete(['summary'])
            path = next(root.glob('completion-*.json'))
            value = json.loads(path.read_text())
            folder = root / value['units']['summary']['path']
            record = folder / 'record.json'
            unrelated = folder / 'unrelated.json'
            unrelated.write_bytes(record.read_bytes())
            descriptor_path = folder / 'complete.json'
            descriptor = json.loads(descriptor_path.read_text())
            descriptor['files'] = {'unrelated.json': descriptor['files']['record.json']}
            descriptor_path.write_text(json.dumps(descriptor))
            value['units']['summary']['descriptor_sha256'] = file_hash(descriptor_path)
            path.write_text(json.dumps(value))
            with self.assertRaisesRegex(ValueError, 'key or scientific record binding'):
                completed(root)

    def test_missing_quadruple_stratum_cannot_be_omitted_from_report(self):
        aggregate = {'query_count': 1130, 'component_count': 3520,
                     'kinds': {'exact': {'by_component_count': {'2': {}, '3': {}}}}}
        clean = {'query_count': 16, 'component_count': 16, 'kinds': {}}
        summary = {'summaries': {str(n): {a: {'combined': aggregate, 'clean': clean}
                                         for a in ('local_alignment_raw', 'local_alignment_reference_selfscore')}
                                 for n in (101, 250, 500, 1000)}}
        with self.assertRaisesRegex(ValueError, 'component-count strata'):
            headline_rows(summary, 'alignment', 'fixture', 'expanded_alignment_retrieval_v1')

    def test_missing_queries_or_changed_catalogues_do_not_become_common_subset(self):
        query_ids = ['q' + str(i) for i in range(1146)]
        original = {'query_ids': query_ids, 'expansion_manifest_sha256': 'manifest', 'retrieval_design_sha256': 'design',
                    'target_aliases': {'a': 'a'}, 'subsets': {str(n): [str(i) for i in range(n)] for n in (101, 250, 500, 1000)}}
        require_common_design([original, copy.deepcopy(original)], query_ids)
        for mutation in ('missing', 'duplicate', 'catalogue', 'mask'):
            changed = copy.deepcopy(original)
            if mutation == 'missing':
                changed['query_ids'].pop()
            elif mutation == 'duplicate':
                changed['query_ids'][-1] = changed['query_ids'][0]
            elif mutation == 'catalogue':
                changed['subsets']['1000'][-1] = 'different'
            else:
                changed['target_aliases']['a'] = 'different'
            with self.assertRaises(ValueError):
                require_common_design([original, changed], query_ids)

    def test_equal_sources_with_non_nested_catalogue_are_rejected(self):
        query_ids = ['q' + str(i) for i in range(1146)]
        source = {'query_ids': query_ids, 'expansion_manifest_sha256': 'm', 'retrieval_design_sha256': 'd', 'target_aliases': {},
                  'subsets': {str(n): [str(i) for i in range(n)] for n in (101, 250, 500, 1000)}}
        source['subsets']['250'][0] = 'not-in-larger-prefix'
        with self.assertRaises(ValueError):
            require_common_design([source, copy.deepcopy(source)], query_ids)


if __name__ == '__main__':
    unittest.main()
