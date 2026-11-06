"""Checks for scientific correctness and data boundaries; no network/model download."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from iri_constructs.t05_data import (build_query, load_catalogue, local_data_path,
                                     read_single_fasta, sha256, validate_protein)
from iri_constructs.schema_check import check_schema_example

try:
    import torch
except ImportError:
    torch = None


class InputTests(unittest.TestCase):
    def test_length_boundary_and_invalid_characters(self):
        validate_protein('A' * 1022)  # Deliberately artificial test fixture, not a source record.
        for sequence in ['', 'A' * 1023, 'A*', '<mask>', 'aa', 'A C']:
            with self.subTest(length=len(sequence)), self.assertRaises(ValueError):
                validate_protein(sequence)

    def test_multirecord_fasta_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            p = Path(directory) / 'fixture.fasta'
            p.write_text('>one\nAC\n>two\nDE\n')
            with self.assertRaises(ValueError):
                read_single_fasta(p, 'one')

    def test_paths_cannot_escape_data_directories(self):
        for path in ['README.md','/tmp/protein.fasta','data/raw/../../outside.fasta']:
            with self.subTest(path=path), self.assertRaises(ValueError):
                local_data_path(ROOT,path)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root/'data/raw').mkdir(parents=True)
            (root/'data/raw/escape').symlink_to(root)
            with self.assertRaises(ValueError):
                local_data_path(root,'data/raw/escape/file.fasta')

    def test_known_boundary_exact_and_unique(self):
        sequences = {'target':'ACD','other':'EF'}  # Artificial toy letters only.
        query, oracle, span = build_query(['other','target'],'target',sequences)
        self.assertEqual(span,(2,5))
        self.assertEqual(query[slice(*span)],oracle)
        self.assertEqual(oracle,sequences['target'])
        with self.assertRaises(ValueError):
            build_query(['target','target'],'target',sequences)

    def test_governance_rejects_pending_before_reading_files(self):
        record = copy.deepcopy(json.loads((ROOT/'data/manifests/t05_benign_catalogue.json').read_text())['records'][0])
        for gate in ['review_status','licensing_status','safety_review']:
            candidate = copy.deepcopy(record); candidate[gate] = 'pending'
            with self.subTest(gate=gate), self.assertRaisesRegex(ValueError,'unapproved'):
                load_catalogue(ROOT,{'manifest_version':'t05-benign-v1','records':[candidate]})

    def test_corrupted_source_checksum_rejected(self):
        record = copy.deepcopy(json.loads((ROOT/'data/manifests/t05_benign_catalogue.json').read_text())['records'][0])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            p=root/record['raw_path']; p.parent.mkdir(parents=True); p.write_bytes(b'changed')
            lic=root/record['license_evidence_path']; lic.write_bytes(b'fixture')
            record['license_evidence_sha256']=sha256(b'fixture')
            with self.assertRaisesRegex(ValueError,'checksum mismatch'):
                load_catalogue(root,{'manifest_version':'t05-benign-v1','records':[record]})

    def test_schema_v02_protein_and_missing_reference(self):
        # Start with existing fictional fixture; never embed a real protein in tests.
        doc=json.loads((ROOT/'data/examples/example_single_cds.json').read_text())
        doc['schema_version']='0.2'
        parent=doc['construct']; parent.update(molecule_type='protein',sequence_length_nt=None,
            sequence_length_aa=1100,coordinate_system='protein_0_based_half_open',primary_archetype='unknown')
        ref={'path':'data/processed/fictional.fasta','format':'fasta','record_id':'fictional',
             'sha256':'0'*64,'length':1100,'alphabet':'protein'}
        parent['sequence_ref']=ref
        doc['features'][0]['coordinate_system']='protein_0_based_half_open'
        doc['features'][0]['location']['strand']=None
        unit=doc['units'][0]; unit.update(model_input_status='ready',sequence_ref=ref,
                                         coordinate_system='protein_0_based_half_open')
        unit['source_segments'][0]['strand']=None
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'example.json'; p.write_text(json.dumps(doc))
            self.assertEqual(check_schema_example(p),[])
            del unit['sequence_ref']; p.write_text(json.dumps(doc))
            self.assertTrue(any('sequence_ref' in error for error in check_schema_example(p)))
            unit['sequence_ref']=ref; parent['sequence_length_nt']=3300
            p.write_text(json.dumps(doc))
            self.assertTrue(any('nucleotide' in error for error in check_schema_example(p)))


@unittest.skipIf(torch is None,'Install requirements-t05.lock.txt for tensor checks')
class TensorTests(unittest.TestCase):
    def test_pool_excludes_special_and_padding_sentinels(self):
        from iri_constructs.t05_embeddings import mean_pool
        hidden=torch.tensor([[[999.,999.],[2.,4.],[4.,8.],[777.,777.],[888.,888.]]])
        pooled,mask=mean_pool(hidden,torch.tensor([[1,1,1,1,0]]),torch.tensor([[1,0,0,1,0]]))
        self.assertTrue(torch.equal(pooled,torch.tensor([[3.,6.]])))
        self.assertEqual(mask.sum().item(),2)
        with self.assertRaises(ValueError):
            mean_pool(hidden,torch.zeros(1,5),torch.ones(1,5))

    def test_rank_correct_target_and_deterministic_ties(self):
        from iri_constructs.t05_embeddings import rank_cosine
        query=torch.tensor([1.,0.]); refs=torch.tensor([[0.,1.],[1.,0.],[1.,0.]])
        rows=rank_cosine(query,refs,['other','B','A'])
        self.assertEqual([r['accession'] for r in rows],['A','B','other'])
        self.assertEqual(rows[0]['exact_score_tie_count'],2)
        with self.assertRaises(ValueError): rank_cosine(torch.zeros(2),refs,['other','B','A'])
        with self.assertRaises(ValueError): rank_cosine(torch.tensor([float('nan'),0]),refs,['other','B','A'])


if __name__ == '__main__':
    unittest.main()
