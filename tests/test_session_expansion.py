import json
from pathlib import Path
import sys
import tempfile
import unittest
from urllib.parse import urlencode
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from iri_constructs.session_expansion import _eligible, _fasta, _response, checked_path
from iri_constructs.session_store import digest


class ExpansionSourceTests(unittest.TestCase):
    def test_source_response_checks_bytes_and_requested_fields(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); path = root / 'source.tsv'
            body = b'Entry\tLength\nP00001\t8\n'
            path.write_bytes(body)
            provenance = {'sha256': digest(body), 'bytes': len(body), 'status': 200,
                          'retrieved_utc': '2026-10-31T00:00:00Z',
                          'url': 'https://rest.uniprot.org/uniprotkb/search?' + urlencode({'fields': 'accession,length', 'format': 'tsv'})}
            path.with_suffix('.tsv.provenance.json').write_text(json.dumps(provenance))
            labels = {'accession': 'Entry', 'length': 'Length'}
            self.assertEqual(_response(root, 'source.tsv', ['accession', 'length'], labels)[1],
                             [{'accession': 'P00001', 'length': '8'}])
            with self.assertRaisesRegex(ValueError, 'requested fields'):
                _response(root, 'source.tsv', ['accession'], labels)
            path.write_bytes(body.replace(b'8', b'9'))
            with self.assertRaisesRegex(ValueError, 'saved bytes'):
                _response(root, 'source.tsv', ['accession', 'length'], labels)

    def test_path_escape_and_noncanonical_fasta_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaisesRegex(ValueError, 'escapes'):
                checked_path(root, '../outside')
            path = root / 'canonical.fasta'; path.write_text('>A\nACDE\n>A\nACDE\n')
            with self.assertRaisesRegex(ValueError, 'identifier'):
                _fasta(path)
            path.write_text('>A\nACDX\n')
            with self.assertRaisesRegex(ValueError, 'Noncanonical'):
                _fasta(path)

    def test_mechanical_source_gate_does_not_admit_unknown_identity(self):
        row = {'accession': 'P00001', 'reviewed': 'reviewed', 'organism_id': '559292',
               'fragment': '', 'version': '1', 'sequence_version': '1', 'length': '8',
               'protein_existence': 'Evidence at protein level'}
        _eligible(row)
        with self.assertRaisesRegex(ValueError, 'Uncertain'):
            _eligible({**row, 'protein_existence': 'Uncertain'})
        with self.assertRaisesRegex(ValueError, 'source identity'):
            _eligible({**row, 'organism_id': '1'})


if __name__ == '__main__': unittest.main()
