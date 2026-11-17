"""Regenerate the fictional GenBank fixtures in examples/synthetic/.

Every sequence is produced by a seeded random generator. They are software
fixtures: they exercise parser behavior and are not biological observations.
EXPECTED.json states, by construction, the outcome each feature must get; the
tests compare the toolkit's output with those declared expectations.

    python examples/make_synthetic_examples.py
"""
from __future__ import annotations

import io
import json
from pathlib import Path
import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

OUT = Path(__file__).resolve().parent / 'synthetic'
RNG = random.Random(20261117)
CODON = {'A': 'GCT', 'C': 'TGC', 'D': 'GAC', 'E': 'GAA', 'F': 'TTC', 'G': 'GGT', 'H': 'CAC', 'I': 'ATC',
         'K': 'AAA', 'L': 'CTG', 'M': 'ATG', 'N': 'AAC', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGT', 'S': 'TCT',
         'T': 'ACC', 'V': 'GTT', 'W': 'TGG', 'Y': 'TAC'}


def protein(n):
    return 'M' + ''.join(RNG.choice('ACDEFGHIKLNPQRSTVWY') for _ in range(n - 1))


def coding(aa, stop='TAA'):
    return ''.join(CODON[a] for a in aa) + stop


def spacer(n):
    return ''.join(RNG.choice('ACGT') for _ in range(n))


def revcomp(dna):
    return str(Seq(dna).reverse_complement())


class Builder:
    def __init__(self, name, topology='linear', definition='Fictional software fixture; not a biological record'):
        self.name, self.topology, self.definition = name, topology, definition
        self.dna, self.features, self.expected = '', [], {}

    def add(self, dna):
        start = len(self.dna)
        self.dna += dna
        return start, len(self.dna)

    def feature(self, location, kind, gene, expect, **qualifiers):
        q = {'gene': [gene], **{k: v if isinstance(v, list) else [v] for k, v in qualifiers.items()}}
        self.features.append(SeqFeature(location, type=kind, qualifiers=q))
        self.expected[gene] = expect

    def cds(self, gene, aa, expect='translation_validated', *, strand=1, translation='same', stop='TAA', **qualifiers):
        dna = coding(aa, stop)
        start, end = self.add(dna if strand == 1 else revcomp(dna))
        if translation == 'same':
            qualifiers['translation'] = aa
        elif translation is not None:
            qualifiers['translation'] = translation
        self.feature(SimpleLocation(start, end, strand=strand), 'CDS', gene, expect,
                     product='fictional protein ' + gene, codon_start='1', transl_table='11', **qualifiers)
        self.add(spacer(30))
        return start, end

    def text(self, record_id=None):
        record = SeqRecord(Seq(self.dna), id=record_id or self.name + '.1', name=self.name, description=self.definition)
        record.annotations = {'molecule_type': 'DNA', 'topology': self.topology, 'data_file_division': 'SYN',
                              'date': '17-NOV-2026', 'organism': 'synthetic construct', 'source': 'synthetic construct',
                              'taxonomy': ['other sequences', 'artificial sequences']}
        record.features = [SeqFeature(SimpleLocation(0, len(self.dna), strand=1), type='source',
                                      qualifiers={'organism': ['synthetic construct'], 'mol_type': ['other DNA']})] + self.features
        handle = io.StringIO()
        SeqIO.write(record, handle, 'genbank')
        return handle.getvalue()


def single():
    b = Builder('FIXSINGLE')
    s, e = b.add(spacer(40))
    b.feature(SimpleLocation(s, e, strand=1), 'promoter', 'fixP', 'noncoding_metadata', note='fictional promoter')
    b.cds('fixA', protein(60))
    s, e = b.add(spacer(35))
    b.feature(SimpleLocation(s, e, strand=1), 'terminator', 'fixT', 'noncoding_metadata')
    return b


def edge_cases():
    b = Builder('FIXEDGE', topology='circular')
    # Origin-spanning CDS: its 3' part is the first bases of the record, its 5' part the last.
    wrap_aa = protein(40)
    wrap = coding(wrap_aa)
    tail = wrap[66:]
    b.add(tail)
    b.add(spacer(30))
    b.cds('fwd', protein(50))
    b.cds('rev', protein(45), strand=-1)
    b.cds('notrans', protein(30), 'missing_translation', translation=None)
    wrong = protein(30)
    b.cds('mismatch', protein(30), 'translation_mismatch', translation=wrong)
    b.cds('nostop', protein(30), 'incomplete_cds', stop='')
    b.cds('pseudo', protein(30), 'exceptional_cds_annotation', pseudo='')
    aa = protein(30)
    b.cds('noncanon', aa, 'noncanonical_translation_residues', translation=aa[:10] + 'X' + aa[11:])
    b.cds('alt_table', 'M' + 'W' + protein(20)[1:], 'translation_validated')  # replaced below with a table-4 codon
    # Table 4 reads TGA as Trp; rewrite the second codon and the qualifier.
    f = b.features[-1]
    s = int(f.location.start)
    b.dna = b.dna[:s + 3] + 'TGA' + b.dna[s + 6:]
    f.qualifiers['transl_table'] = ['4']
    # 5'-partial CDS: fuzzy start, so no exact protein can be asserted.
    partial = protein(31)
    s, e = b.add(coding(partial)[3:])
    b.feature(SimpleLocation(BeforePosition(s), e, strand=1), 'CDS', 'partial5', 'location_fuzzy',
              product='fictional protein partial5', codon_start='1', transl_table='11', translation=partial[1:])
    b.add(spacer(30))
    # Two-segment reverse-strand join (segments listed in biological order).
    split_aa = protein(40)
    split = coding(split_aa)
    first, second = split[:60], split[60:]
    s2, e2 = b.add(revcomp(second))
    b.add(spacer(25))
    s1, e1 = b.add(revcomp(first))
    b.feature(CompoundLocation([SimpleLocation(s1, e1, strand=-1), SimpleLocation(s2, e2, strand=-1)], operator='join'),
              'CDS', 'revjoin', 'translation_validated', product='fictional protein revjoin',
              codon_start='1', transl_table='11', translation=split_aa)
    b.add(spacer(30))
    s, e = b.add(spacer(20))
    for _ in range(2):  # identical duplicates must both survive, numbered
        b.features.append(SeqFeature(SimpleLocation(s, e, strand=1), type='misc_feature', qualifiers={'gene': ['dup'], 'note': ['fictional duplicate']}))
    b.expected['dup'] = 'noncoding_metadata'
    b.feature(SimpleLocation(s + 5, s + 5, strand=1), 'misc_feature', 'site', 'location_zero_length', note='fictional between-bases site')
    b.feature(CompoundLocation([SimpleLocation(s, s + 5, strand=1), SimpleLocation(s + 10, s + 15, strand=1)], operator='order'),
              'misc_feature', 'ordered', 'location_unsupported_operator')
    b.feature(SimpleLocation(0, 12, strand=1, ref='FIXREMOTE.1'), 'misc_feature', 'remote', 'location_remote')
    b.add(spacer(30))
    head_start, head_end = b.add(wrap[:66])
    b.feature(CompoundLocation([SimpleLocation(head_start, head_end, strand=1), SimpleLocation(0, len(tail), strand=1)], operator='join'),
              'CDS', 'wrap', 'translation_validated', product='fictional protein wrap',
              codon_start='1', transl_table='11', translation=wrap_aa)
    return b


def main():
    OUT.mkdir(exist_ok=True)
    expected = {}
    one = single()
    (OUT / 'single_cds.gb').write_text(one.text())
    expected['FIXSINGLE.1'] = one.expected
    edge = edge_cases()
    text = edge.text()
    # A feature past the end of the sequence cannot be written by Biopython; add it as text.
    oob = f'     misc_feature    {len(edge.dna) - 5}..{len(edge.dna) + 40}\n                     /gene="beyond"\n'
    text = text.replace('ORIGIN', oob + 'ORIGIN', 1)
    edge.expected['beyond'] = 'location_out_of_bounds'
    (OUT / 'edge_cases.gb').write_text(text)
    expected['FIXEDGE.1'] = edge.expected
    long = Builder('FIXLONG')
    long.cds('long', protein(1100))
    (OUT / 'long_protein.gb').write_text(long.text())
    expected['FIXLONG.1'] = long.expected
    # Multi-record file: valid, malformed (broken location), valid.
    first, third = Builder('FIXMULTIA'), Builder('FIXMULTIC')
    first.cds('multiA', protein(40))
    third.cds('multiC1', protein(35))
    third.cds('multiC2', protein(42), strand=-1)
    bad = Builder('FIXBROKEN')
    bad.cds('broken', protein(20))
    start = int(bad.features[0].location.start) + 1
    broken = bad.text().replace(f'CDS             {start}..', 'CDS             garbage(', 1)
    assert 'garbage(' in broken
    (OUT / 'multi_record_one_bad.gb').write_text(first.text() + broken + third.text())
    expected['FIXMULTIA.1'], expected['FIXMULTIC.1'] = first.expected, third.expected
    # A tool export without accession/version, as plasmid editors commonly write.
    export = Builder('Exported')
    export.cds('exp', protein(33))
    (OUT / 'no_accession_export.gb').write_text(export.text(record_id='.'))
    expected['Exported'] = export.expected
    (OUT / 'sequence_only.fasta').write_text('>FIXFASTA fictional sequence-only fixture\n' + spacer(120) + '\n')
    (OUT / 'not_a_sequence.gb').write_text('This file is not a sequence record.\n')
    (OUT / 'EXPECTED.json').write_text(json.dumps(expected, indent=2, sort_keys=True) + '\n')


if __name__ == '__main__':
    main()
