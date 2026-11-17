"""Exploratory oracle diagnostic of SPEC_EXTENSION.md (needs the model environment).

For fusion-class units, embed only the segment that aligns to the labeling
anchor and score it against the frozen panel. The cut uses the answer key, so
this is an upper bound on what within-protein splitting could recover for a
mean-pooled embedding. It is not a decomposition method or a toolkit feature.

    HF_HOME=results/component_study_v1/model-cache python experiments/component_study_v1/oracle_segment_diagnostic.py
"""
from __future__ import annotations

import json
from pathlib import Path

from build_panel import Aligner, sha
from iri_constructs.component_eval import unit_retrieval
from iri_constructs.metrics import cosine, mean_pool
from iri_constructs.model_adapters import ESM2_REVISION, LocalAdapter, prepare_device

ROOT = Path(__file__).resolve().parents[2]
V1 = ROOT / 'docs' / 'results' / 'component-study-v1'
RUNS = (('v1', V1, ROOT / 'results' / 'component_study_v1'),
        ('extension', V1 / 'extension-1', ROOT / 'results' / 'component_study_v1_ext1'))
UNIPROT = ROOT / 'data' / 'raw' / 'component_study_v1' / 'uniprot'
K, MIN_SEGMENT = 5, 30


def fasta(path):
    sequences, current = {}, None
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            current = line[1:].split()[0]
            sequences[current] = ''
        else:
            sequences[current] += line
    return sequences


def main():
    import torch
    torch.set_num_threads(2)
    adapter = LocalAdapter('esm2', prepare_device('cpu', 20261117), ESM2_REVISION, 1022, offline=True)
    aligner = Aligner()
    references = json.loads((RUNS[0][2] / 'panel.private.json').read_text())['references']
    family = {r: v['family'] for r, v in references.items()}
    reference_vectors = {r: mean_pool(adapter.residues(v['sequence'])) for r, v in references.items()}
    anchors, items, seen = {}, [], set()
    for name, record, private in RUNS:
        freeze = json.loads((record / 'panel.freeze.json').read_text())
        sequences = fasta(private / 'bundle' / 'proteins.fasta')
        whole = json.loads((private / 'scores.esm2.json').read_text())['scores']
        for u in freeze['units']:
            if not u['supported'] or u['sequence_sha256'] in seen:
                continue
            seen.add(u['sequence_sha256'])
            sequence = sequences[u['unit_id']]
            assert sha(sequence) == u['sequence_sha256']
            for fam in u['supported_families']:
                accession = u['label_evidence'][fam]['anchor']
                if accession not in anchors:
                    anchors[accession] = json.loads((UNIPROT / f'{accession}.json').read_text())['sequence']['value']
                blocks = aligner.aligner.align(sequence, anchors[accession])[0].aligned[0]
                start, end = int(blocks[0][0]), int(blocks[-1][1])
                items.append({'run': name, 'label_class': u['label_class'], 'family': fam, 'length': len(sequence),
                              'segment': sequence[start:end], 'segment_fraction': (end - start) / len(sequence),
                              'whole': unit_retrieval(whole[u['sequence_sha256']], {fam}, family, K)})
    cache = {}
    for item in items:
        segment = item.pop('segment')
        item['segment_length'] = len(segment)
        if len(segment) < MIN_SEGMENT:
            item['oracle'] = None
            continue
        if segment not in cache:
            vector = mean_pool(adapter.residues(segment))
            cache[segment] = {r: cosine(vector, v) for r, v in reference_vectors.items()}
        item['oracle'] = unit_retrieval(cache[segment], {item['family']}, family, K)
    adapter.close()

    def summary(rows):
        rows = [r for r in rows if r['oracle']]
        return {'family_slots': len(rows),
                'whole_protein_hit_at_1': sum(r['whole']['hit_at_1'] for r in rows), 'whole_protein_hit_at_5': sum(r['whole']['hit_at_k'] for r in rows),
                'oracle_segment_hit_at_1': sum(r['oracle']['hit_at_1'] for r in rows), 'oracle_segment_hit_at_5': sum(r['oracle']['hit_at_k'] for r in rows),
                'lost_whole_recovered_by_segment': sum(not r['whole']['hit_at_k'] and r['oracle']['hit_at_k'] for r in rows),
                'found_whole_lost_by_segment': sum(r['whole']['hit_at_k'] and not r['oracle']['hit_at_k'] for r in rows),
                'median_segment_fraction_of_protein': round(sorted(r['segment_fraction'] for r in rows)[len(rows) // 2], 3) if rows else None}
    out = {'status': 'Exploratory oracle diagnostic; the segment boundary comes from the label alignment (the answer key).',
           'interpretation': 'Upper bound on what within-protein splitting could recover for mean-pooled ESM-2 8M; not a method.',
           'scope': 'distinct supported sequences of v1 and extension 1; one row per (sequence, labeled family)',
           'model_calls_for_segments': len(cache),
           'by_label_class': {c: summary([r for r in items if r['label_class'] == c]) for c in ('full', 'embedded', 'multi_family')},
           'fusion_classes_by_run': {run: summary([r for r in items if r['run'] == run and r['label_class'] != 'full']) for run in ('v1', 'extension')},
           'fusion_classes_total': summary([r for r in items if r['label_class'] != 'full']),
           'fusion_slots_still_missed_by_segment': [{'family': r['family'], 'label_class': r['label_class'], 'segment_length': r['segment_length'],
                                                     'first_rank': r['oracle']['first_same_family_rank']}
                                                    for r in items if r['label_class'] != 'full' and r['oracle'] and not r['oracle']['hit_at_k']]}
    (V1 / 'extension-1' / 'exploratory.oracle-segment.json').write_text(json.dumps(out, indent=1, sort_keys=True) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
