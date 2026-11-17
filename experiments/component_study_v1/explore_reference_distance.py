"""Exploratory (not in SPEC.md): how close are the references that make retrieval succeed?

For every distinct supported single-family unit, measure identity to each
same-family reference, then ask what per-unit retrieval looks like if only
*remote* same-family references (identity below a threshold) may count, with
closer same-family references removed from the ranking. Added after the
registered results were seen, to learn whether the near-ceiling per-unit result
reflects easy references. It evaluates the model/scoring method, not preprocessing.

    python experiments/component_study_v1/explore_reference_distance.py
"""
from __future__ import annotations

import json
from pathlib import Path
from statistics import median

from build_panel import Aligner, sha
from iri_constructs.component_eval import ranked

ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / 'results' / 'component_study_v1'
RECORD = ROOT / 'docs' / 'results' / 'component-study-v1'
REMOTE = (0.40, 0.30)


def main():
    freeze = json.loads((RECORD / 'panel.freeze.json').read_text())
    references = json.loads((PRIVATE / 'panel.private.json').read_text())['references']
    sequences, current = {}, None
    for line in (PRIVATE / 'bundle' / 'proteins.fasta').read_text().splitlines():
        if line.startswith('>'):
            current = line[1:].split()[0]
            sequences[current] = ''
        else:
            sequences[current] += line
    units = {}
    for u in freeze['units']:
        if u['supported'] and u['label_class'] == 'full':
            units.setdefault(u['sequence_sha256'], (sequences[u['unit_id']], u['supported_families'][0]))
    aligner, family = Aligner(), freeze['reference_family']
    identity = {}
    for key, (sequence, fam) in units.items():
        assert sha(sequence) == key
        identity[key] = {}
        for r, f in family.items():
            if f == fam:
                alignment = aligner.aligner.align(sequence, references[r]['sequence'])[0]
                a, b = sequence, references[r]['sequence']
                same = sum(x == y for (a0, a1), (b0, b1) in zip(*alignment.aligned) for x, y in zip(a[a0:a1], b[b0:b1]))
                # Identity over the shorter sequence: unaligned residues count against it.
                identity[key][r] = round(same / min(len(a), len(b)), 4)
    nearest = sorted(max(v.values()) for v in identity.values())
    out = {'status': 'Exploratory; added after registered results were seen. Not a registered endpoint.',
           'scope': 'distinct supported single-family (full) units', 'units': len(units),
           'identity_definition': 'identical aligned residues in the best local BLOSUM62 alignment / length of the shorter sequence',
           'nearest_same_family_reference_identity': {
               'min': nearest[0], 'median': round(median(nearest), 4), 'max': nearest[-1],
               'units_by_bin': {'<0.30': sum(x < .3 for x in nearest), '0.30-0.50': sum(.3 <= x < .5 for x in nearest),
                                '0.50-0.70': sum(.5 <= x < .7 for x in nearest), '>=0.70': sum(x >= .7 for x in nearest)}},
           'remote_only': {}}
    for method in ('esm2', 'alignment'):
        scores = json.loads((PRIVATE / f'scores.{method}.json').read_text())['scores']
        out['remote_only'][method] = {}
        for threshold in REMOTE:
            ranks = []
            for key, (_, fam) in units.items():
                remote = {r for r, v in identity[key].items() if v < threshold}
                if not remote:
                    continue
                keep = {r: s for r, s in scores[key].items() if family[r] != fam or r in remote}
                order = ranked(keep)
                ranks.append(next(i + 1 for i, r in enumerate(order) if r in remote))
            out['remote_only'][method][f'identity_below_{threshold}'] = {
                'units_with_a_remote_reference': len(ranks), 'hit_at_1': sum(r == 1 for r in ranks),
                'hit_at_5': sum(r <= 5 for r in ranks), 'median_first_rank': median(ranks) if ranks else None,
                'worst_first_rank': max(ranks) if ranks else None}
    (RECORD / 'exploratory.reference-distance.json').write_text(json.dumps(out, indent=1, sort_keys=True) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
