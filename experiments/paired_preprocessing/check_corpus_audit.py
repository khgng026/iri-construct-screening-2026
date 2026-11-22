"""Check the saved acquisition and audit metadata without loading sequences."""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

from .verify import _covered, _location


def _expected_records(acquisition, findings):
    expected = {}
    files = 0
    for sample, entry in acquisition['samples'].items():
        for batch in entry['files']:
            files += 1
            if batch['byte_identical'] is not True or batch['sha256'] != batch['frozen_sha256']:
                findings['nonidentical_batch_files'] += 1
            if Counter(batch['frozen_versions']) != Counter(batch['returned_versions']):
                findings['batch_version_list_mismatches'] += 1
            for accession in batch['frozen_versions']:
                key = (sample, accession)
                if key in expected:
                    findings['duplicate_expected_records'] += 1
                expected[key] = batch['file']
    if not expected:
        raise ValueError('acquisition metadata contains no expected records')
    return expected, files


def _named_regions(record):
    named = {}
    for feature in record['named_subregions']:
        if feature['location_status'] != 'exact' or feature['strand'] not in ('+', '-'):
            continue
        locus = _location({'strand': {'+': 1, '-': -1}[feature['strand']], 'spans': feature['segments']},
                          record['length'], 'named subregion')
        for family in feature['named_families']:
            named.setdefault(family, []).append(locus)
    return named


def _count_regions(record, unit, named_regions, counts, findings):
    parent = _location({'strand': unit['strand'], 'spans': unit['segments']}, record['length'], 'exported CDS')
    for family, evidence in unit['anchor_evidence'].items():
        counts['anchor_regions'] += 1
        retranslated = evidence['dna_region_retranslates'] is True
        counts['anchor_regions_retranslate'] += retranslated
        if not retranslated:
            findings['failed_region_retranslation_flags'] += 1
        try:
            region = _location({'strand': parent[0], 'spans': evidence['dna_spans']},
                               record['length'], 'anchor region')
        except (ValueError, TypeError, KeyError):
            findings['invalid_region_locations'] += 1
            continue
        if _covered(region, parent) != 1:
            findings['regions_outside_parent_cds'] += 1
            continue
        start, end = evidence['unit_region_aa']
        area = sum(right - left for left, right in region[1])
        if not (type(start) is int and type(end) is int and 0 <= start < end <= unit['length_aa']
                and area == 3 * (end - start)):
            findings['region_length_mismatches'] += 1
            continue
        if retranslated and any(_covered(region, named) == _covered(named, region) == 1
                                for named in named_regions.get(family, ())):
            counts['anchor_regions_matching_a_named_source_feature'] += 1


def _shared_protein_groups(fingerprints):
    """Summarize connected groups formed by sharing any reported protein hash."""
    parents = list(range(len(fingerprints)))

    def root(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    owner = {}
    for index, proteins in enumerate(fingerprints):
        for protein in sorted(proteins):
            other = owner.setdefault(protein, index)
            parents[root(index)] = root(other)
    groups = Counter(root(index) for index in range(len(fingerprints)))
    return {'groups': len(groups), 'largest_group': max(groups.values(), default=0),
            'distinct_protein_hashes': len(owner)}


def check_audit(acquisition, audit):
    """Check persisted report consistency and return aggregate findings."""
    if acquisition['schema'] != 'paired-corpus-acquisition-1' or audit['schema'] != 'paired-corpus-audit-1':
        raise ValueError('unexpected acquisition or audit schema')
    counts, findings, seen = Counter(), Counter(), Counter()
    expected, files = _expected_records(acquisition, findings)
    fingerprints = []
    for record in audit['constructs']:
        key = (record['sample'], record['accession'])
        seen[key] += 1
        counts['records'] += 1
        if record['status'] != 'ok':
            findings['records_not_processed'] += 1
            continue
        if key in expected and record['batch_file'] != expected[key]:
            findings['records_in_unexpected_batch'] += 1
        named_regions = _named_regions(record)
        proteins = set()
        for unit in record['cds']:
            if unit['extraction_state'] != 'exported':
                continue
            counts['exported_units'] += 1
            proteins.add(unit['protein_sha256'])
            for field in ('frozen_unit_match', 'frozen_label_match'):
                if unit.get(field) is True:
                    counts[field] += 1
                else:
                    findings[f'{field}_missing_or_false'] += 1
            _count_regions(record, unit, named_regions, counts, findings)
        fingerprints.append(proteins)

    findings['missing_audit_records'] = len(expected.keys() - seen.keys())
    findings['unexpected_audit_records'] = len(seen.keys() - expected.keys())
    findings['duplicate_audit_records'] = sum(n - 1 for n in seen.values())
    for field in ('records', 'exported_units', 'frozen_unit_match', 'frozen_label_match', 'anchor_regions',
                  'anchor_regions_retranslate'):
        if audit['checks'].get(field, 0) != counts[field]:
            findings['counter_mismatches'] += 1
    return {
        'expected_batch_files': files,
        'expected_records': len(expected),
        'saved_counts': {**dict(counts), 'anchor_regions_matching_a_named_source_feature':
                         counts['anchor_regions_matching_a_named_source_feature']},
        'findings': {key: value for key, value in findings.items() if value},
        'shared_protein_group_diagnostics': _shared_protein_groups(fingerprints),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--acquisition', type=Path, required=True)
    parser.add_argument('--audit', type=Path, required=True)
    args = parser.parse_args()
    try:
        acquisition = json.loads(args.acquisition.read_text(encoding='utf-8'))
        audit = json.loads(args.audit.read_text(encoding='utf-8'))
        report = check_audit(acquisition, audit)
    except (OSError, ValueError, TypeError, KeyError) as exc:
        parser.exit(2, f'Unable to check audit metadata: {exc}\n')
    print(json.dumps(report, indent=2, sort_keys=True))
    raise SystemExit(2 if report['findings'] else 0)


if __name__ == '__main__':
    main()
