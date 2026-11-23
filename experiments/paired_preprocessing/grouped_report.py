"""Summarize paired recall with the evaluator's construct groups."""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from .verify import SCHEMA, _index, _int, _require

GROUP_SCHEMA = 'paired-recall-groups-1'


def _group_manifest(manifest):
    _require(manifest['schema'] == GROUP_SCHEMA, f'expected group schema {GROUP_SCHEMA}')
    records = _index(manifest['constructs'], 'id', 'group manifest')
    partitions = {}
    for row in records.values():
        group, partition = row['group'], row['partition']
        _require(isinstance(group, str) and group.strip() and isinstance(partition, str) and partition.strip(),
                 'each construct needs a group and partition')
        _require(partitions.setdefault(group, partition) == partition, 'a group crosses partitions')
    return records


def _counts(run, cap):
    tier = run[str(cap)] if str(cap) in run else run[cap]
    rows = _index(tier['per_construct'], 'construct_id', 'paired report')
    for row in rows.values():
        n, recovered = row['eligible_components'], row['recovered_components']
        _require(_int(n, 0) and _int(recovered, 0) and recovered <= n, 'invalid component counts')
    return rows


def _compare(base, candidate, manifest):
    _require(base.keys() == candidate.keys(), 'candidate and baseline construct sets differ')
    _require(base.keys() <= manifest.keys(), 'not all scored constructs have a group')
    partitions = {manifest[key]['partition'] for key in base}
    _require(len(partitions) == 1, 'a comparison must use one partition')
    eligible, differences = Counter(), Counter()
    for key, row in base.items():
        other = candidate[key]
        n = row['eligible_components']
        _require(n == other['eligible_components'], 'candidate and baseline eligible counts differ')
        if n:
            group = manifest[key]['group']
            eligible[group] += n
            differences[group] += other['recovered_components'] - row['recovered_components']

    count, difference = sum(eligible.values()), sum(differences.values())
    _require(count > 0, 'no eligible documented components')
    # Leave one group out of the pooled estimate to show sensitivity to individual groups.
    omissions = [(difference - differences[g]) / (count - n) for g, n in eligible.items() if n < count]
    return {'constructs': len(base), 'groups_with_eligible_components': len(eligible), 'eligible_components': count,
            'largest_group_share_of_components': max(eligible.values()) / count,
            'pooled_recall_gain': difference / count,
            'equal_group_recall_gain': sum(differences[g] / n for g, n in eligible.items()) / len(eligible),
            'groups_improved': sum(x > 0 for x in differences.values()),
            'groups_declined': sum(x < 0 for x in differences.values()),
            'groups_tied': sum(x == 0 for x in differences.values()),
            'leave_one_group_out_pooled_range': [min(omissions), max(omissions)] if omissions else None}


def grouped_report(report, groups):
    """Report pooled and equally weighted group effects at the existing query caps."""
    _require(report['schema'] == SCHEMA, f'expected paired report schema {SCHEMA}')
    manifest = _group_manifest(groups)
    baseline = report['baseline']
    reference = report['annotation_reference']
    caps = report['protocol']['query_caps']
    base_counts = {cap: _counts(report['runs'][baseline], cap) for cap in caps}
    first = base_counts[caps[0]]
    _require(all(rows.keys() == first.keys() and all(
        row['eligible_components'] == first[key]['eligible_components'] for key, row in rows.items())
        for rows in base_counts.values()), 'baseline construct set or eligible counts differ across caps')
    candidates = {name: run for name, run in report['runs'].items() if name not in (baseline, reference)}
    _require(candidates, 'no candidate method to compare')
    return {'baseline': baseline, 'primary_query_cap': report['protocol']['primary_query_cap'],
            'candidate_minus_baseline': {
                name: {cap: _compare(base_counts[cap], _counts(run, cap), manifest) for cap in caps}
                for name, run in candidates.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--groups', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    try:
        report_bytes, group_bytes = args.report.read_bytes(), args.groups.read_bytes()
        summary = grouped_report(json.loads(report_bytes), json.loads(group_bytes))
        summary['input_sha256'] = {'report': hashlib.sha256(report_bytes).hexdigest(),
                                   'groups': hashlib.sha256(group_bytes).hexdigest()}
        with args.out.open('x', encoding='utf-8') as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write('\n')
    except (OSError, ValueError, TypeError, KeyError, IndexError) as exc:
        parser.exit(2, f'Cannot summarize paired groups: {exc}\n')


if __name__ == '__main__':
    main()
