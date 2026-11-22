#!/usr/bin/env python3
"""Independent metadata/rank audit of pinned Nishit Kalani public artifacts.

No sequence generation, model inference, retrieval scoring, or source edits.
Run from the session directory: python scripts/audit_nishit_saved.py
"""
from __future__ import annotations
import csv
import hashlib
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

SESSION = Path(__file__).resolve().parents[1]
SOURCE = SESSION / 'sources/nishit-working'
RELEASE = SESSION / 'sources/nishit-benchmark'
OUT = SESSION / 'evidence/nishit-audit'


def rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def measures(pairs, depth=10):
    """All-case recall/MRR; conditional rank with explicit coverage, no zero imputation."""
    n = len(pairs)
    have = [rank for _, rank in pairs if rank is not None]
    if any(rank < 1 for rank in have):
        raise ValueError('Ranks must be positive or None')
    groups = defaultdict(list)
    for payload, rank in pairs:
        groups[payload].append(rank)
    means = [statistics.mean(r for r in rr if r is not None)
             for rr in groups.values() if any(r is not None for r in rr)]
    return {
        'n': n, 'n_scorable': len(have), 'n_unscorable': n-len(have),
        'n_payloads': len(groups), 'n_payloads_with_rank': len(means),
        'recall_at_10': sum(r <= depth for r in have)/n if n else None,
        'mrr': sum(1/r for r in have)/n if n else None,
        'mean_rank_conditional': statistics.mean(have) if have else None,
        'mean_rank_macro_conditional': statistics.mean(means) if means else None,
        'upstream_zero_imputed_mean_rank_macro': sum(means)/len(groups) if groups else None,
        'recall_at_10_macro': statistics.mean(
            sum(r is not None and r <= depth for r in rr)/len(rr) for rr in groups.values()
        ) if groups else None,
        'mrr_macro': statistics.mean(
            sum(1/r for r in rr if r is not None)/len(rr) for rr in groups.values()
        ) if groups else None,
    }


def dump_csv(path, records):
    if not records:
        return
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader(); w.writerows(records)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    bd = SOURCE/'benchmark/v1'
    manifest = json.loads((bd/'MANIFEST.json').read_text())
    checks=[]
    for name, data in manifest['files'].items():
        if not isinstance(data, dict):
            continue
        checks.append({'file': name, 'expected_sha256': data['sha256'],
                       'actual_sha256': sha(bd/name), 'matches':sha(bd/name)==data['sha256']})
    # The independent release repo can have its benchmark at a different relative root.
    release_manifests=list(RELEASE.rglob('MANIFEST.json'))
    release_checks=[]
    for mf in release_manifests:
        for name, data in json.loads(mf.read_text())['files'].items():
            if isinstance(data, dict):
                release_checks.append({'manifest':str(mf.relative_to(RELEASE)), 'file':name,
                    'matches':sha(mf.parent/name)==data['sha256'],
                    'same_as_working_release':sha(mf.parent/name)==sha(bd/name)})
    missing=[]; canonical=[]
    for kind in ['alignment','embedding']:
        baseline=json.loads((bd/f'baseline_{kind}.json').read_text(encoding='utf-8-sig'))
        for item in baseline['provenance']:
            missing.append({'kind':kind, **item, 'exists_in_clone':(SOURCE/item['source']).exists()})
        for r in rows(bd/f'baseline_{kind}.csv'):
            if r['tier']=='synth_real' and r['stratum']=='f < 0.50' and float(r['threshold']) in (1.,.25):
                canonical.append({'evidence_level':'saved_author_summary_only', **r})
    dump_csv(OUT/'canonical_synth_real_summary.csv',canonical)
    # Gated target-availability counts do not require sequence inference.
    targets=defaultdict(set)
    catalogue=set((bd/'catalogue_accessions.txt').read_text().splitlines())
    for r in rows(bd/'homologs.csv'):
        if float(r['identity'])<=.25 and r['homolog_accession'] in catalogue:
            targets[r['payload_accession']].add(r['homolog_accession'])
    items=[r for r in rows(bd/'items.csv') if r['tier']=='synth_real' and r['stratum']=='f < 0.50']
    per=Counter(r['payload_accession'] for r in items)
    absent={a:n for a,n in per.items() if not targets[a]}
    summary_effect=[]
    for r in canonical:
        if float(r['threshold']) != .25:
            continue
        # Here all 24 payload groups have six items and seven groups have no target.
        # Conditional micro=conditional macro exactly before rounding in this balanced stratum.
        mean=float(r['mean_rank']); macro=float(r['mean_rank_macro'])
        summary_effect.append({'kind':r['arm_kind'],'arm':r['arm'],
            'published_conditional_micro_mean_rank':mean,
            'published_zero_imputed_macro_mean_rank':macro,
            'expected_zero_imputed_macro_from_micro':round(mean*(len(per)-len(absent))/len(per),2),
            'n_payloads':len(per),'n_payloads_no_target':len(absent),
            'zero_imputation_shrinkage_fraction':len(absent)/len(per),
            'evidence_level':'algebra_from_saved_summary_and_independent_target_availability'})
    dump_csv(OUT/'canonical_mean_rank_consequence.csv',summary_effect)
    # Recompute every shipped per-item twilight CSV independently. These are older/other runs,
    # explicitly not substitutions for the missing canonical raw prediction files.
    raw_inventory=[]; recomputed=[]; summary_checks=[]
    for f in sorted((SOURCE/'results').rglob('twilight*.csv')):
        rr=rows(f); cells=defaultdict(list)
        for r in rr:
            rank=int(r['best_target_rank']) if r['best_target_rank'] else None
            cells[(float(r['threshold']),r['arm'])].append((r['payload_accession'],rank))
        raw_inventory.append({'path':str(f.relative_to(SOURCE)),'sha256':sha(f),
            'n_prediction_rows':len(rr),'n_constructs':len({r['plasmid_id'] for r in rr}),
            'n_cells':len(cells),'n_payloads':len({r['payload_accession'] for r in rr}),
            'stores_full_ranked_lists':False,'stores_target_rank_only':True})
        summary=f.with_name(f.stem+'_summary.json')
        ss=json.loads(summary.read_text(encoding='utf-8-sig')) if summary.exists() else None
        for (thr,arm),pairs in sorted(cells.items()):
            m=measures(pairs)
            recomputed.append({'source':str(f.relative_to(SOURCE)), 'threshold':thr,'arm':arm,
                'evidence_level':'metrics_recomputed_from_saved_target_ranks',**m})
            if ss is None:
                continue
            saved=ss.get('by_threshold',{}).get(f'{thr:.2f}',{})
            have=[r for _,r in pairs if r is not None]
            # Historical twilight summaries intentionally drop unscorable rows.
            expected={arm+'_top1':sum(r==1 for r in have),
                      arm+'_top10':sum(r<=10 for r in have),arm+'_n':len(have),
                      arm+'_mean_rank':round(statistics.mean(have),2) if have else None,
                      arm+'_mrr':round(sum(1/r for r in have)/len(have),3) if have else None}
            for name,value in expected.items():
                if name in saved:
                    summary_checks.append({'source':str(f.relative_to(SOURCE)),'threshold':thr,
                        'field':name,'saved':saved[name],'recomputed':value,'matches':saved[name]==value})
    dump_csv(OUT/'saved_rank_recomputed.csv',recomputed)
    dump_csv(OUT/'historical_summary_checks.csv',summary_checks)
    # Concrete numerical counterexample: adding one unavailable payload makes old mean rank
    # appear twice as good, even though available ranks are unchanged. Recall/MRR correctly halve.
    minimal={'before':measures([('available',8)]),
             'after':measures([('available',8),('unavailable',None)])}
    payload_rows=[]
    metadata_path=OUT/'payload_metadata.tsv'
    metadata={}
    if metadata_path.exists():
        metadata={r['Entry']:r for r in csv.DictReader(metadata_path.open(),delimiter='\t')}
    for a,n in sorted(per.items()):
        payload_rows.append({'accession':a,'n_f_lt_half_items':n,'n_targets_at_025':len(targets[a]),
            'protein_name':metadata.get(a,{}).get('Protein names',''),
            'organism':metadata.get(a,{}).get('Organism',''),
            'metadata_only_not_safety_admission':True})
    dump_csv(OUT/'synthetic_payload_inventory.csv',payload_rows)
    result={'source_revision':'5dfef4734f07b981fa6305383fee170d6af28a6d',
        'release_revision':'aec12a9488631194892240f7b30ee15e652aed08',
        'working_release_file_hash_checks':checks,'other_release_file_checks':release_checks,
        'canonical_raw_source_availability':missing,
        'canonical_synth_real_strict_target_availability':{'n':len(items),'n_payloads':len(per),
            'items_per_payload':dict(per),'n_unscorable':sum(absent.values()),
            'n_unscorable_payloads':len(absent),'unscorable_payloads':absent},
        'mean_rank_counterexample':minimal,'raw_rank_inventory':raw_inventory,
        'recomputed_prediction_rows':sum(x['n_prediction_rows'] for x in raw_inventory),
        'recomputed_metric_cells':len(recomputed),'historical_summary_checks':len(summary_checks),
        'historical_summary_mismatches':[x for x in summary_checks if not x['matches']],
        'inference_performed':False,'full_score_matrices_found':False,
        'limits':['Canonical per-item maskfix predictions are absent; canonical summary percentages are parsed, not rerun.',
                  'Other saved twilight runs differ in catalogue/protocol; their numerical audit is separate evidence.',
                  'Target-rank-only records cannot reconstruct a fixed total top10 merged list.',
                  'Published benchmark labels one designated payload, not every biological component.',
                  'Sequence-derived labels and catalogue holdout do not establish biological function or pretraining novelty.']}
    (OUT/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'hash_checks_pass':all(x['matches'] for x in checks),
        'canonical_raw_files_present':sum(x['exists_in_clone'] for x in missing),
        'canonical_raw_files_expected':len(missing),'saved_prediction_rows':result['recomputed_prediction_rows'],
        'metric_cells':len(recomputed),'summary_checks':len(summary_checks),
        'summary_mismatches':len(result['historical_summary_mismatches']),
        'strict_n_unscorable':sum(absent.values()),'strict_payloads_no_targets':len(absent)},indent=2))


if __name__=='__main__':
    main()
