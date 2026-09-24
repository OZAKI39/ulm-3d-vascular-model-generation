#!/usr/bin/env python3
"""Separate integration convergence from sampling uncertainty and censoring."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from particle_3d.particle82_diagnostics import basin_rows
from particle_3d.particle82_provenance import atomic_json


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--baseline',required=True);p.add_argument('--refined',required=True);p.add_argument('--output',required=True);a=p.parse_args()
    baseline=Path(a.baseline);refined=Path(a.refined)
    x=basin_rows(baseline);y=basin_rows(refined)
    if len(x)!=len(y):raise ValueError('Point convergence population mismatch')
    for left in sorted((baseline/'shards').glob('*.npz')):
        right=refined/'shards'/left.name
        if not np.array_equal(np.load(left)['positions'],np.load(right)['positions']):raise ValueError('Point convergence changed inlet positions')
    changes=[dict(tracer_id=l['tracer_id'],baseline=l['outlet'],refined=r['outlet'],baseline_reason=l['end_reason'],refined_reason=r['end_reason']) for l,r in zip(x,y) if l['outlet']!=r['outlet']]
    bs=json.loads((baseline/'POINT_TRACER_SUMMARY.json').read_text());rs=json.loads((refined/'POINT_TRACER_SUMMARY.json').read_text())
    bounds={o:dict(observed_lower_fraction=n/rs['count'],unknown_fate_upper_fraction=(n+rs['unresolved'])/rs['count']) for o,n in rs['outlet_counts'].items()}
    atomic_json(a.output,dict(count=len(x),identical_seed_positions=True,outlet_class_changes=len(changes),changes=changes,
        baseline_summary=bs,refined_summary=rs,unresolved_not_assigned_to_any_outlet=True,
        partial_identification_bounds_due_to_unknown_fates=bounds,
        point_tracer_three_basins_observed=rs['three_outlets_observed'],
        sampling_only_flux_match=all(r['sampling_compatibility'] for r in rs['fraction_vs_frozen'].values()),
        flux_sanity_status='PASS' if rs['unresolved']==0 and all(r['sampling_compatibility'] for r in rs['fraction_vs_frozen'].values()) else 'NOT_ESTABLISHED',
        interpretation='No arbitrary percentage threshold. Binomial tests/intervals are explicit; censored point paths invalidate an unqualified all-arrivals outlet-fraction interpretation. No Frozen field correction is applied.'))
