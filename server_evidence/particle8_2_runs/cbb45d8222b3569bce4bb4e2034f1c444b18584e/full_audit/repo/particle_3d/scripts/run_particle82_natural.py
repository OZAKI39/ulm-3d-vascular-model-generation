#!/usr/bin/env python3
"""Stage F: fresh natural flux-weighted MB ledger, only after A–E evidence."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle82_batch import run_batch
from particle_3d.particle82_diagnostics import new_ledger
from particle_3d.particle82_provenance import atomic_json,require_remote


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['output-dir','host-provenance','point-basin-summary','admission-audit','stop-audit','extended-audit','timestep-audit','scaling-benchmark','budget']:p.add_argument('--'+n,required=True)
    a=p.parse_args();read=lambda f:json.loads(Path(f).read_text());host=read(a.host_provenance);require_remote(host,host['hostname'])
    point=read(a.point_basin_summary);admission=read(a.admission_audit);stop=read(a.stop_audit);extended=read(a.extended_audit)
    timestep=read(a.timestep_audit);scaling=read(a.scaling_benchmark);budget=read(a.budget)
    if point['count']<100000 or not point['three_outlets_observed']:raise ValueError('Blocking point basin gate')
    if admission['totals']['scheduled']!=2200 or stop['count']!=720:raise ValueError('B/C accounting gate failed')
    if not scaling['all_exact'] or not all(r['all_original_prefixes_exact'] for r in extended['factors']):raise ValueError('D/parallel parity gate failed')
    if not timestep['rows']:raise ValueError('Missing timestep audit')
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True);ledger_path=out/'birth_ledger.json'
    ledger=read(ledger_path) if ledger_path.exists() else new_ledger(budget['natural_initial_scheduled'])
    if len(ledger['events'])>budget['natural_max_scheduled']:raise ValueError('Predeclared acquisition budget exceeded')
    atomic_json(ledger_path,ledger)
    config=dict(dt_s=.00025,guard_factor=16,dataset_role='NATURAL_FLUX_WEIGHTED_DATASET',
        guard_choice='Predeclared largest completed diagnostic guard factor, preserving original P65 rules; not production convergence',
        scientific_constants=dict(H_D=.45,C_MB_m3=8.5e12,isotropic_V0=True))
    atomic_json(out/'STAGE_GATES.json',dict(point_basin=point,admission_totals=admission['totals'],
        prefix_parity=True,parallel_parity=True,timestep_groups=timestep['groups'],budget=budget,
        flow_flux_sanity='RETAIN_POINT_UNRESOLVED_AND_STATISTICAL_DISCREPANCY_AS_A_SEPARATE_LIMITATION_NO_FROZEN_FLOW_CORRECTION'))
    run_batch(ledger,out,scaling['chosen_workers'],1,len(ledger['events']),True,1,config,host)
