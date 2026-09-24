#!/usr/bin/env python3
"""Read-only all-step certificates, saved contact states and escalation gate."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import argparse,time,json
from collections import Counter
import numpy as np
from run_particle9a2 import R,D,O,ROLES,read,sha
from particle_3d.particle81_simulation import environment,dump
from particle_3d.particle_shapes import Sphere
from particle_3d.nearfield_handoff import wall_handoff_certificate,partitioned_wall_handoff_certificate
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.stationary_radius_audit import state_audit
from particle_3d.particle8_replay import canonical_hash

def audit_one(spec):
    event,folder=spec;pid=event['particle_id'];path=Path(folder)/'trajectories'/f'mb_{pid:06d}.json'
    m=read(path);a=np.load(path.with_suffix('.npz'))['samples'];env=environment();ro=env.wall.roundoff_m
    assert sha(path.with_suffix('.npz'))==m['samples_sha256']
    assert canonical_hash(event)==m['birth_metadata_sha256']
    assert np.array_equal(a[0,1:4],event['birth_center_m']) and m['radius_m']==event['radius_m']
    assert np.isfinite(a).all() and (np.diff(a[:,0])>0).all()
    lower=NearFieldRegularizationV1().lower_handoff_gap(event['radius_m'])['h_lower_m']
    failed=[];proofs=Counter();minimum_point_gap=float('inf')
    for i,(x,y) in enumerate(zip(a[:-1,1:4],a[1:,1:4])):
        start,end=Sphere(x,event['radius_m']),Sphere(y,event['radius_m'])
        safe,proof=wall_handoff_certificate(start,end,env.wall,lower)
        if not safe and proof['proof']=='HANDOFF_SWEEP_NOT_CERTIFIED':safe,proof=partitioned_wall_handoff_certificate(start,end,env.wall,lower)
        proofs[proof['proof']]+=1
        if not safe:failed.append(dict(segment=i,proof=proof))
    for x in a[:,1:4]:minimum_point_gap=min(minimum_point_gap,env.wall.nearest_center_triangle(x)[1]-event['radius_m'])
    final_state=None;supported=False
    stationary=not m['completed'] and any(k in str(m['failure_detail']) for k in ['EXACT_STATIONARY','ROUNDOFF_SCALE_STAGNATION'])
    if stationary:
        try:
            final_state=state_audit(pid,event['radius_m'],a[-1,1:4])
            budget=final_state['solver']['contact_kkt']['velocity_budget_m_s']
            supported=bool(final_state['contact_rank']==3 and final_state['retained_rows_independent']
                and np.max(np.abs(final_state['constrained_velocity'][:3]))<=budget
                and np.all(np.array(final_state['multipliers'])>=0)
                and final_state['downstream']['feasible_downstream_direction'] is False)
        except Exception as exc:final_state=dict(audit_error=str(exc))
    if m['completed']:
        hit=env.classifier.first_event(a[-2,1:4],a[-1,1:4]);assert hit is not None and hit.role==m['exit_outlet']
    outside='CENTER_OUTSIDE' in str(m['failure_detail'])
    row=dict(particle_id=pid,outlet=m['exit_outlet'] or 'NO_EXIT',completed=m['completed'],stationary=stationary,
        supported_stationary=supported,reason=m['end_reason'],failure_detail=m['failure_detail'],
        diameter_m=event['diameter_m'],minimum_saved_gap_m=float(a[:,14].min()),minimum_saved_g_nf_m=float(a[:,15].min()),
        independent_minimum_wall_gap_m=float(minimum_point_gap),continuous_segments=len(a)-1,
        continuous_proof_counts=dict(proofs),continuous_failed=failed,final_state_audit=final_state,
        penetration=bool(min(minimum_point_gap,a[:,14].min()) < -ro),
        handoff_violation=bool(a[:,15].min() < -ro or minimum_point_gap-lower < -ro),
        inlet_escape_or_outside=outside,numerical_failure=bool(not m['completed'] and not supported),
        sample_count=len(a),birth_metadata_sha256=m['birth_metadata_sha256'],trajectory_sha256=m['samples_sha256'])
    return row

def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['smoke30','formal500']);p.add_argument('--workers',type=int,default=6);a=p.parse_args()
    folder=R/'smoke30' if a.stage=='smoke30' else O
    ledger=read(folder/'admission/birth_ledger.json');events=ledger['events']
    assert len(events)==(30 if a.stage=='smoke30' else 500)
    environment();rows=[];start=time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(audit_one,(e,str(folder))) for e in events]):rows.append(f.result())
    rows.sort(key=lambda r:r['particle_id']);dump(D/(a.stage+'_trajectory_audit.json'),rows)
    gate=dict(count=len(rows),outlet_counts={k:sum(r['outlet']==k for r in rows) for k in ROLES},
        stationary_count=sum(r['stationary'] for r in rows),supported_stationary_count=sum(r['supported_stationary'] for r in rows),
        penetration_count=sum(r['penetration'] for r in rows),handoff_violation_count=sum(r['handoff_violation'] for r in rows),
        inlet_escape_or_outside_count=sum(r['inlet_escape_or_outside'] for r in rows),
        numerical_failure_ids=[r['particle_id'] for r in rows if r['numerical_failure']],
        uncertified_segment_count=sum(len(r['continuous_failed']) for r in rows),
        continuous_segments_rechecked=sum(r['continuous_segments'] for r in rows),workers=a.workers,wall_seconds=time.time()-start,
        no_new_systematic_failure=all(not(r['penetration'] or r['handoff_violation'] or r['inlet_escape_or_outside'] or r['numerical_failure'] or r['continuous_failed']) for r in rows),
        stationary_semantics='CURRENT_RIGID_SPHERE_MODEL_AT_EXISTING_HANDOFF; NOT_PHYSIOLOGICAL_TRAPPING')
    dump(D/(a.stage+'_gate.json'),gate);print(json.dumps(gate,indent=2),flush=True)
if __name__=='__main__':main()
