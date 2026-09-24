#!/usr/bin/env python3
"""Verify saved stationary states; no trajectory integration or model change."""
from pathlib import Path
import json,time,hashlib,socket
import numpy as np
from particle_3d.particle81_simulation import environment,dump
from particle_3d.particle_shapes import Sphere
from particle_3d.particle65_motion import assemble_v1
from particle_3d.particle9a_motion import augment_planar_system
from particle_3d.wall_gap import wall_gap
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.nearfield_handoff import handoff_constraints
from particle_3d.resistance_solver import contact_jacobian
from particle_3d.particle9a2_equilibrium_audit import independent_contact_equilibrium
from particle_3d.routing_stationary_audit import downstream_feasibility

ROOT=Path(__file__).resolve().parents[2];D=ROOT/'particle_3d/reports/particle9a2_inlet_sampling/data';O=ROOT/'particle_3d/outputs/particle9a2_2mmps'

def main():
    rows=json.loads((D/'formal500_trajectory_audit.json').read_text());env=environment();policy=NearFieldRegularizationV1();records=[];start=time.time()
    for row in rows:
        if not row['stationary']:continue
        state=row['final_state_audit'];pid=row['particle_id'];x=np.array(state['final_position_m']);radius=state['radius_m'];shapes={pid:Sphere(x,radius)};field=env.field.sample(x)
        base=assemble_v1(shapes,{pid:np.r_[field.velocity_m_s,.5*field.vorticity_s_inv]},env.mu,env.wall,policy=policy)
        system=augment_planar_system(base,shapes,env.mu,env.wall,lambda _:field.velocity_gradient_s_inv,wall_gap)
        cs,_=handoff_constraints(shapes,env.wall,policy);J,_=contact_jacobian(system,shapes,cs)
        R=system.matrix.toarray();b=system.rhs;U=np.array(state['constrained_velocity'])
        cert=independent_contact_equilibrium(R,b,J,U)
        effective_force=b[:3]-R[:3,3:]@np.linalg.solve(R[3:,3:],b[3:])
        drive=downstream_feasibility(np.array([c.normal for c in cs]),effective_force)
        record=dict(particle_id=pid,diameter_m=2*radius,free_FEM_direction_feasible=state['downstream']['feasible_downstream_direction'],
            effective_translation_drive_N=effective_force.tolist(),effective_drive_cone_check=drive,
            original_solver_velocity_budget_m_s=state['solver']['contact_kkt']['velocity_budget_m_s'],certificate=cert,
            verified=bool(cert['recorded_solution_verified'] and cert['independent_translation_zero_within_roundoff']))
        if pid in [371,435]:
            p=D/f'equilibrium_system_{pid}.npz';np.savez_compressed(p,R=R,b=b,J=J,recorded_U=U)
        records.append(record)
    result=dict(role='READ_ONLY_INDEPENDENT_SAVED_STATE_QP_CERTIFICATE',records=records,stationary_count=len(records),
        all_verified=all(r['verified'] for r in records),
        geometrically_FEM_downstream_feasible_ids=[r['particle_id'] for r in records if r['free_FEM_direction_feasible']],
        wall_seconds=time.time()-start,physical_model_or_trajectories_changed=False,
        interpretation='Geometric FEM-direction feasibility and affine-hydrodynamic equilibrium are different tests. KKT + SPD certifies the unique instantaneous model solution, not absence of geometric escape paths.')
    result.update(host=socket.gethostname(),audit_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        independent_solver_source_sha256=hashlib.sha256((ROOT/'particle_3d/src/particle_3d/particle9a2_equilibrium_audit.py').read_bytes()).hexdigest())
    dump(D/'stationary_equilibrium.json',result)
    # Preserve the conservative geometry-only decision. This phase verifies the
    # original mechanical model; it never edits any trajectory or birth metadata.
    archive=D.parent/'reference/initial_stationary_gate';archive.mkdir(parents=True,exist_ok=True)
    for name in ['formal500_gate.json','formal500_trajectory_audit.json']:
        if not (archive/name).exists():(archive/name).write_bytes((D/name).read_bytes())
    indexed={r['particle_id']:r for r in records}
    for row in rows:
        row['geometry_no_FEM_downstream_direction']=bool(row['stationary'] and row['final_state_audit'].get('downstream',{}).get('feasible_downstream_direction') is False)
        if row['stationary']:
            record=indexed[row['particle_id']]
            row['independent_equilibrium_certificate']=record['certificate']
            row['supported_stationary']=record['verified']
            row['stationary_classification']='GEOMETRIC_FEM_DOWNSTREAM_BLOCKED' if row['geometry_no_FEM_downstream_direction'] else 'AFFINE_HYDRODYNAMIC_EQUILIBRIUM_WITH_GEOMETRIC_FEM_DIRECTION'
        row['numerical_failure']=bool(not row['completed'] and not row['supported_stationary'])
    dump(D/'formal500_trajectory_audit.json',rows)
    gate=json.loads((D/'formal500_gate.json').read_text())
    gate.update(supported_stationary_count=sum(r['supported_stationary'] for r in rows),
        geometry_no_FEM_downstream_count=sum(r['geometry_no_FEM_downstream_direction'] for r in rows),
        affine_equilibrium_with_geometric_FEM_direction_ids=result['geometrically_FEM_downstream_feasible_ids'],
        numerical_failure_ids=[r['particle_id'] for r in rows if r['numerical_failure']],
        no_new_systematic_failure=all(not(r['penetration'] or r['handoff_violation'] or r['inlet_escape_or_outside'] or r['numerical_failure'] or r['continuous_failed']) for r in rows),
        audit_semantics_revision='SPD_KKT_MECHANICAL_EQUILIBRIUM_SEPARATE_FROM_GEOMETRIC_FEM_DIRECTION',
        independent_equilibrium_checks=len(records),independent_equilibrium_seconds=result['wall_seconds'],
        original_geometry_only_gate_sha256=hashlib.sha256((archive/'formal500_gate.json').read_bytes()).hexdigest(),
        trajectories_reintegrated_or_changed=False)
    dump(D/'formal500_gate.json',gate)
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))
    for r in records:
        if r['particle_id'] in [371,435]:print(json.dumps(r,indent=2))
if __name__=='__main__':main()
