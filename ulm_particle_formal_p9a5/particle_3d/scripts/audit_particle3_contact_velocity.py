#!/usr/bin/env python3
"""Recompute actual contact residuals from immutable saved trajectories."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json,sys
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.particle_shapes import Sphere,Capsule,Ellipsoid,EPS
from particle_3d.rbc_orientation import rotation_matrix
from particle_3d.wall_geometry import WallGeometry
from particle_3d.wall_gap import touching_contacts
from particle_3d.particle3_cases import write_json,write_rows
DATA=PACKAGE/'reports/particle3/data'


def audit_case(path):
    rows=json.loads(path.read_text());result=[]
    if not rows:return path.stem,0,result
    wall=WallGeometry.from_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    for before,after in zip(rows,rows[1:]):
        if before['contact_state']!='TOUCHING':continue
        mode=before['shape_mode'];center=np.array(before['center_m'])
        if mode=='SPHERE_MB':shape=Sphere(center,before['radius_m'])
        elif mode=='CAPILLARY_DEFORMED':
            shape=Capsule(center,before['capsule_axis_world'],before['R_cap_m'],before['L_cap_m'])
            # Preserve the recorded already-normalized vector bit for bit.
            # A second normalization must not perturb the replayed geometry.
            object.__setattr__(shape,'axis_world',np.array(before['capsule_axis_world']))
        else:shape=Ellipsoid(center,[before[k] for k in ['a_m','b_m','c_m']],rotation_matrix(before['q']))
        contacts=touching_contacts(shape,wall)
        assert contacts,'Saved touching state must have an actual finite contact'
        n=np.array([g.normal_inward for _,g in contacts]);omega=np.array(before['free_omega_s_inv'])
        b=np.zeros(len(contacts)) if mode!='FREE_OBLATE' else np.array([-g.normal_inward@np.cross(omega,g.particle_point_m-center) for _,g in contacts])
        velocity=np.array(after['corrected_velocity_m_s'])
        residual=float(max(0.,np.max(b-n@velocity)))
        free_scale=np.linalg.norm(before['free_velocity_m_s'])
        qp_bound=float(1024*EPS*max(free_scale,np.max(abs(b)),np.finfo(float).tiny))
        # The existing solver coalesces directions within 256 eps. Audit ALL
        # original directions, propagating that pre-existing representation
        # error explicitly; do not pretend the QP tolerance alone covers it.
        # |(n-n_kept) dot V| + |b-b_kept| <= ||n-n_kept|| ||V|| + |b-b_kept|.
        # No solver, wall, or merging threshold is changed by this audit.
        kept=[];representation=[]
        for i,normal in enumerate(n):
            match=next((j for j in kept if np.linalg.norm(normal-n[j])<256*EPS and abs(b[i]-b[j])<=256*EPS*max(abs(b[i]),abs(b[j]),free_scale,np.finfo(float).tiny)),None)
            if match is None:kept.append(i);match=i
            representation.append(np.linalg.norm(normal-n[match])*np.linalg.norm(velocity)+abs(b[i]-b[match]))
        representation_bound=float(max(representation));bound=qp_bound+representation_bound
        result.append(dict(case=path.stem,t0_s=before['time_s'],t1_s=after['time_s'],contact_count=len(contacts),
            actual_normal_residual_m_s=residual,declared_velocity_roundoff_m_s=bound,residual_to_bound_ratio=residual/bound,
            qp_velocity_roundoff_m_s=qp_bound,normal_representation_bound_m_s=representation_bound,qp_only_ratio=residual/qp_bound,
            stored_solver_diagnostic_m_s=after['normal_constraint_error_m_s']))
    return path.stem,len(rows)-1,result


def main():
    paths=sorted(DATA.glob('*_real_*_states.json'));assert len(paths)==18
    rows=[];counts={}
    with ProcessPoolExecutor(max_workers=6) as pool:
        for name,count,records in pool.map(audit_case,paths):
            counts[name]=count;rows.extend(records);print(name,len(records),'contact intervals audited',flush=True)
    summary=dict(status='PASS' if all(r['residual_to_bound_ratio']<=1 for r in rows) else 'FAIL',
        processed_case_count=len(paths),accepted_intervals_by_case=counts,contact_intervals=len(rows),
        max_actual_normal_residual_m_s=max(r['actual_normal_residual_m_s'] for r in rows),
        maximum_residual_to_bound_ratio=max(r['residual_to_bound_ratio'] for r in rows),
        maximum_qp_only_ratio=max(r['qp_only_ratio'] for r in rows),intervals_above_qp_only_bound=sum(r['qp_only_ratio']>1 for r in rows),
        qp_roundoff_factor=1024,normal_merge_direction_factor=256,no_wall_or_velocity_tolerance_change=True,
        rule='Original QP budget 1024*eps*velocity_scale PLUS max(||n-n_kept||*||V_applied||+|b-b_kept|) from the original 256*eps normal-coalescing rule; all contributions stored separately',
        no_state_or_velocity_modified=True,source_role='INDEPENDENT_RECOMPUTATION_FROM_ALL_ACTUAL_CONTACT_NORMALS')
    write_rows(DATA/'04_real_contact_residuals.csv',rows);write_json(DATA/'04_real_contact_velocity_audit.json',summary)
    print(json.dumps(summary),flush=True)
    if summary['status']!='PASS':raise SystemExit(1)


if __name__=='__main__':main()
