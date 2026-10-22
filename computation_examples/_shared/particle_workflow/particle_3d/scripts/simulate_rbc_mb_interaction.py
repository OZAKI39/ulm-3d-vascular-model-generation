"""Small, explicitly kinematic RBC/MB encounter in the new 2 mm/s FEM field.

Reuse P3's volume/area-compatible capsule and P4's finite-body contact solver.
This is a selected local encounter, not a physiological RBC suspension.
"""
from pathlib import Path
from dataclasses import asdict
import argparse, hashlib, json, sys, time
import numpy as np
import pyvista as pv

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.field import FrozenFEMField
from particle_3d.wall_geometry import WallGeometry
from particle_3d.wall_gap import wall_gap
from particle_3d.particle_shapes import Sphere,Capsule,unit
from particle_3d.rbc_distribution import sample_rbc_geometries
from particle_3d.rbc import RBCGeometry
from particle_3d.rbc_capillary_surrogate import area_feasible_interval,capsule_length
from particle_3d.particle4_motion import simulate_world,initial_world,world_trial
from particle_3d.physical_time_refinement import refine_interval
from particle_3d.rbc_mb_encounter import advance_encounter,capsule_sphere_gap
from particle_3d.pair_geometry import pair_gap
from particle_3d.particle3_cases import write_json,write_rows

CASE=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps')
OUT=ROOT/'particle_3d/reports/rbc_mb_interaction'
RBC_ID=101;MB_ID=203;MB_RADIUS=.588015722765549e-6

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def setup():
    a=np.load(CASE/'frozen_flow/flow_arrays_si.npz')
    field=FrozenFEMField(a['points_m'],a['tetra'],a['velocity_m_s'],a['pressure_pa'])
    w=pv.read(CASE/'SV_MESH/mesh-surfaces/WALL.vtp')
    triangles=a['points_m'][np.asarray(w['GlobalNodeID'])[w.faces.reshape(-1,4)[:,1:]]-1]
    wall=WallGeometry(triangles,provenance={'role':'NEW_2MMPS_FEM_WALL','path':str(CASE/'SV_MESH/mesh-surfaces/WALL.vtp')})
    def provider(i,s,t):
        sample=field.sample(s.center_m)
        if not sample.inside_lumen:raise ValueError('CENTER_OUTSIDE_NEW_FEM')
        # P4 capsule axis is retained exactly. No nonspherical torque model.
        return np.array(sample.velocity_m_s),np.zeros(3)
    return field,wall,provider

def initial_shapes(record):
    c=record['rbc'];b=record['mb']
    return {RBC_ID:Capsule(c['center_m'],c['axis_world'],c['radius_m'],c['cylindrical_length_m']),
            MB_ID:Sphere(b['center_m'],b['radius_m'])}

def select(field,wall,provider):
    # Retain exactly the median aspect-ratio RBC used by P3/P4, including D/V.
    population=sample_rbc_geometries(41311,2026092002)
    g=RBCGeometry.from_population(population,41310)
    low,high,budget=area_feasible_interval(g)
    centers=field.points_m[field.tetra].mean(1)
    ids=np.random.default_rng(22).choice(len(centers),6000,replace=False)
    distances=np.array([wall.nearest_center_triangle(p)[1] for p in centers[ids]])
    ids=ids[np.argsort(distances)[-20:][::-1]]
    rows=[]
    for ci in ids:
        c=centers[ci];v=field.sample(c).velocity_m_s;axis=unit(v)
        side=unit(np.cross(axis,np.eye(3)[np.argmin(abs(axis))]));other=np.cross(axis,side)
        for radius in [1.3e-6,1.5e-6,1.7e-6]:
            if not low<=radius<=high:continue
            rbc=Capsule(c,axis,radius,capsule_length(g.volume_m3,radius));wg=wall_gap(rbc,wall)
            if wg.gap_m<.12e-6:continue
            for theta in [25,45,65,85,115,145]:
                for phi in range(0,360,45):
                    normal=np.cos(np.deg2rad(theta))*axis+np.sin(np.deg2rad(theta))*(np.cos(np.deg2rad(phi))*side+np.sin(np.deg2rad(phi))*other)
                    mb=Sphere(rbc.support(normal)+(MB_RADIUS+.035e-6)*normal,MB_RADIUS)
                    s=field.sample(mb.center_m)
                    if not s.inside_lumen:continue
                    gap=wall_gap(mb,wall).gap_m
                    closing=-float(normal@(s.velocity_m_s-v))
                    if gap<.08e-6 or closing<5e-6:continue
                    tangential=np.linalg.norm((s.velocity_m_s-v)+closing*normal)
                    rows.append(dict(rbc=asdict(rbc),mb=asdict(mb),tetra_id=int(ci),theta_deg=theta,phi_deg=phi,
                        rbc_wall_gap_m=wg.gap_m,mb_wall_gap_m=gap,initial_pair_gap_m=pair_gap(rbc,mb).gap_m,
                        closing_speed_m_s=closing,tangential_speed_m_s=tangential,
                        score=float(min(gap,wg.gap_m)*closing),rbc_original_geometry=asdict(g),area_budget_m2=budget))
    rows.sort(key=lambda r:r['score'],reverse=True)
    if not rows:raise RuntimeError('No admissible closing encounter found')
    write_json(OUT/'data/initial_candidates.json',rows)
    print('INITIAL_CANDIDATES',len(rows),flush=True)
    write_json(OUT/'data/initialization.json',rows[0])
    return rows[0]

def run(dt,horizon,label):
    begin=time.time();field,wall,provider=setup()
    p=OUT/'data/initialization.json'
    record=json.loads(p.read_text()) if p.exists() else select(field,wall,provider)
    shapes=initial_shapes(record)
    print('RUN',label,'dt',dt,'horizon',horizon,'initial gap um',pair_gap(*shapes.values()).gap_m*1e6,flush=True)
    end=initial_world(shapes,provider,wall);states=[end.to_dict()];ledger=[]
    from dataclasses import replace
    def trial(old,delta):
        a,b=old.shapes[RBC_ID],old.shapes[MB_ID]
        va,vb=provider(RBC_ID,a,old.time_s)[0],provider(MB_ID,b,old.time_s)[0]
        aa,bb,vel,rec=advance_encounter(a,b,va,vb,delta)
        bound=float(np.sqrt(va@va+vb@vb)*delta)
        gaps={i:wall_gap(s,wall) for i,s in [(RBC_ID,aa),(MB_ID,bb)]}
        # The contact projection is nonexpansive in the joint translation norm;
        # signed wall distance is 1-Lipschitz under translation. This proves the
        # WHOLE curved interval, without substituting a straight endpoint chord.
        if min(g.gap_m for g in old.wall_gaps.values())<=bound:
            from particle_3d.physical_time_refinement import TrialNeedsSubdivision
            raise TrialNeedsSubdivision('CURVED_WALL_MOTION_BOUND')
        rec.update(wall_motion_bound_m=bound,wall_clearance_at_start_m=min(g.gap_m for g in old.wall_gaps.values()),
                   canonical_ids=[('PAIR',RBC_ID,MB_ID)] if rec['contact_duration_s']>0 else [])
        new=replace(old,time_s=old.time_s+delta,shapes={RBC_ID:aa,MB_ID:bb},velocities=vel,
            pair_gaps=[capsule_sphere_gap(aa,bb)],wall_gaps=gaps,projection=rec)
        states.append(new.to_dict())
        if len(states)%100==0:print('ACCEPTED',len(states),'time',new.time_s,'gap',new.pair_gaps[0].gap_m,flush=True)
        return new,'ANALYTIC_P4_HEMISPHERE_CONTACT_AND_CONTINUOUS_WALL_CLEARANCE_BOUND'
    for step in range(1,int(np.ceil(horizon/dt))+1):
        if min(g.gap_m for g in end.wall_gaps.values())<.05e-6:
            raise RuntimeError('LOCAL_ENCOUNTER_ENDS_BEFORE_WALL_CONTACT: shorten horizon; no position correction')
        end,_=refine_interval(end,min(step*dt,horizon),trial,ledger=ledger)
    write_json(OUT/f'data/{label}_states.json',states);write_json(OUT/f'data/{label}_ledger.json',ledger)
    rows=[]
    for s in states:
        a,b=s['particles'];projection=s['projection']
        rows.append(dict(time_s=s['time_s'],rbc_x_m=a['center_m'][0],rbc_y_m=a['center_m'][1],rbc_z_m=a['center_m'][2],
            mb_x_m=b['center_m'][0],mb_y_m=b['center_m'][1],mb_z_m=b['center_m'][2],
            pair_gap_m=s['pair_gaps'][0]['gap_m'],rbc_wall_gap_m=a['wall_gap_m'],mb_wall_gap_m=b['wall_gap_m'],
            contact_count=len(projection.get('canonical_ids',[])),
            mb_speed_m_s=float(np.linalg.norm(b['velocity_m_s'])),rbc_speed_m_s=float(np.linalg.norm(a['velocity_m_s']))))
    write_rows(OUT/f'data/{label}_trajectory.csv',rows)
    result=dict(label=label,dt_s=dt,horizon_s=horizon,final_time_s=end.time_s,accepted_intervals=len(ledger),
        minimum_pair_gap_m=min(r['pair_gap_m'] for r in rows),contact_states=sum(r['contact_count']>0 for r in rows),
        minimum_wall_gap_m=min(min(r['rbc_wall_gap_m'],r['mb_wall_gap_m']) for r in rows),
        physical_time_coverage_error_s=abs(sum(r['dt_s'] for r in ledger)-horizon),elapsed_seconds=time.time()-begin,
        source_npz_sha256=sha(CASE/'frozen_flow/flow_arrays_si.npz'),
        source_wall_sha256=sha(CASE/'SV_MESH/mesh-surfaces/WALL.vtp'),
        model='P4_KINEMATIC_HARD_FRICTIONLESS_CONTACT_WITH_FIXED_P3_CAPSULE',
        RBC_membrane_mechanics=False,nonspherical_lubrication=False,two_way_CFD_coupling=False,
        shape_during_encounter='FIXED_VOLUME_AREA_COMPATIBLE_CAPSULE_FIXED_AXIS',
        population_role='SELECTED_LOCAL_TWO_BODY_ENCOUNTER_NOT_CONCENTRATION_OR_HEMATOCRIT')
    write_json(OUT/f'data/{label}_summary.json',result);print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dt',type=float,default=2e-5);p.add_argument('--horizon',type=float,default=.002)
    p.add_argument('--label',default='trial');a=p.parse_args();OUT.mkdir(exist_ok=True);(OUT/'data').mkdir(exist_ok=True)
    run(a.dt,a.horizon,a.label)
