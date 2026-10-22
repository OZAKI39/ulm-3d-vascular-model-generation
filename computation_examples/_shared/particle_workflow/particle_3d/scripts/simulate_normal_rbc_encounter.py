"""Bounded idealized demonstration with a NORMAL undeformed biconcave RBC."""
from pathlib import Path
import sys,json,hashlib,argparse
import numpy as np
import pyvista as pv
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.normal_rbc_encounter import *
from particle_3d.particle3_cases import write_json,write_rows
OUT=ROOT/'particle_3d/reports/rbc_mb_normal'
TUBE_RADIUS=7.e-6;MEAN_SPEED=.002;HORIZON=.030

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def initial():
    a=np.array([0.,2.e-6,0.]);dy=-2.e-6;distance=RBC_RADIUS+MB_RADIUS+.5e-6
    b=a+np.array([-np.sqrt(distance**2-dy**2),dy,0.])
    return a,b

def run(dt,label):
    a,b=initial();positions=[np.array([a,b])];times=[0.];intervals=[];velocities=[];gaps=[];walls=[]
    for k in range(int(round(HORIZON/dt))):
        vx=poiseuille(a,TUBE_RADIUS,MEAN_SPEED);vy=poiseuille(b,TUBE_RADIUS,MEAN_SPEED)
        aa,bb,v,rec=advance_equatorial(a,b,vx,vy,dt)
        # Enclosing balls bound the rigid body's entire tube clearance.
        lower=min(TUBE_RADIUS-np.linalg.norm(a[1:])-RBC_RADIUS,TUBE_RADIUS-np.linalg.norm(b[1:])-MB_RADIUS)
        motion=float(np.sqrt(vx@vx+vy@vy)*dt)
        if lower<=motion:raise ValueError('Continuous tube-clearance bound failed')
        rec.update(t0_s=k*dt,t1_s=(k+1)*dt,wall_lower_bound_m=lower,maximum_motion_m=motion)
        intervals.append(rec);a,b=aa,bb;positions.append(np.array([a,b]));times.append((k+1)*dt)
        velocities.append(v);gaps.append(rec['gap_m']);walls.append(lower)
    x=np.array(positions);t=np.array(times);v=np.array(velocities)
    np.savez_compressed(OUT/f'data/{label}.npz',times_s=t,positions_m=x,velocities_m_s=v)
    write_json(OUT/f'data/{label}_intervals.json',intervals)
    rows=[]
    for k,time in enumerate(t):
        rows.append(dict(time_s=float(time),rbc_x_m=x[k,0,0],rbc_y_m=x[k,0,1],rbc_z_m=x[k,0,2],
            mb_x_m=x[k,1,0],mb_y_m=x[k,1,1],mb_z_m=x[k,1,2],
            gap_m=float(np.linalg.norm(x[k,0]-x[k,1])-RBC_RADIUS-MB_RADIUS),
            contact_duration_s=0. if k==0 else intervals[k-1]['contact_duration_s']))
    write_rows(OUT/f'data/{label}_trajectory.csv',rows)
    contact=[r for r in intervals if r['contact_duration_s']>0]
    result=dict(label=label,dt_s=dt,physical_duration_s=HORIZON,accepted_intervals=len(intervals),
        first_contact_time_s=contact[0]['t0_s']+contact[0]['hit_time_s'],last_contact_time_s=contact[-1]['t0_s']+contact[-1]['contact_duration_s'],
        contact_intervals=len(contact),minimum_pair_gap_m=min(gaps),minimum_wall_bound_m=min(walls),
        final_gap_m=gaps[-1],final_relative_x_m=float(x[-1,1,0]-x[-1,0,0]),
        rbc_diameter_m=2*RBC_RADIUS,rbc_volume_m3=RBC_VOLUME,mb_radius_m=MB_RADIUS,
        tube_radius_m=TUBE_RADIUS,mean_flow_m_s=MEAN_SPEED,field='IDEALIZED_ANALYTICAL_POISEUILLE_NOT_ORIGINAL_FEM',
        shape='NORMAL_BICONCAVE_RBC_FIXED_REST_GEOMETRY',contact_model='P4_KINEMATIC_FRICTIONLESS_EQUATORIAL_RIM',
        deformation=False,membrane_solver=False,lubrication=False,CFD_feedback=False,position_projection=False,
        shape_source=SHAPE_SOURCE,shape_initial_z_scale=volume_scale(),volume_scale_role='INITIAL_REST_SHAPE_VOLUME_CALIBRATION_NOT_DYNAMIC_DEFORMATION')
    write_json(OUT/f'data/{label}_summary.json',result)
    print(json.dumps(result,indent=2),flush=True)

def geometry():
    points,faces=mesh_arrays();m=pv.PolyData(points,np.c_[np.full(len(faces),3),faces].ravel())
    m=m.compute_normals(auto_orient_normals=True,consistent_normals=True)
    assert m.n_open_edges==0
    error=abs(m.volume-RBC_VOLUME)/RBC_VOLUME
    assert error<.001
    m.save(OUT/'data/normal_rbc_si.vtp')
    theta=np.linspace(0,np.pi/2,20001);p=surface_points(theta,np.zeros_like(theta))
    record=dict(diameter_um=2*RBC_RADIUS*1e6,volume_fL=RBC_VOLUME*1e18,mesh_volume_fL=m.volume*1e18,
        mesh_relative_volume_error=error,center_thickness_um=2*p[0,2]*1e6,maximum_thickness_um=2*p[:,2].max()*1e6,
        double_concavity=bool(p[:,2].max()>p[0,2]*2),closed_manifold=m.n_open_edges==0,
        coefficient_source=SHAPE_SOURCE,coefficients=COEFFICIENTS.tolist(),original_diameter_and_volume_preserved=True)
    write_json(OUT/'data/GEOMETRY.json',record);print(record,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dt',type=float,default=1e-5);p.add_argument('--label',default='fine');args=p.parse_args()
    (OUT/'data').mkdir(exist_ok=True,parents=True);geometry();run(args.dt,args.label)
