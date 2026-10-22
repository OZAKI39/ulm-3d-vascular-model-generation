"""Shared-clock co-visualization of existing independent RBC / MB advection."""
from pathlib import Path
import sys,json,hashlib,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.normal_rbc_encounter import RBC_RADIUS,RBC_VOLUME,MB_RADIUS,poiseuille
from particle_3d.particle3_cases import write_json,write_rows
OUT=ROOT/'particle_3d/reports/rbc_mb_coflow'
RADIUS=7e-6;MEAN=.002;DURATION=.03;FRAMES=432

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def make_data():
    records=[]
    for k,x in enumerate(np.arange(-112,17,16)):
        p=np.array([x,1.6,0.])*1e-6;v=poiseuille(p,RADIUS,MEAN)
        records.append(dict(id=101+k,species='RBC',initial_m=p.tolist(),velocity_m_s=v.tolist(),
            radius_m=RBC_RADIUS,volume_m3=RBC_VOLUME,normal=[0.,0.,1.]))
    for lane,(y,z) in enumerate([(-3.6,-.3),(-.6,2.7)]):
        for k,x in enumerate(np.arange(-104,27,10)):
            p=np.array([x+(4 if lane else 0),y,z])*1e-6;v=poiseuille(p,RADIUS,MEAN)
            records.append(dict(id=201+100*lane+k,species='MB',initial_m=p.tolist(),velocity_m_s=v.tolist(),radius_m=MB_RADIUS))
    t=np.linspace(0,DURATION,FRAMES);initial=np.array([r['initial_m'] for r in records]);velocity=np.array([r['velocity_m_s'] for r in records])
    positions=initial[None,:,:]+t[:,None,None]*velocity[None,:,:]
    np.savez_compressed(OUT/'data/coflow.npz',time_s=t,positions_m=positions,velocity_m_s=velocity)
    source=ROOT/'particle_3d/reports/rbc_mb_normal/data/normal_rbc_si.vtp'
    shutil.copy2(source,OUT/'data/normal_rbc_si.vtp')
    geom=json.loads((source.parent/'GEOMETRY.json').read_text())
    manifest=dict(records=records,counts={'RBC':sum(r['species']=='RBC' for r in records),'MB':sum(r['species']=='MB' for r in records)},
        duration_s=DURATION,frames=FRAMES,tube_radius_m=RADIUS,mean_flow_m_s=MEAN,view_x_bounds_m=[-20e-6,20e-6],
        model='EXISTING_INDEPENDENT_OVERDAMPED_ADVECTION_ON_ONE_SHARED_PHYSICAL_CLOCK',
        solution='EXACT_ADVECTION_IN_X_INVARIANT_POISEUILLE_FIELD',
        RBC_shape='UNCHANGED_NORMAL_BICONCAVE_DISC',RBC_normal=[0,0,1],
        RBC_orientation_note='Jeffery short axis remains z in this symmetry plane; spin about z leaves the axisymmetric shape unchanged.',
        MB_response='EXISTING_NO_OTHER_FORCE_STOKES_EQUILIBRIUM_V_EQUALS_U',
        scene='USER_APPROVED_WIDER_IDEALIZED_VESSEL_NOT_ORIGINAL_FEM',
        pair_contact=False,hydrodynamic_pair_coupling=False,new_coupled_solver_used=False,
        particle_count_role='DISPLAY_SEEDING_NOT_CONCENTRATION_OR_HEMATOCRIT',
        display_window_role='FIXED_LOCAL_OBSERVATION_WINDOW_NOT_PHYSICAL_INLET_OR_OUTLET',
        original_geometry=geom,shape_source_sha256=sha(source),array_sha256=sha(OUT/'data/coflow.npz'))
    write_json(OUT/'data/SCENE.json',manifest)
    rows=[]
    for k,time in enumerate(t):
        for j,r in enumerate(records):
            x=positions[k,j];rows.append(dict(frame=k,physical_time_s=time,particle_id=r['id'],species=r['species'],
                x_m=x[0],y_m=x[1],z_m=x[2],vx_m_s=velocity[j,0],vy_m_s=velocity[j,1],vz_m_s=velocity[j,2]))
    write_rows(OUT/'data/trajectories.csv',rows)
    print(json.dumps(dict(counts=manifest['counts'],duration_s=DURATION,total_positions=len(rows)),indent=2))

if __name__=='__main__':
    (OUT/'data').mkdir(parents=True,exist_ok=True);make_data()
