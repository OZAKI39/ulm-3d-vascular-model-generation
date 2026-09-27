"""Prepare small shared-clock trajectories with hydrodynamic orientation."""
from pathlib import Path
import sys, json, hashlib, shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.coflow_rotation import PoiseuilleField, integrate_orientations
from particle_3d.rbc_orientation import quaternion_from_short_axis, short_axis
from particle_3d.rbc_distribution import sample_rbc_geometries
from particle_3d.rbc import RBCGeometry
from particle_3d.normal_rbc_encounter import RBC_RADIUS, RBC_VOLUME, MB_RADIUS
from particle_3d.particle3_cases import write_json, write_rows

OUT=ROOT/'particle_3d/reports/rbc_mb_flow_rotation'
RADIUS=7e-6; MEAN=.002; DURATION=.03; FRAMES=432
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def seed_scene():
    field=PoiseuilleField(RADIUS,MEAN)
    geometry=RBCGeometry.from_population(sample_rbc_geometries(41311,2026092002),41310)
    records=[]
    axes=[[.3,.7,.64],[.8,.4,.45],[-.5,.6,.62]]
    for k,x in enumerate(np.arange(-112,17,16)):
        records.append(dict(id=101+k,species='RBC',initial_m=(np.array([x,1.6,0.])*1e-6).tolist(),
            radius_m=RBC_RADIUS,volume_m3=RBC_VOLUME,
            initial_q_wxyz=quaternion_from_short_axis(axes[k%3]).tolist()))
    # Bounding spheres remain separated for arbitrary RBC orientations.
    for lane,(y,z) in enumerate([(-3.6,-.3),(-2.,3.)]):
        for k,x in enumerate(np.arange(-104,27,10)):
            records.append(dict(id=201+100*lane+k,species='MB',
                initial_m=(np.array([x+4*lane,y,z])*1e-6).tolist(),radius_m=MB_RADIUS,
                initial_q_wxyz=[1.,0.,0.,0.]))
    for r in records:
        s=field.sample(r['initial_m'])
        r.update(velocity_m_s=s.velocity_m_s.tolist(),gradient_s_inv=s.velocity_gradient_s_inv.tolist(),
                 vorticity_s_inv=s.vorticity_s_inv.tolist())
    return records,field,geometry


def export_field(field):
    import pyvista as pv
    grid=pv.ImageData(dimensions=(81,29,29),spacing=(.5e-6,)*3,origin=(-20e-6,-7e-6,-7e-6))
    p=grid.points; radial=np.sum(p[:,1:]**2,axis=1)
    velocity=np.zeros_like(p);velocity[:,0]=2*MEAN*(1-radial/RADIUS**2)
    g=np.zeros((len(p),3,3));g[:,0,1:]=-4*MEAN*p[:,1:]/RADIUS**2
    curl=np.column_stack([np.zeros(len(p)),g[:,0,2],-g[:,0,1]])
    grid['Inside']=np.asarray(radial<=RADIUS**2,dtype=np.uint8)
    grid['Velocity_m_s']=velocity;grid['Gradient_s_inv']=g.reshape(-1,9)
    grid['Vorticity_s_inv']=curl;grid['Speed_mm_s']=np.linalg.norm(velocity,axis=1)*1000
    volume=grid.threshold([.5,1.5],scalars='Inside',all_scalars=True)
    volume.save(OUT/'data/analytical_flow_si.vtu')
    return dict(points=volume.n_points,cells=volume.n_cells,grid_spacing_m=.5e-6,
                gradient_layout='G_ij = du_i/dx_j; row-major',source='ANALYTICAL_POISEUILLE_NOT_FEM')


def main():
    for name in ['data','figures','animations']:(OUT/name).mkdir(parents=True,exist_ok=True)
    records,field,geometry=seed_scene()
    t=np.linspace(0,DURATION,FRAMES)
    p0=np.array([r['initial_m'] for r in records]);v=np.array([r['velocity_m_s'] for r in records])
    q0=np.array([r['initial_q_wxyz'] for r in records]);positions=p0[None]+t[:,None,None]*v[None]
    q,omega,n=integrate_orientations(t,q0,9,p0,geometry.jeffery_lambda,field,5e-6)
    q_coarse,_,_=integrate_orientations(t,q0,9,p0,geometry.jeffery_lambda,field,10e-6)
    np.savez_compressed(OUT/'data/motion.npz',time_s=t,positions_m=positions,velocity_m_s=v,
        quaternion_wxyz=q,angular_velocity_s_inv=omega,quaternion_coarse_wxyz=q_coarse,
        rbc_short_axis=short_axis(q[:,:9]))
    old=ROOT/'particle_3d/reports/rbc_mb_normal/data'
    shutil.copy2(old/'normal_rbc_si.vtp',OUT/'data/normal_rbc_si.vtp')
    geom=json.loads((old/'GEOMETRY.json').read_text())
    scene=dict(records=records,counts={'RBC':9,'MB':28},duration_s=DURATION,frames=FRAMES,
        tube_radius_m=RADIUS,mean_flow_m_s=MEAN,view_x_bounds_m=[-20e-6,20e-6],
        field_information='EXPLICIT_ANALYTICAL_VELOCITY_GRADIENT_AND_VORTICITY_AVAILABLE',
        model='EXISTING_TRANSLATION_AND_JEFFERY_RBC_HALF_CURL_MB_ROTATION',
        scene='USER_APPROVED_WIDER_IDEALIZED_VESSEL',RBC_shape='RIGID_NORMAL_BICONCAVE',
        rbc_rotation_mobility='EXISTING_VOLUME_MATCHED_OBLATE_JEFFERY_APPROXIMATION',
        rbc_jeffery_lambda=float(geometry.jeffery_lambda),rbc_aspect_ratio=float(geometry.r),
        MB_marker='BODY_FIXED_DISPLAY_TEXTURE_NOT_MEMBRANE_DYNAMICS',
        pair_contact=False,hydrodynamic_pair_coupling=False,new_coupled_solver_used=False,
        concentration_role='DISPLAY_SEEDING_ONLY',rotation_substep_limit_s=5e-6,rotation_substeps=n,
        translation='EXACT_IN_X_INVARIANT_FIELD',original_geometry=geom,
        shape_source_sha256=sha(old/'normal_rbc_si.vtp'),array_sha256=sha(OUT/'data/motion.npz'),
        field_export=export_field(field),
        literature_url='https://doi.org/10.1017/jfm.2021.543')
    write_json(OUT/'data/SCENE.json',scene)
    rows=[]
    for k,time in enumerate(t):
        for j,r in enumerate(records):
            row=dict(frame=k,physical_time_s=time,particle_id=r['id'],species=r['species'])
            for name,val in zip(['x_m','y_m','z_m'],positions[k,j]):row[name]=val
            for name,val in zip(['q_w','q_x','q_y','q_z'],q[k,j]):row[name]=val
            for name,val in zip(['omega_x_s_inv','omega_y_s_inv','omega_z_s_inv'],omega[k,j]):row[name]=val
            for name,val in zip(['vx_m_s','vy_m_s','vz_m_s'],v[j]):row[name]=val
            rows.append(row)
    write_rows(OUT/'data/trajectories_and_rotation.csv',rows)
    print(json.dumps(dict(counts=scene['counts'],rbc_lambda=geometry.jeffery_lambda,substeps=n),indent=2))


if __name__=='__main__':main()
