from pathlib import Path
import json
import sys
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'particle_3d/src'))
sys.path.insert(0,str(ROOT/'particle_3d/scripts'))
from prepare_rbc_mb_flow_rotation import OUT, seed_scene, sha
from particle_3d.coflow_rotation import PoiseuilleField, integrate_orientations
from particle_3d.rbc_orientation import rotation_matrix, short_axis
from particle_3d.rbc import RBCState
from particle_3d.rbc_integrator import advance_single_rbc
from particle_3d.microbubble import MicrobubbleState
from particle_3d.integrator import advance_single_microbubble


def saved():
    return json.loads((OUT/'data/SCENE.json').read_text()),np.load(OUT/'data/motion.npz')


def rotation_distance(a,b):
    return 4*np.arctan2(np.minimum(np.linalg.norm(a-b,axis=-1),np.linalg.norm(a+b,axis=-1)),
                       np.maximum(np.linalg.norm(a-b,axis=-1),np.linalg.norm(a+b,axis=-1)))


def test_field_gradient_curl_no_slip_and_mean_flow():
    f=PoiseuilleField();p=np.array([3.,1.6,-1.1])*1e-6;h=1e-10
    numeric=np.column_stack([(f.sample(p+np.eye(3)[i]*h).velocity_m_s-
                             f.sample(p-np.eye(3)[i]*h).velocity_m_s)/(2*h) for i in range(3)])
    sample=f.sample(p)
    np.testing.assert_allclose(numeric,sample.velocity_gradient_s_inv,atol=1e-8)
    g=numeric;curl=np.array([g[2,1]-g[1,2],g[0,2]-g[2,0],g[1,0]-g[0,1]])
    np.testing.assert_allclose(curl,sample.vorticity_s_inv,atol=1e-8)
    assert np.trace(g)==0
    np.testing.assert_allclose(f.sample([0,f.radius_m,0]).velocity_m_s,0,atol=1e-18)
    r=np.linspace(0,f.radius_m,2001)
    mean=np.trapezoid(2*f.mean_m_s*(1-r*r/f.radius_m**2)*2*r,r)/f.radius_m**2
    assert abs(mean/f.mean_m_s-1)<1e-6


def test_saved_motion_matches_existing_particle_step_apis():
    scene,d=saved();records,field,g=seed_scene()
    q0=np.array([r['initial_q_wxyz'] for r in records]);p0=np.array([r['initial_m'] for r in records])
    t=5e-6
    q,_,_=integrate_orientations([0,t],q0,9,p0,g.jeffery_lambda,field,t)
    for j in [0,1,2]:
        state=RBCState(j,p0[j],q0[j],np.zeros(3),np.zeros(3),g)
        new=advance_single_rbc(state,field,t)
        np.testing.assert_allclose(new.quaternion_wxyz,q[-1,j],atol=1e-14,rtol=0)
        np.testing.assert_allclose(new.position_m,p0[j]+t*np.array(records[j]['velocity_m_s']),atol=1e-18)
    for j in [9,25]:
        state=MicrobubbleState(j,p0[j],records[j]['radius_m'],np.zeros(3),np.zeros(3))
        new=advance_single_microbubble(state,field,.03)
        np.testing.assert_allclose(new.position_m,d['positions_m'][-1,j],atol=1e-18)
        np.testing.assert_allclose(new.angular_velocity_s_inv,d['angular_velocity_s_inv'][-1,j],atol=1e-12)


def test_rbc_axis_against_independent_matrix_exponential_solution():
    scene,d=saved();g=np.array(scene['records'][0]['gradient_s_inv']);lam=scene['rbc_jeffery_lambda']
    e=.5*(g+g.T);w=.5*(g-g.T);p0=d['rbc_short_axis'][0]
    errors=[]
    for k,t in enumerate(d['time_s']):
        exact=(expm((w+lam*e)*t)@p0.T).T
        exact/=np.linalg.norm(exact,axis=1)[:,None]
        errors.append(np.max(np.linalg.norm(d['rbc_short_axis'][k]-exact,axis=1)))
    assert max(errors)<.003


def test_rbc_full_quaternion_against_independent_high_accuracy_ode():
    scene,d=saved();g=np.array(scene['records'][0]['gradient_s_inv']);e=.5*(g+g.T)
    curl=np.array(scene['records'][0]['vorticity_s_inv']);lam=scene['rbc_jeffery_lambda']
    for j in [0,1,2]:
        def rhs(t,q):
            r=Rotation.from_quat(np.r_[q[1:],q[0]])
            p=r.apply([0,0,1]);omega=.5*curl+lam*np.cross(p,e@p)
            return .5*np.r_[-omega@q[1:],q[0]*omega+np.cross(omega,q[1:])]
        ref=solve_ivp(rhs,[0,.03],d['quaternion_wxyz'][0,j],method='DOP853',
                      rtol=1e-11,atol=1e-13,t_eval=d['time_s'])
        assert ref.success
        assert rotation_distance(d['quaternion_wxyz'][:,j],ref.y.T).max()<.002


def test_mb_body_orientation_matches_analytic_rotation_not_cosmetic_spin():
    scene,d=saved()
    for j in [9,23,36]:
        omega=.5*np.array(scene['records'][j]['vorticity_s_inv'])
        exact=Rotation.from_rotvec(d['time_s'][:,None]*omega).as_matrix()
        np.testing.assert_allclose(rotation_matrix(d['quaternion_wxyz'][:,j]),exact,rtol=0,atol=2e-14)
        np.testing.assert_allclose(d['angular_velocity_s_inv'][:,j],np.tile(omega,(432,1)),atol=1e-12)
        assert np.linalg.norm(omega)>.0


def test_refinement_rigid_shapes_shared_clock_and_nonzero_rotation():
    scene,d=saved();q=d['quaternion_wxyz'];qc=d['quaternion_coarse_wxyz']
    assert q.shape==(432,37,4) and d['positions_m'].shape==(432,37,3)
    np.testing.assert_allclose(np.linalg.norm(q,axis=-1),1,rtol=0,atol=5e-15)
    assert rotation_distance(q,qc).max()<.002
    for j in range(9):
        axis=d['rbc_short_axis'][:,j]
        assert np.max(np.linalg.norm(np.cross(axis,axis[0]),axis=1))>.5
    r=rotation_matrix(q)
    np.testing.assert_allclose(r@np.swapaxes(r,-1,-2),np.broadcast_to(np.eye(3),r.shape),atol=3e-15)
    assert sha(OUT/'data/motion.npz')==scene['array_sha256']
    assert sha(OUT/'data/normal_rbc_si.vtp')==scene['shape_source_sha256']


def test_all_time_pair_and_wall_bounds_for_every_orientation():
    scene,d=saved();records=scene['records'];p=d['positions_m'][0];v=d['velocity_m_s']
    for j,r in enumerate(records):
        assert np.linalg.norm(p[j,1:])+r['radius_m']<scene['tube_radius_m']
        for k,other in enumerate(records[:j]):
            delta=p[j]-p[k];dv=v[j]-v[k]
            t=0 if dv@dv==0 else np.clip(-delta@dv/(dv@dv),0,.03)
            assert np.linalg.norm(delta+t*dv)>r['radius_m']+other['radius_m']
    for positions in d['positions_m']:
        seen={'RBC':0,'MB':0}
        for point,record in zip(positions,records):
            if abs(point[0])<20e-6+record['radius_m']:seen[record['species']]+=1
        assert seen['RBC']>=2 and seen['MB']>=4
