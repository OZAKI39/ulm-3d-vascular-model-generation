from pathlib import Path
import sys,json
import numpy as np
import pyvista as pv
from scipy.integrate import quad,solve_ivp
import pytest
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'particle_3d/src'));sys.path.insert(0,str(ROOT/'particle_3d/scripts'))
from particle_3d.normal_rbc_encounter import *
from simulate_normal_rbc_encounter import OUT,TUBE_RADIUS,MEAN_SPEED


def test_normal_biconcave_shape_volume_manifold_and_indentation():
    mesh=pv.read(OUT/'data/normal_rbc_si.vtp');geometry=json.loads((OUT/'data/GEOMETRY.json').read_text())
    assert mesh.n_open_edges==0 and geometry['double_concavity']
    assert geometry['maximum_thickness_um']>2*geometry['center_thickness_um']
    # Independent meridional disk integral, not the beta implementation.
    value=quad(lambda theta:4*np.pi*(RBC_RADIUS*np.sin(theta))*surface_points(theta,0)[2]*RBC_RADIUS*np.cos(theta),
               0,np.pi/2,epsabs=1e-30,epsrel=1e-12)[0]
    assert abs(value/RBC_VOLUME-1)<1e-12
    assert abs(mesh.volume/RBC_VOLUME-1)<.001
    assert abs(np.ptp(mesh.points[:,0])-2*RBC_RADIUS)<1e-18


def test_equatorial_exact_gap_against_independent_render_mesh_distance():
    import vtk
    mesh=pv.read(OUT/'data/normal_rbc_si.vtp');locator=vtk.vtkStaticCellLocator();locator.SetDataSet(mesh);locator.BuildLocator()
    for phi in np.linspace(0,2*np.pi,19):
        pos=(RBC_RADIUS+MB_RADIUS+1e-8)*np.array([np.cos(phi),np.sin(phi),0.])
        point=[0.,0.,0.];cell=vtk.mutable(0);sub=vtk.mutable(0);distance=vtk.mutable(0.)
        locator.FindClosestPoint(pos,point,cell,sub,distance)
        gap=np.sqrt(float(distance))-MB_RADIUS
        # A triangulated mesh is slightly inset from the analytic smooth rim.
        assert -1e-15<=gap-1e-8<.002e-6


def test_contact_solution_matches_independent_projected_velocity_ode():
    radius=RBC_RADIUS+MB_RADIUS;x=np.zeros(3);y=radius*np.array([-.8,-.6,0.])
    vx=np.array([.0035,0,0]);vy=np.array([.004,0,0]);w=vx-vy;dt=.0004
    a,b,v,rec=advance_equatorial(x,y,vx,vy,dt)
    def rhs(t,r):
        n=r/np.linalg.norm(r);return w-min(float(w@n),0.)*n
    ref=solve_ivp(rhs,[0,dt],x-y,rtol=2e-12,atol=1e-19)
    assert np.linalg.norm((a-b)-ref.y[:,-1])<1e-16
    np.testing.assert_allclose(a+b,x+y+dt*(vx+vy),rtol=0,atol=1e-19)
    assert rec['contact_duration_s']==dt and not rec['position_projection']
    assert abs((a-b)@(v[0]-v[1]))<1e-20
    with pytest.raises(ValueError,match='equatorial'):
        advance_equatorial(x,y+[0,0,1e-6],vx,vy,dt)


def test_actual_contact_release_and_continuous_clearance():
    data=np.load(OUT/'data/fine.npz');intervals=json.loads((OUT/'data/fine_intervals.json').read_text());x=data['positions_m'];t=data['times_s']
    assert t[0]==0 and abs(t[-1]-.03)<1e-15 and np.all(np.diff(t)>0)
    assert len(intervals)==3000 and np.all(x[:,:,2]==0)
    gaps=np.linalg.norm(x[:,0]-x[:,1],axis=1)-RBC_RADIUS-MB_RADIUS
    assert np.min(gaps)>-1e-16 and gaps[0]>.4e-6 and gaps[-1]>5e-6
    assert x[0,1,0]<x[0,0,0] and x[-1,1,0]>x[-1,0,0]
    assert any(r['contact_duration_s']>0 for r in intervals)
    assert intervals[0]['contact_duration_s']==intervals[-1]['contact_duration_s']==0
    for k,r in enumerate(intervals):
        assert r['wall_lower_bound_m']>r['maximum_motion_m']>0
        if k%50==0:
            for j in [0,1]:np.testing.assert_allclose(r['free_velocities_m_s'][j],poiseuille(x[k,j],TUBE_RADIUS,MEAN_SPEED),rtol=0,atol=1e-15)
    assert min(r['wall_lower_bound_m'] for r in intervals)>.7e-6
    assert abs(sum(r['dt_s'] for r in intervals)-.03)<1e-14


def test_time_step_refinement_and_replay_follow_actual_curved_contact():
    from render_normal_rbc_encounter import Replay
    a=np.load(OUT/'data/coarse.npz');b=np.load(OUT/'data/fine.npz')
    assert np.linalg.norm(a['positions_m']-b['positions_m'][::2],axis=2).max()<.005e-6
    replay=Replay()
    for k in np.linspace(0,len(replay.times)-1,51,dtype=int):
        pos,_,_,_=replay.at(replay.times[k]);np.testing.assert_allclose(pos,replay.positions[k],rtol=0,atol=2e-18)
    for t in np.linspace(0,.03,432):
        pos,_,rec,_=replay.at(t)
        assert np.linalg.norm(pos[0]-pos[1])-RBC_RADIUS-MB_RADIUS>-1e-16
        assert not rec['position_projection']
