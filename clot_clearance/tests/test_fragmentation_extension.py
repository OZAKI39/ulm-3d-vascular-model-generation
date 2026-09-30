import json
from pathlib import Path
import numpy as np
from pd_clot.geometry import make_cloud
from pd_clot.mechanics import evaluate_cloud
from pd_clot.fragment_mechanics import prepare_supported_shape,evaluate_supported
from pd_clot.fragment_topology import BondDamage,FragmentTracker,surface_particles
from pd_clot.fragment_fluid import FragmentFluid,relaxation_half_step


def fixture():
    c=json.loads((Path(__file__).resolve().parents[1]/'configs/fragmentation/pilot_001.json').read_text())
    c['clot']['cells']=[4,4,4]
    return c,make_cloud(c['clot'])


def test_full_rank_matches_existing_NOSB():
    c,g=fixture();b=np.ones(len(g.pairs));rng=np.random.default_rng(801)
    x=g.X+rng.normal(size=g.X.shape)*1e-7
    shape=prepare_supported_shape(g,b,c['safety'])
    assert np.all(shape[3]==3)
    first=evaluate_cloud(g,x,b,c['material'],c['safety']);second=evaluate_supported(g,x,b,c['material'],c['safety'],shape)
    for a,z in zip(first,second):np.testing.assert_allclose(a,z,atol=2e-10,rtol=1e-8)


def test_rank_reduced_energy_gradient_and_objectivity():
    c,g=fixture();p=g.pairs
    # Manufactured reduced-rank unit test only. This is never used by the demo runner.
    b=(np.abs(g.X[p[:,0],2]-g.X[p[:,1],2])<1e-12).astype(float)
    shape=prepare_supported_shape(g,b,c['safety']);assert np.all(shape[3]==2)
    rng=np.random.default_rng(1);x=g.X+rng.normal(size=g.X.shape)*1e-7
    result=evaluate_supported(g,x,b,c['material'],c['safety'],shape);f=result[0]
    d=rng.normal(size=x.shape);d/=np.linalg.norm(d);epsilon=1e-9
    def E(y):
        z=evaluate_supported(g,y,b,c['material'],c['safety'],shape)
        return np.dot(g.volume,z[3]+z[4])
    np.testing.assert_allclose((E(x+epsilon*d)-E(x-epsilon*d))/(2*epsilon),-np.sum(f*d),rtol=2e-5,atol=1e-11)
    angle=.45;R=np.array([[np.cos(angle),0,np.sin(angle)],[0,1,0],[-np.sin(angle),0,np.cos(angle)]])
    rotated=evaluate_supported(g,x@R.T+.001,b,c['material'],c['safety'],shape)
    np.testing.assert_allclose(rotated[0],f@R.T,atol=1e-14,rtol=1e-7)
    assert np.linalg.norm(f.sum(axis=0))<1e-15
    assert np.linalg.norm(np.cross(x,f).sum(axis=0))<1e-17


def test_cyclic_failure_threshold_and_no_healing():
    c,g=fixture();c['damage'].update(Q_ref=.01,C_damage=.0001,DeltaN=1000,D_break=.2)
    d=BondDamage(3,c['damage']);d.advance(np.array([0.,.02,.03]))
    assert d.active.tolist()==[True,True,False]
    np.testing.assert_allclose(d.g,[1,.9,0])
    old=d.D.copy();d.advance(np.zeros(3));np.testing.assert_array_equal(d.D,old)
    assert d.g[-1]==0


def test_stable_lineage_clearance_no_double_count():
    c,g=fixture();b=np.ones(len(g.pairs));tracker=FragmentTracker(g,.002)
    x=g.X.copy();v=np.zeros_like(x)
    attached,_,_=tracker.update(b,x,v,0);assert attached.all()
    upper=g.X[:,2]>g.X[:,2].mean()
    b[upper[g.pairs[:,0]]!=upper[g.pairs[:,1]]]=0 # graph-only unit test
    attached,components,stats=tracker.update(b,x,v,1000)
    assert stats['detached_components']==1
    free_id=np.unique(tracker.labels[~attached]);assert len(free_id)==1
    first_labels=tracker.labels.copy();x[upper,0]+=.004
    attached,components,stats=tracker.update(b,x,v,2000)
    np.testing.assert_array_equal(tracker.labels,first_labels)
    assert stats['cleared_volume_fraction']==.5
    assert len(tracker.crossings)==1
    tracker.update(b,x,v,3000);assert len(tracker.crossings)==1
    assert np.isclose(stats['attached_volume_fraction']+stats['detached_volume_fraction'],1)


def test_dynamic_surface_excludes_anchors():
    c,g=fixture();b=np.ones(len(g.pairs));surface,ratio=surface_particles(g,b,c['surface'])
    assert not surface[g.fixed].any()
    damaged=b*.4;new,newratio=surface_particles(g,damaged,c['surface'])
    assert np.all(newratio<=ratio);assert new.sum()>=surface.sum()


def test_hydrodynamic_relaxation_exact_and_fixed_particles():
    v=np.zeros((8,3));u=np.tile([.1,0,0],(8,1));mobile=np.arange(8)>1
    for _ in range(10):v=relaxation_half_step(v,u,mobile,1e-5,1e-4)
    np.testing.assert_allclose(v[mobile,0],.1*(1-np.exp(-1)),rtol=1e-14)
    assert np.all(v[~mobile]==0)


def test_prescribed_vortex_divergence():
    c,g=fixture();c['transport']['background_velocity_multiplier']=0;provider=FragmentFluid(c)
    points=g.X.copy();h=1e-8;div=np.zeros(len(points))
    for k in range(3):
        d=np.eye(3)[k]*h
        div+=(provider.velocity(points+d,0)[:,k]-provider.velocity(points-d,0)[:,k])/(2*h)
    assert np.max(np.abs(div))<1e-6


def test_line_fragment_keeps_force_and_isolated_points_have_no_internal_force():
    c,g=fixture();b=np.zeros(len(g.pairs));b[0]=.8
    shape=prepare_supported_shape(g,b,c['safety'])
    i,j=g.pairs[0];assert shape[3][i]==shape[3][j]==1
    x=g.X.copy();x[j]+=.02*(g.X[j]-g.X[i])
    f,*_=evaluate_supported(g,x,b,c['material'],c['safety'],shape)
    assert np.linalg.norm(f[i])>0
    np.testing.assert_allclose(f[i],-f[j],atol=1e-18)
    assert np.all(f[shape[3]==0]==0)
