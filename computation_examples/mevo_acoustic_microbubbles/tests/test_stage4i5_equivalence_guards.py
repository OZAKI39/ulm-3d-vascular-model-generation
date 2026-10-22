"""Physical-decision regressions, including evolving and reordered contact populations."""
import numpy as np
import pytest

from mevo_mb.bubbles.contact import contact_candidates
from mevo_mb.bubbles.neighbor_cache import ContactNeighborCache
from mevo_mb.coupling.termination import SimulationTerminationPolicy


def test_termination_requires_every_physical_and_metric_condition():
    values=dict(objective="terminal_branch_assignment_only",terminal_outlet_proven=True,
                no_return_proven=True,excluded_pair_influence_proven=True,require_downstream_metrics=False)
    query=dict(time_s=1.,observation_limit_s=6.,total_particles=256,terminally_assigned=256)
    assert SimulationTerminationPolicy(**values).reason(**query)
    for key in ("terminal_outlet_proven","no_return_proven","excluded_pair_influence_proven"):
        assert SimulationTerminationPolicy(**(values|{key:False})).reason(**query) is None
    assert SimulationTerminationPolicy(**(values|dict(require_downstream_metrics=True))).reason(**query) is None
    assert SimulationTerminationPolicy(**values).reason(**(query|dict(terminally_assigned=255))) is None
    assert SimulationTerminationPolicy().reason(**query) is None
    assert SimulationTerminationPolicy().reason(**(query|dict(time_s=6.)))=="observation_limit_reached"


@pytest.mark.parametrize("kind",["sparse","dense","wall_near","high_velocity","variable_radius"])
def test_neighbor_cache_zero_miss(kind):
    rng=np.random.default_rng(20261005)
    r=rng.uniform(1e-6,3e-6,128)
    x=rng.normal(size=(128,3))*(2e-5 if kind=="sparse" else 4e-6)
    if kind=="wall_near":
        x[:,2]=r
    cache=ContactNeighborCache(8e-6,1e-11)
    for step in range(150):
        x+=rng.normal(size=x.shape)*(1e-5 if kind=="high_velocity" else 5e-8)
        if kind=="variable_radius" and step%7==0:
            r*=rng.uniform(.9,1.1,len(r))
        if step==60:
            x,r=x[::-1].copy(),r[::-1].copy()
        if step==90:
            x,r=x[:-3],r[:-3]
        gap=1e-11 if step%2 else 0.
        reference=contact_candidates(x,r,gap,"ckdtree")
        assert np.array_equal(cache.candidates(x,r,gap),reference)
    assert cache.rebuilds>0
    if kind!="high_velocity":
        assert cache.reuses>0


def test_neighbor_boundary_and_gap_growth():
    r=np.full(2,1e-6)
    x=np.array([[0.,0.,0.],[5e-6,0.,0.]])
    cache=ContactNeighborCache(2e-6,0.)
    assert len(cache.candidates(x,r))==0
    for distance in [3.000000000000001e-6,2e-6,np.nextafter(2e-6,0),1e-6]:
        x[1,0]=distance
        assert np.array_equal(cache.candidates(x,r),contact_candidates(x,r,0.,"ckdtree"))
    x[1,0]=9e-6
    assert np.array_equal(cache.candidates(x,r,8e-6),contact_candidates(x,r,8e-6,"ckdtree"))


def test_exact_wall_cache_handles_ragged_mutable_batches_and_mesh_rebind():
    import pyvista as pv
    import yaml
    from mevo_mb.coupling.bifurcation import BifurcationGeometry
    from mevo_mb.coupling.wall_backends import ExactCoordinateWallCache
    from scripts.stage4.geometry_certificate import install
    with open("configs/stage4/bifurcation.yaml") as stream:
        config=yaml.safe_load(stream)["geometry"]
    original,candidate=BifurcationGeometry(config),BifurcationGeometry(config)
    sphere=pv.Sphere(radius=.001).triangulate()
    for g in (original,candidate):
        g.bind_mesh_surface(sphere)
    cache=ExactCoordinateWallCache(candidate,capacity=40)
    candidate.level_and_gradient=cache.level_and_gradient
    install(original,per_particle=True)
    install(candidate,per_particle=True)
    rng=np.random.default_rng(42)
    points=rng.normal(size=(50,3))*.001
    for index in (np.arange(50),np.arange(3,11),np.arange(50),np.arange(20,40)):
        a,b=original.level_and_gradient(points[index]),candidate.level_and_gradient(points[index])
        assert all(np.array_equal(x,y) for x,y in zip(a,b))
        pa,ha=original.project(points[index],1.78e-6)
        pb,hb=candidate.project(points[index],1.78e-6)
        assert np.array_equal(pa,pb) and np.array_equal(ha,hb)
        b[0][:]=999.
        b[1][:]=999.
    for g in (original,candidate):
        g.bind_mesh_surface(pv.Sphere(radius=.002).triangulate())
    assert all(np.array_equal(a,b) for a,b in zip(original.level_and_gradient(points),candidate.level_and_gradient(points)))
    assert cache.hits>0 and len(cache.cache)<=40
