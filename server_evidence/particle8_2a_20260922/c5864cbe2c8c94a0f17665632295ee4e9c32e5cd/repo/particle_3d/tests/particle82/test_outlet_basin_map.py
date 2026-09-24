from types import SimpleNamespace
import numpy as np
import pyvista as pv
import pytest
from particle_3d.field import FrozenFEMField
from particle_3d.validation_boundary import ValidationBoundaryClassifier
from particle_3d.particle82_point_native import NativePointTracer


def fixture_field(affine=False):
    points=np.array([[0,0,0],[10e-6,0,0],[0,10e-6,0],[0,0,10e-6]],float)
    velocity=np.tile([1e-3,0,0],(4,1))
    if affine:velocity[:,0]=100*points[:,0]
    field=FrozenFEMField(points,np.array([[0,1,2,3]]),velocity,np.zeros(4))
    cap=pv.PolyData(points[[1,2,3]],np.array([3,0,1,2]));cap.point_data['GlobalNodeID']=np.array([2,3,4])
    boundaries={'OUTLET_01':cap}
    return SimpleNamespace(field=field,boundaries=boundaries,classifier=ValidationBoundaryClassifier(boundaries))


def test_native_sampler_preserves_original_affine_p1():
    env=fixture_field(True);solver=NativePointTracer(env)
    for p in [[1e-6,1e-6,1e-6],[3e-6,2e-6,1e-6]]:
        point=np.array(p);v=np.zeros(3)
        assert solver.lib.p82_sample(point.ctypes.data,0,v.ctypes.data)==0
        np.testing.assert_allclose(v,env.field.sample(point).velocity_m_s,atol=1e-18,rtol=0)


@pytest.mark.parametrize('affine',[False,True])
def test_diagnostic_advection_matches_analytic_exit(affine):
    # Two tiny synthetic unit cases; never part of any natural MB catalog.
    solver=NativePointTracer(fixture_field(affine));r=solver.trace(np.array([1e-6]*3),error=1e-13)
    assert r['outlet']=='OUTLET_01'
    np.testing.assert_allclose(r['path'][-1,1:],[8e-6,1e-6,1e-6],rtol=0,atol=1e-18)
    exact=np.log(8)/100 if affine else .007
    assert abs(r['path'][-1,0]-exact)<2e-8
    assert np.all(np.diff(r['path'][:,0])>0)
