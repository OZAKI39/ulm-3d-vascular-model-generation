import json
import numpy as np
import pyvista as pv
import pytest
from scipy.spatial import Delaunay
from pd_clot.regularization.resolved_flow import ResolvedStreamingFieldProvider

A=np.array([[.2,.3,0],[0,-.1,.4],[.1,0,-.1]])
B=np.array([.001,.002,.003])


def fixture(tmp_path,suffix,planar=False,with_gradient=False):
    xyz=np.indices((2,2,1 if planar else 2)).reshape(3,-1).T*.001
    u=xyz@A.T+B;p=400+xyz@np.array([100.,200.,300.])
    path=tmp_path/('field'+suffix)
    if suffix=='.csv':
        data=np.column_stack((xyz*1000,u*1000,p/1000));header='x,y,z,ux,uy,uz,p'
        if with_gradient:
            data=np.column_stack((data,np.tile(A.ravel(),(len(xyz),1))))
            header+=','+','.join(f'd{v}_d{x}' for v in ['u','v','w'] for x in ['x','y','z'])
        np.savetxt(path,data,delimiter=',',header=header,comments='')
    else:
        if suffix=='.vtp':mesh=pv.PolyData(xyz*1000)
        else:
            cells=Delaunay(xyz).simplices;mesh=pv.UnstructuredGrid(np.column_stack((np.full(len(cells),4),cells)).ravel(),np.full(len(cells),10,np.uint8),xyz*1000)
        mesh['velocity']=u*1000;mesh['pressure']=p/1000
        if with_gradient:mesh['velocity_gradient']=np.tile(A.ravel(),(len(xyz),1))
        mesh.save(path)
    meta=dict(length_unit='mm',velocity_unit='mm/s',pressure_unit='kPa',time_unit='ms',normal_convention='outward_solid',sample_domain='convex_hull')
    q=np.array([[.0002,.0003,0 if planar else .0004]])
    return path,meta,q


@pytest.mark.parametrize('suffix',['.vtu','.vtk','.vtp','.csv'])
def test_independent_affine_velocity_pressure_gradient_and_SI(tmp_path,suffix):
    path,meta,q=fixture(tmp_path,suffix)
    p=ResolvedStreamingFieldProvider(path,meta,.003,q,tmp_path/'report.json')
    s=p.sample(q,2.)
    np.testing.assert_allclose(s.velocity,q@A.T+B,atol=1e-14)
    np.testing.assert_allclose(s.pressure,400+q@np.array([100.,200.,300.]),atol=1e-10)
    np.testing.assert_allclose(s.gradient,A[None],atol=1e-11)
    n=np.array([[1.,0,0]]);expected=(-s.pressure[0]*np.eye(3)+.003*(A+A.T))@n[0]
    np.testing.assert_allclose(p.traction(q,n,0)[0],expected,atol=1e-11)
    assert p.report['surface_coverage_fraction']==1
    with pytest.raises(ValueError,match='outside'):p.sample(np.array([[.01,0,0]]),0)


@pytest.mark.parametrize('suffix',['.vtp','.csv'])
def test_planar_import_requires_full_supplied_gradient(tmp_path,suffix):
    path,meta,q=fixture(tmp_path,suffix,planar=True,with_gradient=False)
    with pytest.raises(ValueError,match='3-D gradient'):ResolvedStreamingFieldProvider(path,meta,.003,q)
    path,meta,q=fixture(tmp_path,suffix,planar=True,with_gradient=True)
    p=ResolvedStreamingFieldProvider(path,meta,.003,q)
    np.testing.assert_allclose(p.sample(q,0).gradient,A[None],atol=1e-12)
    with pytest.raises(ValueError,match='affine plane'):p.sample(q+[0,0,.0001],0)


def test_missing_units_bad_normal_nan_divergence_and_coverage_rejected(tmp_path):
    path,meta,q=fixture(tmp_path,'.csv')
    no_units=dict(meta);no_units.pop('pressure_unit')
    with pytest.raises(ValueError,match='pressure_unit'):ResolvedStreamingFieldProvider(path,no_units,.003,q)
    bad_normal=dict(meta,normal_convention='outward_fluid')
    with pytest.raises(ValueError,match='outward_solid'):ResolvedStreamingFieldProvider(path,bad_normal,.003,q)
    with pytest.raises(ValueError,match='outside'):ResolvedStreamingFieldProvider(path,meta,.003,q+.1,tmp_path/'rejected.json')
    assert json.loads((tmp_path/'rejected.json').read_text())['status']=='REJECTED'
    raw=np.genfromtxt(path,delimiter=',',names=True);raw['ux']+=np.asarray(raw['x'])*10
    np.savetxt(path,np.column_stack([raw[n] for n in raw.dtype.names]),delimiter=',',header=','.join(raw.dtype.names),comments='')
    with pytest.raises(ValueError,match='divergence'):ResolvedStreamingFieldProvider(path,meta,.003,q)
    raw['p'][0]=np.nan
    np.savetxt(path,np.column_stack([raw[n] for n in raw.dtype.names]),delimiter=',',header=','.join(raw.dtype.names),comments='')
    with pytest.raises(ValueError,match='NaN'):ResolvedStreamingFieldProvider(path,meta,.003,q)
