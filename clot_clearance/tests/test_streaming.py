import numpy as np
import pytest
import pyvista as pv
from pd_clot.streaming import FileStreaming,fluid_traction,split_traction,pipe_field,map_forces
from pd_clot.geometry import make_cloud
from test_mechanics import model


def manufactured(points):
    G=np.array([[0,3.,2.],[0,0,-1.],[0,0,0]])
    u=points@G.T+np.array([.01,.02,.03]);p=2+points@np.array([2.,3.,4.])
    return u,p,G


@pytest.mark.parametrize('extension',['csv','vtu','vtk'])
def test_import_velocity_pressure_linear_exact(tmp_path,extension):
    grid=pv.ImageData(dimensions=(4,4,4),spacing=(.001,.001,.001),origin=(-.001,-.001,-.001)).cast_to_unstructured_grid().triangulate()
    u,p,G=manufactured(grid.points);grid['Velocity']=u;grid['Pressure']=p
    path=tmp_path/f'field.{extension}'
    if extension=='csv':np.savetxt(path,np.column_stack((grid.points,u,p)),delimiter=',',header='x,y,z,ux,uy,uz,p',comments='')
    else:grid.save(path)
    provider=FileStreaming(path,.003)
    points=np.array([[.0003,.0006,.0008],[.0001,.0002,.0003]])
    n=np.array([[0.,0.,1.],[1.,0.,0.]])
    expected=fluid_traction(np.tile(G,(len(points),1,1)),manufactured(points)[1],n,.003)
    np.testing.assert_allclose(provider.traction(points,n,0),expected,atol=1e-12)
    with pytest.raises(ValueError,match='outside'):
        provider.traction(np.array([[.1,.1,.1]]),n[:1],0)


@pytest.mark.parametrize('extension',['csv','vtp'])
def test_direct_surface_traction(tmp_path,extension):
    xy=np.array([[x,y,0.] for x in [0.,.001,.002] for y in [0.,.001,.002]])
    traction=np.column_stack((1+2*xy[:,0],3+4*xy[:,1],np.zeros(len(xy))))
    path=tmp_path/f'traction.{extension}'
    if extension=='csv':np.savetxt(path,np.column_stack((xy,traction)),delimiter=',',header='x,y,z,tx,ty,tz',comments='')
    else:
        poly=pv.PolyData(xy);poly['traction_Pa']=traction;poly.save(path)
    provider=FileStreaming(path,.003)
    x=np.array([[.0003,.0007,0.]])
    np.testing.assert_allclose(provider.traction(x,np.array([[0.,0.,1.]]),0),[[1.0006,3.0028,0]],atol=1e-12)


def test_traction_sign_and_geometric_force(model):
    c,g,b=model
    normals=g.face_normal;pressure=np.full(len(normals),12.)
    tr=fluid_traction(np.zeros((len(normals),3,3)),pressure,normals,.003)
    normal,tangent=split_traction(tr,normals)
    np.testing.assert_allclose(np.sum(normal*normals,axis=1),-12)
    np.testing.assert_allclose(tangent,0)
    np.testing.assert_allclose(map_forces(g,tr).sum(axis=0),(tr*g.face_area[:,None]).sum(axis=0),atol=1e-20)


def test_pipe_flux_and_no_slip(model):
    c,_,_=model;p=c['pipe'];R=p['radius_m']
    r=np.linspace(0,R,1001);points=np.column_stack((r*0,r,r*0));u,pressure,grad=pipe_field(points,p)
    Q=np.trapezoid(2*np.pi*r*u[:,0],r)
    assert abs(Q/p['flow_rate_m3_s']-1)<2e-6
    assert abs(u[-1,0])<1e-15
