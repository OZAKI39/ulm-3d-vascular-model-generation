"""Independent NumPy trilinear/RK2 and SciPy reference. Never imports C++ code."""
from reference_flow_sampler import ReferenceField
from scipy.integrate import solve_ivp
import numpy as np

def velocity(field,x):
 u,status=field.sample(x)
 if np.any(status):raise ValueError('INVALID_QUERY '+str(status.tolist()))
 return u[0] if np.asarray(x).ndim==1 else u

def rk2(field,x,dt):
 x=np.asarray(x);k1=velocity(field,x);mid=x+.5*dt*k1;k2=velocity(field,mid);new=x+dt*k2
 return new,velocity(field,new),mid

def dop853(field,x0,times,affine=None):
 if affine is None:fun=lambda t,x:velocity(field,x)
 else:
  b,A=affine;fun=lambda t,x:np.asarray(b)+np.asarray(A)@x
 sol=solve_ivp(fun,[0,float(times[-1])],x0,t_eval=times,method='DOP853',rtol=1e-12,atol=1e-17)
 if not sol.success:raise RuntimeError(sol.message)
 return sol.y.T,{'method':'DOP853','rtol':1e-12,'atol_m':1e-17,'nfev':sol.nfev}

class WallReference:
 def __init__(self,filename):
  import vtk
  self.vtk=vtk;reader=vtk.vtkSTLReader();reader.SetFileName(str(filename));reader.Update();self.mesh=reader.GetOutput();self.distance=vtk.vtkImplicitPolyDataDistance();self.distance.SetInput(self.mesh)
 def gap(self,x,radius=0):return abs(self.distance.EvaluateFunction(x))-radius
 def enclosed(self,xyz):
  from vtk.util.numpy_support import numpy_to_vtk,vtk_to_numpy
  vtk=self.vtk;points=vtk.vtkPoints();points.SetData(numpy_to_vtk(np.asarray(xyz),deep=True));cloud=vtk.vtkPolyData();cloud.SetPoints(points);e=vtk.vtkSelectEnclosedPoints();e.SetInputData(cloud);e.SetSurfaceData(self.mesh);e.SetTolerance(1e-9);e.Update();return vtk_to_numpy(e.GetOutput().GetPointData().GetArray('SelectedPoints'))!=0
 def segment_safe(self,a,b,threshold,depth=0):
  da=self.gap(a);db=self.gap(b);length=float(np.linalg.norm(b-a))
  if min(da,db)<threshold:return False
  if min(da,db)-length/2>=threshold:return True
  mid=(a+b)/2
  if depth>=20:return False # conservative stop if the distance certificate remains unresolved
  return self.segment_safe(a,mid,threshold,depth+1) and self.segment_safe(mid,b,threshold,depth+1)
