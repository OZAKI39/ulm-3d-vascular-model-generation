"""Independent Python reference. Does not load or call production C++.

Same frozen work-conjugate convention; NumPy/SciPy linear algebra and VTK
geometry serve as the independent audit path, never the production timestep.
"""
from pathlib import Path
import numpy as np
import h5py
from scipy.spatial import cKDTree
import vtk
from vtk.util.numpy_support import vtk_to_numpy

class Lookup:
 def __init__(self,path):
  with h5py.File(path) as f:self.e=f['epsilon'][:];self.R=f['R_total_scaled'][:]
 def excess(self,epsilon):
  e=np.asarray(epsilon);shape=e.shape;e=e.ravel();j=np.searchsorted(self.e,np.clip(e,self.e[0],self.e[-1]),side='right');j=np.clip(j,1,len(self.e)-1);i=j-1
  z=np.log(np.clip(e,self.e[0],self.e[-1]));w=(z-np.log(self.e[i]))/(np.log(self.e[j])-np.log(self.e[i]));r=(1-w[:,None,None])*self.R[i]+w[:,None,None]*self.R[j]-np.eye(6)
  r[e>20]=0;return r.reshape(shape+(6,6))
 def global_excess(self,a,mu,h,n):
  q=frame(n);s=np.sqrt([6*np.pi*mu*a]*3+[8*np.pi*mu*a**3]*3);rot=np.zeros((6,6));rot[:3,:3]=q;rot[3:,3:]=q
  return rot@(self.excess(h/a)*np.outer(s,s))@rot.T
def frame(n,previous=None):
 n=np.asarray(n,dtype=float);n/=np.linalg.norm(n);axis=np.argmin(np.abs(n));t=np.eye(3)[axis]-n*n[axis];t/=np.linalg.norm(t)
 if previous is not None and np.dot(t,previous[:,0])<0:t=-t
 return np.column_stack([t,np.cross(n,t),n])

class Geometry:
 def __init__(self,path):
  r=vtk.vtkSTLReader();r.SetFileName(str(path));r.MergingOff();r.Update();self.poly=r.GetOutput();assert self.poly.GetNumberOfCells()>0
  self.locator=vtk.vtkStaticCellLocator();self.locator.SetDataSet(self.poly);self.locator.BuildLocator()
  self.implicit=vtk.vtkImplicitPolyDataDistance();self.implicit.SetInput(self.poly)
  pts=vtk_to_numpy(self.poly.GetPoints().GetData()).astype(float);ids=vtk_to_numpy(self.poly.GetPolys().GetData()).reshape(-1,4)[:,1:];self.tri=pts[ids]
  ab=self.tri[:,1]-self.tri[:,0];ac=self.tri[:,2]-self.tri[:,0];cr=np.cross(ab,ac);norm=np.linalg.norm(cr,axis=1)
  self.normal=cr/norm[:,None];self.area=norm/2
  self.centers=self.tri.mean(axis=1);self.tree=cKDTree(self.centers);self.bound=np.linalg.norm(self.tri-self.centers[:,None],axis=2).max(axis=1);self.maxbound=self.bound.max()
 def closest(self,x):
  p=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);d2=vtk.reference(0.)
  self.locator.FindClosestPoint(np.asarray(x).tolist(),p,cell,sub,d2)
  return np.asarray(p),int(cell),np.sqrt(float(d2))
 def inside(self,x):return self.implicit.EvaluateFunction(np.asarray(x).tolist())<0
 def patch_ids(self,p,radius):
  candidates=self.tree.query_ball_point(p,2*radius+self.maxbound);out=[]
  # VTK EvaluatePosition independently determines exact triangle distance.
  cp=[0.,0.,0.];sub=vtk.reference(0);param=[0.,0.,0.];d2=vtk.reference(0.);weights=[0.,0.,0.]
  for i in sorted(candidates):
   self.poly.GetCell(i).EvaluatePosition(p,cp,sub,param,d2,weights)
   if float(d2)<=4*radius*radius:out.append(i)
  return np.asarray(out,dtype=int)
 def query(self,x,radius,with_patch=True):
  p,i,d=self.closest(x);n=(np.asarray(x)-p)/d;normal=self.normal[i];angle=np.degrees(np.arccos(np.clip(abs(np.dot(n,normal)),0,1)))
  result={'closest':p,'triangle_id':i,'distance':d,'gap':d-radius,'normal':n,'triangle_normal':normal,'nearest_normal_angle_deg':angle}
  if with_patch:
   ids=self.patch_ids(p,radius);xyz=(self.tri[ids]-p)/radius;w=np.repeat(self.area[ids]/3,3);v=xyz.reshape(-1,3);mean=np.average(v,axis=0,weights=w);cov=(v-mean).T@(w[:,None]*(v-mean))/w.sum();eig,Q=np.linalg.eigh(cov);pn=Q[:,0]
   if pn@n<0:pn=-pn
   angles=np.degrees(np.arccos(np.clip(np.abs(self.normal[ids]@pn),0,1)));order=np.argsort(angles,kind='stable');a=angles[order];weights=self.area[ids][order];p95=a[np.searchsorted(np.cumsum(weights),.95*weights.sum())]
   rms=np.sqrt(max(0,eig[0]));result.update(rms_over_a=rms,normal_spread_deg=p95,patch_triangles=len(ids),plane_normal=pn,planar=bool(rms<=.10 and p95<=20 and angle<=20))
  return result
 def swept_safe(self,x,y,a,tol=1e-12,depth=0):
  dx=self.closest(x)[2];dy=self.closest(y)[2];length=np.linalg.norm(np.asarray(y)-x);minimum=a-tol
  if min(dx,dy)<minimum:return False
  if min(dx,dy)-length/2>=minimum:return True
  if depth>=30:return False
  mid=(np.asarray(x)+y)/2
  return self.swept_safe(x,mid,a,tol,depth+1) and self.swept_safe(mid,y,a,tol,depth+1)
