"""Nominal Pa on a frozen material top surface; dead load, not follower pressure.

The same reference material patch remains loaded after deformation/detachment.
No claim of acoustic, streaming or fluid jet realizability is made.
"""
import numpy as np
MODES=('SURFACE_SHEAR','LOCAL_IMPINGEMENT','IDEAL_LOCAL_PULL')
RISE=.002
END=.0025
def envelope(time):
 return .5*(1-np.cos(np.pi*np.clip(time/RISE,0,1)))
def traction(face_position,face_normal,time,parameters):
 x=np.atleast_2d(face_position);n=np.broadcast_to(face_normal,x.shape);mode=parameters['mode']
 assert mode in MODES
 direction=np.broadcast_to([1.,0,0],x.shape) if mode==MODES[0] else n*(-1 if mode==MODES[1] else 1)
 mask=np.ones(len(x),bool) if mode==MODES[0] else x[:,0]<-25e-6+1e-14
 return parameters['peak_pa']*envelope(time)*direction*mask[:,None]
class Load:
 end=END
 def __init__(self,ref,parameters):
  self.params=parameters;X=ref['X'];top=ref['top_ids']-1
  xs=np.unique(X[top,0]);ys=np.unique(X[top,1]);lookup={(round(X[i,0]/50e-6,6),round(X[i,1]/50e-6,6)):i for i in top};q=[]
  for i in range(len(xs)-1):
   for j in range(len(ys)-1):q.append([lookup[(round(x/50e-6,6),round(y/50e-6,6))] for x,y in [(xs[i],ys[j]),(xs[i+1],ys[j]),(xs[i+1],ys[j+1]),(xs[i],ys[j+1])]])
  self.quads=np.array(q);self.positions=X[self.quads].mean(1);self.normals=np.tile([0.,0.,1.],(len(q),1));self.face_area=np.linalg.norm(np.cross(X[self.quads[:,1]]-X[self.quads[:,0]],X[self.quads[:,3]]-X[self.quads[:,0]]),axis=1);self.N=len(X)
  self.mask=np.ones(len(q),bool) if parameters['mode']==MODES[0] else self.positions[:,0]<-25e-6+1e-14
  self.area=float(self.face_area[self.mask].sum());self.peak_force=self.force(END)
 def force(self,t):
  ft=traction(self.positions,self.normals,t,self.params)*self.face_area[:,None];f=np.zeros((self.N,3))
  for j in range(4):np.add.at(f,self.quads[:,j],ft/4)
  return f
 def audit(self,ref):
  X=ref['X'];ff=traction(self.positions,self.normals,END,self.params)*self.face_area[:,None];f=self.peak_force
  vel=np.random.default_rng(20261016).normal(size=X.shape);vf=vel[self.quads].mean(1)
  return dict(area_m2=self.area,reference_nominal_traction=True,shape_weights=[.25]*4,negative_weights=0,max_face_to_node_fraction=.25,max_nodal_force_n=float(np.linalg.norm(f,axis=1).max()),force_error_n=float(np.linalg.norm(f.sum(0)-ff.sum(0))),moment_error_nm=float(np.linalg.norm(np.cross(X,f).sum(0)-np.cross(self.positions,ff).sum(0))),power_error_w=float(abs(np.sum(f*vel)-np.sum(ff*vf))),direction=np.array([1,0,0]) if self.params['mode']==MODES[0] else np.array([0,0,-1 if self.params['mode']==MODES[1] else 1]),total_peak_force_n=f.sum(0))
