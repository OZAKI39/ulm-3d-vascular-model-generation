"""Independent NumPy wall/pair/KKT audit; never calls production code.
Small coupled tests enumerate constraint faces. Large MPI tests use independent
6x6 blocks only after proving that no pair lies in the interaction range.
"""
import itertools
import numpy as np
from reference_wall_hydrodynamics_v0 import Lookup
from reference_pair_frozen import assemble,FlowGradientReference
class RuntimeReference:
 def __init__(self,root,contract):
  self.c=contract;self.cfg=contract['config'];self.a=np.array(contract['radii_m']);self.ids=np.array(contract['ids']);self.table=Lookup(root/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5');self.flow=FlowGradientReference(root/'fields'/self.cfg['field'].split('/')[-1]);self.pairs=list(itertools.combinations(range(len(self.a)),2)) if self.cfg['pair_enabled'] else []
 def system(self,x):
  a=self.a;U,O,E,_=self.flow.query(x);R,b,active=assemble(x,a,U,O,E,self.pairs);drag=np.column_stack([np.repeat((6*np.pi*.001*a)[:,None],3,axis=1),np.repeat((8*np.pi*.001*a**3)[:,None],3,axis=1)]).ravel()
  q0=np.column_stack([U,O]).ravel();drive=np.array([self.cfg['drive_'+s] for s in 'xyz']+[self.cfg['torque_'+s] for s in 'xyz']);d=np.tile(drive,(len(a),1))
  if self.cfg['drive_pattern']:d[self.ids%3!=0]=0
  b+=drag*d.ravel()
  if self.cfg['wall_resistance']:
   for i in range(len(a)):
    W=self.table.global_excess(a[i],.001,x[i,2]-a[i],[0,0,1]);s=slice(6*i,6*i+6);R[s,s]+=W;b[s]+=W@q0[s]
  return R,b,active
 def solve(self,x,base=None,dt=0):
  R,b,active=self.system(x);a=self.a;scale=np.column_stack([np.ones((len(a),3)),np.repeat((1/a)[:,None],3,axis=1)]).ravel();A=R*scale[:,None]*scale[None,:]/1e-8;b=b*scale/1e-8
  J=[];bounds=[]
  if base is not None and dt>0:
   if self.cfg['hard_wall']:
    for i in range(len(a)):
     row=np.zeros(len(b));row[6*i+2]=1;J.append(row);bounds.append(-max(base[i,2]-a[i],0)/dt)
   for i,j in self.pairs:
    d=base[i]-base[j];r=np.linalg.norm(d);h=r-a[i]-a[j]
    if h>=.2*a[i]*a[j]/(a[i]+a[j]):continue
    row=np.zeros(len(b));row[6*i:6*i+3]=d/r;row[6*j:6*j+3]=-d/r;J.append(row);bounds.append(-max(h,0)/dt)
  if not active:
   # Exact block separation, not neglect of any active pair.
   yy=[]
   for i in range(len(a)):
    sl=slice(6*i,6*i+6);aa=A[sl,sl];bb=b[sl];y=np.linalg.solve(aa,bb)
    if base is not None and dt>0 and self.cfg['hard_wall']:
     bound=-max(base[i,2]-a[i],0)/dt
     if y[2]<bound:
      n=np.eye(6)[2];z=np.linalg.solve(aa,n);y+=z*((bound-y[2])/z[2])
    yy.extend(y)
   y=np.array(yy)
  else:
   y=np.linalg.solve(A,b)
   if J:
    J=np.array(J);bounds=np.array(bounds);m=len(J)
    if np.min(J@y-bounds)<-1e-14:
     assert m<=12,'Independent enumerator bounded to small coupled cases';found=False
     for count in range(1,m+1):
      for face in itertools.combinations(range(m),count):
       B=J[list(face)];K=np.block([[A,-B.T],[B,np.zeros((count,count))]])
       try:z=np.linalg.solve(K,np.r_[b,bounds[list(face)]])
       except np.linalg.LinAlgError:continue
       if np.min(z[len(b):])>=-1e-12 and np.min(J@z[:len(b)]-bounds)>=-1e-12:y=z[:len(b)];found=True;break
      if found:break
     if not found:raise RuntimeError('No admissible independent wall/pair KKT solution')
  return (y*scale).reshape(-1,6)
