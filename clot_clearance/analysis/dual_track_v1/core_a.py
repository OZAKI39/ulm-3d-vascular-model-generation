from analysis.dual_track_v1.common import *
import ctypes,math
from scipy.spatial import cKDTree
PERIOD=.0004

def delta_y(x,b):
 d=x[b[:,1]]-x[b[:,0]];d[:,1]-=PERIOD*np.rint(d[:,1]/PERIOD);return d

def model(tip=-.001):
 dims=np.array([.005,PERIOD,.001]);counts=np.array([100,8,20]);axes=[(np.arange(n)+.5)*DX-d/2 for n,d in zip(counts,dims)];X=np.stack(np.meshgrid(*axes,indexing='ij'),axis=-1).reshape(-1,3);n=len(X)
 xp=X+dims/2;b=cKDTree(xp,boxsize=[.5,PERIOD,.1]).query_pairs(DELTA*(1+1e-12),output_type='ndarray');b=b[np.lexsort((b[:,1],b[:,0]))]
 g=dict(X=X,bonds=b,bond_types=np.ones(len(b),int),clot_count=n,anchors=np.zeros(n,bool),base=np.zeros(n,bool),surface=np.zeros(n,bool),areas=np.zeros(n),normals=np.zeros_like(X),info={})
 e=Energy(g);a,c=X[b[:,0]],X[b[:,1]];cross=(a[:,2]<0)!=(c[:,2]<0);frac=np.divide(-a[:,2],c[:,2]-a[:,2],out=np.zeros(len(b)),where=c[:,2]!=a[:,2]);xx=a[:,0]+frac*(c[:,0]-a[:,0]);alive=np.ascontiguousarray(~(cross&(xx<tip-1e-12)),dtype=np.uint8)
 low=np.isclose(X[:,2],X[:,2].min(),rtol=0,atol=1e-12);high=np.isclose(X[:,2],X[:,2].max(),rtol=0,atol=1e-12);mask=np.zeros_like(X,bool);mask[low|high,2]=True
 return g,e,alive,mask,low,high,xx,cross

class Energy:
 def __init__(self,g):
  self.X=g['X'];self.b=np.ascontiguousarray(g['bonds'],dtype=np.int32);self.V=DX**3;self.r=np.linalg.norm(delta_y(self.X,self.b),axis=1);self.h=np.where(abs(self.r-DELTA)<=DX/2,(DELTA+DX/2-self.r)/DX,1.);self.m=np.bincount(self.b.ravel(),weights=np.repeat(self.r*self.V*self.h,2),minlength=len(self.X))
  self.lib=ctypes.CDLL(str(ROOT/'analysis/dual_track_v1/energy_y.so'));self.fun=self.lib.evaluate_periodic;self.fun.restype=ctypes.c_double;self.fun.argtypes=[ctypes.c_int]*2+[ctypes.c_void_p]*6+[ctypes.c_double]*3+[ctypes.c_void_p,ctypes.c_double]
 def eval(self,x,alive,force=False):
  x=np.ascontiguousarray(x);f=np.zeros_like(x) if force else None;arr=[x,self.b,alive,self.r,self.h,self.m];u=self.fun(len(x),len(self.b),*[a.ctypes.data for a in arr],self.V,K,G,f.ctypes.data if force else None,PERIOD);return (u,f) if force else u
 def theta(self,x,alive):
  extension=np.linalg.norm(delta_y(x,self.b),axis=1)-self.r;extension[abs(extension)<2.2204e-16]=0;return 3*self.V/self.m*np.bincount(self.b.ravel(),weights=np.repeat(self.h*extension*alive,2),minlength=len(x))

def contact(x,g,alive):
 xp=x.copy();xp[:,0]+=.05;xp[:,2]+=.05;xp[:,1]=(xp[:,1]+PERIOD/2)%PERIOD;n=len(x);p=cKDTree(xp,boxsize=[.1,PERIOD,.1]).query_pairs(DX,output_type='ndarray');b=g['bonds'];p=p[~np.isin(p[:,0]*n+p[:,1],b[alive.astype(bool),0]*n+b[alive.astype(bool),1])];d=-delta_y(x,p);r=np.linalg.norm(d,axis=1);gap=DX-r;c=15*(18*K/(math.pi*DELTA**4))*DX**6/DELTA;f=np.zeros_like(x)
 if len(p):
  ff=c*gap[:,None]*d/r[:,None];np.add.at(f,p[:,0],ff);np.add.at(f,p[:,1],-ff)
 return float(.5*c*np.sum(gap**2)),f,p
