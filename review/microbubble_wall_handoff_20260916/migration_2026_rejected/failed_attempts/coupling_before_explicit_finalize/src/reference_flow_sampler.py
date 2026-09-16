"""Independent NumPy/HDF5 trilinear reference; no C++ calls or C++ results used."""
import h5py,numpy as np
from itertools import product
class ReferenceField:
 def __init__(self,path):
  with h5py.File(path,'r') as h:
   assert h.attrs['format_version']=='FROZEN_FLOW_FIELD_V0';assert h.attrs['coordinate_units']=='m' and h.attrs['velocity_units']=='m/s'
   self.dims=h['dims'][:];self.origin=h['origin_m'][:];self.dx=float(h['dx_m'][()]);self.ids=h['linear_index'][:];self.fluid=h['fluid_linear_index'][:];self.u=h['Velocity_m_s'][:]
 def xyz(self,ids):
  nx,ny,_=self.dims;return self.origin+np.column_stack((ids%nx,(ids//nx)%ny,ids//(nx*ny)))*self.dx
 def sample(self,positions):
  x=np.atleast_2d(np.asarray(positions,dtype=float));q=(x-self.origin)/self.dx;status=np.zeros(len(x),dtype=int);finite=np.all(np.isfinite(x),axis=1);status[~finite]=4
  valid=finite&np.all((q>=0)&(q<self.dims-1),axis=1);status[finite&~valid]=1
  out=np.full(x.shape,np.nan);ii=np.flatnonzero(valid)
  if not len(ii):return out,status
  base=np.floor(q[ii]).astype(np.int64);frac=q[ii]-base;values=np.zeros((len(ii),3));good=np.ones(len(ii),dtype=bool)
  for corner in product((0,1),repeat=3):
   ijk=base+corner;index=ijk[:,0]+self.dims[0]*(ijk[:,1]+self.dims[1]*ijk[:,2]);fs=np.searchsorted(self.fluid,index);isfluid=(fs<len(self.fluid));isfluid[isfluid]&=self.fluid[fs[isfluid]]==index[isfluid]
   vs=np.searchsorted(self.ids,index);exists=vs<len(self.ids);exists[exists]&=self.ids[vs[exists]]==index[exists]
   status[ii[good&~isfluid]]=2;status[ii[good&isfluid&~exists]]=3;good &= isfluid&exists
   w=np.prod(np.where(np.asarray(corner),frac,1-frac),axis=1);ok=isfluid&exists;values[ok]+=self.u[vs[ok]]*w[ok,None]
  out[ii[good]]=values[good];return out,status
 def interior_lower_ids(self,pad=0):
  nx,ny,nz=self.dims;ids=self.fluid.astype(np.int64);ijk=np.column_stack((ids%nx,(ids//nx)%ny,ids//(nx*ny)))
  good=np.all((ijk>=pad)&(ijk<self.dims-1-pad),axis=1)
  for c in product(range(-pad,2),repeat=3):
   cand=ijk+np.array(c);idx=cand[:,0]+nx*(cand[:,1]+ny*cand[:,2]);slots=np.searchsorted(self.ids,idx);ok=slots<len(self.ids);ok[ok]&=self.ids[slots[ok]]==idx[ok];good &= ok
  return ids[good].astype(np.uint64)
