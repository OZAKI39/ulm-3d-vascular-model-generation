"""Frozen reference-geometry metrics; crossing loss is never itself a crack surface."""
from pathlib import Path
import sys,numpy as np
from scipy.ndimage import label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
sys.path.insert(0,str(Path(__file__).parents[1]/'shear_boundary_notch_v2'))
from audit import regions
DX=50e-6
class Topology:
 def __init__(self,g,mask):
  self.X=g['X'];self.b=g['bonds'];self.mask=mask;self.reg,self.tip,self.grip,self.band=regions(self.X,self.b);self.mid=self.X[self.b].mean(1);self.r=np.linalg.norm(self.X[self.b[:,1]]-self.X[self.b[:,0]],axis=1);self.pp=mask[self.b,2].all(1);self.none=~mask[self.b,2].any(1);a,c=self.X[self.b[:,0]],self.X[self.b[:,1]];self.cross=(a[:,2]<-1e-12)!=(c[:,2]<-1e-12)
  at=a+np.divide(-a[:,2],c[:,2]-a[:,2],out=np.zeros(len(self.b)),where=abs(c[:,2]-a[:,2])>1e-15)[:,None]*(c-a);ij=np.floor((at[:,:2]+[.0001,.0003]+1e-12)/DX).astype(int);ok=self.cross&(ij[:,0]>=0)&(ij[:,0]<8)&(ij[:,1]>=0)&(ij[:,1]<12);self.crossids=np.flatnonzero(ok);self.binids=ij[ok,0]*12+ij[ok,1];self.initial_counts=np.bincount(self.binids,minlength=96).reshape(8,12);self.voxel=np.floor((self.mid+[.0003,.0003,.0005])/DX+1e-8).astype(int)
 def crossing(self,alive):
  remaining=np.bincount(self.binids,weights=alive[self.crossids],minlength=96).reshape(8,12).astype(int);lost=(remaining==0)&(self.initial_counts>0);advance=np.cumprod(lost,axis=0).sum(0)*DX;return remaining,lost,advance
 def cluster(self,alive):
  dead=~alive.astype(bool);grid=np.zeros((13,13,21),bool);grid[tuple(self.voxel[dead].T)]=True;lab,nlab=label(grid);sizes=np.bincount(lab.ravel())[1:];largest=int(sizes.argmax()+1) if len(sizes) else 0;points=np.argwhere(lab==largest) if largest else np.empty((0,3));m=(points+.5)*DX-[.0003,.0003,.0005];meaningful=np.flatnonzero(sizes>=10)+1;kept=dead&np.isin(lab[tuple(self.voxel.T)],meaningful);near=(abs(self.mid[:,0]+.0001)<=DX)&(abs(self.mid[:,2])<=150.5e-6/2)
  pca={}
  if len(m)>=3:
   _,sv,vt=np.linalg.svd(m-m.mean(0),full_matrices=False);pca=dict(normal_reference=vt[-1].tolist(),rms_thickness_m=float(sv[-1]/np.sqrt(len(m))),principal_spreads_m=(sv/np.sqrt(len(m))).tolist(),y_coverage_m=float(np.ptp(m[:,1])),bounds_m=[m.min(0).tolist(),m.max(0).tolist()])
  return dict(occupied_voxels=int(grid.sum()),largest_voxels=int(sizes.max()) if len(sizes) else 0,secondary_voxel_fraction=float((sizes.sum()-sizes.max())/sizes.sum()) if len(sizes) else 0,meaningful_cluster_count=int(len(meaningful)),meaningful_intersects_tip=bool(np.any(kept&near)),largest_pca=pca)
 def components(self,alive):
  bb=self.b[alive.astype(bool)];n=len(self.X);M=coo_matrix((np.ones(2*len(bb)),(np.r_[bb[:,0],bb[:,1]],np.r_[bb[:,1],bb[:,0]])),shape=(n,n));nc,lab=connected_components(M,directed=False);sizes=np.sort(np.bincount(lab))[::-1];return sizes,lab
 def summary(self,alive,components=True):
  dead=~alive.astype(bool);N=int(dead.sum());remain,lost,advance=self.crossing(alive);out=dict(broken=N,region_counts=np.bincount(self.reg[dead],minlength=4).tolist(),target_band_count=int(self.band[dead].sum()),off_plane_count=int((~self.cross[dead]).sum()),both_axial_prescribed_broken=int((dead&self.pp).sum()),fully_free_endpoint_broken=int((dead&self.none).sum()),target_band_fraction=float(self.band[dead].mean()) if N else 0,grip_fraction=float((self.reg[dead]==1).mean()) if N else 0,off_plane_fraction=float((~self.cross[dead]).mean()) if N else 0,front_advance_m=advance.tolist(),central_advancing_strips_2dx=int(np.sum(advance[4:8]>=2*DX-1e-12)),all_advancing_strips_2dx=int(np.sum(advance>=2*DX-1e-12)),crossing_connectivity_area_proxy_m2=float(advance.sum()*DX),surface_status='CROSSING_CONNECTIVITY_PROXY_ONLY')
  if components:sizes,_=self.components(alive);out['component_sizes']=sizes.tolist()
  return out
