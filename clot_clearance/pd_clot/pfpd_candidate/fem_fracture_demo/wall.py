from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from dolfinx import fem
from .native import Coupled
from ..curved_clot_audit.native_fem import Native
from ..curved_clot_audit.common import read,dump,digest,OLD

class Hemisphere:
 def __init__(self,O):
  self.O=Path(O);self.C=read(self.O/'CONFIGURATION.json');n=Native('fine',degree=2,restore=False);self.R=n.R;self.domain=n.domain;self.V=n.V;self.u=n.u;self.coords=n.coords;self.fixed=np.repeat(n.fixed,3);self.wall=n.fixed
  self.boundary=fem.Function(self.V);bc=fem.dirichletbc(self.boundary,np.flatnonzero(n.fixed).astype(np.int32));C=self.C['hemisphere']
  self.s=Coupled(self.domain,self.V,self.u,[bc],self.fixed,quadrature=4,Gc=.96,ell=C['ell_m'],kres=C['k_res'],prefix='wall_fracture_',options={'ksp_type':'preonly','pc_type':'lu','pc_factor_mat_solver_type':'mumps'})
  self.seed=np.array([],int);self.s.phase_problem(self.seed);self.X=self.domain.geometry.x.copy();self.cells=self.domain.geometry.dofmap.copy();self.dc=self.s.D.dofmap.list.copy();self.vertex_u=cKDTree(self.coords).query(self.X)[1];self.vertex_d=cKDTree(self.s.dcoords).query(self.X)[1]
 def prescribe(self,strain):
  self.boundary.x.array[:]=0;self.boundary.x.array.reshape(-1,3)[self.wall,0]=strain*self.coords[self.wall,0];self.u.x.array[self.fixed]=self.boundary.x.array[self.fixed];return self.boundary.x.array.copy()
 def checkpoint_identity(self,out):
  p=Path(out)/'native_identity.npz';np.savez_compressed(p,mesh_X=self.X,mesh_cells=self.cells,u_coords=self.coords,d_coords=self.s.dcoords,u_dofmap=self.V.dofmap.list,d_dofmap=self.dc,vertex_u=self.vertex_u,vertex_d=self.vertex_d,seed_dofs=self.seed,fixed_flat=self.fixed,wall=self.wall)
  Y=self.X[self.cells];h=maxedge=np.max(np.stack([np.linalg.norm(Y[:,i]-Y[:,j],axis=1) for i in range(4) for j in range(i)],axis=1),axis=1)
  dump(Path(out)/'CHECKPOINT_IDENTITY.json',dict(native_identity_sha256=digest(p),mesh_sha256=digest(OLD/'geometry/fine/mesh.npz'),configuration_sha256=digest(self.O/'CONFIGURATION.json'),method_id=self.C['method_id'],material_id='OGDEN1_BOVINE_LITERATURE_CANDIDATE_V1',u_degree=2,damage_degree=1,R_m=self.R,ell_m=self.s.ell,h_over_ell_quantiles=np.quantile(h/self.s.ell,[0,.5,.95,1]),typical_h_le_ell_over_3=bool(np.median(h)<=self.s.ell/3),reaction_definition='Q_lambda=sum(Rx*X); plotted force=Q_lambda/R, conjugate displacement=(lambda-1)*R; not net resultant, which is zero',bulk_only=True,adhesion_debonding=False))
