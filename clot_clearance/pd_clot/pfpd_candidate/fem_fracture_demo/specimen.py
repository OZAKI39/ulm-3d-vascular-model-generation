from pathlib import Path
import numpy as np
import basix.ufl,ufl
from mpi4py import MPI
from dolfinx import mesh,fem
from scipy.spatial import cKDTree
from .native import Coupled
from ..curved_clot_audit.common import read,dump,digest

class Specimen:
 def __init__(self,O,prefix='demo_',options=None):
  self.O=Path(O);self.C=read(self.O/'CONFIGURATION.json');self.G=self.C['geometry'];p=self.O/'mesh/mesh.npz';z=np.load(p)
  coord=basix.ufl.element('Lagrange','tetrahedron',1,shape=(3,));self.domain=mesh.create_mesh(MPI.COMM_SELF,z['cells'].astype(np.int64),ufl.Mesh(coord),z['X'])
  self.V=fem.functionspace(self.domain,('Lagrange',2,(3,)));self.u=fem.Function(self.V);self.coords=self.V.tabulate_dof_coordinates();self.boundary=fem.Function(self.V)
  self.top=np.isclose(self.coords[:,1],.005,atol=1e-12,rtol=0);self.ligament=(abs(self.coords[:,1])<1e-12)&(self.coords[:,0]>=.013-1e-12);self.zsym=abs(self.coords[:,2])<1e-12
  fixed=np.zeros((len(self.coords),3),bool);fixed[self.top]=True;fixed[self.ligament,1]=True;fixed[self.zsym,2]=True;self.fixed=fixed.ravel()
  bcs=[fem.dirichletbc(self.boundary,np.flatnonzero(self.top).astype(np.int32))]
  bcs.extend([fem.dirichletbc(0.,(3*np.flatnonzero(mask)+comp).astype(np.int32),self.V.sub(comp)) for mask,comp in [(self.ligament,1),(self.zsym,2)]])
  C=self.C['fracture'];self.s=Coupled(self.domain,self.V,self.u,bcs,self.fixed,quadrature=C['quadrature'],Gc=self.C['material']['Gc_J_m2'],ell=C['ell_m'],kres=C['k_res'],prefix=prefix,options=options)
  dx=self.s.dcoords;self.seed=np.flatnonzero((abs(dx[:,1])<1e-12)&(dx[:,0]<=.013+1e-12));self.s.phase_problem(self.seed)
  self.X=self.domain.geometry.x.copy();self.cells=self.domain.geometry.dofmap.copy();self.dc=self.s.D.dofmap.list.copy();self.vertex_u=cKDTree(self.coords).query(self.X)[1];self.vertex_d=cKDTree(dx).query(self.X)[1]
 def prescribe(self,strain):
  self.boundary.x.array[:]=0;self.boundary.x.array.reshape(-1,3)[self.top,1]=strain*.01/2;self.u.x.array[self.fixed]=self.boundary.x.array[self.fixed];return self.boundary.x.array.copy()
 def checkpoint_identity(self,out):
  p=Path(out)/'native_identity.npz';np.savez_compressed(p,mesh_X=self.X,mesh_cells=self.cells,u_coords=self.coords,d_coords=self.s.dcoords,u_dofmap=self.V.dofmap.list,d_dofmap=self.dc,vertex_u=self.vertex_u,vertex_d=self.vertex_d,seed_dofs=self.seed,fixed_flat=self.fixed,top=self.top)
  dump(Path(out)/'CHECKPOINT_IDENTITY.json',dict(native_identity_sha256=digest(p),mesh_sha256=digest(self.O/'mesh/mesh.npz'),configuration_sha256=digest(self.O/'CONFIGURATION.json'),method_id=self.C['method_id'],material_id='OGDEN1_BOVINE_LITERATURE_CANDIDATE_V1',u_degree=2,damage_degree=1,full_force_factor=2.,full_energy_factor=4.,symmetry=self.G['symmetry']))
