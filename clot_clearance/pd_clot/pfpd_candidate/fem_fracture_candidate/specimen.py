"""Native 3D specimen identity and actual displacement/symmetry constraints."""
from pathlib import Path
import numpy as np
import basix.ufl,ufl
from mpi4py import MPI
from dolfinx import mesh,fem
from scipy.spatial import cKDTree
from .native import Coupled
from ..curved_clot_audit.common import read,dump,digest
class Specimen:
 def __init__(self,O,prefix='specimen_'):
  self.O=Path(O);self.C=read(self.O/'CONFIGURATION.json');self.G=read(self.O/'FRACTURE_GEOMETRY_CONFIGURATION.json');p=self.O/'shared/specimen/final_mesh.npz';z=np.load(p);coord=basix.ufl.element('Lagrange','tetrahedron',1,shape=(3,));self.domain=mesh.create_mesh(MPI.COMM_SELF,z['cells'].astype(np.int64),ufl.Mesh(coord),z['X']);self.V=fem.functionspace(self.domain,('Lagrange',self.C['fracture']['degree_u'],(3,)));self.u=fem.Function(self.V);self.coords=self.V.tabulate_dof_coordinates();self.top=np.isclose(self.coords[:,1],.005,atol=1e-12,rtol=0);self.bottom=np.isclose(self.coords[:,1],-.005,atol=1e-12,rtol=0);self.grips=self.top|self.bottom;self.boundary=fem.Function(self.V);fixed=np.zeros(len(self.u.x.array),bool);fixed.reshape(-1,3)[self.grips]=True
  plane=mesh.locate_entities_boundary(self.domain,2,lambda x:np.isclose(x[2],0,atol=1e-12,rtol=0));zdofs=fem.locate_dofs_topological(self.V.sub(2),2,plane);fixed[zdofs]=True;self.fixed=fixed
  bcs=[fem.dirichletbc(self.boundary,np.flatnonzero(self.grips).astype(np.int32)),fem.dirichletbc(0.,zdofs,self.V.sub(2))]
  self.s=Coupled(self.domain,self.V,self.u,bcs,fixed,quadrature=2,prefix=prefix);dx=self.s.dcoords;self.seed=np.flatnonzero((abs(dx[:,1])<1e-12)&(dx[:,0]<=.013+1e-12));assert len(self.seed)>20;self.s.phase_problem(self.seed)
  self.X=self.domain.geometry.x.copy();self.cells=self.domain.geometry.dofmap.copy();self.dc=self.s.D.dofmap.list.copy();self.vertex_u=cKDTree(self.coords).query(self.X)[1];self.vertex_d=cKDTree(dx).query(self.X)[1]
 def prescribe(self,strain):
  self.boundary.x.array[:]=0;self.boundary.x.array.reshape(-1,3)[self.top,1]=strain*.01/2;self.boundary.x.array.reshape(-1,3)[self.bottom,1]=-strain*.01/2;self.u.x.array[self.fixed]=self.boundary.x.array[self.fixed];return self.boundary.x.array.copy()
 def static_checkpoint(self,out):
  p=Path(out)/'native_identity.npz';np.savez_compressed(p,mesh_X=self.X,mesh_cells=self.cells,u_coords=self.coords,d_coords=self.s.dcoords,u_dofmap=self.V.dofmap.list,d_dofmap=self.dc,vertex_u=self.vertex_u,vertex_d=self.vertex_d,seed_dofs=self.seed,fixed_flat=self.fixed,top=self.top,bottom=self.bottom,u_owner=self.V.dofmap.index_map.local_to_global(np.arange(self.V.dofmap.index_map.size_local,dtype=np.int32)),d_owner=self.s.D.dofmap.index_map.local_to_global(np.arange(self.s.D.dofmap.index_map.size_local,dtype=np.int32)))
  dump(Path(out)/'CHECKPOINT_IDENTITY.json',dict(native_identity_sha256=digest(p),mesh_sha256=digest(self.O/'shared/specimen/final_mesh.npz'),configuration_sha256=digest(self.O/'CONFIGURATION.json'),method_id=self.C['fracture']['method_id'],material_id='OGDEN1_BOVINE_LITERATURE_CANDIDATE_V1',u_degree=self.C['fracture']['degree_u'],damage_degree=1,seed_actual_xmax_m=float(self.s.dcoords[self.seed,0].max()),full_thickness_factor=2.,no_display_reconstruction=True))
