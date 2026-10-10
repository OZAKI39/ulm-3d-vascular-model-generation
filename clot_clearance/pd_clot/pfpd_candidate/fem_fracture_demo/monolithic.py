"""PETSc reduced-space VI for the same coupled AT2 energy and accepted bounds.

Optional acceleration of staggered minimization. No clipping or altered history.
"""
import numpy as np,time
import ufl,basix.ufl
from dolfinx import fem
from dolfinx.fem.petsc import NonlinearProblem,assemble_residual,assemble_jacobian
from scipy.spatial import cKDTree
from petsc4py import PETSc

class CoupledVI:
 def __init__(self,s,prefix='coupled_vi_'):
  self.s=s;cell=s.domain.basix_cell();U=basix.ufl.element('Lagrange',cell,2,shape=(3,));D=basix.ufl.element('Lagrange',cell,1)
  W=fem.functionspace(s.domain,basix.ufl.mixed_element([U,D]));self.w=fem.Function(W);VU,umap=W.sub(0).collapse();VD,dmap=W.sub(1).collapse();self.umap=np.array(umap);self.dmap=np.array(dmap)
  dist,self.order_u=cKDTree(s.coords).query(VU.tabulate_dof_coordinates());assert dist.max()<1e-12
  dist,self.order_d=cKDTree(s.dcoords).query(VD.tabulate_dof_coordinates());assert dist.max()<1e-12
  self.bc_value=fem.Function(VU);mf=s.fixed_flat.reshape(-1,3)[self.order_u].ravel();child=np.flatnonzero(mf).astype(np.int32)
  self.bc=fem.dirichletbc(self.bc_value,[self.umap[child].astype(np.int32),child],W.sub(0))
  u,d=ufl.split(self.w);v,eta=ufl.TestFunctions(W);du,dd=ufl.TrialFunctions(W);g=(1-s.kres)*(1-d)**2+s.kres;gp=-2*(1-s.kres)*(1-d);gpp=2*(1-s.kres);i,j,k,l=ufl.indices(4)
  # Positive row scaling gives both equation blocks force units; roots/bounds unchanged.
  scale=1/s.ell
  R=(ufl.inner(g*s.P+s.Pv,ufl.grad(v))+scale*(gp*s.Wi*eta+s.Gc*(d*eta/s.ell+s.ell*ufl.dot(ufl.grad(d),ufl.grad(eta)))))*s.dx
  Jac=((g*s.A[i,j,k,l]+s.Av[i,j,k,l])*ufl.grad(du)[k,l]*ufl.grad(v)[i,j]+gp*dd*ufl.inner(s.P,ufl.grad(v))+scale*(gp*eta*ufl.inner(s.P,ufl.grad(du))+gpp*s.Wi*dd*eta+s.Gc*(dd*eta/s.ell+s.ell*ufl.dot(ufl.grad(dd),ufl.grad(eta)))))*s.dx
  if s.external is not None:raise ValueError('External work adapter required for this monolithic traction solve')
  self.problem=NonlinearProblem(R,self.w,J=Jac,bcs=[self.bc],petsc_options_prefix=prefix,petsc_options=dict(snes_type='vinewtonrsls',snes_linesearch_type='bt',snes_rtol=1e-10,snes_atol=1e-13,snes_max_it=50,ksp_type='preonly',pc_type='lu',pc_factor_mat_solver_type='mumps',snes_error_if_not_converged=False))
  self.lower=self.problem.x.duplicate();self.upper=self.problem.x.duplicate()
  self.assembly_seconds={'residual':0.,'jacobian':0.}
  def update(x,tangent=True):
   self.w.x.array[:]=x.array_r;self.w.x.scatter_forward();s.u.x.array.reshape(-1,3)[self.order_u]=self.w.x.array[self.umap].reshape(-1,3);s.d.x.array[self.order_d]=self.w.x.array[self.dmap];s.update(tangent=tangent)
  self.update=update
  def residual(sn,x,b):
   lifted=x.array_r.copy();self.bc.set(lifted)
   try:update(x,tangent=not np.array_equal(lifted,x.array_r))
   except ValueError:b.set(np.inf);return
   tm=time.perf_counter();assemble_residual(self.w,self.problem.F,self.problem.J,[self.bc],sn,x,b);self.assembly_seconds['residual']+=time.perf_counter()-tm
  def jacobian(sn,x,A,P):
   update(x);tm=time.perf_counter();assemble_jacobian(self.w,self.problem.J,self.problem.preconditioner,[self.bc],sn,x,A,P);self.assembly_seconds['jacobian']+=time.perf_counter()-tm
  self.problem.solver.setFunction(residual,self.problem.b);self.problem.solver.setJacobian(jacobian,self.problem.A,self.problem.P_mat)
 def solve(self,prescribed):
  s=self.s;self.w.x.array[self.umap]=s.u.x.array.reshape(-1,3)[self.order_u].ravel();self.w.x.array[self.dmap]=s.d.x.array[self.order_d];self.bc_value.x.array[:]=prescribed.reshape(-1,3)[self.order_u].ravel()
  self.lower.set(-PETSc.INFINITY);self.upper.set(PETSc.INFINITY);lower=s.old.x.array.copy();lower[s.seed]=1.
  self.lower.array[self.dmap]=lower[self.order_d];self.upper.array[self.dmap]=1.;self.problem.solver.setVariableBounds(self.lower,self.upper)
  self.problem.solve();self.update(self.problem.x);return int(self.problem.solver.getConvergedReason())
