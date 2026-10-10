"""DOLFINx quadrature material adapter and PETSc SNES/VI orchestration.

The nonlinear and bound solvers are PETSc; this module supplies variational
residuals and the analytic Ogden energy derivatives at their quadrature points.
"""
import time
import numpy as np
import basix,basix.ufl,ufl
from dolfinx import fem,mesh
from dolfinx.fem.petsc import NonlinearProblem,assemble_residual,assemble_jacobian,assemble_matrix,assemble_vector
from petsc4py import PETSc
from ..materials_literature.ogden import evaluate,MU,BULK

def scalar(form):return float(fem.assemble_scalar(fem.form(form)))
def vector(form):
    v=assemble_vector(fem.form(form));v.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE);a=v.array.copy();v.destroy();return a

def volumetric_tangent(F):
    J=np.linalg.det(F);G=np.linalg.inv(F).swapaxes(-1,-2)
    return BULK*( (J*(2*J-1))[...,None,None,None,None]*np.einsum('...ij,...kl->...ijkl',G,G) -(J*(J-1))[...,None,None,None,None]*np.einsum('...il,...kj->...ijkl',G,G))

class Coupled:
    def __init__(self,domain,V,u,bcs,fixed_flat,quadrature=2,Gc=.96,ell=.0005,kres=1e-6,external=None,prefix='ogden_',options=None):
        self.domain=domain;self.V=V;self.u=u;self.bcs=bcs;self.fixed_flat=fixed_flat;self.coords=V.tabulate_dof_coordinates();self.nc=domain.topology.index_map(3).size_local;self.entities=np.arange(self.nc,dtype=np.int32);self.q=quadrature;self.Gc=Gc;self.ell=ell;self.kres=kres
        self.D=fem.functionspace(domain,('Lagrange',1));self.d=fem.Function(self.D);self.old=fem.Function(self.D);self.dcoords=self.D.tabulate_dof_coordinates()
        self.constant_F=(V.ufl_element().degree==1)
        self.qpoints,self.qweights=basix.make_quadrature(domain.basix_cell(),1 if self.constant_F else quadrature);self.nq=len(self.qpoints);self.F=ufl.Identity(3)+ufl.grad(u);self.Fexpr=fem.Expression(self.F,self.qpoints)
        self.dx=ufl.Measure('dx',domain=domain,metadata={'quadrature_degree':quadrature,'quadrature_rule':'default'})
        def Q(shape):
            element=basix.ufl.element('DG',domain.basix_cell(),0,shape=shape) if self.constant_F else basix.ufl.quadrature_element(domain.basix_cell(),degree=quadrature,value_shape=shape)
            return fem.Function(fem.functionspace(domain,element))
        self.P=Q((3,3));self.Pv=Q((3,3));self.A=Q((3,3,3,3));self.Av=Q((3,3,3,3));self.W=Q(());self.Wi=Q(());self.Wv=Q(());self.qids=self.W.function_space.dofmap.list;self._cached=None;self.material_calls=0;self.material_seconds=0.
        self.g=(1-kres)*(1-self.d)**2+kres;self.Uel=(self.g*self.Wi+self.Wv)*self.dx;self.Ufrac=Gc*(self.d**2/(2*ell)+ell/2*ufl.dot(ufl.grad(self.d),ufl.grad(self.d)))*self.dx;self.Pi=self.Uel+self.Ufrac
        v=ufl.TestFunction(V);du=ufl.TrialFunction(V);i,j,k,l=ufl.indices(4);self.internal=ufl.inner(self.g*self.P+self.Pv,ufl.grad(v))*self.dx
        self.Ru=self.internal if external is None else self.internal-external;self.Ju=(self.g*self.A[i,j,k,l]+self.Av[i,j,k,l])*ufl.grad(du)[k,l]*ufl.grad(v)[i,j]*self.dx
        self.external=external;self.full_extern=np.zeros(len(u.x.array)) if external is None else vector(external)
        opts=dict(snes_type='newtontr',snes_tr_fallback_type='dogleg',snes_tr_delta0=1e-3,snes_rtol=1e-10,snes_atol=1e-13,snes_max_it=80,ksp_type='gmres',ksp_rtol=1e-10,ksp_atol=1e-20,ksp_max_it=1000,ksp_gmres_restart=100,pc_type='gamg',mg_levels_ksp_type='chebyshev',mg_levels_pc_type='jacobi',snes_error_if_not_converged=False,ksp_error_if_not_converged=False)
        if options:opts.update(options)
        self.mechanics=NonlinearProblem(self.Ru,u,J=self.Ju,bcs=bcs,petsc_options_prefix=prefix,petsc_options=opts)
        self.assembly_seconds={'residual':0.,'jacobian':0.}
        def residual(snes,x,b):
            u.x.array[:]=x.array_r;u.x.scatter_forward()
            lifted=u.x.array.copy()
            for bc in bcs:bc.set(lifted)
            try:self.update(tangent=not np.array_equal(lifted,u.x.array))
            except ValueError:
                self._cached=None;b.set(np.inf);return # PETSc BT reduces nonfinite trial merit; binding has no domain-error setter
            tm=time.perf_counter();assemble_residual(u,self.mechanics.F,self.mechanics.J,bcs,snes,x,b);self.assembly_seconds['residual']+=time.perf_counter()-tm
        def jacobian(snes,x,A,P):
            u.x.array[:]=x.array_r;u.x.scatter_forward();self.update();tm=time.perf_counter();assemble_jacobian(u,self.mechanics.J,self.mechanics.preconditioner,bcs,snes,x,A,P);self.assembly_seconds['jacobian']+=time.perf_counter()-tm
        self.mechanics.solver.setFunction(residual,self.mechanics.b);self.mechanics.solver.setJacobian(jacobian,self.mechanics.A,self.mechanics.P_mat)
        def objective(snes,x):
            u.x.array[:]=x.array_r;u.x.scatter_forward()
            try:self.update(tangent=False)
            except ValueError:self._cached=None;return np.inf
            return scalar(self.Uel)-float(self.full_extern@u.x.array)
        self.mechanics.solver.setObjective(objective)
        # Work-conjugate rigid modes, eliminated constrained components.
        xyz=self.coords/max(np.ptp(self.coords,axis=0));B=np.zeros((len(xyz),3,6));B[:,0,0]=1;B[:,1,1]=1;B[:,2,2]=1;B[:,1,3]=-xyz[:,2];B[:,2,3]=xyz[:,1];B[:,0,4]=xyz[:,2];B[:,2,4]=-xyz[:,0];B[:,0,5]=-xyz[:,1];B[:,1,5]=xyz[:,0];bb=B.reshape(-1,6);bb[fixed_flat]=0;Qm,_=np.linalg.qr(bb);self.rigid=[]
        for i in range(6):v=self.mechanics.A.createVecRight();v.array[:]=Qm[:,i];self.rigid.append(v)
        self.ns=PETSc.NullSpace().create(vectors=self.rigid,comm=domain.comm);self.mechanics.A.setNearNullSpace(self.ns)
        eta=ufl.TestFunction(self.D);self.Rd=ufl.derivative(self.Pi,self.d,eta);self.Jd=ufl.derivative(self.Rd,self.d,ufl.TrialFunction(self.D));self.phase=None
    def update(self,force=False,tangent=True):
        if not force and self._cached is not None and np.array_equal(self._cached,self.u.x.array) and (not tangent or getattr(self,'_tangent_current',False)):return
        start=time.perf_counter();minimum=np.inf;maximum=-np.inf;finite=True
        for lo in range(0,self.nc,512):
            ids=self.entities[lo:lo+512];F=self.Fexpr.eval(self.domain,ids).reshape(-1,3,3);z=evaluate(F,tangent);qi=self.qids[ids].ravel();self.P.x.array.reshape(-1,9)[qi]=z['P_iso'].reshape(-1,9);self.Pv.x.array.reshape(-1,9)[qi]=z['P_vol'].reshape(-1,9)
            if tangent:
                Av=volumetric_tangent(F);self.A.x.array.reshape(-1,81)[qi]=(z['A']-Av).reshape(-1,81);self.Av.x.array.reshape(-1,81)[qi]=Av.reshape(-1,81)
            self.W.x.array[qi]=z['W0'];self.Wi.x.array[qi]=z['W_iso'];self.Wv.x.array[qi]=z['W_vol'];minimum=min(minimum,float(z['J'].min()));maximum=max(maximum,float(z['J'].max()));finite &= all(np.isfinite(v).all() for v in z.values())
        if not finite:raise ValueError('Nonfinite constitutive fields')
        self.minJ=minimum;self.maxJ=maximum;self._cached=self.u.x.array.copy();self._tangent_current=tangent;self.material_calls+=1;self.material_seconds+=time.perf_counter()-start
    def phase_problem(self,seed,options=None):
        self.seed=np.asarray(seed,dtype=int);opts=dict(snes_type='vinewtonrsls',snes_rtol=1e-10,snes_atol=1e-13,snes_max_it=50,ksp_type='cg',ksp_rtol=1e-11,ksp_atol=1e-20,ksp_max_it=1000,pc_type='gamg',snes_error_if_not_converged=False,ksp_error_if_not_converged=False)
        if options:opts.update(options)
        self.phase=NonlinearProblem(self.Rd,self.d,J=self.Jd,petsc_options_prefix='phase_',petsc_options=opts);self.lower=self.phase.x.duplicate();self.upper=self.phase.x.duplicate();self.upper.set(1.)
        self.mass=vector(ufl.TestFunction(self.D)*self.dx);self.phase_scale=self.mass*self.Gc/self.ell
    def phase_solve(self):
        self.lower.array[:]=self.old.x.array;self.lower.array[self.seed]=1.;self.phase.solver.setVariableBounds(self.lower,self.upper);self.phase.solve();return int(self.phase.solver.getConvergedReason())
    def kkt(self):
        gradient=vector(self.Rd);scaled=gradient/self.phase_scale;current=self.d.x.array;projected=current-np.minimum(np.maximum(current-scaled,self.lower.array_r),self.upper.array_r)
        # Projection is for a residual measurement only. Never modify saved d.
        return float(np.max(abs(projected))),gradient
    def mechanical_quality(self,force_floor,force_scale,length_scale,prescribed):
        self.update();res=vector(self.Ru);fixed=self.fixed_flat;free=~fixed;reaction=np.where(fixed,res,0.);current=self.coords+self.u.x.array.reshape(-1,3);external=self.full_extern.reshape(-1,3);R=reaction.reshape(-1,3);q=dict(free_residual_N=float(np.linalg.norm(res[free])),free_residual_relative=float(np.linalg.norm(res[free])/max(np.linalg.norm(self.full_extern[free]),np.linalg.norm(reaction),force_floor)),global_force_relative=float(np.linalg.norm((R+external).sum(0))/force_scale),current_moment_relative=float(np.linalg.norm(np.cross(current,R+external).sum(0))/(force_scale*length_scale)),BC_error_over_L=float(np.max(abs(self.u.x.array[fixed]-prescribed[fixed]))/length_scale),J_min=self.minJ,J_max=self.maxJ,finite=bool(np.isfinite(self.u.x.array).all() and np.isfinite(self.d.x.array).all()),mechanical_reason=int(self.mechanics.solver.getConvergedReason()))
        return q,reaction,res
    def energies(self):
        self.update();return dict(U_el_J=scalar(self.Uel),U_frac_J=scalar(self.Ufrac),U_total_J=scalar(self.Pi),residual_stiffness_energy_J=scalar(self.kres*self.Wi*self.dx),intact_W_iso_J=scalar(self.Wi*self.dx),intact_W_vol_J=scalar(self.Wv*self.dx))
