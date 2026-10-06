"""Intact backends for the existing boundary conversion and solver adapter.

PD2 residuals delegate to the sealed variational candidate. Only its exact
directional derivative is assembled here for the existing TRF/LSMR solver.
No energy Hessian matrix, eigenanalysis or new nonlinear iteration is used.
"""
import time
from dataclasses import dataclass
import numpy as np
from scipy.sparse import csr_matrix, diags
from ..higher_order_rk.variational_energy import bond_state
from ..higher_order_rk.adjoint_accumulation import bond_adjoint, nodal_adjoint
from ..svk_spectral_split import svk_split_energy, green_lagrange
from ...regularization import equilibrium as legacy


class PD2Backend:
    def __init__(self, candidate, family):
        self.candidate = candidate
        self.family = family
        self.material = candidate.material

    def evaluate(self, x): return self.candidate.evaluate(x)
    def energy(self, x): return self.candidate.energy(x)
    def force(self, x): return self.evaluate(x).force
    def node_F(self, x): return self.candidate.nodes(x)
    def bond_Ftilde(self, x):
        c = self.candidate
        return bond_state(c.reference, x, c.nodes(x), slice(None), c.material)[0]
    def psi_plus_diagnostic(self, x):
        return svk_split_energy(green_lagrange(self.bond_Ftilde(x)), self.material.lame, self.material.mu)[0]

    def stress_force(self, P):
        """Frozen-stress channel diagnostic using the accepted force pullback."""
        c = self.candidate; r = c.reference
        G = np.zeros((len(r.X),3,3)); grad = np.zeros_like(r.X)
        for k in range(0,len(r.pairs),c.bond_block):
            sl=slice(k,k+c.bond_block)
            stress = np.broadcast_to(P,(len(r.pairs[sl]),3,3)) if np.ndim(P)==2 else P[sl]
            bond_adjoint(r,stress,sl,G,grad)
        return -nodal_adjoint(c.operator,G,grad,c.directed_block)


class C1Backend:
    """Compatibility wrapper only; historical large campaigns are never rerun."""
    def __init__(self,family,material): self.family,self.material=family,material
    def evaluate(self,x):
        from ..baqp import evaluate
        return evaluate(self.family,x,self.material)
    def energy(self,x):return self.evaluate(x).energy
    def force(self,x):return self.evaluate(x).force
    def node_F(self,x):return self.evaluate(x).Fbar
    def bond_Ftilde(self,x):return self.evaluate(x).Ftilde
    def psi_plus_diagnostic(self,x):
        return svk_split_energy(green_lagrange(self.bond_Ftilde(x)),self.material.lame,self.material.mu)[0]


@dataclass
class Linearization:
    candidate_state: object
    Ftilde: np.ndarray
    S: np.ndarray


class PD2SolverCallback:
    def __init__(self, backend, external, wall_limit=3600):
        self.backend=backend; self.external=np.asarray(external)
        self.x=None; self.state=None; self.calls=0; self.tangent_calls=0
        self.started=time.monotonic();self.wall_limit=wall_limit
        self.assembly='PD2 exact directional derivative; CSR linear-map aggregation'
        c=backend.candidate;r=c.reference;op=c.operator;n=len(r.X);m=len(r.pairs)
        self.D=[]
        for a in range(3):
            mat=csr_matrix((op.g[:,a],(op.source,op.target)),shape=(n,n))
            mat=mat-diags(np.bincount(op.source,weights=op.g[:,a],minlength=n))
            self.D.append(mat.tocsr())
        self.DT=[a.T.tocsr() for a in self.D]
        cols=np.arange(m)
        self.inc=[csr_matrix((np.ones(m),(r.pairs[:,k],cols)),shape=(n,m)) for k in (0,1)]

    def __call__(self,x):
        if time.monotonic()-self.started>self.wall_limit:raise TimeoutError('Predeclared per-case mechanics wall limit')
        b=self.backend;f=b.family
        if self.x is None or not np.array_equal(x,self.x):
            st=b.evaluate(x)
            Ft=bond_state(b.candidate.reference,x,st.Fbar,slice(None),b.material)[0]
            E=.5*(Ft.swapaxes(1,2)@Ft-np.eye(3))
            S=b.material.lame*np.trace(E,axis1=1,axis2=2)[:,None,None]*np.eye(3)+2*b.material.mu*E
            self.state=Linearization(st,Ft,S);self.x=x.copy();self.calls+=1
            print('PD2 residual evaluation',self.calls,'tangent actions',self.tangent_calls,flush=True)
        st=self.state.candidate_state;i,j=f.pairs.T
        strain=np.maximum(np.linalg.norm(x[j]-x[i],axis=1)/f.length-1,0.)
        return legacy.ForceState(st.force,self.external,strain,st.energy,float(np.linalg.det(self.state.Ftilde).min()),np.zeros_like(x),np.zeros(len(x),bool))

    def tangent(self,v,state=None):
        """d[-grad U]/dx applied to v; identical accepted graph, sparse sums."""
        self.tangent_calls+=1;st=self.state if state is None else state
        r=self.backend.candidate.reference;i,j=r.pairs.T
        dF=np.stack([D@v for D in self.D],axis=2)
        av=.5*(dF[i]+dF[j]);xi=r.xi;r2=r.r2
        dFt=av+(v[j]-v[i]-np.einsum('bij,bj->bi',av,xi))[:,:,None]*xi[:,None,:]/r2[:,None,None]
        F=st.Ftilde;mu=self.backend.material.mu;lam=self.backend.material.lame
        dE=.5*(dFt.swapaxes(1,2)@F+F.swapaxes(1,2)@dFt)
        dS=2*mu*dE+lam*np.trace(dE,axis1=1,axis2=2)[:,None,None]*np.eye(3)
        dP=dFt@st.S+F@dS
        pxi=np.einsum('bij,bj->bi',dP,xi)
        A=.5*r.m[:,None,None]*(dP-pxi[:,:,None]*xi[:,None,:]/r2[:,None,None])
        G=((self.inc[0]@A.reshape(-1,9))+(self.inc[1]@A.reshape(-1,9))).reshape(-1,3,3)
        t=r.m[:,None]*pxi/r2[:,None]
        grad=self.inc[1]@t-self.inc[0]@t
        for a in range(3):grad+=self.DT[a]@G[:,:,a]
        return -grad
