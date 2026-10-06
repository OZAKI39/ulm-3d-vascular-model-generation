"""Parallel intact mechanics candidate; never imported by production BAQP."""
from dataclasses import dataclass
import numpy as np
from ..energy import Material
from .rk_kinematics import deformation
from .variational_energy import bond_prefactors,energy_from_nodes
from .variational_force import energy_force


def frozen(a):
    a=np.asarray(a);a.setflags(write=False);return a


@dataclass(frozen=True)
class Reference:
    X: np.ndarray
    volume: np.ndarray
    pairs: np.ndarray
    xi: np.ndarray
    r2: np.ndarray
    m: np.ndarray

    @classmethod
    def from_family(cls,f):
        return cls(*(frozen(a) for a in [f.X,f.volume,f.pairs,f.xi,
            np.sum(f.xi*f.xi,axis=1),bond_prefactors(f)]))


@dataclass(frozen=True)
class GradientOperator:
    name: str
    source: np.ndarray
    target: np.ndarray
    g: np.ndarray

    def __post_init__(self):
        if self.name not in ('C1','PD2'):raise ValueError('Only frozen C1 / PD2 operators')
        if self.g.shape!=(len(self.source),3) or self.target.shape!=self.source.shape or not np.isfinite(self.g).all():
            raise ValueError('Invalid physical gradient weights')
        for a in [self.source,self.target,self.g]:a.setflags(write=False)


def c1_operator(f,source,target):
    """Historical C1 weights including Vj; no use of any PD2 tensor in old Z."""
    i,j=f.pairs.T;n=len(f.X);K=np.zeros((n,3,3))
    q=f.omega[:,None,None]*f.xi[:,:,None]*f.xi[:,None,:]
    np.add.at(K,i,q*f.volume[j,None,None]);np.add.at(K,j,q*f.volume[i,None,None])
    solved=np.linalg.solve(K,np.broadcast_to(np.eye(3),K.shape))
    from ..kernel import radial
    xi=f.X[target]-f.X[source]
    g=np.einsum('bij,bj->bi',solved[source],xi)
    g*= (radial(np.linalg.norm(xi,axis=1)/f.horizon)*f.volume[target])[:,None]
    return GradientOperator('C1',source,target,g)


@dataclass
class CandidateState:
    energy: float
    force: np.ndarray
    Fbar: np.ndarray


class VariationalBAQP:
    def __init__(self,reference,gradient_operator,material=None,bond_block=20000,directed_block=100000):
        self.reference=reference;self.operator=gradient_operator
        self.material=Material('svk',1000.,2000./3) if material is None else material
        if self.material.model!='svk' or self.material.mu!=1000. or self.material.lame!=0.:
            raise ValueError('This candidate is the intact lambda=0, mu=1000 Pa SVK model only')
        self.bond_block=bond_block;self.directed_block=directed_block

    def nodes(self,x):
        x=np.asarray(x,float)
        if x.shape!=self.reference.X.shape or not np.isfinite(x).all():raise ValueError('Invalid position')
        return deformation(self.reference.X,x,self.operator.source,self.operator.target,self.operator.g)

    def energy(self,x):
        return energy_from_nodes(self.reference,x,self.nodes(x),self.material,self.bond_block)

    def evaluate(self,x):
        F=self.nodes(x)
        U,force=energy_force(self.reference,self.operator,x,F,self.material,self.bond_block,self.directed_block)
        return CandidateState(U,force,F)
