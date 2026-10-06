"""Force is minus the adjoint of the same intact bond-energy graph."""
import numpy as np
from .variational_energy import bond_state
from .adjoint_accumulation import bond_adjoint,nodal_adjoint


def energy_force(reference,operator,x,F,material,bond_block=20000,directed_block=100000):
    G=np.zeros_like(F);gradU=np.zeros_like(x);energy=0.
    for start in range(0,len(reference.pairs),bond_block):
        sl=slice(start,start+bond_block)
        _,W,P=bond_state(reference,x,F,sl,material)
        energy+=float(np.sum(reference.m[sl]*W))
        bond_adjoint(reference,P,sl,G,gradU)
    nodal_adjoint(operator,G,gradU,directed_block)
    return energy,-gradU
