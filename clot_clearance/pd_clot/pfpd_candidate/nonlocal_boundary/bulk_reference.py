"""Nonaffine input adapter for the previously verified selective BAQP force.

Only the analytic nodal positions differ. All neighbor/volume/K/F/Z/force
methods and closure checks are inherited unchanged from AffineStencil.
"""
import numpy as np
from ..force_stencil_audit import AffineStencil, H
from .nonaffine_manufactured import QuadraticAxial, QuadraticShear, stress_channel


class ManufacturedStencil(AffineStencil):
    def __init__(self, field, n=27, center=(0., 0., 0.)):
        super().__init__(n, epsilon=field.epsilon)
        self.X = self.X+np.asarray(center)
        self.x = field.position(self.X)
        # Match evaluate(): reconstructed gradients use x-X, including rounding.
        self.u = self.x-self.X

    def force(self, pid): return self._force_channel(pid, 'original')['force']


class StressChannelStencil(ManufacturedStencil):
    """Same selective assembly, projecting only the actual bond Piola stress."""
    def __init__(self, field, n=27, center=(0.,0.,0.), channel='combined'):
        super().__init__(field,n,center); self.stress_channel = channel

    def bond_stress(self, low, high):
        return stress_channel(super().bond_stress(low,high),self.stress_channel)


def oracle_job(spec):
    """Small process-worker job; no global bond tensors or physics modification."""
    field = (QuadraticShear(gamma=spec.get('gamma',.003), beta_s=spec.get('beta',.001))
             if spec.get('kind')=='shear' else QuadraticAxial(beta=spec.get('beta', .001)))
    n = spec.get('n', 27); offset = np.asarray(spec.get('offset', [0, 0, 0]), dtype=int)
    point = np.array([spec.get('x', 0.), spec.get('y', 0.), spec['Z']])
    stencil = (StressChannelStencil(field,n,point-offset*H,spec['channel']) if 'channel' in spec
               else ManufacturedStencil(field, n, point-offset*H))
    pid = stencil.id_at_offset(offset); f = stencil.force(pid)
    return dict(**spec, observed_X_m=stencil.X[pid], force_N=f, q_bulk_N_m3=f/H**3,
                nodal_volume_m3=float(stencil.V[pid]), needed_Fbar_nodes=len(stencil._Fbar),
                cached_Z_nodes=len(stencil._Z)//3)
