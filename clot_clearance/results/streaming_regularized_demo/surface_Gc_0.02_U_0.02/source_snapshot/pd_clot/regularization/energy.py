"""SI energy accounting for the existing NOSB kernel and dissipative updates.

Damage dissipation is measured as a potential drop at fixed coordinates, not
invented as Gc times a guessed crack area. Negative drops stop the candidate.
"""
from dataclasses import dataclass
import numpy as np


def potential(cloud, mechanics):
    elastic = float(np.dot(cloud.volume, mechanics[3]))
    stabilization = float(np.dot(cloud.volume, mechanics[4]))
    if not np.isfinite([elastic, stabilization]).all():
        raise FloatingPointError('Nonfinite material energy')
    return elastic, stabilization


def kinetic(mass, velocity):
    return float(.5 * np.sum(mass[:, None] * velocity**2))


def bond_quadrature_volume(cloud):
    """Partition each particle volume among its ORIGINAL weighted family.

    v_ij = Vi wij Vj/Si + Vj wij Vi/Sj, Si=sum_j wij Vj.
    Unordered-bond sum equals total material volume. No fracture Gc formula is
    implied: a separate numerical coupon determines the critical energy density.
    """
    i, j = cloud.pairs.T
    total = np.bincount(i, weights=cloud.weight*cloud.volume[j], minlength=len(cloud.X))
    total += np.bincount(j, weights=cloud.weight*cloud.volume[i], minlength=len(cloud.X))
    if np.any(total <= 0):
        raise ValueError('Initial cloud has unsupported particles')
    return cloud.weight*cloud.volume[i]*cloud.volume[j]*(1/total[i]+1/total[j])


@dataclass
class EnergyLedger:
    initial_energy_J: float
    external_work_J: float = 0.
    damping_dissipation_J: float = 0.
    damage_dissipation_J: float = 0.
    transport_work_J: float = 0.
    negative_damage_roundoff_J: float = 0.

    def damage_drop(self, before_J, after_J, allocation_weights):
        drop = float(before_J-after_J)
        tolerance = 1e-12*max(abs(before_J), abs(after_J), 1e-16)+1e-20
        if drop < -tolerance:
            raise FloatingPointError(f'Damage increased stored energy by {-drop:.9g} J; energy admissibility failed')
        if drop < 0:
            self.negative_damage_roundoff_J += -drop
        released = max(0., drop)
        weights = np.maximum(np.asarray(allocation_weights), 0.)
        if released > tolerance and weights.sum() == 0:
            raise ValueError('Unattributed fracture energy')
        self.damage_dissipation_J += released
        return released*weights/weights.sum() if weights.sum() else np.zeros_like(weights)

    def row(self, elastic_J, stabilization_J, kinetic_J):
        residual = (elastic_J+stabilization_J+kinetic_J-self.initial_energy_J
                    +self.damping_dissipation_J+self.damage_dissipation_J
                    -self.external_work_J-self.transport_work_J)
        scale = max(abs(self.external_work_J)+abs(self.transport_work_J),
                    abs(elastic_J)+abs(stabilization_J)+abs(kinetic_J), 1e-30)
        row = dict(elastic_energy_J=elastic_J, stabilization_energy_J=stabilization_J,
                   kinetic_energy_J=kinetic_J, external_work_J=self.external_work_J,
                   transport_work_J=self.transport_work_J,
                   damping_dissipation_J=self.damping_dissipation_J,
                   damage_dissipation_estimate_J=self.damage_dissipation_J,
                   numerical_energy_residual_J=residual,
                   relative_numerical_energy_residual=residual/scale,
                   negative_damage_roundoff_J=self.negative_damage_roundoff_J)
        if not np.isfinite(list(row.values())).all():
            raise FloatingPointError('Nonfinite energy ledger')
        return row
