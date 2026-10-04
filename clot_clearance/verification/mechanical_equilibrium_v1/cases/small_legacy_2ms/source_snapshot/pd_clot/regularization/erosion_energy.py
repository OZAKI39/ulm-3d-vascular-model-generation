"""Signed energy accounting: the strain-rate erosion law is NOT a Gc law.

Recomputing the damaged NOSB shape may increase stored energy at fixed x.
Retain that sign explicitly; do not clip it into positive dissipation or tune
the material. The old cyclic EnergyLedger and its rejection rule are untouched.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class ErosionEnergyLedger:
    initial_energy_J: float
    external_work_J: float=0.
    damping_dissipation_J: float=0.
    transport_work_J: float=0.
    damage_stored_energy_drop_J: float=0.
    damage_positive_energy_drop_J: float=0.
    damage_energy_increase_J: float=0.
    damage_energy_increase_count: int=0

    def damage_drop(self,before_J,after_J,allocation_weights):
        drop=float(before_J-after_J)
        if not np.isfinite(drop):raise FloatingPointError('Nonfinite damage energy change')
        weights=np.maximum(np.asarray(allocation_weights),0.)
        tolerance=1e-12*max(abs(before_J),abs(after_J),1e-16)+1e-20
        if abs(drop)>tolerance and weights.sum()==0:raise ValueError('Unattributed damage energy change')
        self.damage_stored_energy_drop_J+=drop
        self.damage_positive_energy_drop_J+=max(drop,0.)
        self.damage_energy_increase_J+=max(-drop,0.)
        if drop < -tolerance:self.damage_energy_increase_count+=1
        return drop*weights/weights.sum() if weights.sum() else np.zeros_like(weights)

    def row(self,elastic_J,stabilization_J,kinetic_J):
        residual=(elastic_J+stabilization_J+kinetic_J-self.initial_energy_J+
                  self.damping_dissipation_J+self.damage_stored_energy_drop_J-
                  self.external_work_J-self.transport_work_J)
        scale=max(abs(self.external_work_J)+abs(self.transport_work_J),
                  abs(elastic_J)+abs(stabilization_J)+abs(kinetic_J),1e-30)
        row=dict(elastic_energy_J=elastic_J,stabilization_energy_J=stabilization_J,kinetic_energy_J=kinetic_J,
                 external_work_J=self.external_work_J,transport_work_J=self.transport_work_J,
                 damping_dissipation_J=self.damping_dissipation_J,
                 damage_stored_energy_drop_J=self.damage_stored_energy_drop_J,
                 damage_positive_energy_drop_J=self.damage_positive_energy_drop_J,
                 damage_energy_increase_J=self.damage_energy_increase_J,
                 damage_energy_increase_count=self.damage_energy_increase_count,
                 numerical_energy_residual_J=residual,relative_numerical_energy_residual=residual/scale)
        if not np.isfinite(list(row.values())).all():raise FloatingPointError('Nonfinite erosion energy ledger')
        row['damage_energy_status']=('ENERGY_DISSIPATION_NOT_ESTABLISHED' if self.damage_energy_increase_count else
                                     'NONNEGATIVE_STORED_ENERGY_DROPS_OBSERVED_ONLY')
        return row
