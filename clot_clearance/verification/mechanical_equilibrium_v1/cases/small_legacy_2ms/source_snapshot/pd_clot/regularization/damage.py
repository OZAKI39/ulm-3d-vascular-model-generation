"""Coupon-calibrated cyclic bond-energy scale, irreversible D in [0,1]."""
import json
from pathlib import Path
import warnings
import numpy as np
from ..fragment_topology import BondDamage as LegacyBondDamage
from .coupon import cyclic_increment, young_modulus
from .energy import bond_quadrature_volume


class LegacyDamageWithWarning(LegacyBondDamage):
    def __init__(self, count, config):
        super().__init__(count, config)
        self.forced_failure_verification = config.get('D_break', 1.) < 1
        if self.forced_failure_verification:
            warnings.warn('D_break < 1: forced_failure_verification = true; legacy verification only', RuntimeWarning, stacklevel=2)


class EnergyBondDamage:
    def __init__(self, cloud, config, calibration):
        self.c = config['regularized_damage']
        if self.c.get('D_break', 1.) != 1. or config['damage'].get('D_break', 1.) != 1.:
            raise ValueError('Regularized candidate requires D_break = 1')
        if config['damage']['mode'] != 'energy_regularized_cyclic':
            raise ValueError('Wrong damage mode')
        if isinstance(calibration, (str, Path)):
            calibration = json.loads(Path(calibration).read_text())
        for name, actual in [('spacing_m', cloud.spacing), ('horizon_m', cloud.horizon)]:
            if not np.isclose(calibration[name], actual, rtol=1e-12, atol=1e-15):
                raise ValueError(f'Cannot reuse fracture calibration: mismatched {name}')
        if calibration['constitutive_material'] != config['material']:
            raise ValueError('Cannot reuse calibration for different material')
        for name in ['C_E', 'm_E', 'threshold_ratio']:
            if calibration['fatigue_parameters'].get(name, 0.) != self.c.get(name, 0.):
                raise ValueError(f'Cannot reuse calibration for different {name}')
        if not np.isclose(calibration['Gc_demo_J_m2'], self.c['Gc_demo_J_m2'], rtol=1e-12):
            raise ValueError('Wrong Gc target in calibration')
        self.quadrature_volume = bond_quadrature_volume(cloud)
        self.critical_J = calibration['critical_energy_density_Pa']*self.quadrature_volume
        self.young_Pa = young_modulus(config['material'])
        self.D = np.zeros(len(cloud.pairs))
        self.g = np.ones(len(cloud.pairs))
        self.active = np.ones(len(cloud.pairs), bool)
        self.dissipated_J = np.zeros(len(cloud.pairs))
        self.last_driver_J = np.zeros(len(cloud.pairs))
        self.forced_failure_verification = False

    def advance(self, amplitude, cycles):
        if not np.isfinite(amplitude).all() or np.any(amplitude < 0):
            raise ValueError('Invalid measured cyclic strain amplitude')
        self.last_driver_J = .5*self.young_Pa*np.asarray(amplitude)**2*self.quadrature_volume
        increment = cyclic_increment(self.last_driver_J, self.critical_J, cycles,
                                     self.c['C_E'], self.c['m_E'], self.c.get('threshold_ratio', 0.))
        old, old_active = self.D.copy(), self.active.copy()
        self.D = np.where(self.active, np.minimum(1., self.D+increment), self.D)
        self.active &= self.D < 1.
        self.g = 1.-self.D
        assert np.all(self.D >= old) and not np.any(self.active & ~old_active)
        return np.flatnonzero(old_active & ~self.active), old, increment
