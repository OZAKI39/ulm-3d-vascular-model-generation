"""P5 sphere-only, leading normal resistance; no production cutoff or gap floor."""
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import numpy as np
from .audit import sha256
from .particle_shapes import Sphere, unit

VALIDATION_NEARFIELD_RATIO_MAX = .01
FORMULATION = 'AFFINE_BLOCK_RESISTANCE_BALANCE_V0'


class NonsphericalLubricationNotFrozen(ValueError):
    def __init__(self):
        self.record = dict(status='NONSPHERICAL_LUBRICATION_NOT_FROZEN',
            active=False, reason='NONSPHERICAL_MODEL_NOT_FROZEN',
            gap_m=None, relevant_radius_m=None, gap_ratio=None, eligibility_role='PHYSICAL_API_REJECTION')
        super().__init__('NONSPHERICAL_LUBRICATION_NOT_FROZEN')


def require_sphere(shape):
    if not isinstance(shape, Sphere):
        raise NonsphericalLubricationNotFrozen()


def viscosity_from_frozen(fem):
    """Cross-check two frozen scientific inputs, including rho*nu, before use."""
    root = Path(fem)/'frozen_reference'
    summary, solver = root/'baseline_summary.json', root/'run/solver.xml'
    physics = json.loads(summary.read_text(), parse_float=Decimal)['physics']
    xml = ET.parse(solver)
    rho = Decimal(xml.findtext('.//Add_equation/Density'))
    mu = Decimal(xml.findtext('.//Viscosity/Value'))
    expected = Decimal('0.00345312')
    actual = physics['dynamic_viscosity_pa_s']
    if not (mu == actual == expected and rho == physics['density_kg_m3'] == Decimal('1056')
            and physics['kinematic_viscosity_m2_s'] == Decimal('3.27e-6')
            and rho*physics['kinematic_viscosity_m2_s'] == mu):
        raise ValueError(f'VISCOSITY_PROVENANCE_MISMATCH expected={expected} actual={actual},xml={mu} sources={summary},{solver}')
    return float(mu), dict(status='PASS', dynamic_viscosity_pa_s=float(mu),
        density_kg_m3=float(rho), kinematic_viscosity_m2_s=float(physics['kinematic_viscosity_m2_s']),
        sources=[dict(path=str(p), sha256=sha256(p)) for p in [summary, solver]])


def sphere_self_diagonal(shape, mu):
    require_sphere(shape)
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError('Positive finite dynamic viscosity required')
    a = shape.radius_m
    return np.array([6*np.pi*mu*a]*3+[8*np.pi*mu*a**3]*3)


@dataclass(frozen=True)
class PhysicalNearField:
    """Geometric specification only: the assembler computes its coefficient.

    Arbitrary validation coefficients cannot enter this physical interface.
    VALIDATION_GAPS_ONLY permits prescribed synthetic ratios, never real replay.
    """
    particle_i_id: int
    particle_j_id: int | None
    gap_m: float
    normal: np.ndarray
    roundoff_m: float
    eligibility_role: str = 'VALIDATION_ONLY'


def evaluate_block(spec, particles, mu):
    if type(spec) is not PhysicalNearField:
        raise TypeError('VALIDATION_BLOCK_CANNOT_ENTER_PHYSICAL_SOLVER')
    a = particles[spec.particle_i_id]
    require_sphere(a)
    relevant = a.radius_m
    if spec.particle_j_id is not None:
        if spec.particle_j_id == spec.particle_i_id:
            raise ValueError('Distinct pair IDs required')
        b = particles[spec.particle_j_id]
        require_sphere(b)
        relevant = a.radius_m*b.radius_m/(a.radius_m+b.radius_m)
    if spec.eligibility_role not in ['VALIDATION_ONLY', 'VALIDATION_GAPS_ONLY']:
        raise ValueError('Production lubrication cutoff is NOT_FROZEN')
    h, ro = float(spec.gap_m), float(spec.roundoff_m)
    if not np.isfinite([h, ro]).all() or ro < 0:
        raise ValueError('Finite gap and nonnegative geometry roundoff required')
    ratio = h/relevant
    reason = ('CONTACT_REGIME' if h <= ro else
              'VALIDATION_NEAR_FIELD' if spec.eligibility_role == 'VALIDATION_GAPS_ONLY' or ratio <= VALIDATION_NEARFIELD_RATIO_MAX
              else 'OUTSIDE_VALIDATION_ASYMPTOTIC_RANGE')
    active = reason == 'VALIDATION_NEAR_FIELD'
    # No max(h, floor): contact is a separate geometric constraint.
    coefficient = 6*np.pi*mu*relevant**2/h if active else 0.
    return coefficient, unit(spec.normal), dict(particle_i_id=spec.particle_i_id,
        particle_j_id=spec.particle_j_id, kind='WALL' if spec.particle_j_id is None else 'PAIR',
        gap_m=h, geometry_roundoff_m=ro, relevant_radius_m=relevant, gap_ratio=ratio,
        eligibility_role=spec.eligibility_role, active=active, reason=reason,
        coefficient_kg_s=coefficient, model='LEADING_NORMAL_ASYMPTOTIC_V0')
