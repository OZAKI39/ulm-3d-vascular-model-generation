"""Isolated AT1 constrained minimization; not PFPD v1 Eq.3.25.

Bond densities are Pa. Multiply by the existing BAQP pair quadrature for J.
No erosion state, mechanics, kinematic degradation, or history proxy is used.
"""
from dataclasses import dataclass
import numpy as np
from .normalization import c0

NAME = 'PFPD_AT1_CONSISTENT_CANDIDATE'


def critical_driver(Gc, delta, kernel='quadratic'):
    if not np.isfinite([Gc, delta]).all() or Gc <= 0 or delta <= 0:
        raise ValueError('Positive finite Gc and delta required')
    value = float(Gc / (2 * c0(kernel)['c0'] * delta))
    if not np.isfinite(value) or value <= 0:
        raise ValueError('Unrepresentable critical driver')
    return value


def _inputs(Y, Yc, s_old=0.):
    y, yc, old = (np.array(a, dtype=float, copy=True) for a in
                  np.broadcast_arrays(Y, Yc, s_old))
    if not all(np.isfinite(a).all() for a in (y, yc, old)):
        raise ValueError('Finite phase and drivers required')
    if np.any(y < 0) or np.any(yc <= 0) or np.any((old < 0) | (old > 1)):
        raise ValueError('Require Y >= 0, Yc > 0, and 0 <= s_old <= 1')
    return y, yc, old


def at1_trial_state(Y, Yc):
    y, yc, _ = _inputs(Y, Yc)
    ratio = np.divide(yc, y, out=np.ones_like(y), where=y > yc)
    return np.where(y > yc, 1 - ratio, 0.)


@dataclass(frozen=True)
class AT1Update:
    s_old: np.ndarray
    s_trial: np.ndarray
    s_new: np.ndarray
    active_set_status: np.ndarray
    derivative_at_solution: np.ndarray
    elastic_before: np.ndarray
    elastic_after: np.ndarray
    crack_before: np.ndarray
    crack_after: np.ndarray
    total_before: np.ndarray
    total_after: np.ndarray


class HighDamageGate(RuntimeError):
    """Proposal is available for audit; it has NOT been accepted or clipped."""
    def __init__(self, update, cap):
        super().__init__('HIGH_DAMAGE_KINEMATIC_COUPLING_NOT_IMPLEMENTED')
        self.update = update
        self.cap = cap


def at1_constrained_update(s_old, Y, Yc, s_cap=None):
    y, yc, old = _inputs(Y, Yc, s_old)
    if s_cap is not None and (not np.isfinite(s_cap) or not 0 < s_cap <= 1):
        raise ValueError('Require 0 < s_cap <= 1')
    trial = at1_trial_state(y, yc)
    new = np.maximum(old, trial)
    status = np.where(new == old, 'LOWER_BOUND', 'INTERIOR')
    status = np.where(old == 1, 'FIXED_SINGLETON', status)
    with np.errstate(over='raise', invalid='raise'):
        eb, ea = (1-old)**2*y, (1-new)**2*y
        cb, ca = 2*yc*old, 2*yc*new
        q = -2*(1-new)*y + 2*yc
        update = AT1Update(old, trial, new, status, q, eb, ea, cb, ca, eb+cb, ea+ca)
    if not all(np.isfinite(getattr(update, k)).all() for k in update.__dataclass_fields__
               if k != 'active_set_status'):
        raise FloatingPointError('Nonfinite AT1 update')
    # Reaching the cap also stops, as required by the coupled-run acceptance gate.
    if s_cap is not None and np.any(new >= s_cap):
        raise HighDamageGate(update, s_cap)
    return update


def audit_update(update, Y, Yc, *, kkt_tolerance=5e-12, energy_tolerance=5e-13):
    y, yc, _ = _inputs(Y, Yc, update.s_old)
    scale = np.maximum(y, yc)
    lower = update.active_set_status == 'LOWER_BOUND'
    interior = update.active_set_status == 'INTERIOR'
    scaled_q = update.derivative_at_solution / scale
    energy_change = (update.total_after-update.total_before) / scale
    passed = bool(np.all(update.s_new >= update.s_old) and np.all(update.s_new <= 1)
                  and np.all(scaled_q[lower] >= -kkt_tolerance)
                  and np.all(abs(scaled_q[interior]) <= kkt_tolerance)
                  and np.all(energy_change <= energy_tolerance))
    return dict(passed=passed, maximum_scaled_interior_derivative=float(np.max(
        abs(scaled_q[interior]), initial=0.)), minimum_scaled_lower_derivative=float(
        np.min(scaled_q[lower], initial=0.)), maximum_scaled_energy_increase=float(
        np.max(energy_change, initial=0.)))
