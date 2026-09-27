"""Body-to-world Hamilton quaternions (w,x,y,z); p=R(q)@(0,0,1).

Omega is a world-frame angular velocity: q_new=delta_q_world ⊗ q_old.
Helpers support independent validation batches; no inter-particle coupling.
"""
import numpy as np
from .microbubble import positive_scalar


def normalize_quaternion(value):
    q = np.asarray(value, dtype=np.float64)
    if q.ndim < 1 or q.shape[-1] != 4 or not np.isfinite(q).all():
        raise ValueError("quaternion must be finite (...,4) wxyz")
    norm = np.linalg.norm(q, axis=-1, keepdims=True)
    if not np.isfinite(norm).all() or np.any(norm == 0):
        raise ValueError("quaternion norm must be positive and finite")
    return q / norm


def quaternion_multiply(left, right):
    a, b = np.asarray(left, dtype=np.float64), np.asarray(right, dtype=np.float64)
    scalar = a[..., :1] * b[..., :1] - np.sum(a[..., 1:] * b[..., 1:], axis=-1, keepdims=True)
    vector = a[..., :1] * b[..., 1:] + b[..., :1] * a[..., 1:] + np.cross(a[..., 1:], b[..., 1:])
    return np.concatenate((scalar, vector), axis=-1)


def rotation_matrix(q):
    q = normalize_quaternion(q)
    w, x, y, z = np.moveaxis(q, -1, 0)
    return np.stack((1 - 2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y),
                     2*(x*y+w*z), 1 - 2*(x*x+z*z), 2*(y*z-w*x),
                     2*(x*z-w*y), 2*(y*z+w*x), 1 - 2*(x*x+y*y)), axis=-1).reshape(q.shape[:-1] + (3, 3))


def short_axis(q):
    return rotation_matrix(q)[..., :, 2]


def quaternion_from_short_axis(axis):
    p = np.asarray(axis, dtype=np.float64)
    if p.shape != (3,) or not np.isfinite(p).all() or np.linalg.norm(p) == 0:
        raise ValueError("short axis must be a finite nonzero shape (3,) vector")
    p = p / np.linalg.norm(p)
    if p[2] < -1 + 16*np.finfo(float).eps:
        return np.array([0., 1., 0., 0.])
    return normalize_quaternion(np.array([1 + p[2], -p[1], p[0], 0.]))


def advance_orientation(q, omega_world_s_inv, dt_s):
    """Frozen-old-Omega incremental rotation, normalization, equivalent sign."""
    dt = positive_scalar(dt_s, "dt_s")
    old = normalize_quaternion(q)
    omega = np.asarray(omega_world_s_inv, dtype=np.float64)
    if omega.shape != old.shape[:-1] + (3,) or not np.isfinite(omega).all():
        raise ValueError("Omega must be finite and match quaternion batch shape")
    magnitude = np.linalg.norm(omega, axis=-1)
    half_angle = magnitude * dt / 2
    # sin(theta/2)/|omega|, with its analytic zero limit dt/2.
    factor = (dt / 2) * np.sinc(half_angle / np.pi)
    delta = np.concatenate((np.cos(half_angle)[..., None], omega * factor[..., None]), axis=-1)
    updated = normalize_quaternion(quaternion_multiply(delta, old))
    updated = np.where((np.sum(updated * old, axis=-1) < 0)[..., None], -updated, updated)
    # Zero Omega preserves the rotation, also for a non-unit input representation.
    return np.where((magnitude == 0)[..., None], old, updated)


def decompose_gradient(gradient_s_inv):
    g = np.asarray(gradient_s_inv, dtype=np.float64)
    if g.shape != (3, 3) or not np.isfinite(g).all():
        raise ValueError("G must be finite shape (3,3), in s^-1")
    return .5 * (g + g.T), .5 * (g - g.T)


def angular_velocity(p, jeffery_lambda, gradient_s_inv, vorticity_s_inv=None):
    """Omega_world=curl(u)/2+lambda*(p cross E@p); p is the short axis."""
    g = np.asarray(gradient_s_inv, dtype=np.float64)
    e, _ = decompose_gradient(g)
    curl = np.array([g[2,1]-g[1,2], g[0,2]-g[2,0], g[1,0]-g[0,1]]) if vorticity_s_inv is None else np.asarray(vorticity_s_inv)
    p = np.asarray(p, dtype=np.float64)
    return .5 * curl + np.asarray(jeffery_lambda)[..., None] * np.cross(p, p @ e.T)


def jeffery_axis_derivative(p, jeffery_lambda, gradient_s_inv):
    e, w = decompose_gradient(gradient_s_inv)
    p = np.asarray(p, dtype=np.float64)
    ep = p @ e.T
    return p @ w.T + np.asarray(jeffery_lambda)[..., None] * (ep - np.sum(p * ep, axis=-1, keepdims=True) * p)


def shape_axis_distance(p, other):
    """Unoriented axis distance in radians; p and -p are the same shape."""
    p, other = np.asarray(p), np.asarray(other)
    p = p / np.linalg.norm(p, axis=-1, keepdims=True)
    other = other / np.linalg.norm(other, axis=-1, keepdims=True)
    return np.arccos(np.clip(np.abs(np.sum(p * other, axis=-1)), 0., 1.))
