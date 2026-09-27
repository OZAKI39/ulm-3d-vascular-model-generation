"""Explicit float64 barycentric mathematics; no VTK interpolation."""
import numpy as np

EPS = np.finfo(np.float64).eps


def readonly(value, dtype=np.float64):
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


class TetraGeometry:
    """Barycentric weights tell how strongly each of four vertices contributes.

    A = [x1-x0, x2-x0, x3-x0]; N[1:] = inv(A) @ (x-x0).
    Error budget tau = 64 eps cond_inf(A) (1 + coordinate_scale / ||A||inf).
    The scale ratio accounts for subtraction of untranslated SI coordinates;
    cond_inf accounts for the tetra's shape. 64 covers the short subtraction,
    3x3 inversion and matrix/vector operations, including the N0 sum.
    This is a roundoff allowance, not a chosen physical boundary thickness.
    Ill-resolved tetrahedra (tau >= sqrt(eps)) are rejected, never repaired.
    """

    def __init__(self, points_m, tetra):
        self.points = readonly(points_m)
        raw = np.asarray(tetra)
        if raw.ndim != 2 or raw.shape[1] != 4 or not np.issubdtype(raw.dtype, np.integer):
            raise ValueError("tetra must be an integer array with shape (M,4)")
        self.tetra = readonly(raw, np.int64)
        if self.points.ndim != 2 or self.points.shape[1] != 3 or not np.isfinite(self.points).all():
            raise ValueError("points_m must be finite with shape (N,3)")
        if not len(raw) or raw.min() < 0 or raw.max() >= len(self.points):
            raise ValueError("empty tetra mesh or invalid global point index")
        xyz = self.points[self.tetra]
        self.origins = readonly(xyz[:, 0])
        a = np.swapaxes(xyz[:, 1:] - xyz[:, :1], 1, 2)
        try:
            inverse = np.linalg.inv(a)
        except np.linalg.LinAlgError as exc:
            raise ValueError("degenerate tetrahedron") from exc
        scale = np.linalg.norm(a, ord=np.inf, axis=(1, 2))
        condition = scale * np.linalg.norm(inverse, ord=np.inf, axis=(1, 2))
        coordinate_scale = np.max(np.abs(xyz), axis=(1, 2))
        tau = 64 * EPS * condition * (1 + coordinate_scale / scale)
        if not np.isfinite(tau).all() or np.any(tau >= np.sqrt(EPS)):
            raise ValueError("tetra geometry is not resolvable with the float64 roundoff contract")
        self.inverse = readonly(inverse)
        self.scale_m = readonly(scale)
        self.condition_inf = readonly(condition)
        self.weight_tolerance = readonly(tau)
        # If each N can be negative by tau, any coordinate can escape a cell
        # AABB by at most 4*tau*||A||inf. Bounds padding only gathers candidates;
        # each candidate still must pass its own barycentric test below.
        self.candidate_padding_m = float(np.max(4 * tau * scale))
        self.volume_m3 = readonly(np.abs(np.linalg.det(a)) / 6)
        gradients = np.concatenate((-inverse.sum(axis=1)[:, None, :], inverse), axis=1)
        self.shape_gradients_m_inv = readonly(gradients)

    def weights(self, position_m, tetra_ids):
        ids = np.asarray(tetra_ids, dtype=np.int64)
        tail = np.einsum("...ij,...j->...i", self.inverse[ids],
                         np.asarray(position_m, dtype=np.float64) - self.origins[ids])
        return np.concatenate((1 - tail.sum(axis=-1, keepdims=True), tail), axis=-1)

    def contains(self, weights, tetra_ids):
        tau = self.weight_tolerance[tetra_ids]
        return (np.min(weights, axis=-1) >= -tau) & (np.max(weights, axis=-1) <= 1 + tau)
