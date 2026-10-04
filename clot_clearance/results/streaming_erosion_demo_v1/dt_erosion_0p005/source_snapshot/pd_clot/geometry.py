from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree


@dataclass
class Cloud:
    X: np.ndarray
    volume: np.ndarray
    fixed: np.ndarray
    exposed: np.ndarray
    face_particle: np.ndarray
    face_position: np.ndarray
    face_normal: np.ndarray
    face_area: np.ndarray
    pairs: np.ndarray
    xi: np.ndarray
    length: np.ndarray
    weight: np.ndarray
    spacing: float
    horizon: float


def make_cloud(c):
    h = float(c['particle_spacing_m'])
    shape = np.asarray(c['cells'], dtype=int)
    if h <= 0 or np.any(shape < 2):
        raise ValueError('Positive spacing and at least two cells per axis required')
    ijk = np.indices(shape).reshape(3, -1).T
    X = (ijk + .5) * h + np.asarray(c['origin_m'])
    fixed = X[:, 2] < c['origin_m'][2] + c['fixed_base_thickness_m']
    if not fixed.any() or fixed.all():
        raise ValueError('Fixed base must contain some, but not all particles')
    fp, fx, fn = [], [], []
    for axis in range(3):
        for side in [-1, 1]:
            if axis == 2 and side == -1:
                continue  # support interface is not exposed to fluid
            ids = np.flatnonzero(ijk[:, axis] == (0 if side == -1 else shape[axis] - 1))
            n = np.eye(3)[axis] * side
            fp.extend(ids); fx.extend(X[ids] + .5 * h * n)
            fn.extend(np.tile(n, (len(ids), 1)))
    fp = np.asarray(fp, dtype=np.int64)
    exposed = np.zeros(len(X), dtype=bool); exposed[fp] = True
    horizon = h * c['horizon_ratio']
    pairs = cKDTree(X).query_pairs(horizon, output_type='ndarray').astype(np.int64)
    xi = X[pairs[:, 1]] - X[pairs[:, 0]]
    length = np.linalg.norm(xi, axis=1)
    # Nonnegative compact influence; no dimensional normalization hidden here.
    weight = (1 - length / horizon) ** 2
    return Cloud(X, np.full(len(X), h**3), fixed, exposed, fp,
                 np.asarray(fx), np.asarray(fn), np.full(len(fp), h*h),
                 pairs, xi, length, weight, h, horizon)
