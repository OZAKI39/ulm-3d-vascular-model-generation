"""Explicit reference-domain faces and quadrature; geometric rules only."""
from dataclasses import dataclass
import numpy as np
from scipy.spatial import cKDTree
from ..kernel import Family

FACE_NAMES = ('x-', 'x+', 'y-', 'y+', 'z-', 'z+')


@dataclass(frozen=True)
class BlockBoundary:
    lower: np.ndarray
    upper: np.ndarray
    dirichlet_faces: tuple = ('z-', 'z+')

    def __post_init__(self):
        lo, hi = np.asarray(self.lower, float), np.asarray(self.upper, float)
        if lo.shape != (3,) or hi.shape != (3,) or np.any(hi <= lo):
            raise ValueError('Finite positive Cartesian block required')
        if not np.isfinite(np.r_[lo, hi]).all() or not set(self.dirichlet_faces) <= set(FACE_NAMES):
            raise ValueError('Invalid block or boundary partition')
        object.__setattr__(self, 'lower', lo); object.__setattr__(self, 'upper', hi)

    @property
    def center(self): return .5*(self.lower+self.upper)

    @property
    def lengths(self): return self.upper-self.lower

    @property
    def volume(self): return float(np.prod(self.lengths))

    @property
    def normals(self): return np.array([s*np.eye(3)[a] for a in range(3) for s in (-1, 1)])

    @property
    def neumann_faces(self): return tuple(k for k in FACE_NAMES if k not in self.dirichlet_faces)

    def face_distances(self, X):
        X = np.asarray(X)
        return np.column_stack([X[:, a]-self.lower[a] if s == -1 else self.upper[a]-X[:, a]
                                for a in range(3) for s in (-1, 1)])

    def distance(self, X): return self.face_distances(X).min(axis=1)

    def project(self, X, face=None):
        """Nearest face with deterministic corner ties; return all face normals separately when needed."""
        X = np.asarray(X, float)
        indices = np.argmin(abs(self.face_distances(X)), axis=1) if face is None else np.full(len(X), FACE_NAMES.index(face))
        Y = np.clip(X, self.lower, self.upper)
        for j, (a, s) in enumerate((a, s) for a in range(3) for s in (-1, 1)):
            Y[indices == j, a] = self.lower[a] if s == -1 else self.upper[a]
        return Y, self.normals[indices], indices

    def tangent_basis(self, face):
        n = self.normals[FACE_NAMES.index(face)]
        e = np.eye(3)[np.argmin(abs(n))]
        t1 = np.cross(n, e); return t1, np.cross(n, t1)

    def masks(self, X, delta):
        d = self.distance(X); tol = np.max(self.lengths)*1e-12
        physical = np.all(self.face_distances(X) >= -tol, axis=1)
        return dict(physical=physical, boundary=physical & (d <= tol),
                    direct_horizon_region=physical & (d < delta),
                    geometric_core=physical & (d >= delta),
                    exterior_interaction=(~physical) & (np.linalg.norm(np.asarray(X)-np.clip(X, self.lower, self.upper), axis=1) < delta))

    def face_mask(self, X, face):
        return abs(self.face_distances(X)[:, FACE_NAMES.index(face)]) <= np.max(self.lengths)*1e-12

    def dirichlet_mask(self, X):
        mask = np.zeros(len(X), dtype=bool)
        for face in self.dirichlet_faces:
            mask |= self.face_mask(X, face)
        return mask

    def face_quadrature(self, X):
        """Cartesian nodal trapezoidal areas; edges belong to each incident face."""
        X = np.asarray(X)
        nodal_widths = []
        for axis in range(3):
            levels = np.unique(X[:, axis])
            edges = np.r_[levels[0], .5*(levels[:-1]+levels[1:]), levels[-1]]
            nodal_widths.append(np.diff(edges)[np.searchsorted(levels, X[:, axis])])
        rows = []
        for j, face in enumerate(FACE_NAMES):
            ids = np.flatnonzero(self.face_mask(X, face)); axes = [a for a in range(3) if a != j//2]
            area = nodal_widths[axes[0]][ids]*nodal_widths[axes[1]][ids]
            rows.append(dict(face=face, ids=ids, X=X[ids], normal=self.normals[j], area=area))
        return rows


def block_family(n, lengths, delta, center=(0., 0., 0.)):
    lengths = np.broadcast_to(lengths, (3,)).astype(float); center = np.asarray(center)
    shape = np.broadcast_to(n, (3,)).astype(int)
    geometry = BlockBoundary(center-lengths/2, center+lengths/2)
    ijk = np.indices(tuple(shape)).reshape(3, -1).T
    X = np.column_stack([np.linspace(geometry.lower[a], geometry.upper[a], shape[a])[ijk[:, a]] for a in range(3)])
    h = lengths/(shape-1); V = np.prod(h)*np.prod(np.where((ijk == 0) | (ijk == shape-1), .5, 1.), axis=1)
    pairs = cKDTree(X).query_pairs(delta, output_type='ndarray')
    pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
    return geometry, Family.from_arrays(X, V, pairs, delta)


def influence_support(geometry, family):
    """BAQP dependency envelope, NOT a chosen physical BC thickness."""
    seed = geometry.distance(family.X) < family.horizon
    mask = seed.copy(); i, j = family.pairs.T
    for _ in range(2):
        old = mask.copy(); np.logical_or.at(mask, i, old[j]); np.logical_or.at(mask, j, old[i])
    return mask, seed
