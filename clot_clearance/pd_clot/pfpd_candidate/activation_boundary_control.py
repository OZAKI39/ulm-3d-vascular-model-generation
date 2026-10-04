"""Isolated component-wise clamp control; the baseline model is not changed.

Geometry, BAQP, material, phase law and PETSc options come from the existing
candidate. This adapter changes only prescribed displacement components.
"""
import ast
import resource
import time
from functools import lru_cache
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from scripts import c7_early_branch as historical
from scripts import c7_reduced_vi as reduced
from scripts import petsc_vi_audit as prior
from . import petsc_vi_candidate as vi

FULL = 'FULL_CLAMP_BASELINE'
AXIAL = 'AXIAL_ONLY_CANDIDATE'


class BoundaryProblem(historical.SavedProblem):
    def __init__(self, n, strain, mode=FULL, lower=None):
        super().__init__(n, strain, lower)
        if mode not in (FULL, AXIAL):
            raise ValueError(mode)
        self.mode = mode
        X = self.family.X
        self.bottom = X[:, 2] == X[:, 2].min()
        fixed = np.zeros_like(X, dtype=bool)
        fixed[self.bottom | self.top, 2] = True
        self.gauge = []
        if mode == FULL:
            fixed[self.bottom | self.top, :] = True
        else:
            # Both loaded planes remove axial translation and rotations about
            # x/y. Remove the two lateral translations and remaining z rotation.
            bottom = np.flatnonzero(self.bottom)
            center = bottom[np.argmin(np.linalg.norm(X[bottom, :2], axis=1))]
            assert np.linalg.norm(X[center, :2]) < 1e-15
            candidates = bottom[np.isclose(X[bottom, 1], 0., atol=1e-15)]
            reference = candidates[np.argmax(X[candidates, 0])]
            fixed[center, :2] = True
            fixed[reference, 1] = True
            self.gauge = [(int(center), 0), (int(center), 1), (int(reference), 1)]
        self.fixed_dof = fixed
        # Boolean component indexing remains compatible with the existing
        # force and raw residual formulas; only pack/unpack need an adapter.
        self.free = ~fixed
        self.nx = int(self.free.sum())

    def prescribe(self, x):
        x = np.array(x, dtype=float, copy=True)
        x[self.fixed_dof] = self.boundary[self.fixed_dof]
        return x

    def unpack(self, z):
        z = np.asarray(z, float)
        if z.shape != (self.nx + self.ns,) or not np.isfinite(z).all():
            raise ValueError('Invalid joint state')
        x = self.boundary.copy()
        x[self.free] = self.family.X[self.free] + self.length_scale*z[:self.nx]
        return x, z[self.nx:].copy()

    def affine(self):
        x = self.family.X.copy()
        x[:, 2] += self.strain*(x[:, 2] + historical.old.v2.L/2)
        return x

    def gauge_rank(self):
        X = self.family.X/historical.old.v2.L
        R = np.zeros((len(X), 3, 6))
        R[:, :, :3] = np.eye(3)
        for k in range(3):
            R[:, :, 3+k] = np.cross(np.eye(3)[k], X)
        return int(np.linalg.matrix_rank(R[self.fixed_dof]))


@lru_cache(maxsize=8)
def _worker_problem(n, strain, mode):
    return BoundaryProblem(n, strain, mode)


def column_chunk(task):
    n, strain, mode, lower, ids, z, step, columns = task
    full = _worker_problem(n, strain, mode)
    # This task's mechanical solves keep s exactly zero. Phase support is
    # deliberately absent until a separately justified early-branch stage.
    assert not np.any(lower)
    p = reduced.ReducedProblem(full, ids)
    nx, ns = p.nx, p.ns
    class Columns:
        def __init__(self, rows):
            self.data = np.empty((rows, len(columns)))
            self.lookup = {j:i for i,j in enumerate(columns)}
        def __setitem__(self, key, value):
            self.data[:, self.lookup[key[1]]] = value
    obj = SimpleNamespace(M=Columns(nx), Yx=Columns(ns))
    _, loop = vi._source_nodes()
    loop.iter = ast.Name(id='columns', ctx=ast.Load())
    scope = dict(vars(vi.verified))
    scope.update(self=obj, p=p, z=z, step=step, nx=nx, ns=ns,
                 floor=p.config['mechanical_equilibrium']['force_floor_N'], columns=columns)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[loop], type_ignores=[])),
                 '<unchanged force/driver FD columns with component BC>', 'exec'), scope)
    return columns, obj.M.data, obj.Yx.data, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024


def parallel_columns(obj, p, z, step, pool, workers):
    tasks = [(p.n, p.strain, p.full.mode, p.full.lower, p.ids, z, step, c.tolist())
             for c in np.array_split(np.arange(p.nx), workers) if len(c)]
    chunks = pool.map(column_chunk, tasks) if pool else map(column_chunk, tasks)
    obj.worker_peak_RSS = 0
    for ids, M, Yx, rss in chunks:
        obj.M[:, ids] = M
        obj.Yx[:, ids] = Yx
        obj.worker_peak_RSS = max(obj.worker_peak_RSS, rss)


class Adapter(vi.PETScAdapter):
    def jacobian(self, snes, v, J, P):
        start = time.perf_counter()
        p, z = self.p, self.gather(v)
        scope = dict(vars(vi.verified))
        scope.update(_parallel_columns=parallel_columns, pool=self.pool, workers=self.workers)
        exec(compile(vi.adapted_constructor(), '<verified raw residual blocks>', 'exec'), scope)
        block = SimpleNamespace()
        scope['build'](block, p, z)
        assert p.ns == 0, 'This driver is restricted to intact mechanics.'
        # Retain structural zeros exactly as the historical adapter does. A
        # later Newton iterate may fill a zero entry in this mechanical block.
        nx = p.nx
        matrix = sparse.csr_matrix((block.M.ravel(), np.tile(np.arange(nx), nx),
                                    np.arange(nx+1)*nx), shape=(nx, nx))
        J.setValuesCSR(matrix.indptr.astype(self.PETSc.IntType), matrix.indices.astype(self.PETSc.IntType), matrix.data)
        J.assemble()
        self.jacobian_seconds.append(time.perf_counter()-start)
        self.worker_peak_RSS = max(getattr(self, 'worker_peak_RSS', 0), block.worker_peak_RSS)


def solve_intact(p, x, pool, workers):
    from petsc4py import PETSc
    rp = reduced.ReducedProblem(p, [])
    adapter = Adapter(rp, PETSc.COMM_SELF, pool, workers)
    adapter.assign(adapter.v, rp.compress(p.pack(p.prescribe(x), p.lower)))
    snes = PETSc.SNES().create(PETSc.COMM_SELF)
    cfg = prior.config()['primary']
    snes.setType(cfg['type'])
    snes.setFunction(adapter.function, adapter.f)
    snes.setJacobian(adapter.jacobian, adapter.J)
    snes.setVariableBounds(adapter.lo, adapter.hi)
    snes.setTolerances(atol=cfg['atol'], rtol=cfg['rtol'], stol=cfg['stol'], max_it=cfg['max_it'])
    snes.getLineSearch().setType(cfg['line_search'])
    ksp = snes.getKSP(); ksp.setType(cfg['ksp'])
    pc = ksp.getPC(); pc.setType(cfg['pc']); pc.setFactorSolverType(cfg['factor'])
    history = []
    def monitor(sn, it, norm):
        history.append(dict(iteration=int(it), scaled_force_norm=float(norm)))
        print('EQUILIBRIUM', p.n, p.mode, p.strain, int(it), float(norm), flush=True)
    snes.setMonitor(monitor)
    start = time.perf_counter()
    try:
        snes.solve(None, adapter.v)
        z = rp.expand(adapter.gather(adapter.v))
        x, s = p.unpack(z)
        info = dict(reason=int(snes.getConvergedReason()), iterations=int(snes.getIterationNumber()),
                    wall_seconds=time.perf_counter()-start, settings=cfg, history=history,
                    jacobian_seconds=adapter.jacobian_seconds,
                    worker_peak_RSS_bytes=getattr(adapter, 'worker_peak_RSS', 0))
    finally:
        snes.destroy(); adapter.destroy()
    return x, s, info
