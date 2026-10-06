"""Reuse the existing SciPy mechanics solver with a matrix-free BAQP tangent.

No new nonlinear iteration is defined. The historical _scipy function code is
executed unchanged; only its final dense condition-number diagnostic is disabled
for a LinearOperator. The default C1 residual calls unchanged baqp.evaluate;
an explicit PD2 backend delegates to the accepted variational force candidate.
"""
import time
import types
from types import SimpleNamespace
import numpy as np
from scipy.sparse.linalg import LinearOperator
from scipy.sparse import csr_matrix
from ...regularization import equilibrium as legacy
from ..baqp import evaluate


class IntactForce:
    def __init__(self, family, material, external, wall_limit=1800, assembly='sparse'):
        self.family, self.material = family, material
        self.external = np.asarray(external); self.x = None; self.state = None
        self.calls = 0; self.tangent_calls = 0; self.started = time.monotonic(); self.wall_limit = wall_limit
        self.assembly = assembly; self.incidence = None

    def scatter(self, values, endpoint):
        """Same ordered bond-to-node sum; CSR only speeds the aggregation."""
        f = self.family; n = len(f.X)
        if self.assembly == 'numpy':
            result = np.zeros((n, *values.shape[1:]))
            np.add.at(result, f.pairs[:, endpoint], values)
            return result
        if self.assembly != 'sparse': raise ValueError('Unknown tangent assembly')
        if self.incidence is None:
            m = len(f.pairs); columns = np.arange(m)
            self.incidence = [csr_matrix((np.ones(m), (f.pairs[:, k], columns)), shape=(n, m)) for k in (0, 1)]
        return (self.incidence[endpoint]@values.reshape(len(values), -1)).reshape((n, *values.shape[1:]))

    def __call__(self, x):
        if time.monotonic()-self.started > self.wall_limit:
            raise TimeoutError('Predeclared per-case mechanics wall limit')
        if self.x is None or not np.array_equal(x, self.x):
            self.state = evaluate(self.family, x, self.material)
            self.x = x.copy(); self.calls += 1
        s = self.state; f = self.family; i, j = f.pairs.T
        strain = np.maximum(np.linalg.norm(x[j]-x[i], axis=1)/f.length-1, 0.)
        return legacy.ForceState(s.force, self.external, strain, s.energy,
                                 float(np.linalg.det(s.Ftilde).min()), np.zeros_like(x),
                                 np.zeros(len(x), dtype=bool))

    def tangent(self, v, state=None):
        """Exact directional derivative of the frozen-topology intact SVK law."""
        self.tangent_calls += 1
        f, st = self.family, self.state if state is None else state; i, j = f.pairs.T; V = f.volume
        dv = v[j]-v[i]
        dF = self.scatter(dv[:, :, None]*st.shape_i[:, None, :]*V[j, None, None], 0)
        dF -= self.scatter(dv[:, :, None]*st.shape_j[:, None, :]*V[i, None, None], 1)
        av = .5*(dF[i]+dF[j]); r2 = f.length**2
        residual = dv-np.einsum('bij,bj->bi', av, f.xi)
        dFt = av+residual[:, :, None]*f.xi[:, None, :]/r2[:, None, None]
        F = st.Ftilde
        E = .5*(F.swapaxes(1, 2)@F-np.eye(3))
        S = self.material.lame*np.trace(E, axis1=1, axis2=2)[:, None, None]*np.eye(3)+2*self.material.mu*E
        dE = .5*(dFt.swapaxes(1, 2)@F+F.swapaxes(1, 2)@dFt)
        dS = self.material.lame*np.trace(dE, axis1=1, axis2=2)[:, None, None]*np.eye(3)+2*self.material.mu*dE
        dP = dFt@S+F@dS
        del dFt, dE, dS, E, S, av
        projection = np.eye(3)-f.xi[:, :, None]*f.xi[:, None, :]/r2[:, None, None]
        term = f.wbar[:, None, None]*(dP@projection)
        dZ = self.scatter(term*V[j, None, None], 0)+self.scatter(term*V[i, None, None], 1)
        direct = f.wbar[:, None]*np.einsum('bij,bj->bi', dP, f.xi)/r2[:, None]
        ti = direct+np.einsum('bij,bj->bi', dZ[i], st.shape_i)
        tj = -direct+np.einsum('bij,bj->bi', dZ[j], st.shape_j)
        pair = (ti-tj)*V[i, None]*V[j, None]
        return self.scatter(pair, 0)-self.scatter(pair, 1)


class _LinalgReportProxy:
    def __getattr__(self, key): return getattr(np.linalg, key)
    def cond(self, matrix):
        # _scipy casts this diagnostic to float; the summary converts the
        # unavailable value to JSON null after the unchanged solver returns.
        return float('nan') if isinstance(matrix, LinearOperator) else np.linalg.cond(matrix)


class _NumpyReportProxy:
    linalg = _LinalgReportProxy()
    def __getattr__(self, key): return getattr(np, key)


def solve(family, material, prescribed, fixed, initial, external, h, force_scale, config, *, mechanics=None):
    free = ~fixed
    view = SimpleNamespace(X=np.where(fixed[:, None], prescribed, family.X), fixed=fixed,
                           spacing=h, pairs=family.pairs, volume=family.volume)
    if mechanics is None:
        callback = IntactForce(family, material, external, config['solver_max_wall_seconds_per_case'])
    else:
        from .mechanics_backend import PD2SolverCallback
        callback = PD2SolverCallback(mechanics, external, config['solver_max_wall_seconds_per_case'])
    x0 = initial.copy(); x0[fixed] = prescribed[fixed]
    origin = view.X[free]
    scale = max(np.linalg.norm(external[free]), force_scale)
    def jacobian(q):
        x = view.X.copy(); x[free] = origin+h*q.reshape((-1, 3)); callback(x)
        count = 3*int(free.sum()); linearization_state = callback.state
        def action(v):
            full = np.zeros_like(x); full[free] = h*np.asarray(v).reshape((-1, 3))
            return callback.tangent(full, linearization_state)[free].ravel()/scale
        # The internal force Jacobian is the negative symmetric energy Hessian.
        return LinearOperator((count, count), matvec=action, rmatvec=action, dtype=float)
    ec = dict(method='scipy_reference', maximum_iterations=config['solver_max_nfev'],
              relative_force_residual_tolerance=config['mechanical_residual_gate'],
              force_floor_N=force_scale, driver_drift_tolerance=.01, driver_absolute_floor=1e-12,
              kinetic_to_stored_energy_tolerance=1e-6, energy_floor_J=1e-24,
              minimum_consecutive_converged_checks=1, check_interval=1,
              scipy_reference_max_free_dofs=3*len(family.X),
              scipy=dict(method='trf', jac=jacobian, xtol=1e-12, ftol=1e-12, gtol=1e-12,
                         max_nfev=config['solver_max_nfev'], x_scale=1., loss='linear',
                         tr_solver='lsmr', tr_options=dict(atol=1e-10, btol=1e-10,
                                                          maxiter=config['solver_lsmr_maxiter'])))
    c = dict(mechanical_equilibrium=ec, safety=dict(maximum_displacement_m=.01, maximum_radius_ratio=1.),
             pipe=dict(radius_m=1.), clot=dict(density_kg_m3=1000.))
    legacy.validate_options(view, x0, np.zeros_like(x0), np.ones(len(family.pairs)), c, 'scipy_reference', True)
    # The sole reporting adapter avoids attempting dense SVD of a matrix-free
    # Jacobian. Solver code, tolerances, residuals and acceptance remain visible.
    scope = dict(legacy._scipy.__globals__); scope['np'] = _NumpyReportProxy()
    solver = types.FunctionType(legacy._scipy.__code__, scope, name=legacy._scipy.__name__)
    started = time.monotonic()
    try:
        result = solver(view, x0, callback, c)
        x = result.x; info = result.summary
    except TimeoutError as exc:
        x = callback.x.copy()
        callback.wall_limit = float('inf')  # Read cached final diagnostics only.
        info = dict(converged=False, status='PREDECLARED_WALL_LIMIT', reason=str(exc),
                    final=legacy.diagnostics(view, x, callback(x), c))
    info.update(wall_seconds=time.monotonic()-started, BAQP_evaluations=callback.calls,
                tangent_actions=callback.tangent_calls, final_jacobian_condition=None,
                condition_report='NOT_ESTIMATED_MATRIX_FREE', no_new_solver=True,
                tangent_assembly=callback.assembly)
    print('SOLVE', info['status'], info['BAQP_evaluations'], info['tangent_actions'], info['wall_seconds'], flush=True)
    return x, info
