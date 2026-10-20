"""Independent endpoint-star variational LPS with nonuniform nodal quadrature.

A_i=V_i(K-5G/3), L_ib=3 V_j h_b/m_i,
D_ib=15G V_i V_j h_b/(r_b m_i), k_ib=9K V_i V_j h_b/(r_b m_i).
These reduce to the frozen V3 equations for equal nodal volumes.
"""

import numpy as np
from scipy.sparse import coo_matrix
from .config import K as DEFAULT_K, G as DEFAULT_G


class SofteningLPS:
    def __init__(self, ref, K=DEFAULT_K, G=DEFAULT_G):
        self.X = ref["X"].copy()
        self.b = ref["bonds"].copy()
        self.r = ref["r0"].copy()
        self.h = ref["taper"].copy()
        self.m = ref["m"].copy()
        self.V = ref["V"].copy()
        self.K = K
        self.G = G
        assert K > 5 * G / 3 and (self.V > 0).all() and (self.m > 0).all()
        self.n = len(self.X)
        self.nb = len(self.b)
        self.i, self.j = self.b.T
        self.center = np.r_[self.i, self.j]
        self.other = np.r_[self.j, self.i]
        self.edge = np.tile(np.arange(self.nb), 2)
        self.l = 3 * self.V[self.other] * self.h[self.edge] / self.m[self.center]
        self.Ddir = (
            15
            * G
            * self.V[self.center]
            * self.V[self.other]
            * self.h[self.edge]
            / (self.r[self.edge] * self.m[self.center])
        )
        self.krefdir = (
            9
            * K
            * self.V[self.center]
            * self.V[self.other]
            * self.h[self.edge]
            / (self.r[self.edge] * self.m[self.center])
        )
        self.A = self.V * (K - 5 * G / 3)
        self.D = self.Ddir[: self.nb] + self.Ddir[self.nb :]
        self.kref = self.krefdir[: self.nb] + self.krefdir[self.nb :]
        self.L = coo_matrix(
            (self.l, (self.center, self.edge)), shape=(self.n, self.nb)
        ).tocsr()

    def geometry(self, x):
        dy = x[self.j] - x[self.i]
        length = np.linalg.norm(dy, axis=1)
        assert (length > 0).all()
        return length - self.r, dy / length[:, None]

    def assemble(self, t, di):
        ef = t[:, None] * di
        f = np.zeros_like(self.X)
        for k in range(3):
            f[:, k] = np.bincount(
                self.i, weights=ef[:, k], minlength=self.n
            ) - np.bincount(self.j, weights=ef[:, k], minlength=self.n)
        return f

    def baseline(self, x, alive=None):
        alive = np.ones(self.nb) if alive is None else np.asarray(alive, float)
        e, di = self.geometry(x)
        theta = self.L @ (e * alive)
        t = alive * (self.D * e + np.asarray(self.L.T @ (self.A * theta)))
        U = 0.5 * np.dot(self.D * alive, e * e) + 0.5 * np.dot(self.A, theta * theta)
        return dict(U=float(U), force=self.assemble(t, di), theta=theta, t=t, e=e)

    def evaluate(self, x, d=None, mode="CANDIDATE_SOFTENING"):
        if mode == "BASELINE_INTACT_LPS":
            return self.baseline(x)
        if mode != "CANDIDATE_SOFTENING":
            raise ValueError(mode)
        d = np.zeros(self.nb) if d is None else np.asarray(d, float)
        assert ((d >= 0) & (d <= 1)).all()
        e, di = self.geometry(x)
        dd = d[self.edge]
        q = 1 - dd
        ee = e[self.edge]
        den = self.Ddir * dd + self.krefdir * q
        w = self.krefdir * q / den
        h = dd / den
        numerator = np.bincount(self.center, weights=self.l * w * ee, minlength=self.n)
        denominator = 1 / self.A + np.bincount(
            self.center, weights=self.l * self.l * h, minlength=self.n
        )
        p = numerator / denominator
        coupling = self.l * p[self.center]
        v = w * ee - h * coupling
        tdir = w * (self.Ddir * ee + coupling)
        t = tdir[: self.nb] + tdir[self.nb :]
        theta = np.bincount(self.center, weights=self.l * v, minlength=self.n)
        compliance = 0.5 * np.sum(
            self.krefdir * q * dd * ((self.Ddir * ee + coupling) / den) ** 2
        )
        local = 0.5 * np.dot(self.Ddir, v * v)
        volume = 0.5 * np.dot(self.A, theta * theta)
        U = local + volume + compliance
        Ydir = 0.5 * self.krefdir * ((self.Ddir * ee + coupling) / den) ** 2
        Y = Ydir[: self.nb] + Ydir[self.nb :]
        return dict(
            U=float(U),
            force=self.assemble(t, di),
            theta=theta,
            t=t,
            e=e,
            v=v,
            Y=Y,
            U_local=float(local),
            U_volume=float(volume),
            U_compliance=float(compliance),
            dual_energy_error_J=float(U - 0.5 * np.dot(e, t)),
            d=d,
        )
