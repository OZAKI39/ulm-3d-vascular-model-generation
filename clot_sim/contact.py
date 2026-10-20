"""Physical N/J particle repulsion, excluding every surviving reference bond."""

import numpy as np
from .config import DX, K, HORIZON, CONFIG

CONTACT = dict(
    effective_radius_m=CONFIG["wall"]["effective_radius"],
    stiffness_N_m=CONFIG["wall"]["stiffness"],
    damping_ratio=CONFIG["wall"]["damping_ratio"],
)


class Contact:

    def __init__(self, ref):
        from scipy.spatial import cKDTree

        self.tree = cKDTree
        self.b = ref["bonds"]
        self.N = len(ref["X"])
        self.codes = self.b[:, 0] * self.N + self.b[:, 1]
        self.vol = ref["V"]
        self.c = 15 * (18 * K / (np.pi * HORIZON**4)) / HORIZON

    def coefficients(self, pairs):
        return self.c * ((self.vol[pairs[:, 0]] + self.vol[pairs[:, 1]]) / 2) ** 2

    def energy_for_pairs(self, x, pairs):
        return float(
            0.5
            * np.sum(
                self.coefficients(pairs)
                * np.minimum(
                    np.linalg.norm(x[pairs[:, 1]] - x[pairs[:, 0]], axis=1) - DX, 0
                )
                ** 2
            )
        )

    def evaluate(self, x, d):
        pairs = self.tree(x).query_pairs(DX, output_type="ndarray")
        pairs = pairs[~np.isin(pairs[:, 0] * self.N + pairs[:, 1], self.codes[d < 1])]
        f = np.zeros_like(x)
        if not len(pairs):
            return (0.0, f, pairs)
        delta = x[pairs[:, 1]] - x[pairs[:, 0]]
        r = np.linalg.norm(delta, axis=1)
        dr = r - DX
        c = self.coefficients(pairs)
        ff = c[:, None] * dr[:, None] * delta / r[:, None]
        for k in range(3):
            f[:, k] = np.bincount(
                pairs[:, 0], weights=ff[:, k], minlength=self.N
            ) - np.bincount(pairs[:, 1], weights=ff[:, k], minlength=self.N)
        return (float(0.5 * np.dot(c, dr * dr)), f, pairs)
