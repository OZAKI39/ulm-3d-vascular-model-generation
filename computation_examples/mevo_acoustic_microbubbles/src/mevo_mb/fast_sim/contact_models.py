"""NumPy reference for symmetric local-mobility contact, with exact net balance."""

import numpy as np


def mobility_tensors(x, physics):
    p = physics

    def phi(e):
        return 0.5 / e - 0.45 * np.log(e) - 3 / 56 * e * np.log(e)

    out = []
    for i in range(len(x)):
        resistance = np.eye(3)
        for j in range(len(x)):
            if i == j:
                continue
            delta = x[j] - x[i]
            d = np.linalg.norm(delta)
            if d <= 0:
                raise ValueError("Coincident bubbles")
            e = (d - 2 * p.radius_m) / p.radius_m
            if e < p.lubrication_activation:
                n = delta / d
                c = 0.5 * (phi(max(e, p.lubrication_floor)) - phi(p.lubrication_activation))
                resistance += c * np.outer(n, n)
        out.append(np.linalg.inv(resistance))
    return np.asarray(out)


def contact_pass(x, mobility, physics):
    dx = np.zeros_like(x)
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            dvec = x[i] - x[j]
            d = np.linalg.norm(dvec)
            overlap = 2 * physics.radius_m - d
            if overlap > 0:
                n = dvec / d
                h = 0.5 * (mobility[i] + mobility[j])
                v = h @ n
                correction = physics.contact_relaxation * 0.5 * overlap * v / (n @ v)
                dx[i] += correction
                dx[j] -= correction
    return x + dx
