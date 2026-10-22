"""Free-space monodisperse no-slip translational RPY, including overlap branch.

No wall images, lubrication, stresslets, Brownian forcing or bubble-slip model.
Near-contact branch regularizes pair mobility, not the separate acoustic force.
Wajnryb et al., JFM 731 R3 (2013), DOI 10.1017/jfm.2013.402.
"""

import numpy as np
import warp as wp


def matrix_numpy(x, radius, viscosity):
    x = np.asarray(x, float)
    if radius <= 0 or viscosity <= 0 or not np.isfinite(x).all():
        raise ValueError("Positive radius/viscosity and finite positions required")
    n = len(x)
    mobility = np.zeros((n, 3, n, 3))
    self_m = 1 / (6 * np.pi * viscosity * radius)
    for i in range(n):
        mobility[i, :, i, :] = self_m * np.eye(3)
        for j in range(i):
            dx = x[i] - x[j]
            d = np.linalg.norm(dx)
            if d == 0:
                raise ValueError("Coincident centers: strictly positive definiteness unavailable")
            nn = np.outer(dx, dx) / d**2
            if d >= 2 * radius:
                block = (
                    (1 + 2 * radius**2 / (3 * d * d)) * np.eye(3) + (1 - 2 * radius**2 / (d * d)) * nn
                ) / (8 * np.pi * viscosity * d)
            else:
                block = self_m * ((1 - 9 * d / (32 * radius)) * np.eye(3) + 3 * d / (32 * radius) * nn)
            mobility[i, :, j, :] = block
            mobility[j, :, i, :] = block.T
    return mobility.reshape(3 * n, 3 * n)


@wp.kernel
def velocity(
    x: wp.array(dtype=wp.vec3d),
    force: wp.array(dtype=wp.vec3d),
    flow: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=int),
    out: wp.array(dtype=wp.vec3d),
    n: int,
    radius: wp.float64,
    viscosity: wp.float64,
    flow_enabled: int,
    dt: wp.float64,
    free_path: wp.array(dtype=wp.float64),
    error: wp.array(dtype=int),
):
    i = wp.tid()
    u = wp.vec3d(wp.float64(0.0))
    mu0 = wp.float64(1.0) / (wp.float64(18.84955592153876) * viscosity * radius)
    if active[i] == 1:
        u = mu0 * force[i]
        start = (i // n) * n
        for j in range(start, start + n):
            if j != i and active[j] == 1:
                dx = x[i] - x[j]
                d = wp.length(dx)
                if d <= wp.float64(0.0):
                    wp.atomic_or(error, 0, 4)
                else:
                    normal = dx / d
                    aa = mu0 * (wp.float64(1.0) - wp.float64(9.0) * d / (wp.float64(32.0) * radius))
                    bb = mu0 * wp.float64(3.0) * d / (wp.float64(32.0) * radius)
                    if d >= wp.float64(2.0) * radius:
                        coeff = wp.float64(1.0) / (wp.float64(25.132741228718345) * viscosity * d)
                        aa = coeff * (
                            wp.float64(1.0) + wp.float64(2.0) * radius * radius / (wp.float64(3.0) * d * d)
                        )
                        bb = coeff * (wp.float64(1.0) - wp.float64(2.0) * radius * radius / (d * d))
                    u = u + aa * force[j] + bb * wp.dot(normal, force[j]) * normal
        if flow_enabled == 1:
            u = u + flow[i]
    out[i] = u
    free_path[i] = free_path[i] + wp.length(u) * dt


@wp.kernel
def copy_velocity(source: wp.array(dtype=wp.vec3d), target: wp.array(dtype=wp.vec3d)):
    i = wp.tid()
    target[i] = source[i]
