"""Resident strict tetrahedral P1 velocity/gradient interpolation in SI.

The float32 BVH is a conservative broad phase only. Membership and interpolation
use original coordinates promoted to float64, with the CPU sampler's 1e-10
barycentric roundoff tolerance. No nearest-neighbour extrapolation is allowed.
"""

import hashlib
from pathlib import Path
from time import perf_counter
import numpy as np
import warp as wp

CFD_SHA256 = "131f7cb9b3e4ac9887bb0aa344d0eb2c0a2c8ac1d3b19fafb093ecd9a1045620"
FLOW_DOMAIN_ERROR = 16384


@wp.func
def barycentric(x: wp.vec3d, origin: wp.vec3d, inverse: wp.mat33d):
    q = inverse * (x - origin)
    return wp.vec4d(wp.float64(1.0) - q[0] - q[1] - q[2], q[0], q[1], q[2])


@wp.func
def contains(w: wp.vec4d):
    valid = True
    for j in range(4):
        if w[j] < wp.float64(-1e-10) or w[j] > wp.float64(1.0 + 1e-10):
            valid = False
    return valid


@wp.kernel
def tetra_query(
    bvh: wp.uint64,
    origins: wp.array(dtype=wp.vec3d),
    inverses: wp.array(dtype=wp.mat33d),
    values: wp.array2d(dtype=wp.vec3d),
    gradients: wp.array(dtype=wp.mat33d),
    center: wp.vec3d,
    padding: wp.float64,
    x: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=int),
    pending: wp.array(dtype=int),
    allow_terminal: int,
    cache: wp.array(dtype=int),
    velocity: wp.array(dtype=wp.vec3d),
    gradient: wp.array(dtype=wp.mat33d),
    valid: wp.array(dtype=int),
    stats: wp.array2d(dtype=int),
    error: wp.array(dtype=int),
    failed_x: wp.array(dtype=wp.vec3d),
):
    i = wp.tid()
    velocity[i] = wp.vec3d(wp.float64(0.0))
    gradient[i] = wp.mat33d(wp.float64(0.0))
    valid[i] = 0
    if active[i] == 0 or (allow_terminal == 1 and pending[i] >= 2 and pending[i] <= 4):
        return
    if error[0] != 0:
        return
    stats[i, 0] = stats[i, 0] + 1
    found = int(-1)
    w = wp.vec4d(wp.float64(0.0))
    previous = cache[i]
    if previous >= 0:
        w = barycentric(x[i], origins[previous], inverses[previous])
        if contains(w):
            found = previous
            stats[i, 1] = stats[i, 1] + 1
    if found < 0:
        p = x[i] - center
        low = p - wp.vec3d(padding)
        high = p + wp.vec3d(padding)
        query = wp.bvh_query_aabb(
            bvh,
            wp.vec3(float(low[0]), float(low[1]), float(low[2])),
            wp.vec3(float(high[0]), float(high[1]), float(high[2])),
        )
        cell = int(0)
        while wp.bvh_query_next(query, cell):
            stats[i, 2] = stats[i, 2] + 1
            weights = barycentric(x[i], origins[cell], inverses[cell])
            if contains(weights) and (found < 0 or cell < found):
                found = cell
                w = weights
    cache[i] = found
    if found < 0:
        stats[i, 3] = stats[i, 3] + 1
        failed_x[i] = x[i]
        wp.atomic_or(error, 0, 16384)
        return
    u = wp.vec3d(wp.float64(0.0))
    for a in range(4):
        u = u + w[a] * values[found, a]
    velocity[i] = u
    gradient[i] = gradients[found]
    valid[i] = 1


@wp.kernel
def affine_query(
    x: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=int),
    u0: wp.vec3d,
    center: wp.vec3d,
    tensor: wp.mat33d,
    velocity: wp.array(dtype=wp.vec3d),
    gradient: wp.array(dtype=wp.mat33d),
    valid: wp.array(dtype=int),
):
    i = wp.tid()
    velocity[i] = wp.vec3d(wp.float64(0.0))
    gradient[i] = tensor
    valid[i] = active[i]
    if active[i] == 1:
        velocity[i] = u0 + tensor * (x[i] - center)


class QueryBuffers:
    def __init__(self, n, device):
        self.cache = wp.full(n, -1, dtype=wp.int32, device=device)
        self.velocity = wp.zeros(n, dtype=wp.vec3d, device=device)
        self.gradient = wp.zeros(n, dtype=wp.mat33d, device=device)
        self.valid = wp.zeros(n, dtype=wp.int32, device=device)
        self.stats = wp.zeros((n, 4), dtype=wp.int32, device=device)
        self.failed_x = wp.zeros(n, dtype=wp.vec3d, device=device)

    def reset(self):
        self.cache.fill_(-1)
        self.stats.zero_()
        self.valid.zero_()
        self.failed_x.zero_()


class GPUFlow:
    """Load the original tetrahedra once; queries/caches remain on the GPU."""

    def __init__(self, path, device="cuda:0"):
        import pyvista as pv

        begin = perf_counter()
        self.path, self.device = Path(path), device
        digest = hashlib.sha256(self.path.read_bytes()).hexdigest()
        if digest != CFD_SHA256:
            raise ValueError("CFD snapshot identity mismatch")
        grid = pv.read(path)
        if set(grid.cells_dict) != {pv.CellType.TETRA}:
            raise ValueError("Only original linear tetrahedra are supported")
        cells = np.asarray(grid.cells_dict[pv.CellType.TETRA], np.int32)
        points = np.asarray(grid.points, np.float64)
        values = np.asarray(grid.point_data["Velocity"], np.float64)
        if not np.isfinite(points).all() or not np.isfinite(values).all() or values.shape != points.shape:
            raise ValueError("Finite SI positions and node velocities required")
        tet = points[cells]
        inverses = np.linalg.inv(np.swapaxes(tet[:, 1:] - tet[:, :1], 1, 2))
        grad = np.swapaxes(values[cells[:, 1:]] - values[cells[:, :1]], 1, 2) @ inverses
        center = points.mean(0)
        quant = np.abs((points - center).astype(np.float32).astype(float) - (points - center)).max()
        padding = max(2e-8, float(quant) * 16)
        lo = (tet.min(1) - center - padding).astype(np.float32)
        hi = (tet.max(1) - center + padding).astype(np.float32)
        self.lower = wp.array(lo, dtype=wp.vec3, device=device)
        self.upper = wp.array(hi, dtype=wp.vec3, device=device)
        self.bvh = wp.Bvh(self.lower, self.upper)
        self.origins = wp.array(tet[:, 0], dtype=wp.vec3d, device=device)
        self.inverses = wp.array(inverses, dtype=wp.mat33d, device=device)
        self.values = wp.array(values[cells], dtype=wp.vec3d, device=device)
        self.gradients = wp.array(grad, dtype=wp.mat33d, device=device)
        self.center = wp.vec3d(*center)
        self.padding = padding
        wp.load_module(module=tetra_query.module, device=device)
        wp.synchronize_device(device)
        self.metadata = dict(
            source=str(path),
            sha256=digest,
            array="Velocity",
            units="m, m/s, 1/s",
            node_count=len(points),
            tetra_count=len(cells),
            coordinate_dtype=str(grid.points.dtype),
            velocity_max_m_s=float(np.linalg.norm(values, axis=1).max()),
            bounds_m=np.column_stack([points.min(0), points.max(0)]).tolist(),
            barycentric_roundoff=1e-10,
            broadphase_padding_m=padding,
            interpolation="original tetrahedral P1; float64 membership/values; padded float32 BVH",
            static_build_s=perf_counter() - begin,
            static_BVH=True,
            gradient="piecewise-constant exact derivative of nodal P1 velocity; discontinuous at cell faces",
            outside="explicit error bit 16384; no nearest-neighbour or zero-velocity continuation",
        )

    def launch(self, x, active, pending, buffers, error, allow_terminal=False):
        wp.launch(
            tetra_query,
            dim=len(x),
            inputs=[
                self.bvh.id,
                self.origins,
                self.inverses,
                self.values,
                self.gradients,
                self.center,
                self.padding,
                x,
                active,
                pending,
                int(allow_terminal),
                buffers.cache,
                buffers.velocity,
                buffers.gradient,
                buffers.valid,
                buffers.stats,
                error,
                buffers.failed_x,
            ],
            device=self.device,
            block_dim=32,
        )

    def sample(self, x):
        """Validation/save-point adapter only. Never called from physical steps."""
        x = np.atleast_2d(x)
        gpu = wp.array(x, dtype=wp.vec3d, device=self.device)
        active = wp.ones(len(x), dtype=int, device=self.device)
        pending = wp.zeros_like(active)
        error = wp.zeros(1, dtype=int, device=self.device)
        b = QueryBuffers(len(x), self.device)
        self.launch(gpu, active, pending, b, error)
        u = b.velocity.numpy()
        valid = b.valid.numpy().astype(bool)
        u[~valid] = np.nan
        return dict(
            velocity_m_s=u,
            gradient_s_inv=b.gradient.numpy(),
            valid=valid,
            tetra_id=b.cache.numpy(),
            error=int(error.numpy()[0]),
            stats=b.stats.numpy(),
        )


class AffineFlow:
    """Mathematical check only; never labelled Y-vessel CFD."""

    def __init__(self, velocity, gradient=None, center=None, device="cuda:0"):
        self.device = device
        self.u0 = wp.vec3d(*velocity)
        self.tensor = wp.mat33d(np.zeros((3, 3)) if gradient is None else np.asarray(gradient))
        self.center = wp.vec3d(*(np.zeros(3) if center is None else center))
        self.metadata = dict(
            scope="SOFTWARE_TEST_ONLY_AFFINE",
            velocity_m_s=list(velocity),
            gradient_s_inv=np.asarray(self.tensor).reshape(3, 3).tolist(),
        )

    def launch(self, x, active, pending, buffers, error, allow_terminal=False):
        wp.launch(
            affine_query,
            dim=len(x),
            inputs=[
                x,
                active,
                self.u0,
                self.center,
                self.tensor,
                buffers.velocity,
                buffers.gradient,
                buffers.valid,
            ],
            device=self.device,
            block_dim=32,
        )
