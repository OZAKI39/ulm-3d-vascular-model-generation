"""Reuse M0 main-bore motion; reject unvalidated junction dynamics instead of false walls."""

import numpy as np
import warp as wp
from .straight_tube_engine import StraightTubeSimulator
from .engine import FastSimulator
from . import kernels as k


@wp.kernel
def reject_port_neighborhood(
    x: wp.array(dtype=wp.vec3d),
    center: wp.vec3d,
    rotation: wp.mat33d,
    port_x: wp.float64,
    half_width: wp.float64,
    error: wp.array(dtype=int),
):
    i = wp.tid()
    q = wp.transpose(rotation) * (x[i] - center)
    if wp.abs(q[0] - port_x) <= half_width:
        wp.atomic_or(error, 0, 65536)


class PDMSMainBoreSimulator(StraightTubeSimulator):
    def __init__(self, *args, top_port, **kwargs):
        self.top_port = top_port
        super().__init__(*args, **kwargs)
        x = np.asarray(args[0])
        q = self.tube.local(x)
        if np.any(
            abs(q[..., 0] - top_port["center_xy_m"][0]) <= top_port["diameter_m"] / 2 + self.p.radius_m
        ):
            raise ValueError(
                "Port junction particle dynamics not admitted; do not apply a false cylindrical wall"
            )
        wp.load_module(module=reject_port_neighborhood.module, device=self.device)

    def launch(self, kernel, dim, args):
        result = super().launch(kernel, dim, args)
        if self._straight_ready and kernel in (k.predict, k.contact):
            buffer = args[2] if kernel is k.predict else args[1]
            t = self.tube
            FastSimulator.launch(
                self,
                reject_port_neighborhood,
                self.total,
                [
                    buffer,
                    wp.vec3d(*t.center),
                    wp.mat33d(*t.rotation.ravel()),
                    self.top_port["center_xy_m"][0],
                    self.top_port["diameter_m"] / 2 + self.p.radius_m,
                    self.error,
                ],
            )
        return result
