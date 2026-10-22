"""Finite-radius straight tube and apparatus; no CFD or triangulated Y geometry."""

from dataclasses import dataclass
import numpy as np


@dataclass
class StraightTube:
    radius_m: float
    length_m: float
    thickness_m: float
    center: np.ndarray
    rotation: np.ndarray
    ends: tuple = ("open_fail_stop", "open_fail_stop")

    @classmethod
    def from_config(cls, g):
        tube = cls(
            g["inner_diameter_m"] / 2,
            g["length_m"],
            g["wall_thickness_m"],
            np.asarray(g["center_m"], float),
            np.asarray(g["rotation_local_to_world"], float),
            tuple(g["end_conditions"]),
        )
        if not (tube.radius_m > 0 and tube.length_m > 0 and tube.thickness_m > 0):
            raise ValueError("Demo geometry requires explicit positive length, lumen radius and thickness")
        if not np.allclose(tube.rotation.T @ tube.rotation, np.eye(3), atol=1e-12) or not np.isclose(
            np.linalg.det(tube.rotation), 1
        ):
            raise ValueError("A proper orthonormal local-to-world rotation is required")
        if len(tube.ends) != 2 or any(e not in ("open_fail_stop", "closed") for e in tube.ends):
            raise ValueError("End conditions must be closed or open_fail_stop (no reflection)")
        return tube

    @property
    def outer_radius_m(self):
        return self.radius_m + self.thickness_m

    def world(self, x):
        return np.asarray(x) @ self.rotation.T + self.center

    def local(self, x):
        return (np.asarray(x) - self.center) @ self.rotation

    def gaps(self, x, radius=0.0):
        q = self.local(x)
        return np.stack(
            [
                (self.radius_m - radius) - np.linalg.norm(q[..., 1:], axis=-1),
                q[..., 0] + self.length_m / 2 - (radius if self.ends[0] == "closed" else 0),
                self.length_m / 2 - q[..., 0] - (radius if self.ends[1] == "closed" else 0),
            ],
            axis=-1,
        )

    def contains(self, x, radius=0.0):
        return np.isfinite(x).all(axis=-1) & (self.gaps(x, radius) >= 0).all(axis=-1)

    def surface(self, radius=None, n=129):
        angle = np.linspace(0, 2 * np.pi, n)
        r = self.outer_radius_m if radius is None else radius
        x, a = np.meshgrid([-self.length_m / 2, self.length_m / 2], angle)
        return self.world(np.stack([x, r * np.cos(a), r * np.sin(a)], axis=-1))


def disc_points(element, n=257, rings=1):
    c, normal = np.array(element["position_m"]), np.array(element["normal"])
    normal = normal / np.linalg.norm(normal)
    u = np.cross(normal, np.eye(3)[np.argmin(abs(normal))])
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    a = np.linspace(0, 2 * np.pi, n)
    rr = np.linspace(0, 1, rings + 1)[1:] * element["size_m"][0] / 2
    return c + (rr[:, None, None] * (np.cos(a)[:, None] * u + np.sin(a)[:, None] * v)).reshape(-1, 3)


def apparatus(g, frequency):
    tube = StraightTube.from_config(g)
    gap, sign = g["gap_m"], g["end_sign"]
    if sign not in (-1, 1) or gap <= 0 or g["pzt_diameter_m"] != 0.009:
        raise ValueError("F10L requires a positive gap, explicit end sign, and 9 mm faces")
    d = tube.outer_radius_m + gap
    centers = [[0, d, 0], [0, -d, 0], [sign * (tube.length_m / 2 + gap), 0, 0]]
    normals = [[0, -1, 0], [0, 1, 0], [-sign, 0, 0]]
    elements = []
    for i, name in enumerate(["side_y_plus", "side_y_minus", "end_x"]):
        c, n = tube.world(centers[i]), np.array(normals[i]) @ tube.rotation.T
        if g.get("pzt_overrides") is not None:
            c, n = [np.asarray(g["pzt_overrides"][i][key], float) for key in ("position_m", "normal")]
        if not np.isfinite(c).all() or not np.isfinite(n).all() or np.linalg.norm(n) == 0:
            raise ValueError("Invalid PZT coordinate or normal")
        n = n / np.linalg.norm(n)
        elements.append(
            dict(
                id=name,
                position_m=c.tolist(),
                normal=n.tolist(),
                geometry="disc",
                size_m=[g["pzt_diameter_m"]],
                amplitude_Pa=0.0,
                phase_rad=0.0,
                enabled=True,
            )
        )
    return tube, elements


def geometry_audit(tube, elements, expected_gap):
    """Exact convex distances via constrained minimization, not display meshes."""
    from scipy.optimize import minimize

    centers = np.array([e["position_m"] for e in elements])
    if len(elements) != 3:
        raise ValueError("Exactly three source faces required")
    local_centers = tube.local(centers)
    no_z_sources = bool(
        np.all(abs(local_centers[:2, 1]) > abs(local_centers[:2, 2]))
        and abs(local_centers[2, 0]) > np.linalg.norm(local_centers[2, 1:])
    )
    if not no_z_sources:
        raise ValueError("Expected two local side-Y sources and one axial source, no Z source")
    normals = np.array([e["normal"] for e in elements])
    radii = np.array([e["size_m"][0] / 2 for e in elements])
    rows, pairs = [], []
    for i, e in enumerate(elements):
        c, n, r = centers[i], normals[i], radii[i]

        def objective(v):
            return np.sum((v[:3] - v[3:]) ** 2)

        def tube_radial(v):
            q = tube.local(v[3:] * 1e-3)
            return (tube.outer_radius_m * 1e3) ** 2 - np.sum((q[1:] * 1e3) ** 2)

        # Optimizer uses mm to avoid tiny SI finite-difference steps.
        cons = [
            {"type": "eq", "fun": lambda v: np.dot(v[:3] - c * 1e3, n)},
            {"type": "ineq", "fun": lambda v: (r * 1e3) ** 2 - np.sum((v[:3] - c * 1e3) ** 2)},
            {"type": "ineq", "fun": tube_radial},
            {"type": "ineq", "fun": lambda v: tube.length_m * 500 - abs(tube.local(v[3:] * 1e-3)[0] * 1e3)},
        ]
        sol = minimize(
            objective,
            np.r_[c, tube.center] * 1e3,
            constraints=cons,
            method="SLSQP",
            options={"ftol": 1e-12, "maxiter": 500},
        )
        if not sol.success:
            raise ValueError("Tube/disc distance optimization failed: " + sol.message)
        distance = np.sqrt(sol.fun) * 1e-3
        target = tube.center if i < 2 else tube.world([np.sign(tube.local(c)[0]) * tube.length_m / 2, 0, 0])
        toward = np.dot(n, target - c) > 0
        if distance <= 0 or not toward:
            raise ValueError("PZT invades tube or faces away from lumen")
        rows.append(
            dict(
                **e,
                closest_outer_envelope_gap_m=float(distance),
                faces_lumen=bool(toward),
                gap_matches_registered=bool(np.isclose(distance, expected_gap, atol=1e-8)),
            )
        )
    for i in range(3):
        for j in range(i + 1, 3):
            cons = []
            for sl, index in [(slice(0, 3), i), (slice(3, 6), j)]:
                cons += [
                    {
                        "type": "eq",
                        "fun": lambda v, sl=sl, index=index: np.dot(
                            v[sl] - centers[index] * 1e3, normals[index]
                        ),
                    },
                    {
                        "type": "ineq",
                        "fun": lambda v, sl=sl, index=index: (radii[index] * 1e3) ** 2
                        - np.sum((v[sl] - centers[index] * 1e3) ** 2),
                    },
                ]
            sol = minimize(
                lambda v: np.sum((v[:3] - v[3:]) ** 2),
                np.r_[centers[i], centers[j]] * 1e3,
                constraints=cons,
                method="SLSQP",
                options={"ftol": 1e-12, "maxiter": 500},
            )
            if not sol.success or sol.fun < 1e-12:
                raise ValueError("Intersecting source faces or unresolved distance")
            pairs.append(
                dict(ids=[elements[i]["id"], elements[j]["id"]], minimum_gap_m=float(np.sqrt(sol.fun) * 1e-3))
            )
    return dict(
        status="PASS",
        pzt=rows,
        pair_gaps=pairs,
        source_count=3,
        no_local_z_sources=no_z_sources,
        injection_assembly_interference="UNKNOWN: apparatus not supplied",
    )


class ZeroFlow:
    mode = "zero"
    metadata = dict(mode="zero", cfd_files=[], analytic=True)

    def query(self, positions, time_s=0.0):
        x = np.asarray(positions)
        if not np.isfinite(x).all():
            raise ValueError("Finite coordinates required")
        return np.zeros_like(x, dtype=float), np.zeros(x.shape[:-1] + (3, 3))
