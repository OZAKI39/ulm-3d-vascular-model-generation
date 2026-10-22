"""Configurable smooth Y lumen and first-crossing branch bookkeeping, in SI."""
from enum import StrEnum
import hashlib
import json

import numpy as np
from scipy.special import logsumexp


class Branch(StrEnum):
    PARENT = "PARENT"
    BRANCH_A = "BRANCH_A"
    BRANCH_B = "BRANCH_B"
    UNASSIGNED = "UNASSIGNED"
    WALL_LOSS = "WALL_LOSS"


class BifurcationGeometry:
    """Smooth union of three tubular branches, with flat physical terminal caps.

    Implicit geometry supports construction/software checks. Production CFD uses
    its actual meshed surface for geometric exclusion, not an analytic wall proxy.
    """
    def __init__(self, config):
        self.config = dict(config)
        self.parent_length = float(config["parent_length_m"])
        self.parent_radius = float(config["parent_radius_m"])
        self.daughter_length = float(config["daughter_length_m"])
        self.daughter_radius = float(config["daughter_radius_m"])
        self.blend = float(config["blending_length_m"])
        angle = np.deg2rad(config["branch_angle_deg"])
        if not 0 < angle < np.pi / 2 or min(self.parent_length, self.parent_radius,
                                           self.daughter_length, self.daughter_radius, self.blend) <= 0:
            raise ValueError("Positive dimensions and a branch angle between 0 and 90 degrees required")
        distance = float(config["classification_distance_m"])
        if not 3 * self.daughter_radius < distance < self.daughter_length:
            raise ValueError("Classification planes must be in separate daughters, upstream of physical outlets")
        self.directions = {
            Branch.BRANCH_A: np.array([np.cos(angle), np.sin(angle), 0.0]),
            Branch.BRANCH_B: np.array([np.cos(angle), -np.sin(angle), 0.0]),
        }
        self._mesh_distance = None
        self._distance_cache = None
        self.distance_cache_hits = 0
        self.distance_cache_misses = 0

    @property
    def identity(self):
        return hashlib.sha256(json.dumps(self.config, sort_keys=True).encode()).hexdigest()

    def planes(self):
        result = {"INLET": dict(origin_m=[-self.parent_length, 0, 0], normal=[-1, 0, 0],
                                radius_m=self.parent_radius),
                  "BIFURCATION": dict(origin_m=[self.config["bifurcation_plane_x_m"], 0, 0],
                                      normal=[1, 0, 0])}
        for name, direction in self.directions.items():
            result["OUTLET_" + name[-1]] = dict(origin_m=(direction * self.daughter_length).tolist(),
                                                normal=direction.tolist(), radius_m=self.daughter_radius)
            result["CLASSIFY_" + name[-1]] = dict(
                origin_m=(direction * self.config["classification_distance_m"]).tolist(),
                normal=direction.tolist(), radius_m=self.daughter_radius)
        return result

    def _tube_values(self, positions):
        x = np.atleast_2d(np.asarray(positions, float))
        directions = [np.array([-1., 0, 0]), *self.directions.values()]
        lengths = [self.parent_length, self.daughter_length, self.daughter_length]
        radii = [self.parent_radius, self.daughter_radius, self.daughter_radius]
        values, gradients = [], []
        for direction, length, radius in zip(directions, lengths, radii):
            along = np.clip(x @ direction, 0, length + 3 * radius)
            delta = x - along[:, None] * direction
            distance = np.linalg.norm(delta, axis=1)
            values.append(distance - radius)
            gradients.append(delta / np.maximum(distance[:, None], 1e-30))
        values = np.stack(values, axis=1)
        weights = np.exp(-values / self.blend - logsumexp(-values / self.blend, axis=1)[:, None])
        return -self.blend * logsumexp(-values / self.blend, axis=1), np.einsum(
            "ni,nij->nj", weights, np.stack(gradients, axis=1))

    def level_and_gradient(self, positions):
        x = np.atleast_2d(np.asarray(positions, float))
        if self._mesh_distance is not None:
            previous = self._distance_cache
            if previous is not None and previous[0].shape == x.shape:
                changed = np.any(x != previous[0],axis=1)
                values,gradients = previous[1].copy(),previous[2].copy()
                self.distance_cache_hits += int((~changed).sum())
            else:
                changed = np.ones(len(x),bool)
                values,gradients = np.empty(len(x)),np.empty_like(x)
            self.distance_cache_misses += int(changed.sum())
            for row in np.flatnonzero(changed):
                p = x[row]
                values[row] = self._mesh_distance.EvaluateFunction(p)
                normal = [0., 0., 0.]
                self._mesh_distance.EvaluateGradient(p, normal)
                gradients[row] = normal
            self._distance_cache = (x.copy(),values.copy(),gradients.copy())
            return values,gradients
        value, gradient = self._tube_values(x)
        for name in ["INLET", "OUTLET_A", "OUTLET_B"]:
            plane = self.planes()[name]
            cap = (x - plane["origin_m"]) @ np.array(plane["normal"])
            take = cap > value
            value[take] = cap[take]
            gradient[take] = plane["normal"]
        return value, gradient

    def bind_mesh_surface(self, surface):
        import vtk
        self._distance_cache = None
        self._mesh_distance = vtk.vtkImplicitPolyDataDistance()
        self._mesh_distance.SetInput(surface)

    def wall_distance(self, positions):
        value, gradient = self.level_and_gradient(positions)
        return -value / np.maximum(np.linalg.norm(gradient, axis=1), 1e-12)

    def wall_gap_and_normal(self, positions, radii):
        """Gap and outward normal from the same exact cached VTK level/gradient."""
        value, gradient = self.level_and_gradient(positions)
        norm = np.linalg.norm(gradient, axis=1)
        valid = np.isfinite(value) & np.isfinite(gradient).all(axis=1) & (norm > 1e-12)
        normal = gradient / np.maximum(norm[:, None], 1e-12)
        gap = -value / np.maximum(norm, 1e-12) - np.asarray(radii)
        return gap, normal, valid

    def contains(self, positions, radii=0.0):
        return self.wall_distance(positions) >= np.asarray(radii)

    def concentration_roi_volume(self):
        """Nominal daughter cylinder segment intersected with its y half-space.

        Both mirrored ROIs have equal volume. The cylinder may differ slightly
        from the meshed lumen; this is an explicitly nominal concentration proxy.
        """
        from scipy.integrate import quad
        radius = self.daughter_radius
        tangent = np.tan(np.deg2rad(self.config["branch_angle_deg"]))

        def area(s):
            h = min(s*tangent, radius)
            return radius**2*np.arccos(-h/radius)+h*np.sqrt(max(radius**2-h**2,0.))

        return quad(area,self.config["concentration_roi_start_m"],
                    self.config["concentration_roi_end_m"],epsabs=1e-18)[0]

    def project(self, positions, radii):
        x = np.asarray(positions).copy()
        radii = np.broadcast_to(radii, (len(x),))
        hit = np.zeros(len(x), bool)
        for _ in range(30):
            value, gradient = self.level_and_gradient(x)
            norm = np.linalg.norm(gradient, axis=1)
            excess = value + radii * norm
            active = excess > 1e-12
            if not active.any():
                return x, hit
            hit |= active
            x[active] -= ((excess[active] + 0.5e-12) / np.maximum(norm[active] ** 2, 1e-24))[:, None] * gradient[active]
        raise RuntimeError("Y lumen projection did not converge")

    def constrain_velocity(self, positions, radii, velocity):
        value, gradient = self.level_and_gradient(positions)
        norm = np.linalg.norm(gradient, axis=1)
        normal = gradient / np.maximum(norm[:, None], 1e-30)
        outward = np.einsum("ij,ij->i", velocity, normal)
        active = (value + np.asarray(radii) * norm >= -1e-12) & (outward > 0)
        result = np.array(velocity, copy=True)
        result[active] -= outward[active, None] * normal[active]
        return result


class BranchClassifier:
    """Only the first outward crossing of a finite classification disk assigns a branch."""
    def __init__(self, geometry, particle_ids, tolerance_m=None):
        self.geometry = geometry
        self.tolerance = geometry.config["classification_tolerance_m"] if tolerance_m is None else tolerance_m
        if self.tolerance < 0:
            raise ValueError("Nonnegative classification tolerance required")
        ids = np.asarray(particle_ids, int)
        if len(set(ids)) != len(ids):
            raise ValueError("Particle IDs must be unique")
        self.assignment = {int(i): Branch.UNASSIGNED for i in ids}
        self.events = []
        self.bifurcation_events = {}

    def observe(self, ids, old, new, old_velocity, new_velocity, t0, t1):
        if not t1 > t0:
            raise ValueError("Crossing segment must advance time")
        new_events = []
        planes = self.geometry.planes()
        for row, pid in enumerate(ids):
            pid = int(pid)
            if pid not in self.assignment:
                raise ValueError("Unknown particle ID")
            if self.assignment[pid] != Branch.UNASSIGNED:
                continue
            candidates = []
            for key, branch in [("BIFURCATION", None), ("CLASSIFY_A", Branch.BRANCH_A),
                                ("CLASSIFY_B", Branch.BRANCH_B)]:
                plane = planes[key]
                normal = np.array(plane["normal"], float)
                start = np.dot(old[row] - plane["origin_m"], normal)
                end = np.dot(new[row] - plane["origin_m"], normal)
                if not start <= self.tolerance < end:
                    continue
                alpha = (self.tolerance - start) / (end - start)
                position = old[row] + alpha * (new[row] - old[row])
                offset = position - plane["origin_m"]
                transverse = offset - np.dot(offset, normal) * normal
                if "radius_m" in plane and np.linalg.norm(transverse) > plane["radius_m"] + self.tolerance:
                    continue
                velocity = old_velocity[row] + alpha * (new_velocity[row] - old_velocity[row])
                event = dict(particle_id=pid, branch=branch, crossing_time_s=float(t0 + alpha * (t1-t0)),
                             position_m=position.tolist(), velocity_m_s=velocity.tolist(),
                             transverse_velocity_m_s=(velocity - np.dot(velocity, normal) * normal).tolist(),
                             plane=key, classification_tolerance_m=float(self.tolerance))
                if branch is None:
                    self.bifurcation_events.setdefault(pid, event)
                else:
                    candidates.append(event)
            if candidates:
                first = min(candidates, key=lambda event: event["crossing_time_s"])
                self.assignment[pid] = first["branch"]
                self.events.append(first)
                new_events.append(first)
        return new_events

    def wall_loss(self, pid, reason):
        pid = int(pid)
        if self.assignment[pid] == Branch.UNASSIGNED:
            self.assignment[pid] = Branch.WALL_LOSS
            self.events.append(dict(particle_id=pid, branch=Branch.WALL_LOSS, reason=reason))

    def summary(self):
        counts = {branch: sum(value == branch for value in self.assignment.values()) for branch in Branch}
        assigned = counts[Branch.BRANCH_A] + counts[Branch.BRANCH_B]
        return dict(N_A=counts[Branch.BRANCH_A], N_B=counts[Branch.BRANCH_B],
                    N_unassigned=counts[Branch.UNASSIGNED], N_wall_loss=counts[Branch.WALL_LOSS],
                    N_assigned=assigned, f_A=counts[Branch.BRANCH_A]/assigned if assigned else None,
                    f_B=counts[Branch.BRANCH_B]/assigned if assigned else None,
                    assignment={str(k): str(v) for k, v in self.assignment.items()},
                    crossings=self.events, bifurcation_crossings=list(self.bifurcation_events.values()))
