"""Topology-common flow conservation with dimension-specific physical models."""

from __future__ import annotations

import math


class FlowScalingMixin:
    """
    Assign terminal shares and update lumen sizes through the selected model.

    ``self.haemodynamics`` is either a true planar strip-flow model or a
    volumetric circular-tube model. The two models have different input
    quantities, cross-section formulae, parameter names, and output metadata.
    """

    def _topological_order(self):
        """Return vessel IDs with every parent before its children."""
        if not self.vessels:
            return []

        roots = [v.vid for v in self.vessels if v.parent_id < 0]
        if len(roots) != 1:
            raise ValueError(
                "The direct inlet-flow model requires exactly one authoritative "
                f"root; found {len(roots)}."
            )

        order = []
        seen = set()
        stack = list(reversed(roots))
        while stack:
            vid = stack.pop()
            if vid in seen:
                raise ValueError("The vessel graph must be acyclic.")
            if vid < 0 or vid >= len(self.vessels):
                raise ValueError(f"Invalid vessel ID {vid} in the tree topology.")
            seen.add(vid)
            order.append(vid)
            stack.extend(reversed(self.vessels[vid].children))

        if len(order) != len(self.vessels):
            raise ValueError("The vessel graph must be one connected rooted tree.")
        return order

    def _allocate_terminal_flows(self, order):
        """Give every current leaf an equal share of the fixed inlet flow."""
        terminals = [vid for vid in order if not self.vessels[vid].children]
        if not terminals:
            raise ValueError("The vessel tree has no terminal leaves.")

        inlet_flow = self.haemodynamics.inlet_value
        terminal_flow = inlet_flow / float(len(terminals))
        if not math.isfinite(terminal_flow) or terminal_flow <= 0.0:
            raise ValueError("Equal terminal flow allocation produced an invalid flow.")

        terminal_set = set(terminals)
        for vessel in self.vessels:
            vessel.prescribed_outflow = (
                terminal_flow if vessel.vid in terminal_set else 0.0
            )
        return terminals, terminal_flow

    def _accumulated_flows(self, order):
        """Accumulate terminal target flows from leaves toward the inlet."""
        flows = {}
        for vid in reversed(order):
            vessel = self.vessels[vid]
            if vessel.children:
                flows[vid] = math.fsum(flows[child] for child in vessel.children)
            else:
                flows[vid] = float(vessel.prescribed_outflow)
        return flows

    def _radius_for_flow(self, flow_rate):
        """Return the fixed-root Murray radius for a positive segment flow."""
        inlet_flow = self.haemodynamics.inlet_value
        return self.haemodynamics.root_size_um * (
            float(flow_rate) / inlet_flow
        ) ** (1.0 / self.haemodynamics.murray_gamma)

    def _apply_radii(self):
        """Reallocate terminal flows, conserve flow, and update all radii directly."""
        if not self.vessels:
            return

        order = self._topological_order()
        self._allocate_terminal_flows(order)
        flows = self._accumulated_flows(order)

        root_id = order[0]
        inlet_flow = self.haemodynamics.inlet_value
        root_error = abs(flows[root_id] - inlet_flow)
        if root_error > max(1.0e-9 * inlet_flow, 1.0e-12):
            raise RuntimeError(
                "Terminal allocation does not sum to the configured "
                f"{self.haemodynamics.flow_quantity} inlet value."
            )

        radii = {vid: self._radius_for_flow(flows[vid]) for vid in order}
        too_small = [
            vid
            for vid, radius in radii.items()
            if radius + 1.0e-12 < self.haemodynamics.minimum_size_um
        ]
        if too_small:
            smallest = min(radii[vid] for vid in too_small)
            raise RuntimeError(
                "The dimensional Murray model produced a lumen "
                f"{self.haemodynamics.size_name} below its configured minimum; "
                "increase the dimensional root size, reduce the final terminal "
                f"count, or lower the minimum. smallest={smallest:.9g} um, "
                f"minimum={self.haemodynamics.minimum_size_um:.9g} um."
            )

        gamma = self.haemodynamics.murray_gamma
        velocity_min = self.haemodynamics.velocity_min_mm_s
        velocity_max = self.haemodynamics.velocity_max_mm_s
        invalid_velocity = []

        for vid in order:
            vessel = self.vessels[vid]
            vessel.flow_rate = flows[vid]
            vessel.radius = radii[vid]
            vessel.mean_velocity = self.haemodynamics.mean_velocity_um_s(
                vessel.flow_rate,
                vessel.radius,
            )

            if vessel.children:
                child_flow = math.fsum(flows[child] for child in vessel.children)
                vessel.flow_conservation_residual = abs(
                    vessel.flow_rate - child_flow
                ) / max(vessel.flow_rate, 1.0e-30)
                child_radius_power = math.fsum(
                    radii[child] ** gamma for child in vessel.children
                )
                vessel.murray_residual = abs(
                    vessel.radius ** gamma - child_radius_power
                ) / max(vessel.radius ** gamma, 1.0e-30)
            else:
                vessel.flow_conservation_residual = 0.0
                vessel.murray_residual = 0.0

            velocity_mm_s = vessel.mean_velocity / 1000.0
            if (
                (velocity_min > 0.0 and velocity_mm_s < velocity_min)
                or (velocity_max > 0.0 and velocity_mm_s > velocity_max)
                or not math.isfinite(velocity_mm_s)
            ):
                invalid_velocity.append((vid, velocity_mm_s))

        if invalid_velocity:
            preview = ", ".join(
                f"id={vid}: {velocity:.6g} mm/s"
                for vid, velocity in invalid_velocity[:8]
            )
            raise RuntimeError(
                f"{self.haemodynamics.geometry_mode} mean velocity is outside "
                "the configured validation "
                f"range [{velocity_min:g}, {velocity_max:g}] mm/s: {preview}"
            )

    def _lumen_measure_cost(self):
        """Return dimensional lumen area (2-D) or volume (3-D)."""
        return self.haemodynamics.lumen_measure_cost(self.vessels)
