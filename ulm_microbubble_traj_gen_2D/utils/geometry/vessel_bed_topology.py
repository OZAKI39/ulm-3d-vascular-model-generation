"""Flow-oriented vessel-bed topology used by molecular-target candidates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel


@dataclass(frozen=True)
class VesselBedUnit:
    """One maximal vessel chain between natural inlet, outlet, or branch boundaries."""

    unit_id: int
    segment_ids: tuple[int, ...]
    parent_unit_id: int
    child_unit_ids: tuple[int, ...]
    root_unit_id: int
    flow_rate_um3_s: float
    length_um: float
    volume_um3: float
    endothelial_wall_area_um2: float
    wall_area_centroid_x_um: float
    wall_area_centroid_z_um: float
    wall_area_second_moment_um4: float
    topology_depth: int
    perfused: bool


@dataclass(frozen=True)
class VesselBedTopology:
    """Directed basic units and the original-segment-to-unit lookup."""

    units: tuple[VesselBedUnit, ...]
    segment_id_to_unit_id: dict[int, int]
    root_unit_ids: tuple[int, ...]

    def descendants(self, unit_id: int) -> tuple[int, ...]:
        """Return all downstream unit IDs in deterministic depth-first order."""

        result: list[int] = []
        stack = list(reversed(self.units[int(unit_id)].child_unit_ids))
        while stack:
            current = int(stack.pop())
            result.append(current)
            stack.extend(reversed(self.units[current].child_unit_ids))
        return tuple(result)


class _DisjointSet:
    """Join segment endpoints that represent the same vascular junction."""

    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        parent = self.parent[value]
        while parent != self.parent[parent]:
            parent = self.parent[parent]
        while value != parent:
            next_value = self.parent[value]
            self.parent[value] = parent
            value = next_value
        return parent

    def union(self, first: int, second: int) -> None:
        root_first = self.find(first)
        root_second = self.find(second)
        if root_first != root_second:
            self.parent[root_second] = root_first


def build_vessel_bed_topology(vessels: list[Vessel]) -> VesselBedTopology:
    """Build flow-oriented maximal chains without changing the input vessels."""

    if not vessels:
        raise ValueError("Cannot build molecular-target candidates without vessels.")

    ordered = sorted(vessels, key=lambda vessel: int(vessel.vid))
    vessel_by_id = {int(vessel.vid): vessel for vessel in ordered}
    if len(vessel_by_id) != len(ordered):
        raise ValueError("Vessel IDs must be unique before candidate topology is built.")

    index_by_id = {int(vessel.vid): index for index, vessel in enumerate(ordered)}
    endpoints = _DisjointSet(2 * len(ordered))
    for vessel in ordered:
        parent_id = int(vessel.parent_id)
        if parent_id < 0:
            continue
        if parent_id not in index_by_id:
            raise ValueError(
                f"Vessel {int(vessel.vid)} references missing parent {parent_id}."
            )
        parent_index = index_by_id[parent_id]
        child_index = index_by_id[int(vessel.vid)]
        endpoints.union(2 * parent_index + 1, 2 * child_index)

    source_node: dict[int, int] = {}
    target_node: dict[int, int] = {}
    incoming: dict[int, list[int]] = {}
    outgoing: dict[int, list[int]] = {}
    for index, vessel in enumerate(ordered):
        segment_id = int(vessel.vid)
        flow = float(vessel.flow_rate)
        if not np.isfinite(flow):
            raise ValueError(f"Vessel {segment_id} has a non-finite flow rate.")
        proximal = endpoints.find(2 * index)
        distal = endpoints.find(2 * index + 1)
        source, target = (proximal, distal) if flow >= 0.0 else (distal, proximal)
        source_node[segment_id] = source
        target_node[segment_id] = target
        outgoing.setdefault(source, []).append(segment_id)
        incoming.setdefault(target, []).append(segment_id)
        incoming.setdefault(source, [])
        outgoing.setdefault(target, [])

    for values in incoming.values():
        values.sort()
    for values in outgoing.values():
        values.sort()

    _validate_directed_forest(tuple(vessel_by_id), source_node, target_node, outgoing)

    chains: list[tuple[int, ...]] = []
    visited: set[int] = set()
    for segment_id in sorted(vessel_by_id):
        source = source_node[segment_id]
        if len(incoming[source]) == 1 and len(outgoing[source]) == 1:
            continue
        chains.append(
            _walk_chain(
                segment_id,
                source_node,
                target_node,
                incoming,
                outgoing,
                visited,
            )
        )
    for segment_id in sorted(vessel_by_id):
        if segment_id not in visited:
            chains.append(
                _walk_chain(
                    segment_id,
                    source_node,
                    target_node,
                    incoming,
                    outgoing,
                    visited,
                )
            )

    chains.sort(key=lambda chain: chain[0])
    segment_to_unit = {
        segment_id: unit_id
        for unit_id, chain in enumerate(chains)
        for segment_id in chain
    }
    parent_ids: list[int] = []
    child_ids: list[tuple[int, ...]] = []
    for unit_id, chain in enumerate(chains):
        first_source = source_node[chain[0]]
        external_parents = sorted(
            {
                segment_to_unit[segment_id]
                for segment_id in incoming[first_source]
                if segment_to_unit[segment_id] != unit_id
            }
        )
        if len(external_parents) > 1:
            raise ValueError("Candidate topology must be a tree and cannot contain merging flows.")
        parent_ids.append(external_parents[0] if external_parents else -1)

        last_target = target_node[chain[-1]]
        children = sorted(
            {
                segment_to_unit[segment_id]
                for segment_id in outgoing[last_target]
                if segment_to_unit[segment_id] != unit_id
            }
        )
        child_ids.append(tuple(children))

    root_ids = tuple(index for index, parent_id in enumerate(parent_ids) if parent_id < 0)
    if not root_ids:
        raise ValueError("Candidate topology has no flow inlet.")

    maximum_flow = max(abs(float(vessel.flow_rate)) for vessel in ordered)
    zero_flow_tolerance = max(
        64.0 * np.finfo(float).eps * maximum_flow,
        np.finfo(float).tiny,
    )
    roots = [_root_unit_id(index, parent_ids) for index in range(len(chains))]
    topology_depths = [_unit_depth(index, parent_ids) for index in range(len(chains))]
    units: list[VesselBedUnit] = []
    for unit_id, chain in enumerate(chains):
        unit_vessels = [vessel_by_id[segment_id] for segment_id in chain]
        length_um = float(sum(float(vessel.length()) for vessel in unit_vessels))
        volume_um3 = float(
            sum(
                np.pi * float(vessel.radius) ** 2 * float(vessel.length())
                for vessel in unit_vessels
            )
        )
        segment_wall_areas = np.asarray(
            [
                2.0 * np.pi * float(vessel.radius) * float(vessel.length())
                for vessel in unit_vessels
            ],
            dtype=np.float64,
        )
        segment_centres_xz = np.asarray(
            [
                0.5
                * (
                    np.asarray(vessel.x_p, dtype=np.float64)[[0, 2]]
                    + np.asarray(vessel.x_d, dtype=np.float64)[[0, 2]]
                )
                for vessel in unit_vessels
            ],
            dtype=np.float64,
        )
        wall_area_um2 = float(np.sum(segment_wall_areas))
        if not np.isfinite(wall_area_um2) or wall_area_um2 <= 0.0:
            raise ValueError(
                f"Vessel unit {unit_id} has non-positive endothelial wall area."
            )
        wall_centroid = np.sum(
            segment_wall_areas[:, None] * segment_centres_xz,
            axis=0,
        ) / wall_area_um2
        wall_second_moment = float(
            np.sum(
                segment_wall_areas
                * np.sum(segment_centres_xz * segment_centres_xz, axis=1)
            )
        )
        inlet_flow = abs(float(unit_vessels[0].flow_rate))
        units.append(
            VesselBedUnit(
                unit_id=unit_id,
                segment_ids=chain,
                parent_unit_id=parent_ids[unit_id],
                child_unit_ids=child_ids[unit_id],
                root_unit_id=roots[unit_id],
                flow_rate_um3_s=inlet_flow,
                length_um=length_um,
                volume_um3=volume_um3,
                endothelial_wall_area_um2=wall_area_um2,
                wall_area_centroid_x_um=float(wall_centroid[0]),
                wall_area_centroid_z_um=float(wall_centroid[1]),
                wall_area_second_moment_um4=wall_second_moment,
                topology_depth=topology_depths[unit_id],
                perfused=inlet_flow > zero_flow_tolerance,
            )
        )

    return VesselBedTopology(
        units=tuple(units),
        segment_id_to_unit_id=segment_to_unit,
        root_unit_ids=root_ids,
    )


def _walk_chain(
    first_segment_id: int,
    source_node: dict[int, int],
    target_node: dict[int, int],
    incoming: dict[int, list[int]],
    outgoing: dict[int, list[int]],
    visited: set[int],
) -> tuple[int, ...]:
    chain: list[int] = []
    current = int(first_segment_id)
    while current not in visited:
        visited.add(current)
        chain.append(current)
        target = target_node[current]
        if len(incoming[target]) != 1 or len(outgoing[target]) != 1:
            break
        current = int(outgoing[target][0])
    return tuple(chain)


def _root_unit_id(unit_id: int, parents: list[int]) -> int:
    current = int(unit_id)
    seen: set[int] = set()
    while parents[current] >= 0:
        if current in seen:
            raise ValueError("Candidate unit hierarchy contains a cycle.")
        seen.add(current)
        current = int(parents[current])
    return current


def _unit_depth(unit_id: int, parents: list[int]) -> int:
    depth = 0
    current = int(unit_id)
    while parents[current] >= 0:
        depth += 1
        current = int(parents[current])
    return depth


def _validate_directed_forest(
    segment_ids: tuple[int, ...],
    source_node: dict[int, int],
    target_node: dict[int, int],
    outgoing: dict[int, list[int]],
) -> None:
    state = {segment_id: 0 for segment_id in segment_ids}

    def visit(segment_id: int) -> None:
        if state[segment_id] == 1:
            raise ValueError("Flow-oriented vessel topology contains a cycle.")
        if state[segment_id] == 2:
            return
        state[segment_id] = 1
        for child_id in outgoing[target_node[segment_id]]:
            if source_node[child_id] == target_node[segment_id]:
                visit(child_id)
        state[segment_id] = 2

    for segment_id in sorted(segment_ids):
        visit(segment_id)
