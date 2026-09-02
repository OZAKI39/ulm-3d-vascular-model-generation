"""Persistence for candidate vessel beds and selected target masks."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .molecular_target_candidates import (
    MolecularTargetCandidate,
    MolecularTargetCandidateCatalog,
    simplify_candidate_selection,
)
from .molecular_target_auto_selection import (
    AutomaticInfluenceAnchorResult,
    AutomaticTargetCandidateEvaluation,
)
from .molecular_target_spatial_heterogeneity import (
    SpatiallyHeterogeneousTargetResult,
)
from ..geometry.vessel_bed_topology import VesselBedTopology, VesselBedUnit


CANDIDATE_SCHEMA_VERSION = "v3"
SUPPORTED_CANDIDATE_SCHEMA_VERSIONS = frozenset(
    {"v1", "v2", CANDIDATE_SCHEMA_VERSION}
)


def save_candidate_catalog(
    npz_path: Path,
    json_path: Path,
    catalog: MolecularTargetCandidateCatalog,
) -> None:
    """Save compact numerical data plus a human-readable candidate table."""

    npz_path = Path(npz_path)
    json_path = Path(json_path)
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    unit_segment_offsets, unit_segment_ids = _pack_rows(
        [unit.segment_ids for unit in catalog.topology.units]
    )
    unit_child_offsets, unit_child_ids = _pack_rows(
        [unit.child_unit_ids for unit in catalog.topology.units]
    )
    candidate_unit_offsets, candidate_unit_ids = _pack_rows(
        [candidate.member_unit_ids for candidate in catalog.candidates]
    )
    np.savez_compressed(
        npz_path,
        schema_version=np.asarray(CANDIDATE_SCHEMA_VERSION),
        x_um=np.asarray(catalog.x_coordinates_um, dtype=np.float64),
        z_um=np.asarray(catalog.z_coordinates_um, dtype=np.float64),
        unit_id_grid=np.asarray(catalog.unit_id_grid, dtype=np.int32),
        candidate_support_mask=np.asarray(catalog.candidate_support_mask, dtype=bool),
        lumen_mask=np.asarray(catalog.lumen_mask, dtype=bool),
        solid_wall_mask=np.asarray(catalog.solid_wall_mask, dtype=bool),
        open_boundary_mask=np.asarray(catalog.open_boundary_mask, dtype=bool),
        unresolved_junction_cells=np.asarray(catalog.unresolved_junction_cells, dtype=np.int64),
        unit_id=np.asarray([unit.unit_id for unit in catalog.topology.units], dtype=np.int32),
        unit_parent_id=np.asarray(
            [unit.parent_unit_id for unit in catalog.topology.units], dtype=np.int32
        ),
        unit_root_id=np.asarray(
            [unit.root_unit_id for unit in catalog.topology.units], dtype=np.int32
        ),
        unit_flow_rate_um3_s=np.asarray(
            [unit.flow_rate_um3_s for unit in catalog.topology.units], dtype=np.float64
        ),
        unit_length_um=np.asarray(
            [unit.length_um for unit in catalog.topology.units], dtype=np.float64
        ),
        unit_volume_um3=np.asarray(
            [unit.volume_um3 for unit in catalog.topology.units], dtype=np.float64
        ),
        unit_endothelial_wall_area_um2=np.asarray(
            [unit.endothelial_wall_area_um2 for unit in catalog.topology.units],
            dtype=np.float64,
        ),
        unit_wall_area_centroid_x_um=np.asarray(
            [unit.wall_area_centroid_x_um for unit in catalog.topology.units],
            dtype=np.float64,
        ),
        unit_wall_area_centroid_z_um=np.asarray(
            [unit.wall_area_centroid_z_um for unit in catalog.topology.units],
            dtype=np.float64,
        ),
        unit_wall_area_second_moment_um4=np.asarray(
            [unit.wall_area_second_moment_um4 for unit in catalog.topology.units],
            dtype=np.float64,
        ),
        unit_topology_depth=np.asarray(
            [unit.topology_depth for unit in catalog.topology.units], dtype=np.int32
        ),
        unit_perfused=np.asarray(
            [unit.perfused for unit in catalog.topology.units], dtype=bool
        ),
        unit_segment_offsets=unit_segment_offsets,
        unit_segment_ids=unit_segment_ids,
        unit_child_offsets=unit_child_offsets,
        unit_child_ids=unit_child_ids,
        candidate_id=np.asarray(
            [candidate.candidate_id for candidate in catalog.candidates], dtype=str
        ),
        candidate_kind=np.asarray(
            [candidate.kind for candidate in catalog.candidates], dtype=str
        ),
        candidate_label=np.asarray(
            [candidate.label for candidate in catalog.candidates], dtype=str
        ),
        candidate_root_unit_id=np.asarray(
            [candidate.root_unit_id for candidate in catalog.candidates], dtype=np.int32
        ),
        candidate_parent_id=np.asarray(
            [candidate.parent_candidate_id or "" for candidate in catalog.candidates],
            dtype=str,
        ),
        candidate_depth=np.asarray(
            [candidate.depth for candidate in catalog.candidates], dtype=np.int32
        ),
        candidate_topology_depth=np.asarray(
            [candidate.topology_depth for candidate in catalog.candidates], dtype=np.int32
        ),
        candidate_inlet_flow_um3_s=np.asarray(
            [candidate.inlet_flow_um3_s for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_root_flow_fraction=np.asarray(
            [candidate.root_flow_fraction for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_network_flow_fraction=np.asarray(
            [candidate.network_flow_fraction for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_volume_um3=np.asarray(
            [candidate.volume_um3 for candidate in catalog.candidates], dtype=np.float64
        ),
        candidate_residence_time_s=np.asarray(
            [candidate.residence_time_s for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_mean_wall_shear_pa=np.asarray(
            [candidate.mean_wall_shear_pa for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_endothelial_wall_area_um2=np.asarray(
            [candidate.endothelial_wall_area_um2 for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_endothelial_wall_area_fraction=np.asarray(
            [candidate.endothelial_wall_area_fraction for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_wall_area_centroid_x_um=np.asarray(
            [candidate.wall_area_centroid_x_um for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_wall_area_centroid_z_um=np.asarray(
            [candidate.wall_area_centroid_z_um for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_radius_of_gyration_um=np.asarray(
            [candidate.radius_of_gyration_um for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        candidate_expected_bubble_visits=np.asarray(
            [candidate.expected_bubble_visits for candidate in catalog.candidates],
            dtype=np.float64,
        ),
        network_endothelial_wall_area_um2=np.asarray(
            catalog.network_endothelial_wall_area_um2, dtype=np.float64
        ),
        network_inlet_flow_um3_s=np.asarray(
            catalog.network_inlet_flow_um3_s, dtype=np.float64
        ),
        injection_rate_per_s=np.asarray(catalog.injection_rate_per_s, dtype=np.float64),
        observation_time_s=np.asarray(catalog.observation_time_s, dtype=np.float64),
        wall_area_weight_um2=np.asarray(catalog.wall_area_weight_um2, dtype=np.float64),
        wall_segment_id_grid=np.asarray(catalog.wall_segment_id_grid, dtype=np.int32),
        accessible_wall_mask=np.asarray(catalog.accessible_wall_mask, dtype=bool),
        expected_bubble_visits_by_unit=np.asarray(
            catalog.expected_bubble_visits_by_unit,
            dtype=np.float64,
        ),
        mapped_endothelial_wall_area_um2=np.asarray(
            catalog.mapped_endothelial_wall_area_um2,
            dtype=np.float64,
        ),
        unmapped_endothelial_wall_area_um2=np.asarray(
            catalog.unmapped_endothelial_wall_area_um2,
            dtype=np.float64,
        ),
        candidate_unit_offsets=candidate_unit_offsets,
        candidate_unit_ids=candidate_unit_ids,
    )

    document = {
        "schema_version": CANDIDATE_SCHEMA_VERSION,
        "meaning": "Selectable candidate vessel beds, not automatically predicted tumour regions.",
        "grid_shape": list(catalog.shape),
        "candidate_count": len(catalog.candidates),
        "network_endothelial_wall_area_um2": _finite_json(
            catalog.network_endothelial_wall_area_um2
        ),
        "network_inlet_flow_um3_s": _finite_json(catalog.network_inlet_flow_um3_s),
        "injection_rate_per_s": _finite_json(catalog.injection_rate_per_s),
        "observation_time_s": _finite_json(catalog.observation_time_s),
        "mapped_endothelial_wall_area_um2": _finite_json(
            catalog.mapped_endothelial_wall_area_um2
        ),
        "unmapped_endothelial_wall_area_um2": _finite_json(
            catalog.unmapped_endothelial_wall_area_um2
        ),
        "automatic_selection_metrics_available": catalog.automatic_metrics_available,
        "unresolved_junction_cells_filled_from_nearest_flow_resolved_basin": int(
            catalog.unresolved_junction_cells
        ),
        "candidates": [_candidate_json(candidate) for candidate in catalog.candidates],
    }
    json_path.write_text(
        json.dumps(document, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def load_candidate_catalog(path: Path) -> MolecularTargetCandidateCatalog:
    """Load a candidate catalog without requiring the original simulation objects."""

    path = Path(path)
    with np.load(path, allow_pickle=False) as data:
        schema = str(np.asarray(data["schema_version"]).item())
        if schema not in SUPPORTED_CANDIDATE_SCHEMA_VERSIONS:
            raise ValueError(
                f"Unsupported candidate catalog schema {schema!r}; expected one of "
                f"{sorted(SUPPORTED_CANDIDATE_SCHEMA_VERSIONS)!r}."
            )
        unit_segment_rows = _unpack_rows(data["unit_segment_offsets"], data["unit_segment_ids"])
        unit_child_rows = _unpack_rows(data["unit_child_offsets"], data["unit_child_ids"])
        unit_depths = _loaded_unit_depths(data, schema)
        units = tuple(
            VesselBedUnit(
                unit_id=int(data["unit_id"][index]),
                segment_ids=unit_segment_rows[index],
                parent_unit_id=int(data["unit_parent_id"][index]),
                child_unit_ids=unit_child_rows[index],
                root_unit_id=int(data["unit_root_id"][index]),
                flow_rate_um3_s=float(data["unit_flow_rate_um3_s"][index]),
                length_um=float(data["unit_length_um"][index]),
                volume_um3=float(data["unit_volume_um3"][index]),
                endothelial_wall_area_um2=_loaded_float(
                    data, "unit_endothelial_wall_area_um2", index
                ),
                wall_area_centroid_x_um=_loaded_float(
                    data, "unit_wall_area_centroid_x_um", index
                ),
                wall_area_centroid_z_um=_loaded_float(
                    data, "unit_wall_area_centroid_z_um", index
                ),
                wall_area_second_moment_um4=_loaded_float(
                    data, "unit_wall_area_second_moment_um4", index
                ),
                topology_depth=int(unit_depths[index]),
                perfused=bool(data["unit_perfused"][index]),
            )
            for index in range(np.asarray(data["unit_id"]).size)
        )
        segment_to_unit = {
            segment_id: unit.unit_id
            for unit in units
            for segment_id in unit.segment_ids
        }
        topology = VesselBedTopology(
            units=units,
            segment_id_to_unit_id=segment_to_unit,
            root_unit_ids=tuple(unit.unit_id for unit in units if unit.parent_unit_id < 0),
        )
        candidate_unit_rows = _unpack_rows(
            data["candidate_unit_offsets"], data["candidate_unit_ids"]
        )
        candidate_ids = np.asarray(data["candidate_id"], dtype=str)
        candidates = tuple(
            MolecularTargetCandidate(
                candidate_id=str(candidate_ids[index]),
                kind=str(data["candidate_kind"][index]),
                label=str(data["candidate_label"][index]),
                root_unit_id=int(data["candidate_root_unit_id"][index]),
                member_unit_ids=candidate_unit_rows[index],
                parent_candidate_id=(
                    str(data["candidate_parent_id"][index]) or None
                ),
                depth=int(data["candidate_depth"][index]),
                topology_depth=int(
                    data["candidate_topology_depth"][index]
                    if "candidate_topology_depth" in data
                    else units[int(data["candidate_root_unit_id"][index])].topology_depth
                ),
                inlet_flow_um3_s=float(data["candidate_inlet_flow_um3_s"][index]),
                root_flow_fraction=float(data["candidate_root_flow_fraction"][index]),
                network_flow_fraction=_loaded_float(
                    data, "candidate_network_flow_fraction", index
                ),
                volume_um3=float(data["candidate_volume_um3"][index]),
                residence_time_s=float(data["candidate_residence_time_s"][index]),
                mean_wall_shear_pa=float(data["candidate_mean_wall_shear_pa"][index]),
                endothelial_wall_area_um2=_loaded_float(
                    data, "candidate_endothelial_wall_area_um2", index
                ),
                endothelial_wall_area_fraction=_loaded_float(
                    data, "candidate_endothelial_wall_area_fraction", index
                ),
                wall_area_centroid_x_um=_loaded_float(
                    data, "candidate_wall_area_centroid_x_um", index
                ),
                wall_area_centroid_z_um=_loaded_float(
                    data, "candidate_wall_area_centroid_z_um", index
                ),
                radius_of_gyration_um=_loaded_float(
                    data, "candidate_radius_of_gyration_um", index
                ),
                expected_bubble_visits=_loaded_float(
                    data, "candidate_expected_bubble_visits", index
                ),
            )
            for index in range(candidate_ids.size)
        )
        return MolecularTargetCandidateCatalog(
            x_coordinates_um=np.asarray(data["x_um"], dtype=np.float64).copy(),
            z_coordinates_um=np.asarray(data["z_um"], dtype=np.float64).copy(),
            unit_id_grid=np.asarray(data["unit_id_grid"], dtype=np.int32).copy(),
            candidate_support_mask=np.asarray(
                data["candidate_support_mask"], dtype=bool
            ).copy(),
            lumen_mask=np.asarray(data["lumen_mask"], dtype=bool).copy(),
            solid_wall_mask=np.asarray(data["solid_wall_mask"], dtype=bool).copy(),
            open_boundary_mask=np.asarray(
                data["open_boundary_mask"], dtype=bool
            ).copy(),
            topology=topology,
            candidates=candidates,
            unresolved_junction_cells=int(
                np.asarray(data["unresolved_junction_cells"]).item()
            ),
            network_endothelial_wall_area_um2=_loaded_scalar(
                data, "network_endothelial_wall_area_um2"
            ),
            network_inlet_flow_um3_s=_loaded_scalar(data, "network_inlet_flow_um3_s"),
            injection_rate_per_s=_loaded_scalar(data, "injection_rate_per_s"),
            observation_time_s=_loaded_scalar(data, "observation_time_s"),
            wall_area_weight_um2=(
                np.asarray(data["wall_area_weight_um2"], dtype=np.float64).copy()
                if "wall_area_weight_um2" in data
                else None
            ),
            wall_segment_id_grid=(
                np.asarray(data["wall_segment_id_grid"], dtype=np.int32).copy()
                if "wall_segment_id_grid" in data
                else None
            ),
            accessible_wall_mask=(
                np.asarray(data["accessible_wall_mask"], dtype=bool).copy()
                if "accessible_wall_mask" in data
                else None
            ),
            expected_bubble_visits_by_unit=(
                np.asarray(
                    data["expected_bubble_visits_by_unit"],
                    dtype=np.float64,
                ).copy()
                if "expected_bubble_visits_by_unit" in data
                else None
            ),
            mapped_endothelial_wall_area_um2=_loaded_scalar(
                data,
                "mapped_endothelial_wall_area_um2",
            ),
            unmapped_endothelial_wall_area_um2=_loaded_scalar(
                data,
                "unmapped_endothelial_wall_area_um2",
            ),
            schema_version=schema,
        )


def save_selected_target_mask(
    path: Path,
    catalog: MolecularTargetCandidateCatalog,
    candidate_ids: list[str] | tuple[str, ...],
    *,
    selection_mode: str = "manual",
    requested_wall_area_fraction: float = math.nan,
    achieved_wall_area_fraction: float = math.nan,
    manually_modified_after_automatic_preview: bool = False,
) -> tuple[str, ...]:
    """Save the canonical coordinate-aware Boolean mask consumed by mask_npz mode."""

    selected = simplify_candidate_selection(catalog, candidate_ids)
    if not selected:
        raise ValueError("Select at least one candidate before saving a molecular target mask.")
    mask = catalog.mask_for_candidate_ids(selected)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        x_um=np.asarray(catalog.x_coordinates_um, dtype=np.float64),
        z_um=np.asarray(catalog.z_coordinates_um, dtype=np.float64),
        target_mask=np.asarray(mask, dtype=bool),
        selected_candidate_ids=np.asarray(selected, dtype=str),
        candidate_schema_version=np.asarray(catalog.schema_version),
        selection_mode=np.asarray(str(selection_mode)),
        requested_wall_area_fraction=np.asarray(
            requested_wall_area_fraction, dtype=np.float64
        ),
        achieved_wall_area_fraction=np.asarray(
            achieved_wall_area_fraction, dtype=np.float64
        ),
        manually_modified_after_automatic_preview=np.asarray(
            manually_modified_after_automatic_preview, dtype=bool
        ),
        target_mask_semantics=np.asarray(
            "Selected candidate vessel-bed union; molecular targets are applied only on eligible solid walls."
        ),
    )
    return selected


def save_spatially_heterogeneous_target_mask(
    path: Path,
    catalog: MolecularTargetCandidateCatalog,
    result: SpatiallyHeterogeneousTargetResult,
) -> None:
    """Save the final positive wall patches and their reproducible spatial realization."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        x_um=np.asarray(catalog.x_coordinates_um, dtype=np.float64),
        z_um=np.asarray(catalog.z_coordinates_um, dtype=np.float64),
        target_mask=np.asarray(result.target_positive_wall_mask, dtype=bool),
        influence_region_mask=np.asarray(result.influence_region_mask, dtype=bool),
        influence_wall_mask=np.asarray(result.influence_wall_mask, dtype=bool),
        anchor_candidate_id=np.asarray(result.anchor_candidate_id),
        influence_center_xz_um=np.asarray(
            result.influence_center_xz_um,
            dtype=np.float64,
        ),
        influence_radius_um=np.asarray(result.influence_radius_um, dtype=np.float64),
        requested_influence_wall_area_fraction=np.asarray(
            result.requested_influence_wall_area_fraction,
            dtype=np.float64,
        ),
        achieved_influence_wall_area_fraction=np.asarray(
            result.achieved_influence_wall_area_fraction,
            dtype=np.float64,
        ),
        requested_positive_wall_fraction_within_influence=np.asarray(
            result.requested_positive_wall_fraction_within_influence,
            dtype=np.float64,
        ),
        achieved_positive_wall_fraction_within_influence=np.asarray(
            result.achieved_positive_wall_fraction_within_influence,
            dtype=np.float64,
        ),
        target_positive_network_wall_area_fraction=np.asarray(
            result.target_positive_network_wall_area_fraction,
            dtype=np.float64,
        ),
        target_correlation_length_um=np.asarray(
            result.correlation_length_um,
            dtype=np.float64,
        ),
        random_seed=np.asarray(result.random_seed, dtype=np.int64),
        random_field_modes=np.asarray(result.random_field_modes, dtype=np.int64),
        random_field_algorithm=np.asarray(result.random_field_algorithm),
        random_wavevectors_um_inv=np.asarray(
            result.random_wavevectors_um_inv,
            dtype=np.float64,
        ),
        random_phases_rad=np.asarray(result.random_phases_rad, dtype=np.float64),
        random_field_threshold=np.asarray(
            result.random_field_threshold,
            dtype=np.float64,
        ),
        influence_wall_flat_indices=np.asarray(
            result.influence_wall_flat_indices,
            dtype=np.int64,
        ),
        influence_wall_random_field=np.asarray(
            result.influence_wall_random_field,
            dtype=np.float64,
        ),
        patch_count=np.asarray(result.patch_count, dtype=np.int64),
        candidate_schema_version=np.asarray(catalog.schema_version),
        selection_mode=np.asarray("automatic_spatial_heterogeneity"),
        target_mask_semantics=np.asarray(
            "Spatially correlated target-positive closed-wall patches inside a compact physical influence region."
        ),
    )


def save_spatially_heterogeneous_target_report(
    path: Path,
    catalog: MolecularTargetCandidateCatalog,
    anchor: AutomaticInfluenceAnchorResult,
    result: SpatiallyHeterogeneousTargetResult,
) -> None:
    """Save biological, numerical, accessibility, and anchor provenance as JSON."""

    anchor_candidate = catalog.candidate_by_id(anchor.anchor_candidate_id)
    document = {
        "schema_version": "v2",
        "meaning": (
            "Reproducible spatially heterogeneous synthetic GBM-style molecular target; "
            "this is not tumour detection or biological calibration."
        ),
        "selection_mode": "automatic_spatial_heterogeneity",
        "anchor": {
            "candidate_id": anchor.anchor_candidate_id,
            "center_x_um": result.influence_center_xz_um[0],
            "center_z_um": result.influence_center_xz_um[1],
            "eligible_candidate_count": anchor.eligible_candidate_count,
            "topology_depth": anchor_candidate.topology_depth,
            "radius_of_gyration_um": anchor_candidate.radius_of_gyration_um,
            "network_flow_fraction": anchor_candidate.network_flow_fraction,
            "expected_bubble_visits": anchor_candidate.expected_bubble_visits,
            "residence_time_s": anchor_candidate.residence_time_s,
            "mean_wall_shear_pa": anchor_candidate.mean_wall_shear_pa,
        },
        "influence_region": {
            "shape": "physical_xz_circle",
            "radius_um": result.influence_radius_um,
            "requested_network_wall_area_fraction": (
                result.requested_influence_wall_area_fraction
            ),
            "achieved_accessible_network_wall_area_fraction": (
                result.achieved_influence_wall_area_fraction
            ),
            "accessible_wall_area_um2": result.influence_accessible_wall_area_um2,
            "inaccessible_wall_area_inside_region_um2": (
                result.inaccessible_wall_area_inside_influence_um2
            ),
        },
        "target_positive_patches": {
            "requested_positive_fraction_within_influence": (
                result.requested_positive_wall_fraction_within_influence
            ),
            "achieved_positive_fraction_within_influence": (
                result.achieved_positive_wall_fraction_within_influence
            ),
            "positive_wall_area_um2": result.target_positive_wall_area_um2,
            "positive_network_wall_area_fraction": (
                result.target_positive_network_wall_area_fraction
            ),
            "grid_connected_patch_count": result.patch_count,
        },
        "random_field": {
            "algorithm": result.random_field_algorithm,
            "correlation_length_um": result.correlation_length_um,
            "correlation_length_grid_cells": result.correlation_length_grid_cells,
            "seed": result.random_seed,
            "mode_count": result.random_field_modes,
            "threshold": result.random_field_threshold,
            "coefficient_arrays_saved_in_target_npz": True,
        },
        "accessibility": {
            "network_inlet_flow_um3_s": catalog.network_inlet_flow_um3_s,
            "injection_rate_per_s": catalog.injection_rate_per_s,
            "observation_time_s": catalog.observation_time_s,
            "minimum_expected_bubble_visits": 1.0,
        },
        "wall_measure": {
            "theoretical_network_area_um2": (
                catalog.network_endothelial_wall_area_um2
            ),
            "mapped_area_um2": catalog.mapped_endothelial_wall_area_um2,
            "unmapped_area_um2": catalog.unmapped_endothelial_wall_area_um2,
        },
        "anchor_candidates": [
            _automatic_evaluation_json(catalog, evaluation)
            for evaluation in anchor.evaluations
        ],
        "default_probability_excludes": [
            "speed",
            "wall shear stress",
            "residence time",
        ],
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, ensure_ascii=False), encoding="utf-8")


def _pack_rows(rows: list[tuple[int, ...]]) -> tuple[np.ndarray, np.ndarray]:
    offsets = np.zeros(len(rows) + 1, dtype=np.int64)
    values: list[int] = []
    for index, row in enumerate(rows):
        values.extend(int(value) for value in row)
        offsets[index + 1] = len(values)
    return offsets, np.asarray(values, dtype=np.int32)


def _unpack_rows(offsets: np.ndarray, values: np.ndarray) -> tuple[tuple[int, ...], ...]:
    packed_offsets = np.asarray(offsets, dtype=np.int64)
    packed_values = np.asarray(values, dtype=np.int32)
    return tuple(
        tuple(int(value) for value in packed_values[packed_offsets[index] : packed_offsets[index + 1]])
        for index in range(packed_offsets.size - 1)
    )


def _candidate_json(candidate: MolecularTargetCandidate) -> dict[str, object]:
    return {
        "candidate_id": candidate.candidate_id,
        "kind": candidate.kind,
        "label": candidate.label,
        "parent_candidate_id": candidate.parent_candidate_id,
        "depth": candidate.depth,
        "topology_depth": candidate.topology_depth,
        "root_unit_id": candidate.root_unit_id,
        "member_unit_ids": list(candidate.member_unit_ids),
        "inlet_flow_um3_s": _finite_json(candidate.inlet_flow_um3_s),
        "root_flow_fraction": _finite_json(candidate.root_flow_fraction),
        "network_flow_fraction": _finite_json(candidate.network_flow_fraction),
        "volume_um3": _finite_json(candidate.volume_um3),
        "residence_time_s": _finite_json(candidate.residence_time_s),
        "mean_wall_shear_pa": _finite_json(candidate.mean_wall_shear_pa),
        "endothelial_wall_area_um2": _finite_json(candidate.endothelial_wall_area_um2),
        "endothelial_wall_area_fraction": _finite_json(
            candidate.endothelial_wall_area_fraction
        ),
        "wall_area_centroid_x_um": _finite_json(candidate.wall_area_centroid_x_um),
        "wall_area_centroid_z_um": _finite_json(candidate.wall_area_centroid_z_um),
        "radius_of_gyration_um": _finite_json(candidate.radius_of_gyration_um),
        "expected_bubble_visits": _finite_json(candidate.expected_bubble_visits),
    }


def _automatic_evaluation_json(
    catalog: MolecularTargetCandidateCatalog,
    evaluation: AutomaticTargetCandidateEvaluation,
) -> dict[str, object]:
    candidate = catalog.candidate_by_id(evaluation.candidate_id)
    return {
        "candidate_id": evaluation.candidate_id,
        "eligible": evaluation.eligible,
        "rejection_reason": evaluation.rejection_reason,
        "endothelial_wall_area_um2": _finite_json(
            candidate.endothelial_wall_area_um2
        ),
        "endothelial_wall_area_fraction": _finite_json(
            evaluation.endothelial_wall_area_fraction
        ),
        "area_log_error": _finite_json(evaluation.area_log_error),
        "topology_depth": candidate.topology_depth,
        "radius_of_gyration_um": _finite_json(candidate.radius_of_gyration_um),
        "network_flow_fraction": _finite_json(candidate.network_flow_fraction),
        "expected_bubble_visits": _finite_json(evaluation.expected_bubble_visits),
        "residence_time_s": _finite_json(candidate.residence_time_s),
        "mean_wall_shear_pa": _finite_json(candidate.mean_wall_shear_pa),
    }


def _loaded_float(
    data: np.lib.npyio.NpzFile,
    key: str,
    index: int,
) -> float:
    return float(data[key][index]) if key in data else math.nan


def _loaded_scalar(data: np.lib.npyio.NpzFile, key: str) -> float:
    return float(np.asarray(data[key]).item()) if key in data else math.nan


def _loaded_unit_depths(data: np.lib.npyio.NpzFile, schema: str) -> np.ndarray:
    if schema in {"v2", CANDIDATE_SCHEMA_VERSION}:
        return np.asarray(data["unit_topology_depth"], dtype=np.int32)
    parents = np.asarray(data["unit_parent_id"], dtype=np.int32)
    depths = np.zeros(parents.size, dtype=np.int32)
    for unit_id in range(parents.size):
        current = int(unit_id)
        while parents[current] >= 0:
            depths[unit_id] += 1
            current = int(parents[current])
    return depths


def _finite_json(value: float) -> float | None:
    return float(value) if math.isfinite(value) else None
