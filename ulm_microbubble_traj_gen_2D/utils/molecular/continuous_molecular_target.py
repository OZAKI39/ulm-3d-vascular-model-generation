"""Migration of a legacy grid target realization onto the v16 wall arclength."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import numpy as np

from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from .molecular_target_spatial_heterogeneity import (
    evaluate_random_fourier_field,
)


TARGET_SCHEMA = "v16_continuous_wall_arclength_target"


def migrate_legacy_grid_target_to_continuous_wall(
    legacy_npz_path: str | Path,
    output_npz_path: str | Path,
    geometry: ContinuousVesselGeometry,
) -> Path:
    """Re-evaluate one saved random field directly on ``Gamma_w``.

    The legacy Boolean pixels are not projected onto the new wall.  Only the
    physical influence circle and deterministic random-Fourier realization are
    retained, so the output is independent of CFD grid spacing.
    """

    source = Path(legacy_npz_path).expanduser().resolve()
    destination = Path(output_npz_path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Legacy molecular target does not exist: {source}")
    required = {
        "influence_center_xz_um",
        "influence_radius_um",
        "requested_positive_wall_fraction_within_influence",
        "target_correlation_length_um",
        "random_seed",
        "random_field_modes",
        "random_field_algorithm",
        "random_wavevectors_um_inv",
        "random_phases_rad",
    }
    with np.load(source, allow_pickle=False) as data:
        missing = sorted(required.difference(data.files))
        if missing:
            raise ValueError(
                "Legacy target NPZ is missing arrays: " + ", ".join(missing)
            )
        center = np.asarray(data["influence_center_xz_um"], dtype=np.float64)
        influence_radius = float(np.asarray(data["influence_radius_um"]).item())
        requested_positive_fraction = float(
            np.asarray(
                data["requested_positive_wall_fraction_within_influence"]
            ).item()
        )
        correlation_length = float(
            np.asarray(data["target_correlation_length_um"]).item()
        )
        random_seed = int(np.asarray(data["random_seed"]).item())
        random_modes = int(np.asarray(data["random_field_modes"]).item())
        random_algorithm = str(
            np.asarray(data["random_field_algorithm"]).item()
        )
        wavevectors = np.asarray(
            data["random_wavevectors_um_inv"], dtype=np.float64
        )
        phases = np.asarray(data["random_phases_rad"], dtype=np.float64)
        anchor = (
            str(np.asarray(data["anchor_candidate_id"]).item())
            if "anchor_candidate_id" in data.files
            else ""
        )
        requested_influence_fraction = (
            float(
                np.asarray(
                    data["requested_influence_wall_area_fraction"]
                ).item()
            )
            if "requested_influence_wall_area_fraction" in data.files
            else np.nan
        )

    if center.shape != (2,) or not np.all(np.isfinite(center)):
        raise ValueError("influence_center_xz_um must be one finite X-Z point.")
    if not np.isfinite(influence_radius) or influence_radius <= 0.0:
        raise ValueError("influence_radius_um must be finite and positive.")
    if not 0.0 <= requested_positive_fraction <= 1.0:
        raise ValueError(
            "requested_positive_wall_fraction_within_influence must be in [0, 1]."
        )

    coordinates = np.asarray(
        geometry.solid_face_center_xz_um, dtype=np.float64
    )
    weights = np.asarray(geometry.solid_face_length_um, dtype=np.float64)
    distance = np.linalg.norm(coordinates - center[None, :], axis=1)
    scale = max(influence_radius, float(np.max(np.abs(center))), 1.0)
    tolerance = 128.0 * np.finfo(np.float64).eps * scale
    influence = distance <= influence_radius + tolerance
    if not np.any(influence):
        raise ValueError("The v16 wall has no elements inside the saved influence circle.")

    field_values = evaluate_random_fourier_field(
        coordinates[influence], center, wavevectors, phases
    )
    local_selected, threshold = _weighted_positive_flags(
        field_values,
        weights[influence],
        np.asarray(geometry.solid_face_ring_index, dtype=np.int64)[influence],
        requested_positive_fraction,
    )
    positive = np.zeros(coordinates.shape[0], dtype=bool)
    positive[np.flatnonzero(influence)[local_selected]] = True
    influence_length = float(np.sum(weights[influence]))
    positive_length = float(np.sum(weights[positive]))
    total_length = float(np.sum(weights))
    patch_count = _continuous_patch_count(
        positive,
        np.asarray(geometry.solid_face_ring_index, dtype=np.int64),
        int(geometry.full_boundary_segment_count),
    )

    payload = {
        "target_geometry_schema": np.asarray(TARGET_SCHEMA),
        "continuous_geometry_hash_sha256": np.asarray(
            geometry.geometry_hash_sha256
        ),
        "legacy_source_npz_path": np.asarray(str(source)),
        "legacy_source_sha256": np.asarray(_file_sha256(source)),
        "migration_semantics": np.asarray(
            "Saved physical influence circle and random Fourier realization "
            "re-evaluated on continuous v16 solid-wall arclength."
        ),
        "anchor_candidate_id": np.asarray(anchor),
        "influence_center_xz_um": center,
        "influence_radius_um": np.asarray(influence_radius, dtype=np.float64),
        "requested_influence_wall_area_fraction": np.asarray(
            requested_influence_fraction, dtype=np.float64
        ),
        "achieved_influence_wall_length_fraction": np.asarray(
            influence_length / total_length, dtype=np.float64
        ),
        "requested_positive_wall_fraction_within_influence": np.asarray(
            requested_positive_fraction, dtype=np.float64
        ),
        "achieved_positive_wall_fraction_within_influence": np.asarray(
            positive_length / influence_length, dtype=np.float64
        ),
        "target_positive_network_wall_length_fraction": np.asarray(
            positive_length / total_length, dtype=np.float64
        ),
        "target_correlation_length_um": np.asarray(
            correlation_length, dtype=np.float64
        ),
        "random_seed": np.asarray(random_seed, dtype=np.int64),
        "random_field_modes": np.asarray(random_modes, dtype=np.int64),
        "random_field_algorithm": np.asarray(random_algorithm),
        "random_wavevectors_um_inv": np.ascontiguousarray(wavevectors),
        "random_phases_rad": np.ascontiguousarray(phases),
        "random_field_threshold": np.asarray(threshold, dtype=np.float64),
        "wall_ring_index": np.asarray(
            geometry.solid_face_ring_index, dtype=np.int64
        ),
        "wall_arclength_start_um": np.asarray(
            geometry.solid_face_arclength_start_um, dtype=np.float64
        ),
        "wall_arclength_end_um": np.asarray(
            geometry.solid_face_arclength_end_um, dtype=np.float64
        ),
        "boundary_ring_length_um": np.asarray(
            geometry.boundary_ring_length_um, dtype=np.float64
        ),
        "wall_start_xz_um": np.asarray(
            geometry.solid_face_start_xz_um, dtype=np.float64
        ),
        "wall_end_xz_um": np.asarray(
            geometry.solid_face_end_xz_um, dtype=np.float64
        ),
        "wall_center_xz_um": coordinates,
        "wall_inward_normal_xz": -np.asarray(
            geometry.solid_face_outward_normal_xz, dtype=np.float64
        ),
        "wall_length_um": weights,
        "wall_inside_influence": np.ascontiguousarray(influence),
        "wall_target_positive": np.ascontiguousarray(positive),
        "influence_wall_random_field": np.ascontiguousarray(field_values),
        "influence_wall_ring_index": np.asarray(
            geometry.solid_face_ring_index, dtype=np.int64
        )[influence],
        "patch_count": np.asarray(patch_count, dtype=np.int64),
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, **payload)
    return destination


def _weighted_positive_flags(
    values: np.ndarray,
    weights: np.ndarray,
    stable_indices: np.ndarray,
    requested_fraction: float,
) -> tuple[np.ndarray, float]:
    if requested_fraction == 0.0:
        return np.zeros(values.size, dtype=bool), float("inf")
    if requested_fraction == 1.0:
        return np.ones(values.size, dtype=bool), float("-inf")
    order = np.lexsort((stable_indices, -values))
    cumulative = np.cumsum(weights[order])
    requested = requested_fraction * float(np.sum(weights))
    count = int(np.argmin(np.abs(cumulative - requested))) + 1
    selected = np.zeros(values.size, dtype=bool)
    selected[order[:count]] = True
    return selected, float(np.min(values[selected]))


def _continuous_patch_count(
    positive: np.ndarray,
    ring_indices: np.ndarray,
    ring_size: int,
) -> int:
    count = 0
    for index, selected in enumerate(positive):
        if not bool(selected):
            continue
        previous = index - 1
        joined = (
            previous >= 0
            and bool(positive[previous])
            and int(ring_indices[index]) - int(ring_indices[previous]) == 1
        )
        if not joined:
            count += 1
    if (
        count > 1
        and positive.size
        and bool(positive[0])
        and bool(positive[-1])
        and (
            int(ring_indices[0]) == 0
            and int(ring_indices[-1]) == int(ring_size) - 1
        )
    ):
        count -= 1
    return count


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()
