"""
Input loading for the field-based microbubble pipeline.
"""

from __future__ import annotations

import math
from pathlib import Path
import re
from typing import Literal

from ulm_vascular_model_generator.utils.io.vessel_transport_export import (
    read_vessel_transport_metadata,
    read_vessel_transport_npz,
)
from ..core.types import PhysicsInput

VascularGeometryMode = Literal["planar_2d", "volumetric_3d"]
_GEOMETRY_MODES: tuple[VascularGeometryMode, ...] = (
    "planar_2d",
    "volumetric_3d",
)


def geometry_mode_from_model_directory_name(
    model_dir: str | Path,
) -> VascularGeometryMode:
    """Read the solver dimension from one generator run-directory name."""

    directory = Path(model_dir).expanduser().resolve()
    name = directory.name.lower()
    matches = [
        mode
        for mode in _GEOMETRY_MODES
        if re.search(rf"(?:^|_){re.escape(mode)}(?:_|$)", name)
    ]
    if len(matches) != 1:
        raise ValueError(
            "The vascular result directory name must contain exactly one "
            "geometry marker, '_planar_2d_' or '_volumetric_3d_': "
            f"{directory}"
        )
    return matches[0]


def validate_physics_input_geometry(
    model_dir: str | Path,
    physics_input: PhysicsInput,
    *,
    expected_mode: VascularGeometryMode | None = None,
) -> VascularGeometryMode:
    """Cross-check directory marker, transport metadata, and coordinates."""

    directory = Path(model_dir).expanduser().resolve()
    directory_mode = geometry_mode_from_model_directory_name(directory)
    metadata_mode = str(
        physics_input.vessel_metadata.get("geometry_mode", "")
    ).strip().lower()
    if metadata_mode not in _GEOMETRY_MODES:
        raise ValueError(
            "The vessel transport metadata must declare geometry_mode as "
            f"'planar_2d' or 'volumetric_3d': {physics_input.vessel_data_path}"
        )
    if metadata_mode != directory_mode:
        raise ValueError(
            "Vascular geometry conflict: directory name declares "
            f"{directory_mode!r}, but transport metadata declares "
            f"{metadata_mode!r}: {directory}"
        )
    if expected_mode is not None and directory_mode != expected_mode:
        raise ValueError(
            f"Expected a {expected_mode} vascular result, got "
            f"{directory_mode}: {directory}"
        )
    if not physics_input.vessels:
        raise ValueError("The vessel transport export contains no vessel segments.")

    y_coordinates: list[float] = []
    coordinate_scale = 1.0
    for vessel in physics_input.vessels:
        proximal = tuple(float(value) for value in vessel.x_p)
        distal = tuple(float(value) for value in vessel.x_d)
        if len(proximal) != 3 or len(distal) != 3:
            raise ValueError("Every vessel endpoint must have exactly three coordinates.")
        if not all(math.isfinite(value) for value in (*proximal, *distal)):
            raise ValueError("Vessel endpoints must contain only finite coordinates.")
        coordinate_scale = max(
            coordinate_scale,
            *(abs(value) for value in (*proximal, *distal)),
        )
        y_coordinates.extend((proximal[1], distal[1]))
        if directory_mode == "planar_2d":
            segment_length = math.hypot(
                distal[0] - proximal[0],
                distal[2] - proximal[2],
            )
            length_description = "projected X-Z"
        else:
            segment_length = math.sqrt(
                sum((distal[index] - proximal[index]) ** 2 for index in range(3))
            )
            length_description = "three-dimensional"
        if segment_length <= 1.0e-12 * coordinate_scale:
            raise ValueError(
                f"Vessel {int(vessel.vid)} has zero {length_description} length."
            )

    if directory_mode == "planar_2d":
        planar_tolerance = 1.0e-9 * coordinate_scale
        y_span = max(y_coordinates) - min(y_coordinates)
        if y_span > planar_tolerance:
            raise ValueError(
                "A directory marked planar_2d contains non-coplanar vessel "
                f"centrelines: Y span is {y_span:.9g} um."
            )
    return directory_mode


def planar_inlet_area_equivalent_thickness_um(vessels) -> float:
    """Return the constant 2-D depth whose inlet area equals circular inlet area.

    For one root of radius ``R`` this is ``pi*R/2``.  Using this depth in
    ``q_2D = Q_3D / depth`` preserves the exported total volumetric flow after
    extrusion and makes the planar inlet mean velocity equal to the source
    circular-tube mean velocity.  A single constant depth also preserves flow
    conservation throughout the planar network.
    """

    roots = [vessel for vessel in vessels if int(vessel.parent_id) < 0]
    if not roots:
        raise ValueError(
            "Cannot derive a planar equivalent thickness without a root inlet."
        )
    radii = [float(vessel.radius) for vessel in roots]
    if any(not math.isfinite(radius) or radius <= 0.0 for radius in radii):
        raise ValueError("Every root inlet radius must be finite and positive.")
    circular_area_um2 = sum(math.pi * radius * radius for radius in radii)
    planar_width_um = sum(2.0 * radius for radius in radii)
    thickness_um = circular_area_um2 / planar_width_um
    if not math.isfinite(thickness_um) or thickness_um <= 0.0:
        raise ValueError("Derived planar equivalent thickness is invalid.")
    return float(thickness_um)


def load_physics_input(
    model_dir: str | Path,
    *,
    planar_extrusion_depth_um: float | None = None,
) -> PhysicsInput:
    """
    Load the matching vessel files from one model result directory.
    """

    model_dir = Path(model_dir).resolve()
    if not model_dir.is_dir():
        raise FileNotFoundError(f"Cannot find the vascular model directory: {model_dir}")

    swc_files = list(model_dir.glob("*.swc"))
    vessel_data_files = list(model_dir.glob("*.vessels.npz"))
    if len(swc_files) != 1 or len(vessel_data_files) != 1:
        raise ValueError(f"{model_dir} must contain exactly one .swc file and one .vessels.npz file.")

    swc_path = swc_files[0]
    vessel_data_path = vessel_data_files[0]
    if vessel_data_path.name.removesuffix(".vessels.npz") != swc_path.stem:
        raise ValueError(f"The .swc and .vessels.npz files in {model_dir} do not have the same name.")
    
    physics_input = PhysicsInput(
        swc_path=swc_path,
        vessel_data_path=vessel_data_path,
        vessels=read_vessel_transport_npz(vessel_data_path),
        vessel_metadata=read_vessel_transport_metadata(vessel_data_path),
    )
    if planar_extrusion_depth_um is not None:
        _adapt_planar_flux_for_extruded_downstream_model(
            physics_input,
            planar_extrusion_depth_um,
        )
    return physics_input


def _adapt_planar_flux_for_extruded_downstream_model(
    physics_input: PhysicsInput,
    extrusion_depth_um: float,
) -> None:
    """Adapt true 2-D generator flux to the downstream extruded representation.

    The vascular generator remains strictly two-dimensional and exports
    ``q_2D`` in um^2/s. The current trajectory pipeline represents particle
    concentration in an explicitly extruded planar domain, so its existing
    internal storage expects ``Q = q_2D * depth`` and divides by the same depth
    when imposing the CFD boundary flux. This adapter is the sole dimensional
    boundary between those independent models.
    """

    mode = str(physics_input.vessel_metadata.get("geometry_mode", "")).lower()
    quantity = str(
        physics_input.vessel_metadata.get("flow_quantity", "")
    ).lower()
    if mode != "planar_2d" or quantity != "planar_flux_per_unit_depth":
        return
    depth = float(extrusion_depth_um)
    if not math.isfinite(depth) or depth <= 0.0:
        raise ValueError(
            "planar_extrusion_depth_um must be finite and positive when loading "
            "a true planar-flux vascular model."
        )
    for vessel in physics_input.vessels:
        vessel.flow_rate *= depth
        vessel.prescribed_outflow *= depth
    physics_input.vessel_metadata[
        "source_flow_quantity"
    ] = "planar_flux_per_unit_depth"
    physics_input.vessel_metadata[
        "runtime_flow_quantity"
    ] = "extruded_volume_flow"
    physics_input.vessel_metadata[
        "runtime_planar_extrusion_depth_um"
    ] = f"{depth:.17g}"
