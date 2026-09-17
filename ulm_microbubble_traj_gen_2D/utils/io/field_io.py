"""Saving field-based simulation outputs."""

from __future__ import annotations

from dataclasses import asdict, dataclass, is_dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from ..core.config import PhysicsConfig
from ..core.types import (
    FiniteElementVelocityField,
    FlowField,
    GridDomain,
    HybridVelocityField,
    ParticleTrajectories,
    RasterizedVessels,
)
from ..flow.hybrid_velocity import hybrid_region_map, validate_hybrid_velocity_field
from ..flow.hybrid_velocity import sample_finite_element_velocity
from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from ..molecular.molecular_target_field import MolecularTargetField
from ..particles.red_blood_cell_transport import RedBloodCellNetworkState

_FIELD_REUSE_CONTRACT_KEY = "_resolved_field_reuse_contract"
_FIELD_REUSE_CONTRACT_SCHEMA = "strict_flow_field_reuse_v2"
_REUSABLE_FIELD_SCHEMA = "v19_hybrid_fem_grid_local_shear"
_LEGACY_REUSABLE_FIELD_SCHEMA = "v18_hybrid_fem_grid_continuous_wall"
_DOLFINX_ACCEPTANCE_SCHEMAS = {
    "dolfinx_fixed_flow_petsc_stokes_v3",
    "dolfinx_fixed_flow_scipy_stokes_v3",
}


class FieldReuseValidationError(ValueError):
    """Raised when a saved flow field cannot be proven compatible."""


@dataclass(frozen=True)
class ReusableFlowField:
    """A fully validated saved flow field and its provenance."""

    flow: FlowField
    field_npz_path: Path
    run_config_path: Path


def build_field_reuse_config_contract(cfg: PhysicsConfig) -> dict[str, Any]:
    """Return the effective input/domain/field contract for strict reuse.

    The normal YAML sections alone cannot represent command-line quick-test
    overrides.  Persisting this resolved contract makes those runs unambiguous
    while the loader can still validate older, non-quick-test result folders.
    """

    return {
        "schema_version": _FIELD_REUSE_CONTRACT_SCHEMA,
        "input": {
            "model_dir": str(Path(cfg.model_dir).resolve()),
        },
        "domain": _config_object_mapping(cfg.domain, "domain"),
        "field": _config_object_mapping(cfg.field, "field"),
        "quick_test": bool(cfg.quick_test),
    }


def load_reusable_flow_field(
    source: str | Path,
    *,
    cfg: PhysicsConfig,
    domain: GridDomain,
    raster: RasterizedVessels,
    continuous_geometry: ContinuousVesselGeometry,
) -> ReusableFlowField:
    """Load a saved accepted field only after strict compatibility checks.

    A result directory or its ``velocity_and_wall_shear_field.npz`` may be
    supplied.  Configuration, grid, rasterized material fields, continuous
    geometry identity, field schema, convergence state, and array shapes are
    all checked before a :class:`FlowField` is returned.  There is deliberately
    no permissive or best-effort fallback.
    """

    field_path, run_config_path = _resolve_reusable_field_paths(source)
    current_contract = build_field_reuse_config_contract(cfg)
    source_contract = _load_source_reuse_contract(run_config_path)
    _require_matching_reuse_contract(current_contract, source_contract)

    try:
        archive = np.load(field_path, allow_pickle=False)
    except (OSError, ValueError) as exc:
        raise FieldReuseValidationError(
            f"Could not read reusable field {field_path}: {exc}"
        ) from exc

    with archive:
        schema = _required_scalar_string(archive, "field_schema_version")
        if schema not in (
            _REUSABLE_FIELD_SCHEMA,
            _LEGACY_REUSABLE_FIELD_SCHEMA,
        ):
            raise FieldReuseValidationError(
                "Reusable field schema mismatch: expected one of "
                f"{(_REUSABLE_FIELD_SCHEMA, _LEGACY_REUSABLE_FIELD_SCHEMA)!r}, "
                f"found {schema!r}."
            )

        _validate_saved_domain(archive, domain)
        _validate_saved_raster(archive, raster)

        geometry_schema = _required_scalar_string(archive, "continuous_geometry_schema")
        current_geometry_schema = str(continuous_geometry.geometry_schema)
        if geometry_schema != current_geometry_schema:
            raise FieldReuseValidationError(
                "Continuous-geometry schema mismatch: saved "
                f"{geometry_schema!r}, current {current_geometry_schema!r}."
            )
        geometry_hash = _required_scalar_string(
            archive, "continuous_geometry_hash_sha256"
        )
        current_geometry_hash = str(continuous_geometry.geometry_hash_sha256)
        if geometry_hash != current_geometry_hash:
            raise FieldReuseValidationError(
                "Continuous-geometry hash mismatch: saved "
                f"{geometry_hash!r}, current {current_geometry_hash!r}."
            )

        solver_metadata = _load_solver_metadata(archive)
        _validate_saved_solver_metadata(solver_metadata)
        flow_arrays = _load_and_validate_flow_arrays(
            archive,
            domain.shape,
            require_local_shear=(schema == _REUSABLE_FIELD_SCHEMA),
        )
        hybrid_velocity = _load_hybrid_velocity_field(archive)
        if "local_shear_stress_pa" not in flow_arrays:
            flow_arrays["local_shear_stress_pa"] = (
                _reconstruct_local_shear_stress_from_fem(
                    domain=domain,
                    raster=raster,
                    hybrid_velocity=hybrid_velocity,
                )
            )
            solver_metadata = dict(solver_metadata)
            solver_metadata["local_shear_definition"] = (
                "mu*sqrt(2*D:D)_newtonian_viscous_stress_magnitude"
            )
            solver_metadata["local_shear_gradient_sampling"] = (
                "reconstructed_from_exported_fem_velocity_at_cartesian_lumen_centres"
            )

    return ReusableFlowField(
        flow=FlowField(
            solver_metadata=solver_metadata,
            hybrid_velocity=hybrid_velocity,
            **flow_arrays,
        ),
        field_npz_path=field_path,
        run_config_path=run_config_path,
    )


def _config_object_mapping(value: object, section_name: str) -> dict[str, Any]:
    if is_dataclass(value):
        mapping = asdict(value)
    elif hasattr(value, "__dict__"):
        mapping = dict(vars(value))
    else:
        raise FieldReuseValidationError(
            f"Current {section_name} configuration is not serializable."
        )
    return {str(key): _plain_config_value(item) for key, item in mapping.items()}


def _plain_config_value(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value.resolve())
    if isinstance(value, tuple):
        return [_plain_config_value(item) for item in value]
    if isinstance(value, list):
        return [_plain_config_value(item) for item in value]
    if isinstance(value, Mapping):
        return {str(key): _plain_config_value(item) for key, item in value.items()}
    if isinstance(value, np.generic):
        return value.item()
    return value


def _resolve_reusable_field_paths(source: str | Path) -> tuple[Path, Path]:
    source_path = Path(source).expanduser().resolve()
    field_path = (
        source_path / "velocity_and_wall_shear_field.npz"
        if source_path.is_dir()
        else source_path
    )
    if not field_path.is_file():
        raise FieldReuseValidationError(f"Reusable field does not exist: {field_path}")
    if field_path.suffix.lower() != ".npz":
        raise FieldReuseValidationError(
            f"Reusable field must be an NPZ file, not {field_path}."
        )
    run_config_path = field_path.parent / "run_config.yaml"
    if not run_config_path.is_file():
        raise FieldReuseValidationError(
            "Strict field reuse requires the sibling source run_config.yaml: "
            f"{run_config_path}"
        )
    return field_path, run_config_path


def _load_source_reuse_contract(run_config_path: Path) -> dict[str, Any]:
    try:
        raw = yaml.safe_load(run_config_path.read_text(encoding="utf-8-sig"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise FieldReuseValidationError(
            f"Could not read source run configuration {run_config_path}: {exc}"
        ) from exc
    if not isinstance(raw, dict):
        raise FieldReuseValidationError(
            f"Source run configuration {run_config_path} must be a mapping."
        )
    for section in ("input", "domain", "field"):
        if not isinstance(raw.get(section), dict):
            raise FieldReuseValidationError(
                f"Source run configuration is missing mapping section {section!r}."
            )

    embedded = raw.get(_FIELD_REUSE_CONTRACT_KEY)
    if not isinstance(embedded, dict):
        raise FieldReuseValidationError(
            f"{_FIELD_REUSE_CONTRACT_KEY} must be present as a mapping."
        )
    return _normalize_reuse_contract(embedded, run_config_path)


def _normalize_reuse_contract(
    contract: Mapping[str, Any], source_path: Path
) -> dict[str, Any]:
    if contract.get("schema_version") != _FIELD_REUSE_CONTRACT_SCHEMA:
        raise FieldReuseValidationError(
            "Unsupported resolved field-reuse contract schema in "
            f"{source_path}: {contract.get('schema_version')!r}."
        )
    input_mapping = contract.get("input")
    domain_mapping = contract.get("domain")
    field_mapping = contract.get("field")
    if not isinstance(input_mapping, Mapping):
        raise FieldReuseValidationError(
            "Resolved reuse contract input must be a mapping."
        )
    if not isinstance(domain_mapping, Mapping):
        raise FieldReuseValidationError(
            "Resolved reuse contract domain must be a mapping."
        )
    if not isinstance(field_mapping, Mapping):
        raise FieldReuseValidationError(
            "Resolved reuse contract field must be a mapping."
        )
    try:
        model_dir = str(Path(input_mapping["model_dir"]).expanduser().resolve())
    except (KeyError, TypeError, ValueError) as exc:
        raise FieldReuseValidationError(
            "Resolved reuse contract model directory is missing or invalid."
        ) from exc
    return {
        "schema_version": _FIELD_REUSE_CONTRACT_SCHEMA,
        "input": {"model_dir": model_dir},
        "domain": {
            str(key): _plain_config_value(value)
            for key, value in domain_mapping.items()
        },
        "field": {
            str(key): _plain_config_value(value) for key, value in field_mapping.items()
        },
        "quick_test": bool(contract.get("quick_test", False)),
    }


def _require_matching_reuse_contract(
    current: Mapping[str, Any], source: Mapping[str, Any]
) -> None:
    for section in ("input", "domain", "field", "quick_test"):
        if source.get(section) != current.get(section):
            raise FieldReuseValidationError(
                "Reusable field configuration mismatch in "
                f"{section!r}: source={source.get(section)!r}, "
                f"current={current.get(section)!r}."
            )


def _required_scalar_string(archive: Any, key: str) -> str:
    if key not in archive.files:
        raise FieldReuseValidationError(
            f"Reusable field is missing required array {key!r}."
        )
    value = np.asarray(archive[key])
    if value.size != 1:
        raise FieldReuseValidationError(
            f"Reusable field array {key!r} must contain one value."
        )
    return str(value.reshape(-1)[0])


def _validate_saved_domain(archive: Any, domain: GridDomain) -> None:
    expected = {
        "origin_um": np.asarray(domain.origin_um),
        "spacing_um": np.asarray([domain.spacing_um], dtype=float),
        "fixed_y_um": np.asarray([domain.fixed_y_um], dtype=float),
        "shape": np.asarray(domain.shape, dtype=int),
        "x_coordinates_um": np.asarray(domain.x_coordinates_um),
        "z_coordinates_um": np.asarray(domain.z_coordinates_um),
    }
    for key, value in expected.items():
        _require_exact_array(archive, key, value, category="domain")


def _validate_saved_raster(archive: Any, raster: RasterizedVessels) -> None:
    expected = {
        "lumen_mask": raster.lumen_mask,
        "lumen_fraction": raster.lumen_fraction,
        "wall_mask": raster.wall_mask,
        "vessel_id": raster.vessel_id,
        "radius_um": raster.radius_um,
        "flow_rate_um3_s": raster.flow_rate_um3_s,
        "q2d_flow_um2_s": raster.q2d_flow_um2_s,
        "viscosity_mpas": raster.viscosity_mpas,
        "distance_to_wall_um": raster.distance_to_wall_um,
        "wall_normal_xz": raster.wall_normal_xz,
    }
    for key, value in expected.items():
        _require_exact_array(archive, key, np.asarray(value), category="raster")

    current_junction = raster.junction_core_mask
    saved_has_junction = "junction_core_mask" in archive.files
    if (current_junction is None) != (not saved_has_junction):
        raise FieldReuseValidationError(
            "Reusable field raster mismatch: junction_core_mask presence differs."
        )
    if current_junction is not None:
        _require_exact_array(
            archive,
            "junction_core_mask",
            np.asarray(current_junction),
            category="raster",
        )


def _require_exact_array(
    archive: Any, key: str, expected: np.ndarray, *, category: str
) -> None:
    if key not in archive.files:
        raise FieldReuseValidationError(
            f"Reusable field is missing required {category} array {key!r}."
        )
    actual = np.asarray(archive[key])
    expected = np.asarray(expected)
    if actual.shape != expected.shape:
        raise FieldReuseValidationError(
            f"Reusable field {category} array {key!r} shape mismatch: "
            f"saved={actual.shape}, current={expected.shape}."
        )
    if actual.dtype != expected.dtype:
        raise FieldReuseValidationError(
            f"Reusable field {category} array {key!r} dtype mismatch: "
            f"saved={actual.dtype}, current={expected.dtype}."
        )
    if np.issubdtype(actual.dtype, np.inexact):
        equal = np.array_equal(actual, expected, equal_nan=True)
    else:
        equal = np.array_equal(actual, expected)
    if not equal:
        raise FieldReuseValidationError(
            f"Reusable field {category} array {key!r} differs from the current run."
        )


def _load_solver_metadata(archive: Any) -> dict[str, object]:
    metadata: dict[str, object] = {}
    for key in archive.files:
        if not key.startswith("solver_"):
            continue
        raw = np.asarray(archive[key])
        if raw.ndim == 0:
            value: object = raw.item()
        elif raw.shape[0] == 1:
            first = raw[0]
            value = (
                first.item()
                if np.asarray(first).ndim == 0
                else np.asarray(first).tolist()
            )
        else:
            value = raw.tolist()
        metadata[key[len("solver_") :]] = value
    return metadata


def _validate_saved_solver_metadata(metadata: Mapping[str, object]) -> None:
    if metadata.get("physical_converged") is not True:
        raise FieldReuseValidationError(
            "Reusable field was not certified as physically converged."
        )
    if metadata.get("solver_mode") != "dolfinx_stokes_gmsh_2d":
        raise FieldReuseValidationError(
            "Reusable field was not generated by the DOLFINx Stokes solver."
        )
    saved_acceptance_schema = metadata.get("physical_acceptance_schema")
    if saved_acceptance_schema not in _DOLFINX_ACCEPTANCE_SCHEMAS:
        raise FieldReuseValidationError(
            "Reusable field has an unsupported DOLFINx acceptance schema."
        )
    try:
        converged_reason = int(
            metadata.get(
                "linear_solver_converged_reason",
                metadata.get("dolfinx_ksp_converged_reason", 0),
            )
        )
        mesh_cell_count = int(metadata["mesh_cell_count"])
    except (KeyError, TypeError, ValueError) as exc:
        raise FieldReuseValidationError(
            "Reusable DOLFINx field is missing valid linear-solver or mesh "
            "acceptance metadata."
        ) from exc
    if converged_reason <= 0 or mesh_cell_count <= 0:
        raise FieldReuseValidationError(
            "Reusable DOLFINx field was not certified by a converged "
            "linear solve on a non-empty boundary-fitted mesh."
        )


def _load_and_validate_flow_arrays(
    archive: Any,
    domain_shape: tuple[int, int],
    *,
    require_local_shear: bool = True,
) -> dict[str, np.ndarray]:
    nx, nz = (int(domain_shape[0]), int(domain_shape[1]))
    expected_shapes: dict[str, tuple[int, ...] | None] = {
        "velocity_xz_um_s": (nx, nz, 2),
        "speed_um_s": (nx, nz),
        "wall_shear_stress_pa": (nx, nz),
        "local_shear_stress_pa": (nx, nz),
        "initial_velocity_xz_um_s": (nx, nz, 2),
        "initial_speed_um_s": (nx, nz),
        "divergence_s_inv": (nx, nz),
        "wall_penetration_um_s": (nx, nz),
        "pressure": (nx, nz),
        "inlet_label": (nx, nz),
        "outlet_label": (nx, nz),
        "boundary_velocity_xz_um_s": (nx, nz, 2),
        "boundary_normal_xz": (nx, nz, 2),
        "boundary_weight": (nx, nz),
        "boundary_edge_length_um": (nx, nz),
        "open_boundary_flux_um2_s": (nx, nz),
        "face_flux_x_um2_s": (nx + 1, nz),
        "face_flux_z_um2_s": (nx, nz + 1),
        "inlet_target_by_label_um2_s": None,
        "outlet_target_by_label_um2_s": None,
        "inlet_actual_by_label_um2_s": None,
        "outlet_actual_by_label_um2_s": None,
        "open_face_cell_ij": None,
        "open_face_index_ij": None,
        "open_face_axis": None,
        "open_face_normal_xz": None,
        "open_face_center_xz_um": None,
        "open_face_length_um": None,
        "open_face_label": None,
        "open_face_kind": None,
        "open_section_point_xz_um": None,
        "open_section_outward_normal_xz": None,
        "open_section_tangent_xz": None,
        "open_section_half_width_um": None,
        "open_section_label": None,
        "open_section_kind": None,
    }
    arrays: dict[str, np.ndarray] = {}
    for key, expected_shape in expected_shapes.items():
        if key not in archive.files:
            if key == "local_shear_stress_pa" and not require_local_shear:
                continue
            raise FieldReuseValidationError(
                f"Reusable field is missing required flow array {key!r}."
            )
        value = np.asarray(archive[key])
        if expected_shape is not None and value.shape != expected_shape:
            raise FieldReuseValidationError(
                f"Reusable flow array {key!r} has shape {value.shape}; "
                f"expected {expected_shape}."
            )
        if value.dtype.hasobject:
            raise FieldReuseValidationError(
                f"Reusable flow array {key!r} has an object dtype."
            )
        if np.issubdtype(value.dtype, np.inexact) and not np.isfinite(value).all():
            raise FieldReuseValidationError(
                f"Reusable flow array {key!r} contains non-finite values."
            )
        arrays[key] = value

    if (
        "local_shear_stress_pa" in arrays
        and np.any(np.asarray(arrays["local_shear_stress_pa"], dtype=float) < 0.0)
    ):
        raise FieldReuseValidationError(
            "Reusable local_shear_stress_pa contains negative values."
        )
    _require_one_dimensional_label_arrays(arrays)
    _require_consistent_open_boundary_arrays(arrays)
    return arrays


def _reconstruct_local_shear_stress_from_fem(
    *,
    domain: GridDomain,
    raster: RasterizedVessels,
    hybrid_velocity: HybridVelocityField,
) -> np.ndarray:
    """Upgrade a v18 cache using its exact exported FEM velocity polynomial."""

    lumen = np.asarray(raster.lumen_mask, dtype=bool)
    indices = np.argwhere(lumen)
    result = np.zeros(domain.shape, dtype=np.float32)
    if indices.size == 0:
        return result
    points = np.column_stack(
        (
            np.asarray(domain.x_coordinates_um, dtype=float)[indices[:, 0]],
            np.asarray(domain.z_coordinates_um, dtype=float)[indices[:, 1]],
        )
    )
    try:
        _, gradient, _ = sample_finite_element_velocity(
            hybrid_velocity.finite_element,
            points,
            np.ones(points.shape[0], dtype=bool),
            use_numba=True,
        )
    except ValueError as exc:
        raise FieldReuseValidationError(
            "Could not reconstruct local shear stress from the legacy FEM cache: "
            f"{exc}"
        ) from exc

    strain_rate = 0.5 * (gradient + np.swapaxes(gradient, 1, 2))
    strain_invariant = np.einsum("nij,nij->n", strain_rate, strain_rate)
    equivalent_shear_rate = np.sqrt(np.maximum(2.0 * strain_invariant, 0.0))
    viscosity_mpas = np.asarray(raster.viscosity_mpas, dtype=float)[lumen]
    positive = viscosity_mpas[
        np.isfinite(viscosity_mpas) & (viscosity_mpas > 0.0)
    ]
    default_viscosity_mpas = float(np.median(positive)) if positive.size else 3.0
    viscosity_pa_s = (
        np.where(
            np.isfinite(viscosity_mpas) & (viscosity_mpas > 0.0),
            viscosity_mpas,
            default_viscosity_mpas,
        )
        * 1.0e-3
    )
    local_shear = viscosity_pa_s * equivalent_shear_rate
    if np.any(~np.isfinite(local_shear)) or np.any(local_shear < 0.0):
        raise FieldReuseValidationError(
            "Reconstructed local shear stress contains invalid values."
        )
    result[lumen] = local_shear.astype(np.float32)
    return result


def _load_hybrid_velocity_field(archive: Any) -> HybridVelocityField:
    """Load the mandatory finite-element half of the hybrid field."""

    required = (
        "hybrid_finite_element_distance_um",
        "hybrid_regular_grid_distance_um",
        "fem_velocity_degree",
        "fem_cell_vertices_xz_um",
        "fem_polynomial_exponents",
        "fem_velocity_coefficients_um_s",
        "fem_bin_origin_xz_um",
        "fem_bin_size_um",
        "fem_bin_shape",
        "fem_bin_offsets",
        "fem_bin_cell_indices",
    )
    missing = [key for key in required if key not in archive.files]
    if missing:
        raise FieldReuseValidationError(
            "Reusable field is missing mandatory hybrid arrays: "
            + ", ".join(missing)
        )

    def scalar(key):
        value = np.asarray(archive[key]).reshape(-1)
        if value.size != 1:
            raise FieldReuseValidationError(
                f"Reusable hybrid field {key!r} must be scalar."
            )
        return value[0]

    bin_shape_array = np.asarray(archive["fem_bin_shape"], dtype=np.int64)
    if bin_shape_array.shape != (2,):
        raise FieldReuseValidationError("fem_bin_shape must have shape (2,).")
    finite_element = FiniteElementVelocityField(
        degree=int(scalar("fem_velocity_degree")),
        cell_vertices_xz_um=np.ascontiguousarray(
            archive["fem_cell_vertices_xz_um"], dtype=np.float64
        ),
        polynomial_exponents=np.ascontiguousarray(
            archive["fem_polynomial_exponents"], dtype=np.int16
        ),
        velocity_coefficients_um_s=np.ascontiguousarray(
            archive["fem_velocity_coefficients_um_s"], dtype=np.float64
        ),
        bin_origin_xz_um=np.ascontiguousarray(
            archive["fem_bin_origin_xz_um"], dtype=np.float64
        ),
        bin_size_um=float(scalar("fem_bin_size_um")),
        bin_shape=(int(bin_shape_array[0]), int(bin_shape_array[1])),
        bin_offsets=np.ascontiguousarray(
            archive["fem_bin_offsets"], dtype=np.int64
        ),
        bin_cell_indices=np.ascontiguousarray(
            archive["fem_bin_cell_indices"], dtype=np.int32
        ),
    )
    hybrid = HybridVelocityField(
        finite_element=finite_element,
        finite_element_distance_um=float(
            scalar("hybrid_finite_element_distance_um")
        ),
        regular_grid_distance_um=float(
            scalar("hybrid_regular_grid_distance_um")
        ),
    )
    try:
        validate_hybrid_velocity_field(hybrid)
    except ValueError as exc:
        raise FieldReuseValidationError(
            f"Reusable hybrid velocity field is invalid: {exc}"
        ) from exc
    return hybrid


def _require_one_dimensional_label_arrays(arrays: Mapping[str, np.ndarray]) -> None:
    for target_name, actual_name in (
        ("inlet_target_by_label_um2_s", "inlet_actual_by_label_um2_s"),
        ("outlet_target_by_label_um2_s", "outlet_actual_by_label_um2_s"),
    ):
        target = arrays[target_name]
        actual = arrays[actual_name]
        if target.ndim != 1 or actual.shape != target.shape:
            raise FieldReuseValidationError(
                f"Reusable flow label arrays {target_name!r}/{actual_name!r} "
                "must be equal-length one-dimensional arrays."
            )


def _require_consistent_open_boundary_arrays(arrays: Mapping[str, np.ndarray]) -> None:
    open_face_count = int(arrays["open_face_axis"].shape[0])
    open_face_shapes = {
        "open_face_cell_ij": (open_face_count, 2),
        "open_face_index_ij": (open_face_count, 2),
        "open_face_axis": (open_face_count,),
        "open_face_normal_xz": (open_face_count, 2),
        "open_face_center_xz_um": (open_face_count, 2),
        "open_face_length_um": (open_face_count,),
        "open_face_label": (open_face_count,),
        "open_face_kind": (open_face_count,),
    }
    open_section_count = int(arrays["open_section_kind"].shape[0])
    open_section_shapes = {
        "open_section_point_xz_um": (open_section_count, 2),
        "open_section_outward_normal_xz": (open_section_count, 2),
        "open_section_tangent_xz": (open_section_count, 2),
        "open_section_half_width_um": (open_section_count,),
        "open_section_label": (open_section_count,),
        "open_section_kind": (open_section_count,),
    }
    for key, expected_shape in {**open_face_shapes, **open_section_shapes}.items():
        if arrays[key].shape != expected_shape:
            raise FieldReuseValidationError(
                f"Reusable flow array {key!r} has shape {arrays[key].shape}; "
                f"expected {expected_shape}."
            )


def save_run_config(path: Path, raw_config: dict[str, Any]) -> None:
    """Write the input YAML content into the result directory."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        yaml.safe_dump(raw_config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )


def save_domain_metadata(path: Path, domain: GridDomain, extra: dict[str, Any]) -> None:
    """Save human-readable metadata describing the grid and solver run."""

    data = {
        "origin_um": [float(x) for x in domain.origin_um],
        "spacing_um": float(domain.spacing_um),
        "shape": [int(domain.shape[0]), int(domain.shape[1])],
        "fixed_y_um": float(domain.fixed_y_um),
        **extra,
    }
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


def particle_boundary_metadata(geometry: ContinuousVesselGeometry) -> dict[str, Any]:
    """Return audit fields for the sole continuous particle-wall geometry."""

    return {
        "particle_boundary_geometry_schema": geometry.geometry_schema,
        "particle_true_gap_definition": (
            "g_R=distance_to_continuous_closed_vessel_boundary-radius_um"
        ),
        "particle_true_gap_discretization": (
            "pre_raster_continuous_boundary_brep_with_recorded_curve_tessellation"
        ),
        "particle_hydrodynamic_gap_definition": (
            "max(g_R,radius_um*xi_min)_mobility_coefficients_only"
        ),
        "particle_molecular_capture_gap_definition": (
            "same_unregularized_continuous_g_R_used_by_position_contact"
        ),
        "molecular_target_surface_support": (
            "subset_of_continuous_closed_wall_excluding_anatomical_open_sections"
        ),
        "particle_inlet_admission_geometry": (
            "continuous_wall_clearance_and_finite_radius_anatomical_inlet_section"
        ),
        "particle_continuous_wall_element_count": int(
            geometry.solid_face_length_um.size
        ),
        "particle_open_section_count": int(geometry.open_section_kind.size),
        "particle_inlet_open_section_count": int(
            np.count_nonzero(geometry.open_section_kind < 0)
        ),
        "particle_outlet_open_section_count": int(
            np.count_nonzero(geometry.open_section_kind > 0)
        ),
        "particle_boundary_curve_quad_segs": int(geometry.curve_quad_segs),
        "particle_boundary_maximum_element_length_um": float(
            geometry.maximum_boundary_element_length_um
        ),
        "particle_continuous_geometry_hash_sha256": str(
            geometry.geometry_hash_sha256
        ),
    }


def _required_local_shear_stress(
    flow: FlowField,
    domain_shape: tuple[int, int],
) -> np.ndarray:
    """Validate the authoritative whole-lumen local shear field before saving."""

    if flow.local_shear_stress_pa is None:
        raise ValueError(
            "FlowField.local_shear_stress_pa is required by the v19 field schema."
        )
    values = np.asarray(flow.local_shear_stress_pa, dtype=np.float32)
    if values.shape != tuple(domain_shape):
        raise ValueError(
            "local_shear_stress_pa must match the domain shape; "
            f"expected {tuple(domain_shape)}, got {values.shape}."
        )
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError(
            "local_shear_stress_pa must contain finite, non-negative values."
        )
    return values


def save_field_npz(
    path: Path,
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    *,
    continuous_geometry: ContinuousVesselGeometry,
) -> None:
    """Save grid fields needed to inspect or reuse the background flow."""

    path.parent.mkdir(parents=True, exist_ok=True)
    authoritative_boundary_arrays = (
        "open_face_cell_ij",
        "open_face_index_ij",
        "open_face_axis",
        "open_face_normal_xz",
        "open_face_center_xz_um",
        "open_face_length_um",
        "open_face_label",
        "open_face_kind",
        "open_section_point_xz_um",
        "open_section_outward_normal_xz",
        "open_section_tangent_xz",
        "open_section_half_width_um",
        "open_section_label",
        "open_section_kind",
    )
    hybrid_region, hybrid_weight = hybrid_region_map(
        raster.distance_to_wall_um,
        raster.lumen_mask,
        flow.hybrid_velocity,
    )
    finite_element = flow.hybrid_velocity.finite_element
    data = {
        "field_schema_version": np.asarray(_REUSABLE_FIELD_SCHEMA),
        "origin_um": domain.origin_um,
        "spacing_um": np.asarray([domain.spacing_um], dtype=float),
        "fixed_y_um": np.asarray([domain.fixed_y_um], dtype=float),
        "shape": np.asarray(domain.shape, dtype=int),
        "x_coordinates_um": domain.x_coordinates_um,
        "z_coordinates_um": domain.z_coordinates_um,
        "lumen_mask": raster.lumen_mask,
        "lumen_fraction": raster.lumen_fraction,
        "wall_mask": raster.wall_mask,
        "vessel_id": raster.vessel_id,
        "radius_um": raster.radius_um,
        "flow_rate_um3_s": raster.flow_rate_um3_s,
        "q2d_flow_um2_s": raster.q2d_flow_um2_s,
        "viscosity_mpas": raster.viscosity_mpas,
        "direction_xz": raster.direction_xz,
        "distance_to_centerline_um": raster.distance_to_centerline_um,
        "distance_to_wall_um": raster.distance_to_wall_um,
        "wall_normal_xz": raster.wall_normal_xz,
        "velocity_xz_um_s": flow.velocity_xz_um_s,
        "speed_um_s": flow.speed_um_s,
        "wall_shear_stress_pa": flow.wall_shear_stress_pa,
        "local_shear_stress_pa": _required_local_shear_stress(flow, domain.shape),
        "hybrid_velocity_region": hybrid_region,
        "hybrid_finite_element_weight": hybrid_weight,
        "hybrid_finite_element_distance_um": np.asarray(
            flow.hybrid_velocity.finite_element_distance_um,
            dtype=np.float64,
        ),
        "hybrid_regular_grid_distance_um": np.asarray(
            flow.hybrid_velocity.regular_grid_distance_um,
            dtype=np.float64,
        ),
        "fem_velocity_degree": np.asarray(
            finite_element.degree, dtype=np.int32
        ),
        "fem_cell_vertices_xz_um": finite_element.cell_vertices_xz_um,
        "fem_polynomial_exponents": finite_element.polynomial_exponents,
        "fem_velocity_coefficients_um_s": (
            finite_element.velocity_coefficients_um_s
        ),
        "fem_bin_origin_xz_um": finite_element.bin_origin_xz_um,
        "fem_bin_size_um": np.asarray(
            finite_element.bin_size_um, dtype=np.float64
        ),
        "fem_bin_shape": np.asarray(finite_element.bin_shape, dtype=np.int32),
        "fem_bin_offsets": finite_element.bin_offsets,
        "fem_bin_cell_indices": finite_element.bin_cell_indices,
    }
    if raster.junction_core_mask is not None:
        data["junction_core_mask"] = raster.junction_core_mask
    if flow.divergence_s_inv is not None:
        data["divergence_s_inv"] = flow.divergence_s_inv
    if flow.wall_penetration_um_s is not None:
        data["wall_penetration_um_s"] = flow.wall_penetration_um_s
    if flow.initial_velocity_xz_um_s is not None:
        data["initial_velocity_xz_um_s"] = flow.initial_velocity_xz_um_s
    if flow.initial_speed_um_s is not None:
        data["initial_speed_um_s"] = flow.initial_speed_um_s
    if flow.pressure is not None:
        data["pressure"] = flow.pressure
    if flow.inlet_label is not None:
        data["inlet_label"] = flow.inlet_label
    if flow.outlet_label is not None:
        data["outlet_label"] = flow.outlet_label
    if flow.boundary_velocity_xz_um_s is not None:
        data["boundary_velocity_xz_um_s"] = flow.boundary_velocity_xz_um_s
    if flow.boundary_normal_xz is not None:
        data["boundary_normal_xz"] = flow.boundary_normal_xz
    if flow.boundary_weight is not None:
        data["boundary_weight"] = flow.boundary_weight
    if flow.boundary_edge_length_um is not None:
        data["boundary_edge_length_um"] = flow.boundary_edge_length_um
    if flow.open_boundary_flux_um2_s is not None:
        data["open_boundary_flux_um2_s"] = flow.open_boundary_flux_um2_s
    if flow.face_flux_x_um2_s is not None:
        data["face_flux_x_um2_s"] = flow.face_flux_x_um2_s
    if flow.face_flux_z_um2_s is not None:
        data["face_flux_z_um2_s"] = flow.face_flux_z_um2_s
    if flow.inlet_target_by_label_um2_s is not None:
        data["inlet_target_by_label_um2_s"] = flow.inlet_target_by_label_um2_s
    if flow.outlet_target_by_label_um2_s is not None:
        data["outlet_target_by_label_um2_s"] = flow.outlet_target_by_label_um2_s
    if flow.inlet_actual_by_label_um2_s is not None:
        data["inlet_actual_by_label_um2_s"] = flow.inlet_actual_by_label_um2_s
    if flow.outlet_actual_by_label_um2_s is not None:
        data["outlet_actual_by_label_um2_s"] = flow.outlet_actual_by_label_um2_s
    for name in authoritative_boundary_arrays:
        value = getattr(flow, name)
        if value is not None:
            data[name] = value
    data.update(
        {
            "continuous_geometry_schema": np.asarray(
                continuous_geometry.geometry_schema
            ),
            "continuous_geometry_hash_sha256": np.asarray(
                continuous_geometry.geometry_hash_sha256
            ),
            "continuous_wall_start_xz_um": continuous_geometry.solid_face_start_xz_um,
            "continuous_wall_end_xz_um": continuous_geometry.solid_face_end_xz_um,
            "continuous_wall_inward_normal_xz": continuous_geometry.solid_face_inward_normal_xz,
            "continuous_wall_length_um": continuous_geometry.solid_face_length_um,
            "continuous_wall_ring_index": continuous_geometry.solid_face_ring_index,
            "continuous_wall_arclength_start_um": continuous_geometry.solid_face_arclength_start_um,
            "continuous_wall_arclength_end_um": continuous_geometry.solid_face_arclength_end_um,
            "continuous_boundary_ring_length_um": np.asarray(
                continuous_geometry.boundary_ring_length_um,
                dtype=np.float64,
            ),
            "continuous_open_section_point_xz_um": continuous_geometry.open_section_point_xz_um,
            "continuous_open_section_outward_normal_xz": continuous_geometry.open_section_outward_normal_xz,
            "continuous_open_section_tangent_xz": continuous_geometry.open_section_tangent_xz,
            "continuous_open_section_half_width_um": continuous_geometry.open_section_half_width_um,
            "continuous_open_section_label": continuous_geometry.open_section_label,
            "continuous_open_section_kind": continuous_geometry.open_section_kind,
            "continuous_open_section_vessel_id": continuous_geometry.open_section_vessel_id,
        }
    )
    for key, value in flow.solver_metadata.items():
        # Iteration histories are already one-dimensional numeric sequences.
        # Store them without an extra leading dimension so reuse restores a
        # simple list instead of a nested ``[[...]]`` structure.  Scalar
        # metadata keeps the established one-element-array representation.
        if key.startswith("residual_history_") and key != "residual_history_schema":
            data[f"solver_{key}"] = np.asarray(value)
        else:
            data[f"solver_{key}"] = np.asarray([value])
    np.savez_compressed(path, **data)


def save_molecular_target_npz(path: Path, target_field: MolecularTargetField) -> None:
    """Save the fixed disease ROI intersected with eligible solid vessel walls."""

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **target_field.to_npz_payload())


def save_red_blood_cell_transport_npz(
    path: Path, network: RedBloodCellNetworkState
) -> None:
    """Save per-vessel haematocrit and RBC-flux branch audit state."""

    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **network.to_npz_payload())


def save_trajectories_npz(path: Path, trajectories: ParticleTrajectories) -> None:
    """Save append-only unique-bubble records and sampled field values."""

    path.parent.mkdir(parents=True, exist_ok=True)
    data = {
        "frame_offsets": trajectories.frame_offsets,
        "record_bubble_id": trajectories.bubble_id,
        "record_positions_um": trajectories.positions_um,
        "record_velocities_um_s": trajectories.velocities_um_s,
        "record_wall_shear_stress_pa": trajectories.wall_shear_stress_pa,
        "record_vessel_id": trajectories.vessel_id,
        "record_active": trajectories.active,
        "record_diameter_um": trajectories.diameter_um,
        "record_wall_gap_um": trajectories.wall_gap_um,
        "record_wall_contact": trajectories.wall_contact,
        "record_wall_normal_xz": trajectories.wall_normal_xz,
        "registry_bubble_id": trajectories.registry_bubble_id,
        "registry_diameter_um": trajectories.registry_diameter_um,
        "birth_frame": trajectories.birth_frame,
        "death_frame": trajectories.death_frame,
        "termination_reason": trajectories.termination_reason,
        "active_count_per_frame": trajectories.active_count_per_frame,
        "injected_count_per_frame": trajectories.injected_count_per_frame,
        "terminated_count_per_frame": trajectories.terminated_count_per_frame,
        "metadata_keys": np.asarray(list(trajectories.metadata.keys()), dtype=str),
        "metadata_values": np.asarray(
            [str(v) for v in trajectories.metadata.values()], dtype=str
        ),
    }
    optional_arrays = {
        "record_realized_velocities_um_s": trajectories.realized_velocities_um_s,
        "record_fluid_velocities_um_s": trajectories.fluid_velocities_um_s,
        "record_angular_velocity_rad_s": trajectories.angular_velocity_rad_s,
        "record_rotation_angle_rad": trajectories.rotation_angle_rad,
        "record_collision_force_xz_pn": trajectories.collision_force_xz_pn,
        "record_collision_neighbor_count": trajectories.collision_neighbor_count,
        "record_gap_ratio": trajectories.gap_ratio,
        "record_near_wall_weight": trajectories.near_wall_weight,
        "record_opposite_wall_hydrodynamic_validity_warning": trajectories.two_wall_warning,
        "record_contact_constraint_active": trajectories.contact_constraint_active,
        "record_contact_reaction_force_pn": trajectories.contact_reaction_force_pn,
        "record_contact_free_normal_velocity_um_s": (
            trajectories.contact_free_normal_velocity_um_s
        ),
        "record_contact_constrained_normal_velocity_um_s": (
            trajectories.contact_constrained_normal_velocity_um_s
        ),
        "registry_scheduled_injection_time_s": (
            trajectories.registry_scheduled_injection_time_s
        ),
        "registry_admission_time_s": trajectories.registry_admission_time_s,
        "registry_exit_time_s": trajectories.registry_exit_time_s,
        "registry_inlet_wait_time_s": trajectories.registry_inlet_wait_time_s,
        "record_cardiac_multiplier": trajectories.cardiac_multiplier,
        "cardiac_waveform_time_s": trajectories.cardiac_waveform_time_s,
        "cardiac_waveform_multiplier": trajectories.cardiac_waveform_multiplier,
        "cardiac_path_distance_um": trajectories.cardiac_path_distance_um,
        "cardiac_delay_s": trajectories.cardiac_delay_s,
        "record_bond_count_expected": trajectories.bond_count_expected,
        "record_bond_total_tangential_extension_um": (
            trajectories.bond_total_tangential_extension_um
        ),
        "record_bond_mean_tangential_extension_um": (
            trajectories.bond_mean_tangential_extension_um
        ),
        "record_bond_force_xz_pn": trajectories.bond_force_xz_pn,
        "record_bond_force_tangent_pn": trajectories.bond_force_tangent_pn,
        "record_bond_force_normal_pn": trajectories.bond_force_normal_pn,
        "record_bond_torque_pn_um": trajectories.bond_torque_pn_um,
        "record_single_bond_tension_pn": trajectories.single_bond_tension_pn,
        "record_bond_formation_rate_bonds_s": trajectories.bond_formation_rate_bonds_s,
        "record_bond_dissociation_rate_s_inv": trajectories.bond_dissociation_rate_s_inv,
        "record_target_reaction_area_um2": trajectories.target_reaction_area_um2,
        "record_available_ligand_count": trajectories.available_ligand_count,
        "record_available_target_count": trajectories.available_target_count,
        "record_target_overlap_fraction": trajectories.target_overlap_fraction,
        "registry_final_bond_count_expected": trajectories.registry_final_bond_count_expected,
        "registry_final_bond_total_tangential_extension_um": (
            trajectories.registry_final_bond_total_tangential_extension_um
        ),
        "registry_target_exposure_time_s": (
            trajectories.registry_target_exposure_time_s
        ),
        "registry_target_exposure_event_count": (
            trajectories.registry_target_exposure_event_count
        ),
        "registry_target_reaction_area_time_um2_s": (
            trajectories.registry_target_reaction_area_time_um2_s
        ),
        "registry_target_exposure_right_censored": (
            trajectories.registry_target_exposure_right_censored
        ),
        "registry_target_exposure_quantitative_applicability_fraction": (
            trajectories.registry_target_exposure_quantitative_applicability_fraction
        ),
        "record_red_blood_cell_velocity_xz_um_s": (
            trajectories.red_blood_cell_velocity_xz_um_s
        ),
        "record_red_blood_cell_drift_velocity_xz_um_s": (
            trajectories.red_blood_cell_drift_velocity_xz_um_s
        ),
        "record_red_blood_cell_fick_velocity_xz_um_s": (
            trajectories.red_blood_cell_fick_velocity_xz_um_s
        ),
        "record_red_blood_cell_local_vessel_diameter_um": (
            trajectories.red_blood_cell_local_vessel_diameter_um
        ),
        "record_red_blood_cell_discharge_hematocrit": (
            trajectories.red_blood_cell_discharge_hematocrit
        ),
        "record_red_blood_cell_tube_hematocrit": (
            trajectories.red_blood_cell_tube_hematocrit
        ),
        "record_red_blood_cell_shear_rate_s_inv": (
            trajectories.red_blood_cell_shear_rate_s_inv
        ),
        "record_red_blood_cell_cfl_width_um": trajectories.red_blood_cell_cfl_width_um,
        "record_red_blood_cell_target_gap_um": trajectories.red_blood_cell_target_gap_um,
        "record_red_blood_cell_transverse_diffusivity_um2_s": (
            trajectories.red_blood_cell_transverse_diffusivity_um2_s
        ),
        "record_red_blood_cell_margination_length_um": (
            trajectories.red_blood_cell_margination_length_um
        ),
        "record_red_blood_cell_margination_time_s": (
            trajectories.red_blood_cell_margination_time_s
        ),
        "record_red_blood_cell_scale_activation": (
            trajectories.red_blood_cell_scale_activation
        ),
        "record_red_blood_cell_nearest_wall_unique": (
            trajectories.red_blood_cell_nearest_wall_unique
        ),
        "record_red_blood_cell_hematocrit_in_quantitative_range": (
            trajectories.red_blood_cell_hematocrit_in_quantitative_range
        ),
        "record_red_blood_cell_shear_rate_in_quantitative_range": (
            trajectories.red_blood_cell_shear_rate_in_quantitative_range
        ),
        "record_red_blood_cell_quantitative_applicability": (
            trajectories.red_blood_cell_quantitative_applicability
        ),
        "record_red_blood_cell_transverse_space_valid": (
            trajectories.red_blood_cell_transverse_space_valid
        ),
        "registry_final_vessel_id": trajectories.registry_final_vessel_id,
        "topological_commitment_parent_vessel_id": (
            trajectories.topological_commitment_parent_vessel_id
        ),
        "topological_commitment_child_vessel_id": (
            trajectories.topological_commitment_child_vessel_id
        ),
        "topological_commitment_point_xz_um": (
            trajectories.topological_commitment_point_xz_um
        ),
        "topological_commitment_downstream_normal_xz": (
            trajectories.topological_commitment_downstream_normal_xz
        ),
        "topological_commitment_tangent_xz": (
            trajectories.topological_commitment_tangent_xz
        ),
        "topological_commitment_half_width_um": (
            trajectories.topological_commitment_half_width_um
        ),
        "topological_commitment_transition_end_distance_um": (
            trajectories.topological_commitment_transition_end_distance_um
        ),
        "topological_commitment_distance_um": (
            trajectories.topological_commitment_distance_um
        ),
        "topological_event_bubble_id": trajectories.topological_event_bubble_id,
        "topological_event_time_s": trajectories.topological_event_time_s,
        "topological_event_from_vessel_id": (
            trajectories.topological_event_from_vessel_id
        ),
        "topological_event_to_vessel_id": (trajectories.topological_event_to_vessel_id),
        "topological_event_section_index": (
            trajectories.topological_event_section_index
        ),
        "topological_event_position_xz_um": (
            trajectories.topological_event_position_xz_um
        ),
    }
    data.update(
        {key: value for key, value in optional_arrays.items() if value is not None}
    )
    np.savez_compressed(path, **data)
