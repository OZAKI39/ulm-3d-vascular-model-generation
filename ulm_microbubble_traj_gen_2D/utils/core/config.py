"""
YAML configuration for field-based microbubble trajectory generation.
"""

import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import yaml

PACKAGE_ROOT        = Path(__file__).resolve().parents[2]
PROJECT_ROOT        = PACKAGE_ROOT.parent
CONFIG_DIR          = PACKAGE_ROOT / "configs"
DEFAULT_CONFIG_PATH = CONFIG_DIR / "physics_flow_config.yaml"


@dataclass(frozen=True)
class DomainConfig:
    """Grid settings used to rasterize the vessel network."""

    grid_spacing_um: float
    padding_um: float
    min_lumen_radius_cells: float
    min_resolved_diameter_cells: float
    max_grid_cells: int
    continuous_boundary_maximum_element_length_um: float = 1.0


@dataclass(frozen=True)
class FieldConfig:
    """Velocity-field and wall-shear-field controls."""

    effective_thickness_um: float
    boundary_depth_cells: float
    flux_tolerance: float
    kinematic_viscosity_um2_s: float
    hybrid_finite_element_distance_um: float
    hybrid_transition_width_um: float
    blood_density_kg_m3: float = 1060.0
    gmsh_bulk_mesh_size_um: float = 0.0
    gmsh_wall_mesh_size_um: float = 0.0
    gmsh_wall_refinement_distance_um: float = 0.0
    gmsh_element_order: int = 1
    dolfinx_velocity_degree: int = 2
    dolfinx_pressure_degree: int = 1
    dolfinx_ksp_rtol: float = 1.0e-8


@dataclass(frozen=True)
class ParticleConfig:
    """Continuous-perfusion physics, recording, and geometry controls."""

    inlet_number_concentration_mb_per_ml: float
    n_steps: int
    dt_s: float
    bubble_diameter_min_um: float
    bubble_diameter_max_um: float
    max_unique_bubbles: int = 0
    max_particle_frame_records: int = 5_000_000
    acceleration_backend: str = "auto"
    contact_geometry_tolerance_um: float = 0.001
    contact_max_time_refinements: int = 12
    wall_contact_threshold_um: float = 0.05


@dataclass(frozen=True)
class ParticleDynamicsConfig:
    """Hydrodynamic mobility, collision, and time-integration controls."""

    time_integrator: str = "euler"
    integration_substeps: int = 1
    near_wall_enabled: bool = True
    collisions_enabled: bool = True
    xi_min: float = 1.0e-3
    xi_near: float = 0.1
    xi_far: float = 1.0
    collision_layer_um: float = 0.05
    collision_relaxation_time_s: float = 0.01
    neighbor_search: str = "all_pairs"
    store_full_diagnostics: bool = True
    two_wall_warning_gap_ratio: float = 1.0


@dataclass(frozen=True)
class CardiacPulsatilityConfig:
    """Periodic carrier-flow modulation and pulse-propagation controls."""

    enabled: bool = False
    waveform: str = "synthetic_ecg_envelope"
    bpm: float = 300.0
    pulse_propagation_velocity_um_s: float = 25_000.0
    initial_phase_fraction: float = 0.0
    waveform_samples_per_cycle: int = 2048
    preserve_cycle_mean_flow: bool = True
    modulation_strength: float = 1.0


@dataclass(frozen=True)
class MolecularTargetSelectionConfig:
    """Manual or deterministic synthetic selection of candidate vessel beds."""

    default_mode: str = "manual"
    influence_region_endothelial_wall_area_fraction: float | None = None
    target_positive_wall_fraction_within_influence: float | None = None
    target_correlation_length_um: float | None = None
    random_seed: int = 42
    random_field_modes: int = 512


@dataclass(frozen=True)
class MolecularTargetConfig:
    """Fixed disease-region target field on the physical vessel wall.

    ``target_density_molecules_per_m2`` is a three-dimensional surface density,
    even though bubble centres are advanced in the x-z plane.  The region mask
    only identifies which closed-wall locations express that density; inlet and
    outlet openings are removed later when the wall field is constructed.
    """

    enabled: bool = False
    region_mode: str = "none"
    mask_npz_path: Path | None = None
    mask_array_key: str = "target_mask"
    x_coordinates_key: str = "x_um"
    z_coordinates_key: str = "z_um"
    target_density_molecules_per_m2: float = 0.0


@dataclass(frozen=True)
class MolecularBindingConfig:
    """Deterministic mean-field parameters for one ligand-target pair.

    Unit-bearing names are intentional.  In particular, the association rate
    is an effective two-surface rate in ``m^2 molecule^-1 s^-1``.  A forward
    rate reported in ``s^-1`` by the reference paper is not interchangeable
    with this quantity and is therefore not accepted under an ambiguous alias.
    """

    enabled: bool = False
    model: str = "deterministic_mean_field"
    ligand_density_molecules_per_m2: float = 0.0
    capture_distance_um: float = 0.0
    rest_length_um: float = 0.0
    association_rate_m2_per_molecule_s: float = 0.0
    zero_force_dissociation_rate_s: float = 0.0
    bond_stiffness_pn_per_um: float = 0.0
    reactive_compliance_nm: float = 0.0
    temperature_k: float = 310.0
    mean_field_warning_count: float = 10.0
    bell_exponent_limit: float = 80.0


@dataclass(frozen=True)
class BindingScenarioSweepConfig:
    """Dimensionless no-experiment study fixed before particle transport.

    The contact pilot measures transport exposure only.  Effective association
    rates are defined independently from a requested Damkoehler number and the
    explicit reference time below, so they cannot change with the contact
    outcome of the same run.  Empty defaults keep ordinary generation unchanged.
    """

    enabled: bool = False
    da_on_reference_time_s: float = 0.0
    da_on_levels: tuple[float, ...] = ()
    capture_distance_to_rest_length_ratios: tuple[float, ...] = ()
    target_density_molecules_per_um2_levels: tuple[float, ...] = ()
    ligand_density_molecules_per_um2_levels: tuple[float, ...] = ()


@dataclass(frozen=True)
class PhysicsConfig:
    """Complete configuration object for the field-based pipeline."""

    source_path: Path
    raw: dict
    random_seed: int
    model_dir: Path
    output_dir: Path
    save_run_config: bool
    save_npz: bool
    domain: DomainConfig
    field: FieldConfig
    particles: ParticleConfig
    particle_dynamics: ParticleDynamicsConfig
    cardiac_pulsatility: CardiacPulsatilityConfig
    red_blood_cell_root_discharge_hematocrit: float
    molecular_target_selection: MolecularTargetSelectionConfig
    molecular_target: MolecularTargetConfig
    molecular_binding: MolecularBindingConfig
    binding_scenario_sweep: BindingScenarioSweepConfig
    quick_test: bool


def _float_tuple(value, name):
    """Convert one YAML sequence to an immutable finite float sequence."""

    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} must be a sequence.")
    converted = []
    for index, item in enumerate(value):
        if (
            isinstance(item, bool)
            or not isinstance(item, (int, float))
            or not math.isfinite(item)
        ):
            raise ValueError(f"{name}[{index}] must be a finite number.")
        converted.append(float(item))
    return tuple(converted)


def load_config(path: str | Path = DEFAULT_CONFIG_PATH, quick_test: bool = False) -> PhysicsConfig:
    """
    Load YAML, apply quick-test overrides, and resolve all filesystem paths.
    """
    # =====================================================================
    # ======== Load YAML ========
    # =====================================================================
    source_path = Path(path).resolve()
    raw         = yaml.safe_load(source_path.read_text(encoding="utf-8-sig"))
    if raw is None:
        raw = {}
    elif not isinstance(raw, dict):
        raise ValueError(f"The YAML root in {source_path} must be a mapping.")
    
    global_random_seed_raw = raw.get("random_seed", 42)
    if (isinstance(global_random_seed_raw, bool) 
        or not isinstance(global_random_seed_raw, int)
        or global_random_seed_raw < 0 
        or global_random_seed_raw > (1 << 63) - 1):
        raise ValueError("random_seed must be a non-negative integer representable as int64.")

    section_names = (
        "input", "output", "domain", "field",
        "particles", "particle_dynamics",
        "cardiac_pulsatility", "red_blood_cell_transport",
        "molecular_target_selection", "molecular_target", "molecular_binding",
        "binding_scenario_sweep", "quick_test")
    sections = []
    for name in section_names:
        value = raw.get(name, {})
        if not isinstance(value, dict):
            raise ValueError(f"The YAML section {name!r} must be a mapping.")
        sections.append(dict(value))

    (input_raw, output_raw, domain_raw, field_raw,
     particles_raw, dynamics_raw,
     cardiac_raw, red_blood_cell_raw,
     target_selection_raw, target_raw, binding_raw, 
     sweep_raw, quick_raw) = sections

    removed_field_options = {
        "solver_mode",
        "max_iterations",
        "min_iterations",
        "divergence_tolerance",
        "max_wall_penetration_um_s",
        "momentum_residual_tolerance",
        "viscous_relaxation_dt_s",
        "pressure_solver_tolerance",
        "pressure_solver_max_iterations",
    }
    configured_removed_options = sorted(
        removed_field_options.intersection(field_raw)
        | removed_field_options.intersection(quick_raw)
    )
    if configured_removed_options:
        raise ValueError(
            "Unsupported field options left over from the removed PhiFlow solver: "
            + ", ".join(configured_removed_options)
        )

    # =====================================================================
    # ======== Apply quick-test overrides ========
    # =====================================================================
    quick_enabled = bool(quick_test or quick_raw.get("enabled", False))
    if quick_enabled:
        for key in ("grid_spacing_um",):
            if key in quick_raw:
                domain_raw[key] = quick_raw[key]
        for key in ("n_steps", "dt_s"):
            if key in quick_raw:
                particles_raw[key] = quick_raw[key]
        if "integration_substeps" in quick_raw:
            dynamics_raw["integration_substeps"] = quick_raw["integration_substeps"]

    # =====================================================================
    # ======== Resolve paths and construct configuration dataclasses ========
    # =====================================================================
    model_dir = (PROJECT_ROOT / Path(input_raw["model_dir"])).resolve()
    if not model_dir.is_dir():
        raise FileNotFoundError(f"Cannot find the vascular model directory: {model_dir}")

    timestamp       = datetime.now().strftime(str(output_raw.get("timestamp_format", "%Y%m%d_%H%M%S")))
    results_root    = (PACKAGE_ROOT / Path(output_raw.get("results_dir", "results"))).resolve()
    output_dir      = results_root / timestamp

    integration_substeps_raw = dynamics_raw.get("integration_substeps", 1)
    if (isinstance(integration_substeps_raw, bool) or not isinstance(integration_substeps_raw, int)):
        raise ValueError("particle_dynamics.integration_substeps must be an integer.")

    dynamics = ParticleDynamicsConfig(
        time_integrator=str(dynamics_raw.get("time_integrator", "euler")).lower(),
        integration_substeps=integration_substeps_raw,
        near_wall_enabled=bool(dynamics_raw.get("near_wall_enabled", True)),
        collisions_enabled=bool(dynamics_raw.get("collisions_enabled", True)),
        xi_min=float(dynamics_raw.get("xi_min", 1.0e-3)),
        xi_near=float(dynamics_raw.get("xi_near", 0.1)),
        xi_far=float(dynamics_raw.get("xi_far", 1.0)),
        collision_layer_um=float(dynamics_raw.get("collision_layer_um", 0.05)),
        collision_relaxation_time_s=float(dynamics_raw.get("collision_relaxation_time_s", 0.01)),
        neighbor_search=str(dynamics_raw.get("neighbor_search", "all_pairs")).lower(),
        store_full_diagnostics=bool(dynamics_raw.get("store_full_diagnostics", True)),
        two_wall_warning_gap_ratio=float(dynamics_raw.get("two_wall_warning_gap_ratio", 1.0)),
    )
    _validate_particle_dynamics(dynamics)

    waveform_samples_raw = cardiac_raw.get("waveform_samples_per_cycle", 2048)
    if isinstance(waveform_samples_raw, bool) or not isinstance(waveform_samples_raw, int):
        raise ValueError("cardiac_pulsatility.waveform_samples_per_cycle must be an integer.")
    cardiac = CardiacPulsatilityConfig(
        enabled=bool(cardiac_raw.get("enabled", False)),
        waveform=str(cardiac_raw.get("waveform", "synthetic_ecg_envelope")).lower(),
        bpm=float(cardiac_raw.get("bpm", 300.0)),
        pulse_propagation_velocity_um_s=float(
            cardiac_raw.get("pulse_propagation_velocity_um_s", 25_000.0)
        ),
        initial_phase_fraction=float(cardiac_raw.get("initial_phase_fraction", 0.0)),
        waveform_samples_per_cycle=waveform_samples_raw,
        preserve_cycle_mean_flow=bool(cardiac_raw.get("preserve_cycle_mean_flow", True)),
        modulation_strength=float(cardiac_raw.get("modulation_strength", 1.0)),
    )
    _validate_cardiac_pulsatility(cardiac)

    red_blood_cell_root_discharge_hematocrit = float(
        red_blood_cell_raw.get("root_discharge_hematocrit", 0.0)
    )
    if (
        not math.isfinite(red_blood_cell_root_discharge_hematocrit)
        or not 0.0 <= red_blood_cell_root_discharge_hematocrit <= 1.0
    ):
        raise ValueError(
            "red_blood_cell_transport.root_discharge_hematocrit must be finite "
            "and in [0, 1]."
        )

    target_selection_random_seed_raw = target_selection_raw.get("random_seed", 42)
    if (
        isinstance(target_selection_random_seed_raw, bool)
        or not isinstance(target_selection_random_seed_raw, int)
    ):
        raise ValueError(
            "molecular_target_selection.random_seed must be an integer."
        )
    random_modes_raw = target_selection_raw.get("random_field_modes", 512)
    if isinstance(random_modes_raw, bool) or not isinstance(random_modes_raw, int):
        raise ValueError(
            "molecular_target_selection.random_field_modes must be an integer."
        )
    influence_fraction_raw = target_selection_raw.get(
        "influence_region_endothelial_wall_area_fraction"
    )
    positive_fraction_raw = target_selection_raw.get(
        "target_positive_wall_fraction_within_influence"
    )
    correlation_length_raw = target_selection_raw.get("target_correlation_length_um")
    target_selection = MolecularTargetSelectionConfig(
        default_mode=str(target_selection_raw.get("default_mode", "manual")).lower(),
        influence_region_endothelial_wall_area_fraction=(
            None if influence_fraction_raw is None else float(influence_fraction_raw)
        ),
        target_positive_wall_fraction_within_influence=(
            None if positive_fraction_raw is None else float(positive_fraction_raw)
        ),
        target_correlation_length_um=(
            None if correlation_length_raw is None else float(correlation_length_raw)
        ),
        random_seed=target_selection_random_seed_raw,
        random_field_modes=random_modes_raw,
    )
    _validate_molecular_target_selection(target_selection)

    mask_npz_value = target_raw.get("mask_npz_path")
    mask_npz_path = None
    if mask_npz_value is not None and str(mask_npz_value).strip():
        mask_npz_path = (source_path.parent / Path(mask_npz_value)).resolve()
    target = MolecularTargetConfig(
        enabled=bool(target_raw.get("enabled", False)),
        region_mode=str(target_raw.get("region_mode", "none")).lower(),
        mask_npz_path=mask_npz_path,
        mask_array_key=str(target_raw.get("mask_array_key", "target_mask")),
        x_coordinates_key=str(target_raw.get("x_coordinates_key", "x_um")),
        z_coordinates_key=str(target_raw.get("z_coordinates_key", "z_um")),
        target_density_molecules_per_m2=float(
            target_raw.get("target_density_molecules_per_m2", 0.0)
        ),
    )
    _validate_molecular_target(target)

    ambiguous_binding_fields = {
        "association_rate_s_inv",
        "forward_rate_s_inv",
        "paper_forward_rate_s_inv",
        "encounter_radius_nm",
        "paper_encounter_radius_nm",
    }
    present_ambiguous_fields = sorted(ambiguous_binding_fields.intersection(binding_raw))
    if present_ambiguous_fields:
        fields = ", ".join(present_ambiguous_fields)
        raise ValueError(
            "molecular_binding contains dimensionally incompatible field(s): "
            f"{fields}. A paper forward rate in s^-1 cannot be used as the effective "
            "2D association_rate_m2_per_molecule_s, and the paper's encounter radius "
            "is not molecular_binding.capture_distance_um."
        )
    binding = MolecularBindingConfig(
        enabled=bool(binding_raw.get("enabled", False)),
        model=str(binding_raw.get("model", "deterministic_mean_field")).lower(),
        ligand_density_molecules_per_m2=float(
            binding_raw.get("ligand_density_molecules_per_m2", 0.0)
        ),
        capture_distance_um=float(binding_raw.get("capture_distance_um", 0.0)),
        rest_length_um=float(binding_raw.get("rest_length_um", 0.0)),
        association_rate_m2_per_molecule_s=float(
            binding_raw.get("association_rate_m2_per_molecule_s", 0.0)
        ),
        zero_force_dissociation_rate_s=float(
            binding_raw.get("zero_force_dissociation_rate_s", 0.0)
        ),
        bond_stiffness_pn_per_um=float(binding_raw.get("bond_stiffness_pn_per_um", 0.0)),
        reactive_compliance_nm=float(binding_raw.get("reactive_compliance_nm", 0.0)),
        temperature_k=float(binding_raw.get("temperature_k", 310.0)),
        mean_field_warning_count=float(binding_raw.get("mean_field_warning_count", 10.0)),
        bell_exponent_limit=float(binding_raw.get("bell_exponent_limit", 80.0)),
    )
    _validate_molecular_binding(binding, target)

    sweep = BindingScenarioSweepConfig(
        enabled=bool(sweep_raw.get("enabled", False)),
        da_on_reference_time_s=float(
            sweep_raw.get("da_on_reference_time_s", 0.0)
        ),
        da_on_levels=_float_tuple(
            sweep_raw.get("da_on_levels"), "binding_scenario_sweep.da_on_levels"
        ),
        capture_distance_to_rest_length_ratios=_float_tuple(
            sweep_raw.get("capture_distance_to_rest_length_ratios"),
            "binding_scenario_sweep.capture_distance_to_rest_length_ratios",
        ),
        target_density_molecules_per_um2_levels=_float_tuple(
            sweep_raw.get("target_density_molecules_per_um2_levels"),
            "binding_scenario_sweep.target_density_molecules_per_um2_levels",
        ),
        ligand_density_molecules_per_um2_levels=_float_tuple(
            sweep_raw.get("ligand_density_molecules_per_um2_levels"),
            "binding_scenario_sweep.ligand_density_molecules_per_um2_levels",
        ),
    )
    _validate_binding_scenario_sweep(sweep, target, binding, dynamics)

    contact_max_time_refinements_raw = particles_raw.get(
        "contact_max_time_refinements", 12
    )
    if (
        isinstance(contact_max_time_refinements_raw, bool)
        or not isinstance(contact_max_time_refinements_raw, int)
    ):
        raise ValueError(
            "particles.contact_max_time_refinements must be an integer."
        )

    config = PhysicsConfig(
        source_path=source_path,
        raw=raw,
        random_seed=int(global_random_seed_raw),
        model_dir=model_dir,
        output_dir=output_dir,
        save_run_config=bool(output_raw.get("save_run_config", True)),
        save_npz=bool(output_raw.get("save_npz", True)),
        domain=DomainConfig(
            grid_spacing_um=float(domain_raw["grid_spacing_um"]),
            padding_um=float(domain_raw["padding_um"]),
            min_lumen_radius_cells=float(domain_raw["min_lumen_radius_cells"]),
            min_resolved_diameter_cells=float(domain_raw.get("min_resolved_diameter_cells", 8.0)),
            max_grid_cells=int(domain_raw["max_grid_cells"]),
            continuous_boundary_maximum_element_length_um=float(
                domain_raw.get(
                    "continuous_boundary_maximum_element_length_um",
                    1.0,
                )
            ),
        ),
        field=FieldConfig(
            effective_thickness_um=float(field_raw.get("effective_thickness_um", 10.0)),
            boundary_depth_cells=float(field_raw.get("boundary_depth_cells", 1.5)),
            flux_tolerance=float(field_raw.get("flux_tolerance", 1.0e-4)),
            kinematic_viscosity_um2_s=float(field_raw.get("kinematic_viscosity_um2_s", 3.0e6)),
            hybrid_finite_element_distance_um=float(
                field_raw["hybrid_finite_element_distance_um"]
            ),
            hybrid_transition_width_um=float(
                field_raw["hybrid_transition_width_um"]
            ),
            blood_density_kg_m3=float(field_raw.get("blood_density_kg_m3", 1060.0)),
            gmsh_bulk_mesh_size_um=float(field_raw.get("gmsh_bulk_mesh_size_um", 0.0)),
            gmsh_wall_mesh_size_um=float(field_raw.get("gmsh_wall_mesh_size_um", 0.0)),
            gmsh_wall_refinement_distance_um=float(
                field_raw.get("gmsh_wall_refinement_distance_um", 0.0)
            ),
            gmsh_element_order=int(field_raw.get("gmsh_element_order", 1)),
            dolfinx_velocity_degree=int(field_raw.get("dolfinx_velocity_degree", 2)),
            dolfinx_pressure_degree=int(field_raw.get("dolfinx_pressure_degree", 1)),
            dolfinx_ksp_rtol=float(field_raw.get("dolfinx_ksp_rtol", 1.0e-8)),
        ),
        particles=ParticleConfig(
            inlet_number_concentration_mb_per_ml=float(
                particles_raw["inlet_number_concentration_mb_per_ml"]
            ),
            n_steps=int(particles_raw["n_steps"]),
            dt_s=float(particles_raw["dt_s"]),
            bubble_diameter_min_um=float(particles_raw["bubble_diameter_min_um"]),
            bubble_diameter_max_um=float(particles_raw["bubble_diameter_max_um"]),
            max_unique_bubbles=int(particles_raw.get("max_unique_bubbles", 0)),
            max_particle_frame_records=int(particles_raw.get("max_particle_frame_records", 5_000_000)),
            acceleration_backend=str(particles_raw.get("acceleration_backend", "auto")).lower(),
            contact_geometry_tolerance_um=float(
                particles_raw.get("contact_geometry_tolerance_um", 0.001)
            ),
            contact_max_time_refinements=contact_max_time_refinements_raw,
            wall_contact_threshold_um=float(particles_raw.get("wall_contact_threshold_um", 0.05)),
        ),
        particle_dynamics=dynamics,
        cardiac_pulsatility=cardiac,
        red_blood_cell_root_discharge_hematocrit=(
            red_blood_cell_root_discharge_hematocrit
        ),
        molecular_target_selection=target_selection,
        molecular_target=target,
        molecular_binding=binding,
        binding_scenario_sweep=sweep,
        quick_test=quick_enabled,
    )
    _validate_domain_config(config.domain)
    _validate_field_config(config.field)
    _validate_particle_config(config.particles)
    return config


def _validate_domain_config(cfg):
    """Validate continuous-boundary controls owned by the domain geometry."""

    maximum_length = float(
        cfg.continuous_boundary_maximum_element_length_um
    )
    if not math.isfinite(maximum_length) or maximum_length <= 0.0:
        raise ValueError(
            "domain.continuous_boundary_maximum_element_length_um must be "
            "finite and positive."
        )


def _validate_field_config(cfg):
    """Validate the boundary-fitted DOLFINx solver controls."""

    viscosity = float(cfg.kinematic_viscosity_um2_s)
    if not math.isfinite(viscosity) or viscosity <= 0.0:
        raise ValueError("field.kinematic_viscosity_um2_s must be finite and positive.")
    density = float(cfg.blood_density_kg_m3)
    if not math.isfinite(density) or density <= 0.0:
        raise ValueError("field.blood_density_kg_m3 must be finite and positive.")
    for name in (
        "gmsh_bulk_mesh_size_um",
        "gmsh_wall_mesh_size_um",
        "gmsh_wall_refinement_distance_um",
    ):
        value = float(getattr(cfg, name))
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"field.{name} must be finite and non-negative.")
    finite_element_distance = float(cfg.hybrid_finite_element_distance_um)
    transition_width = float(cfg.hybrid_transition_width_um)
    if not math.isfinite(finite_element_distance) or finite_element_distance <= 0.0:
        raise ValueError(
            "field.hybrid_finite_element_distance_um must be finite and positive."
        )
    if not math.isfinite(transition_width) or transition_width <= 0.0:
        raise ValueError(
            "field.hybrid_transition_width_um must be finite and positive."
        )
    if int(cfg.gmsh_element_order) != 1:
        raise ValueError(
            "field.gmsh_element_order must be 1. The vessel boundary is "
            "piecewise linear; higher-order Gmsh geometry is unreliable at "
            "merged vascular junction points."
        )
    velocity_degree = int(cfg.dolfinx_velocity_degree)
    pressure_degree = int(cfg.dolfinx_pressure_degree)
    if pressure_degree < 1 or velocity_degree <= pressure_degree:
        raise ValueError(
            "field.dolfinx_velocity_degree must be greater than "
            "field.dolfinx_pressure_degree, which must be at least one."
        )
    if not math.isfinite(float(cfg.dolfinx_ksp_rtol)) or float(cfg.dolfinx_ksp_rtol) <= 0.0:
        raise ValueError("field.dolfinx_ksp_rtol must be finite and positive.")


def _validate_particle_config(cfg):
    """Validate particle configuration once at the YAML boundary."""

    concentration = float(cfg.inlet_number_concentration_mb_per_ml)
    if not math.isfinite(concentration) or concentration <= 0.0:
        raise ValueError(
            "particles.inlet_number_concentration_mb_per_ml must be finite and positive."
        )
    if cfg.n_steps < 1:
        raise ValueError("particles.n_steps must be at least one.")
    if not math.isfinite(cfg.dt_s) or cfg.dt_s <= 0.0:
        raise ValueError("particles.dt_s must be finite and positive.")
    if (
        not math.isfinite(cfg.bubble_diameter_min_um)
        or cfg.bubble_diameter_min_um <= 0.0
        or not math.isfinite(cfg.bubble_diameter_max_um)
        or cfg.bubble_diameter_max_um < cfg.bubble_diameter_min_um
    ):
        raise ValueError(
            "particles requires 0 < bubble_diameter_min_um <= bubble_diameter_max_um."
        )
    if cfg.max_unique_bubbles < 0 or cfg.max_particle_frame_records < 0:
        raise ValueError("Particle capacity limits must be non-negative.")
    if cfg.acceleration_backend not in {"auto", "numba_cpu", "python"}:
        raise ValueError(
            "particles.acceleration_backend must be one of: auto, numba_cpu, python."
        )
    if (
        not math.isfinite(cfg.contact_geometry_tolerance_um)
        or cfg.contact_geometry_tolerance_um <= 0.0
    ):
        raise ValueError(
            "particles.contact_geometry_tolerance_um must be finite and positive."
        )
    if cfg.contact_max_time_refinements < 0:
        raise ValueError(
            "particles.contact_max_time_refinements must be non-negative."
        )
    if not math.isfinite(cfg.wall_contact_threshold_um) or cfg.wall_contact_threshold_um < 0.0:
        raise ValueError("particles.wall_contact_threshold_um must be non-negative.")


def _validate_particle_dynamics(cfg):
    if cfg.time_integrator not in {"euler", "heun"}:
        raise ValueError("particle_dynamics.time_integrator must be 'euler' or 'heun'.")
    if cfg.integration_substeps < 1:
        raise ValueError("particle_dynamics.integration_substeps must be at least one.")
    if not (0.0 < cfg.xi_min <= cfg.xi_near < cfg.xi_far):
        raise ValueError("particle_dynamics requires 0 < xi_min <= xi_near < xi_far.")
    if not math.isfinite(cfg.collision_layer_um) or cfg.collision_layer_um < 0.0:
        raise ValueError("particle_dynamics.collision_layer_um must be non-negative.")
    if not math.isfinite(cfg.collision_relaxation_time_s) or cfg.collision_relaxation_time_s <= 0.0:
        raise ValueError("particle_dynamics.collision_relaxation_time_s must be positive.")
    if cfg.neighbor_search not in {"all_pairs", "cell_list", "auto"}:
        raise ValueError(
            "particle_dynamics.neighbor_search must be 'all_pairs', 'cell_list', or 'auto'."
        )
    if not math.isfinite(cfg.two_wall_warning_gap_ratio) or cfg.two_wall_warning_gap_ratio < 0.0:
        raise ValueError("particle_dynamics.two_wall_warning_gap_ratio must be non-negative.")


def _validate_cardiac_pulsatility(cfg):
    if cfg.waveform not in {"synthetic_ecg_envelope"}:
        raise ValueError(
            "cardiac_pulsatility.waveform must be 'synthetic_ecg_envelope'."
        )
    if not math.isfinite(cfg.bpm) or cfg.bpm <= 0.0:
        raise ValueError("cardiac_pulsatility.bpm must be finite and positive.")
    if (
        not math.isfinite(cfg.pulse_propagation_velocity_um_s)
        or cfg.pulse_propagation_velocity_um_s <= 0.0
    ):
        raise ValueError(
            "cardiac_pulsatility.pulse_propagation_velocity_um_s must be finite and positive."
        )
    if (
        not math.isfinite(cfg.initial_phase_fraction)
        or not 0.0 <= cfg.initial_phase_fraction < 1.0
    ):
        raise ValueError(
            "cardiac_pulsatility.initial_phase_fraction must be in [0, 1)."
        )
    if cfg.waveform_samples_per_cycle < 32:
        raise ValueError(
            "cardiac_pulsatility.waveform_samples_per_cycle must be at least 32."
        )
    if not math.isfinite(cfg.modulation_strength) or not 0.0 <= cfg.modulation_strength <= 1.0:
        raise ValueError(
            "cardiac_pulsatility.modulation_strength must be in [0, 1]."
        )


def _validate_molecular_target(cfg):
    if cfg.region_mode not in {"none", "mask_npz", "continuous_wall_npz"}:
        raise ValueError(
            "molecular_target.region_mode must be 'none', 'mask_npz', or "
            "'continuous_wall_npz'."
        )
    if not math.isfinite(cfg.target_density_molecules_per_m2) or (
        cfg.target_density_molecules_per_m2 < 0.0
    ):
        raise ValueError(
            "molecular_target.target_density_molecules_per_m2 must be finite and non-negative."
        )
    for name, value in (
        ("mask_array_key", cfg.mask_array_key),
        ("x_coordinates_key", cfg.x_coordinates_key),
        ("z_coordinates_key", cfg.z_coordinates_key),
    ):
        if not value.strip():
            raise ValueError(f"molecular_target.{name} must not be empty.")
    if not cfg.enabled:
        return
    if cfg.region_mode == "none":
        raise ValueError(
            "molecular_target.region_mode must identify a disease ROI when molecular_target.enabled is true."
        )
    if cfg.region_mode in {"mask_npz", "continuous_wall_npz"} and cfg.mask_npz_path is None:
        raise ValueError(
            "molecular_target.mask_npz_path is required for NPZ target modes."
        )


def _validate_molecular_target_selection(cfg):
    if cfg.default_mode not in {"manual", "automatic"}:
        raise ValueError(
            "molecular_target_selection.default_mode must be 'manual' or 'automatic'."
        )
    if cfg.random_seed < 0:
        raise ValueError(
            "molecular_target_selection.random_seed must be non-negative."
        )
    if cfg.random_field_modes < 64:
        raise ValueError(
            "molecular_target_selection.random_field_modes must be at least 64."
        )
    influence_fraction = cfg.influence_region_endothelial_wall_area_fraction
    if influence_fraction is not None and (
        not math.isfinite(influence_fraction) or not 0.0 < influence_fraction < 1.0
    ):
        raise ValueError(
            "molecular_target_selection.influence_region_endothelial_wall_area_fraction must be "
            "strictly between 0 and 1."
        )
    positive_fraction = cfg.target_positive_wall_fraction_within_influence
    if positive_fraction is not None and (
        not math.isfinite(positive_fraction) or not 0.0 < positive_fraction <= 1.0
    ):
        raise ValueError(
            "molecular_target_selection.target_positive_wall_fraction_within_influence "
            "must be greater than zero and no greater than one."
        )
    correlation_length = cfg.target_correlation_length_um
    if correlation_length is not None and (
        not math.isfinite(correlation_length) or correlation_length <= 0.0
    ):
        raise ValueError(
            "molecular_target_selection.target_correlation_length_um must be finite and positive."
        )
    if cfg.default_mode != "automatic":
        return
    missing_fields = [
        name
        for name, value in (
            (
                "influence_region_endothelial_wall_area_fraction",
                influence_fraction,
            ),
            (
                "target_positive_wall_fraction_within_influence",
                positive_fraction,
            ),
            ("target_correlation_length_um", correlation_length),
        )
        if value is None
    ]
    if missing_fields:
        raise ValueError(
            "Automatic molecular-target selection requires: "
            f"{', '.join(missing_fields)}."
        )


def _validate_molecular_binding(cfg, target):
    if cfg.model != "deterministic_mean_field":
        raise ValueError(
            "molecular_binding.model must be 'deterministic_mean_field'; discrete stochastic bonds are not supported."
        )
    non_negative_fields = (
        ("ligand_density_molecules_per_m2", cfg.ligand_density_molecules_per_m2),
        ("capture_distance_um", cfg.capture_distance_um),
        ("rest_length_um", cfg.rest_length_um),
        ("association_rate_m2_per_molecule_s", cfg.association_rate_m2_per_molecule_s),
        ("zero_force_dissociation_rate_s", cfg.zero_force_dissociation_rate_s),
        ("bond_stiffness_pn_per_um", cfg.bond_stiffness_pn_per_um),
        ("reactive_compliance_nm", cfg.reactive_compliance_nm),
    )
    for name, value in non_negative_fields:
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"molecular_binding.{name} must be finite and non-negative.")
    positive_fields = (
        ("temperature_k", cfg.temperature_k),
        ("mean_field_warning_count", cfg.mean_field_warning_count),
        ("bell_exponent_limit", cfg.bell_exponent_limit),
    )
    for name, value in positive_fields:
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"molecular_binding.{name} must be finite and positive.")
    if not cfg.enabled:
        return
    if not target.enabled:
        raise ValueError(
            "molecular_binding.enabled requires molecular_target.enabled so bonds cannot form on an undefined wall region."
        )
    if target.target_density_molecules_per_m2 <= 0.0:
        raise ValueError(
            "molecular_binding.enabled requires a positive molecular_target.target_density_molecules_per_m2."
        )
    strictly_positive_fields = (
        ("ligand_density_molecules_per_m2", cfg.ligand_density_molecules_per_m2),
        ("capture_distance_um", cfg.capture_distance_um),
        ("rest_length_um", cfg.rest_length_um),
        (
            "association_rate_m2_per_molecule_s",
            cfg.association_rate_m2_per_molecule_s,
        ),
        ("bond_stiffness_pn_per_um", cfg.bond_stiffness_pn_per_um),
        ("reactive_compliance_nm", cfg.reactive_compliance_nm),
    )
    for name, value in strictly_positive_fields:
        if value <= 0.0:
            raise ValueError(
                f"molecular_binding.{name} must be positive when molecular binding is enabled."
            )


def _validate_binding_scenario_sweep(cfg, target, binding, dynamics):
    reference_time_s = float(cfg.da_on_reference_time_s)
    if not math.isfinite(reference_time_s) or reference_time_s < 0.0:
        raise ValueError(
            "binding_scenario_sweep.da_on_reference_time_s must be finite and "
            "non-negative."
        )
    level_groups = (
        ("da_on_levels", cfg.da_on_levels),
        (
            "capture_distance_to_rest_length_ratios",
            cfg.capture_distance_to_rest_length_ratios,
        ),
        (
            "target_density_molecules_per_um2_levels",
            cfg.target_density_molecules_per_um2_levels,
        ),
        (
            "ligand_density_molecules_per_um2_levels",
            cfg.ligand_density_molecules_per_um2_levels,
        ),
    )
    for name, levels in level_groups:
        if any(level <= 0.0 for level in levels):
            raise ValueError(
                f"binding_scenario_sweep.{name} must contain only finite positive values."
            )
    if not cfg.enabled:
        return
    if reference_time_s <= 0.0:
        raise ValueError(
            "binding_scenario_sweep.enabled requires a positive "
            "da_on_reference_time_s fixed before the transport run."
        )
    if not target.enabled:
        raise ValueError(
            "binding_scenario_sweep.enabled requires molecular_target.enabled for a physical contact pilot."
        )
    if binding.enabled:
        raise ValueError(
            "binding_scenario_sweep requires molecular_binding.enabled=false because the "
            "contact pilot must measure unbound hydrodynamic contact opportunity."
        )
    if binding.rest_length_um <= 0.0:
        raise ValueError(
            "binding_scenario_sweep requires a positive molecular_binding.rest_length_um."
        )
    if not dynamics.store_full_diagnostics:
        raise ValueError(
            "binding_scenario_sweep requires particle_dynamics.store_full_diagnostics=true "
            "so the no-bond contact pilot can reconstruct rotational surface slip."
        )
    for name, levels in level_groups:
        if not levels:
            raise ValueError(f"binding_scenario_sweep.{name} must not be empty when enabled.")
