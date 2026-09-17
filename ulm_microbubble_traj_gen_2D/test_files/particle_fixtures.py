"""Shared production-path fixtures for particle transport tests."""

from __future__ import annotations

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.cardiac.cardiac_pulsatility import (
    build_cardiac_pulsatility,
)
from ulm_microbubble_traj_gen_2D.utils.core.config import (
    CardiacPulsatilityConfig,
    FieldConfig,
    MolecularBindingConfig,
    MolecularTargetConfig,
    ParticleConfig,
    ParticleDynamicsConfig,
)
from ulm_microbubble_traj_gen_2D.utils.core.types import (
    FlowField,
    GridDomain,
    RasterizedVessels,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_field import (
    build_molecular_target_field,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_hydrodynamic_fields import (
    build_particle_hydrodynamic_fields,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_perfusion_transport import (
    advect_particles_with_continuous_perfusion,
)
from ulm_microbubble_traj_gen_2D.test_files.hybrid_velocity_fixture import (
    rectangular_hybrid_velocity,
)


def particle_config(
    *,
    n_steps: int,
    dt_s: float = 0.01,
    inlet_number_concentration_mb_per_ml: float = 5.0e10,
    max_particle_frame_records: int = 5_000_000,
    bubble_diameter_min_um: float = 1.0,
    bubble_diameter_max_um: float = 1.0,
    acceleration_backend: str = "python",
) -> ParticleConfig:
    return ParticleConfig(
        inlet_number_concentration_mb_per_ml=inlet_number_concentration_mb_per_ml,
        n_steps=n_steps,
        dt_s=dt_s,
        bubble_diameter_min_um=bubble_diameter_min_um,
        bubble_diameter_max_um=bubble_diameter_max_um,
        max_particle_frame_records=max_particle_frame_records,
        acceleration_backend=acceleration_backend,
    )


def advect_test_particles(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    vessels: list[object],
    particle_cfg: ParticleConfig,
    dynamics_cfg: ParticleDynamicsConfig,
    *,
    cardiac_cfg: CardiacPulsatilityConfig | None = None,
    molecular_target_cfg: MolecularTargetConfig | None = None,
    molecular_binding_cfg: MolecularBindingConfig | None = None,
    red_blood_cell_network: object | None = None,
    random_seed: int = 42,
):
    """Prepare the same derived inputs that the production runner supplies."""

    field_cfg = FieldConfig(
        effective_thickness_um=1.0,
        boundary_depth_cells=1.5,
        flux_tolerance=1.0e-4,
        kinematic_viscosity_um2_s=1.0e6,
        hybrid_finite_element_distance_um=1.0,
        hybrid_transition_width_um=1.0,
    )
    cardiac_cfg = cardiac_cfg or CardiacPulsatilityConfig()
    molecular_target_cfg = molecular_target_cfg or MolecularTargetConfig()
    molecular_binding_cfg = molecular_binding_cfg or MolecularBindingConfig()
    continuous_geometry = build_continuous_vessel_geometry(vessels, domain)
    hydrodynamic_fields = build_particle_hydrodynamic_fields(
        domain,
        raster,
        flow,
        continuous_geometry=continuous_geometry,
    )
    molecular_target_field = (
        build_molecular_target_field(
            domain,
            hydrodynamic_fields,
            molecular_target_cfg,
        )
        if molecular_target_cfg.enabled
        else None
    )
    cardiac = build_cardiac_pulsatility(domain, raster, vessels, cardiac_cfg)
    return advect_particles_with_continuous_perfusion(
        domain,
        raster,
        flow,
        vessels,
        particle_cfg,
        dynamics_cfg,
        hydrodynamic_fields,
        effective_thickness_um=float(field_cfg.effective_thickness_um),
        boundary_depth_cells=float(field_cfg.boundary_depth_cells),
        random_seed=int(random_seed),
        cardiac=cardiac,
        molecular_binding_cfg=molecular_binding_cfg,
        molecular_target_field=molecular_target_field,
        red_blood_cell_network=red_blood_cell_network,
    )


def straight_channel_case() -> tuple[GridDomain, RasterizedVessels, FlowField]:
    nx, nz = 20, 7
    shape = (nx, nz)
    lumen = np.zeros(shape, dtype=bool)
    lumen[:, 1:6] = True
    sdf = np.zeros(shape, dtype=np.float32)
    sdf[:, 1:6] = np.asarray([1.0, 2.0, 3.0, 2.0, 1.0], dtype=np.float32)
    velocity = np.zeros((*shape, 2), dtype=np.float32)
    velocity[lumen, 0] = 10.0
    inlet = np.zeros(shape, dtype=np.int32)
    inlet[0, 2:5] = 1
    outlet = np.zeros(shape, dtype=np.int32)
    outlet[-1, 2:5] = 1

    domain = GridDomain(
        origin_um=np.asarray([0.0, 0.0, 0.0]),
        spacing_um=1.0,
        shape=shape,
        fixed_y_um=0.0,
        x_coordinates_um=np.arange(nx, dtype=float),
        z_coordinates_um=np.arange(nz, dtype=float),
    )
    raster = RasterizedVessels(
        lumen_mask=lumen,
        wall_mask=lumen.copy(),
        vessel_id=np.zeros(shape, dtype=np.int32),
        radius_um=np.full(shape, 3.0, dtype=np.float32),
        flow_rate_um3_s=np.ones(shape, dtype=np.float32),
        q2d_flow_um2_s=np.ones(shape, dtype=np.float32),
        viscosity_mpas=np.ones(shape, dtype=np.float32),
        direction_xz=np.zeros((*shape, 2), dtype=np.float32),
        distance_to_centerline_um=np.zeros(shape, dtype=np.float32),
        distance_to_wall_um=sdf.copy(),
        wall_normal_xz=np.zeros((*shape, 2), dtype=np.float32),
        lumen_fraction=lumen.astype(np.float32),
    )
    flow = FlowField(
        velocity_xz_um_s=velocity,
        speed_um_s=np.linalg.norm(velocity, axis=-1),
        wall_shear_stress_pa=np.zeros(shape, dtype=np.float32),
        local_shear_stress_pa=np.zeros(shape, dtype=np.float32),
        hybrid_velocity=rectangular_hybrid_velocity(
            domain, velocity_xz=(10.0, 0.0)
        ),
        inlet_label=inlet,
        outlet_label=outlet,
        solver_metadata={"physical_converged": True},
    )
    return domain, raster, flow
