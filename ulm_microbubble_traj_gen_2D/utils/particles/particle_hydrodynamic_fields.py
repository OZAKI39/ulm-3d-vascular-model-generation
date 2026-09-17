"""
Precompute Eulerian fields required by the finite-size particle model.

The flow solver produces a steady finite-element velocity together with a
Cartesian far-field cache, and the vessel rasterizer provides viscosity.
Particle mobility also needs the velocity gradient, while wall distance and
inward normal come directly from the continuous vessel geometry at each live
particle position.

This module prepares the hybrid velocity field, reusable grid-carried material
fields, and lightweight Boolean overlays. It deliberately does not use a
raster wall-distance or wall-normal field for particle mechanics.

Coordinate and unit conventions
--------------------------------
The physical simulation plane is X-Z.  
Velocity components are ordered as (u_x, u_z), 
and the last two axes of the gradient tensor are ordered as [velocity component, spatial coordinate].  
Length and velocity are stored in µm and µm/s, respectively, so differentiating velocity with respect to position directly yields inverse seconds.  
Viscosity arrives from the vascular model in mPa s and is converted here to Pa s because the mobility and force laws use SI dynamic viscosity.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy import ndimage
from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from ..flow.hybrid_velocity import validate_hybrid_velocity_field
from .particle_boundary_support import sample_continuous_boundary_support_to_grid
from ..core.types import FlowField, GridDomain, HybridVelocityField, RasterizedVessels


@dataclass(frozen=True)
class ParticleHydrodynamicFields:
    """
    Grid-carried flow inputs plus the authoritative continuous wall geometry.
    """
    velocity_gradient_s_inv: np.ndarray
    dynamic_viscosity_pa_s: np.ndarray
    solid_site_mask: np.ndarray
    open_boundary_mask: np.ndarray
    boundary_geometry: ContinuousVesselGeometry
    hybrid_velocity: HybridVelocityField


def build_particle_hydrodynamic_fields(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow_field: FlowField,
    *,
    continuous_geometry: ContinuousVesselGeometry,
) -> ParticleHydrodynamicFields:
    """
    Build grid flow fields and attach the required continuous wall geometry.
    """
    # Validate the grid before any derivative or distance-transform work.  
    # The size of grid must be at least 2 x 2 to support finite differences, 
    # and the spacing must be positive and finite to avoid division by zero or NaN propagation.
    shape       = tuple(int(value) for value in domain.shape)
    spacing_um  = float(domain.spacing_um)
    if len(shape) != 2 or min(shape) < 2:
        raise ValueError("Particle hydrodynamic fields require a two-dimensional grid of at least 2 x 2 cells.")
    if not np.isfinite(spacing_um) or spacing_um <= 0.0:
        raise ValueError("domain.spacing_um must be finite and positive.")

    # Gradients are computed in float64 to reduce the loss of precision during floating-point calculations..
    lumen           = np.asarray(raster.lumen_mask, dtype=bool)
    velocity        = np.asarray(flow_field.velocity_xz_um_s, dtype=np.float64)
    viscosity_mpas  = np.asarray(raster.viscosity_mpas, dtype=np.float64)

    if lumen.shape != shape:
        raise ValueError(f"raster.lumen_mask shape {lumen.shape} does not match domain shape {shape}.")
    if velocity.shape != (*shape, 2):
        raise ValueError(f"flow_field.velocity_xz_um_s shape {velocity.shape} does not match expected shape {(*shape, 2)}.")
    if viscosity_mpas.shape != shape:
        raise ValueError(f"raster.viscosity_mpas shape {viscosity_mpas.shape} does not match domain shape {shape}.")
    if not np.all(np.isfinite(velocity)):
        raise ValueError("flow_field.velocity_xz_um_s must contain only finite values.")
    validate_hybrid_velocity_field(flow_field.hybrid_velocity)

    # Calculate and store the complete 2 x 2 velocity-gradient tensor.
    ux = velocity[..., 0]
    uz = velocity[..., 1]
    dux_dx, dux_dz = np.gradient(ux, spacing_um, spacing_um, edge_order=1)
    duz_dx, duz_dz = np.gradient(uz, spacing_um, spacing_um, edge_order=1)

    gradient = np.empty((*shape, 2, 2), dtype=np.float32)
    gradient[..., 0, 0] = dux_dx
    gradient[..., 0, 1] = dux_dz
    gradient[..., 1, 0] = duz_dx
    gradient[..., 1, 1] = duz_dz

    # Extend viscosity into non-lumen storage cells before particle sampling.
    # A near-wall bilinear stencil may include one of those cells; 
    # leaving its raster default at zero would artificially lower drag and 
    # amplify particle velocity exactly where the wall correction is most sensitive.
    dynamic_viscosity_pa_s = _nearest_fluid_dynamic_viscosity_pa_s(viscosity_mpas, lumen)

    # These Boolean overlays support target selection and visualization only.
    # Particle clearance, normals, and contact are queried from the continuous
    # geometry at the live particle position and never interpolated from them.
    solid_sites, open_boundary = sample_continuous_boundary_support_to_grid(
        domain,
        continuous_geometry,
    )

    return ParticleHydrodynamicFields(
        velocity_gradient_s_inv=np.ascontiguousarray(gradient),
        dynamic_viscosity_pa_s=np.ascontiguousarray(dynamic_viscosity_pa_s),
        solid_site_mask=solid_sites,
        open_boundary_mask=open_boundary,
        boundary_geometry=continuous_geometry,
        hybrid_velocity=flow_field.hybrid_velocity,
    )


def _nearest_fluid_dynamic_viscosity_pa_s(viscosity_mpas: np.ndarray, lumen_mask: np.ndarray) -> np.ndarray:
    """
    Extend valid lumen viscosity to every cell, then convert mPa s to Pa s.
    """

    viscosity   = np.asarray(viscosity_mpas, dtype=np.float64)
    lumen       = np.asarray(lumen_mask, dtype=bool)
    valid       = lumen & np.isfinite(viscosity) & (viscosity > 0.0)
    if not np.any(valid):
        raise ValueError("Cannot build particle viscosity field without a positive finite lumen viscosity value.")

    nearest_indices = ndimage.distance_transform_edt(~valid, return_distances=False, return_indices=True)
    extended_mpas   = viscosity[tuple(nearest_indices)]

    # 1 mPa s = 1e-3 Pa s.
    return np.asarray(extended_mpas * 1.0e-3, dtype=np.float32)
