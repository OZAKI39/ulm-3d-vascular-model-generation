# Gmsh / DOLFINx boundary-fitted flow backend

## Environment

DOLFINx is distributed through conda-forge rather than a normal Windows `pip`
wheel. Create the environment on Windows, Linux, or WSL with:

```bash
conda env create -f ulm_microbubble_traj_gen/environment-dolfinx.yml
conda activate ulm-dolfinx
```

The local packages `ulm_microbubble_traj_gen` and
`ulm_vascular_model_generator` must remain importable from the repository root.
On Linux, the backend uses PETSc through `petsc4py`. The conda-forge Windows
build does not include PETSc, so the same DOLFINx variational form and assembled
matrix are solved serially with SciPy. Both paths record their concrete linear
solver and algebraic residual in the saved metadata.

## Configuration

Select the backend in the `field` section:

```yaml
domain:
  continuous_boundary_maximum_element_length_um: 4.0

field:
  gmsh_bulk_mesh_size_um: 5.0
  gmsh_wall_mesh_size_um: 4.0
  gmsh_wall_refinement_distance_um: 1.0
  gmsh_element_order: 1
  dolfinx_velocity_degree: 2
  dolfinx_pressure_degree: 1
  dolfinx_ksp_rtol: 1.0e-8
  blood_density_kg_m3: 1060.0
```

A Gmsh size value of `0` selects an automatic value based on
`domain.grid_spacing_um`. The velocity degree must be greater than the pressure
degree; the default `P2-P1` pair is the Taylor--Hood element. The Gmsh
geometry order is fixed at `1`: the exported lumen boundary is piecewise
linear, and higher-order geometry is unreliable where coincident junction
points are merged.

`continuous_boundary_maximum_element_length_um` controls subdivision of the
already-constructed piecewise-linear lumen ring before it is passed to Gmsh.
It does not move the polygon or change its shape. Keeping this value much
smaller than `gmsh_wall_mesh_size_um` can nevertheless force a much finer mesh
than the requested Gmsh sizes, especially in long, narrow vessels.

## Numerical model

The backend solves steady incompressible Stokes flow on the continuous X-Z
lumen:

\[
-\nabla\cdot\left(2\nu\varepsilon(\mathbf u)\right)+\nabla p=0,
\qquad
\nabla\cdot\mathbf u=0.
\]

Solid vessel walls use no-slip velocity. For a `planar_2d` generator run,
`generate_microbubble_trajectories.py` derives one area-equivalent extrusion
depth from all root inlets,

\[
h_\mathrm{eq}=\frac{\sum_i \pi R_i^2}{\sum_i 2R_i}.
\]

Every root inlet uses a parabolic Dirichlet profile whose integral equals the
exported 3D vessel flow divided by this resolved depth. For a single root this
makes the planar inlet mean velocity exactly equal to the exported circular-
tube mean velocity while retaining a single network-wide depth for flux
conservation. Every terminal outlet uses the same kind of
parabolic velocity Dirichlet profile, directed outward and normalized to its
exported terminal flow. This matches the vascular generator's
`fixed_total_inlet_equal_terminal_shares` model. The solve stops before
assembly if the exported reference inlet and terminal flows are inconsistent.
The solved inlet and outlet fluxes are integrated directly on Gmsh facets;
their totals and every labeled opening must match the exported targets within
`field.flux_tolerance`.

The finite-element pressure uses kinematic units internally. It is converted
with `blood_density_kg_m3` and saved on the Cartesian output grid in mmHg,
as gauge pressure relative to one pinned pressure degree of freedom. The pin
only removes the constant-pressure nullspace; it is not an anatomical pressure
boundary condition.

Gmsh physical tags are stable:

- `1`: fluid surface;
- `10`: solid wall;
- `1000 + label`: inlet;
- `2000 + label`: outlet.

The finite-element velocity and pressure are evaluated at the existing
Cartesian lumen centres after the solve, so particle advection retains the
existing `FlowField` interface. WSS uses the finite-element velocity gradient
evaluated just inside the nearest physical wall and the corresponding normal
from `ContinuousVesselGeometry`; it does not use the raster
distance-transform normal.

`divergence_s_inv` and `face_flux_*` in the exported `FlowField` are
post-sampling Cartesian diagnostics. They are not the DOLFINx acceptance
criterion. Linear-system convergence on the boundary-fitted mixed
finite-element system is recorded under the `linear_solver_*` metadata fields.
