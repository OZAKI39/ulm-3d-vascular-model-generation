# Validated production CFD flow

The production method is
`VALIDATED_TAU1_ADAPTIVE_FLUX_CONTINUOUS_Q_STEADY_LBM`. It consumes the
accepted CFD surface and uses the validated dimensionless Seeder geometry/q
kernel, D3Q19 BGK Musubi, continuous wall q values, `wall_libb`, one
`adaptive_flux_pressure` inlet, and three `pressure_eq` outlets. The former
`MUSUBI_ONLY_RECOVERY`, `mfr_eq`, reference-pair timestep, fixed pressure
offset, and boundary-write-only conservation contracts are retired from the
production default.

## Numerical contract

The Base recommendation remains `dx = 0.20 µm`. With physical kinematic
viscosity `ν = 3.27e-6 m²/s`, production sets

```text
dt = dx² / (6 ν) = 2.038735983690112e-9 s
ν_lattice = 1/6
tau = 1
omega = 1
```

The pressure reference is recomputed for the selected spacing:

```text
P_ref = rho0 * cs² * (dx/dt)²
      = 3387510.7199999993 Pa at Base
```

`P_ref` is an LBM numerical pressure offset. It is not physiological absolute
blood pressure. Scientific reports and default color maps use gauge pressure
`p_gauge = p_solver - P_ref` and pressure differences. The three physical
outlet gauges remain `14.544978101274268`, `132.20454922317552`, and
`-13.700626673311461 Pa`; absolute solver boundary values are derived as
`P_ref + p_gauge`, never hard-coded as resolution-independent physics.

The inlet requirement is `GLOBAL_TARGET_VOLUMETRIC_FLOW`, implemented by
`adaptive_flux_pressure`. The physical targets are
`2.890180380479642e-12 kg/s` and `2.7369132390905703e-15 m³/s`. The controller
converts them dynamically with `targetMassFlow/rho0 * dt/dx³`. A historical
upstream `PARABOLIC` label is retained only in source provenance and no longer
selects a production Musubi boundary.

The pinned corrected Musubi controller forms its denominator from
`MPI_SUM(active local boundary counts)` after solid removal. It also enforces
point-index, finite-sample, active-bitmask, and MPI-reduction invariants with
fail-fast behavior. Both Seeder and Musubi paths are explicit in schema v2;
the executable SHA-256 is recomputed before launch and a mismatch is fatal.

## Steady and conservation acceptance

Future fresh steady runs use physical-time windows:

- short: `0.0002441406727828746 s` (about 119751 Base steps)
- long: `0.0004882813455657492 s` (about 239502 Base steps)

The candidate and confirmation gates include short/long mass residuals,
physical aperture-flow closure, velocity and gauge-pressure residuals, inlet
target error, outlet-fraction drift, density sanity, positive PDFs, maximum
lattice speed, controller target/control error, finite values, and no
significant averaged backflow. A built-in “steady state” message is not by
itself acceptance.

Final solver conservation uses
`MUSUBI_ONE_STEP_DISCRETE_MASS_IDENTITY_V2` with gate `1e-8`. Volumetric flow
is measured only as
`PHYSICAL_INTERIOR_CROSS_SECTION_VELOCITY_FLUX = integral_A(u·n)dA`, using
`STANDARDIZED_INTERIOR_PHYSICAL_PORT_PLANES_V3` and
`CONTINUOUS_APERTURE_GAUSS_MLS_QUADRATIC_V2`. Boundary PDF accounting is not
reported as physical Q.

## Configurations and commands

`configs/cfd_flow.yaml` is schema v2 and represents the future
`FRESH_STEADY` production contract. Because that mode is a multi-hour compute
operation, the entry point deliberately refuses to launch it silently in the
resource-bounded promotion session; a separately authorized long-compute
session must enable the fresh runner rather than presenting a replay as a
fresh solve.

The reproducible low-cost promotion regression is:

```powershell
D:\anaconda3\envs\pmp\python.exe cfd_flow.py configs\cfd_flow_promotion_regression.yaml
```

It performs exactly one fresh-equilibrium 5000-step Base smoke with the
production-generated Lua and corrected binary, then stops. It subsequently
replays the accepted Base restart through production decoding, physical port
flux, QC, VTU export, and visualization. It never calls Seeder and never uses
the smoke transient as the steady solution.

The accepted steady source is explicitly classified as
`VALIDATED_RESEARCH_BASE_ACCEPTED_RESTART`:

- iteration: `598755`
- physical time: `0.001220703363914373 s`
- restart SHA-256:
  `ffcd98b2dc684d1569d937d915b603805809c581d5341e71b17afac2ac64c39f`
- `fresh_full_production_steady_solve = false`

## Outputs and visual review

Each promotion regression writes a new
`outputs/cfd_flow/production_tau1_base_promotion_anchor003274_<timestamp>/`
tree containing input provenance, smoke evidence, steady replay evidence,
machine-readable QC, figures, and an offline HTML report. Primary outputs are:

- `qc/run_summary.json`
- `qc/production_primary_metrics.json`
- `qc/production_primary_metrics.csv`
- `flow/production_steady_flow_field.vtu`
- `visualization/visual_manifest.json`
- `visualization/01_velocity_overview.png` through
  `07_steady_qc_summary.png`
- `visualization/production_review.html`

The VTU contains `velocity_phy`, `velocity_magnitude_m_s`,
`velocity_magnitude_mm_s`, `pressure_gauge_pa`,
`pressure_absolute_solver_pa`, and `rho_lattice`. The absolute-pressure array
is retained for provenance; visualizations use gauge pressure.

Open `production_review.html` directly in a browser; it has no internet, CDN,
or remote-JavaScript dependency.

## Current scientific scope

Base and Coarse steady results are accepted. Their comparison is
`TWO-GRID RESOLUTION SENSITIVITY`, not formal grid convergence. Formal
three-grid GCI was not completed because Fine steady computation was
terminated under resource-budget constraints. No “grid independent proven”
claim is made.

WSS is `DEFERRED_TO_POST_GRID_PRODUCTION_VALIDATION` and is not visualized or
reported as validated in this package.
