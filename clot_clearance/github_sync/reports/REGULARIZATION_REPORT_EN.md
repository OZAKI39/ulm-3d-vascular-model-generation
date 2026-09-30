<!-- Navigation copy only. Original report remains unchanged at ../../REGULARIZATION_REPORT_EN.md -->

# Energy-regularization development report

The existing clot_clearance project was extended with an isolated namespace. Eighteen study runs, verification modules, reports and a new movie are complete. **No main run reached bond failure or detachment; the model is NOT CONVERGED.** Parameters were not weakened to produce a fragmentation movie.

## Status and acceptance

SOFTWARE VERIFIED for the tested implementation: 48 tests pass, all 18 run audits pass, and every saved legacy forced-failure array is bitwise reproducible. ENERGY REGULARIZED applies to the specified numerical opening coupon and damage-energy scale only. It does not establish fracture-energy objectivity under general loading. Main mean damage changes by −56.37%/−67.44% and dissipation by −65.19%/−76.08% on medium/fine versus coarse spacing, hence NOT CONVERGED. New fragment resolution is NOT_ASSESSED_NO_DETACHMENT; 0/0 pulverization ratios are null, not zero. The independently reclassified legacy case is FRAGMENT_RESOLUTION_INADEQUATE (singleton 78.899%, low rank 82.569%). Different old/new loads prevent a matched comparison.

The requested first-detachment acceptance criterion, a time-step window containing first failure, fragmentation localization, finite-fragment dominance and transport convergence are **not satisfied**. Absence of events was retained rather than changing material parameters to pass.

## Preservation and provenance

No Git repository or applicable AGENTS.md was found. All 627 protected pre-existing source/config/result/movie files match their SHA-256 values; compatibility links and unrelated stopped BraVa jobs remain untouched. Every new run retains its configuration, coupon identity, source snapshot and hashes. [Source/output index](../../verification/regularization/SOURCE_OUTPUT_INDEX.md), [preservation](../../verification/regularization/FINAL_PRESERVATION.json), [legacy replay](../../verification/regularization/LEGACY_REPRODUCIBILITY.json), [tests](../../verification/regularization/FINAL_TESTS_SCOPED.log). The initial unrestricted pytest discovery also collected preserved test copies and failed module-name collection; the explicit `tests` suite then passed. VTK emits 40 NumPy shape deprecation warnings. No failed log or source snapshot was removed.

The frozen legacy entry is unchanged. The new `pd_clot.regularization.legacy` wrapper warns and labels D_break<1 as forced_failure_verification. New regularized runs reject any D_break other than one. Inherited Q_ref and other legacy-law fields are not used by the new damage mode.

## Configuration and consistent field

Pipe radius 1.25 mm, length 6 mm, Q=18 mL/min, μ=0.00345312 Pa·s, density 1056 kg/m³. Clot dimensions 1.25×0.75×0.75 mm, volume 0.703125 mm³; original G=1000 Pa, K=50000 Pa and NOSB stabilization retained. The [provided paper](https://doi.org/10.1016/j.apples.2026.100347) and [prior reading notes](../../references/PAPER_READING.md) remain references, not a claim of exact paper reproduction.

The new center is **(−0.825,0,−0.725) mm**, 0.2 mm upstream of the left face at x=−0.625 mm, vertically centered. It is an actual new field input. A hollow diamond marks the localization center, not a downward point load or a resolved bubble surface.

ManufacturedStreamingProvider takes curl(A), A_y=−UL exp(−r²/(2L²)); local u=U exp(−r²/(2L²))(−dz/L,0,dx/L). U=0.01 m/s, L=0.35 mm, modulation 1+0.9 sin(2π·500t), local pressure amplitude zero. The analytic divergence-free vortex is added to the existing pipe u,p,grad. Both t=[−pI+μ(grad u+grad uᵀ)]n and transport use this same field. No Navier–Stokes or coated-bubble problem is solved. Legacy 1 MHz/2 µm fields do not imply resolved ultrasound/bubble physics.

## A. Fracture-energy calibration

The user authorized a demonstration value: **Gc_demo=0.01 J/m²**, not an experimental value. Each spacing uses a separate 0.5 mm cubic, displacement-controlled two-rigid-half opening coupon. Only the symmetric weak-plane crossing bonds accumulate cyclic damage; no bond is manually cut. Projected crack area is 2.5e−7 m², counted once. A scalar root solve calibrates ψcrit against actual NOSB potential-energy loss divided by that area; independent six-point Gauss force integration measures external work. Low-C_E sweeps may extend the same opening increment to three times the initial step budget until the weak plane fails, without changing the baseline completed path.

| Case | h (mm) | Particles | Initial bonds | ψcrit (Pa) | Gc measured (J/m²) | Gc error (%) | First failure / detach N |
|---|---|---|---|---|---|---|---|
| coarse | 0.125 | 360 | 12592 | 0.012520992 | 0.009999937912 | 0.000620878 | Not observed / not observed |
| medium | 0.083333333 | 1215 | 51911 | 0.0345811901 | 0.009999840214 | 0.001597862 | Not observed / not observed |
| fine | 0.0625 | 2880 | 135180 | 0.0479243093 | 0.009999737837 | 0.002621631 | Not observed / not observed |


Independent work residuals are 2.15e−15–2.85e−14, with about 2.5e−9 J dissipated. The configured Gc tolerance is 5%. W_cycle,ij=0.5 E_young a_ij² V_ij, where a is half the stretch range and original weighted-family partition volumes sum to material volume. Wcrit,ij=ψcrit V_ij. ΔD=ΔN C_E max(W_cycle/Wcrit−threshold,0)^m_E, baseline C_E=.001, m_E=1, threshold=0 (no nonzero fatigue threshold). Damage is irreversible and complete at D=1. Calibration reuse rejects inconsistent spacing, horizon, material, fatigue parameters or Gc.

This is an interpretable directional strain-energy proxy, not an exact NOSB per-bond free-energy decomposition. Coupon success is protocol-specific. At fixed coordinates the actual potential decrease is accumulated; per-bond attribution proportional to ΔD×driver remains an estimate. Failure and fracture-energy JSONL ledgers store event cycle, previous damage, driver, allocated dissipation, spacing, horizon and ranks. They are empty in these main runs because no failures occurred.

## B. Spacing sensitivity

All spacings keep physical geometry, material, source/field, Gc, total N=25,000 and horizon/h=3.01. All particles remain attached and rank three. Null event N means not observed by the exposure endpoint.

| Case | Attached (mm³) | Detached (mm³) | Resolved (mm³) | Debris (mm³) | Singleton (mm³) | Rank 0/1 detached (mm³) | Damage dissipation (J) |
|---|---|---|---|---|---|---|---|
| coarse | 0.703125 | 0 | 0 | 0 | 0 | 0 | 2.112127e-17 |
| medium | 0.703125 | 0 | 0 | 0 | 0 | 0 | 7.35266331e-18 |
| fine | 0.703125 | 0 | 0 | 0 | 0 | 0 | 5.05228765e-18 |


The fixed slab uses z<origin_z+0.125 mm, excluding centers exactly at its upper boundary with a 1e−10 h roundoff tolerance. Fixed volume is 1/6,1/9,1/6 on coarse/medium/fine. Thus the boundary discretization is a disclosed confound even with identical continuous geometry. The 10% comparison tolerance is descriptive, not a mathematical convergence criterion. Zero detached volume does not validate fragmentation convergence or size distributions.

## C. Cycle-jump and time-step sensitivity

| Case | ΔN | dt (µs) | Mean D | Δ mean D (%) | Δ dissipation (%) | Max displacement (µm) | Final signed energy residual (%) |
|---|---|---|---|---|---|---|---|
| coarse | 1000 | 1 | 0.000104395759 | 0 | 0 | 0.26126723 | -0.01266413 |
| cycle_500 | 500 | 1 | 0.000104110926 | -0.272839 | -0.143144 | 0.26127049 | -0.01174284 |
| cycle_250 | 250 | 1 | 0.000103968606 | -0.409167 | -0.214566 | 0.26127212 | -0.01025135 |
| dt_half | 1000 | 0.5 | 0.000104397221 | 0.00140048 | 0.00125719 | 0.26126722 | -0.003166068 |
| medium | 1000 | 1 | 4.55492604e-05 | -56.3687 | -65.1883 | 0.28245005 | -0.01935095 |
| fine | 1000 | 1 | 3.39874662e-05 | -67.4436 | -76.0796 | 0.26347362 | -0.02623106 |


No material recalibration occurs between ΔN=1000/500/250. Damage is monotonic and histories are compared at common represented N. Failure/detachment shifts are unobservable, not zero. dt/2 spans the full exposure; it supports sensitivity of this unfractured response only. Mechanical proxy durations are 0.052/0.102/0.202 s, including warmup. Changing ΔN also changes simulated mechanical time; N/f=50 s is not the integrated transport time. External work is never multiplied by ΔN.

## D. Localization

The reference-coordinate sphere has R=0.525 mm; bonds use reference midpoints. Distant bonds are neither deleted nor frozen.

| Case | Damage inside R (%) | Material initially inside R (%) | Damage-weighted distance (mm) | Fixed volume (%) |
|---|---|---|---|---|
| coarse | 50.0597 | 17.7778 | 0.5917913 | 16.6667 |
| medium | 50.7 | 16.7901 | 0.5751262 | 11.1111 |
| fine | 52.7891 | 17.5 | 0.570932 | 16.6667 |


At coarse spacing, 50.06% of damage lies in 17.78% of the material volume. This is preferential localization, not confinement of all damage. Failure and detached-origin localization remain undefined.

## E. Energy

| Final coarse quantity | J |
|---|---|
| elastic_energy_J | 9.32368459013e-14 |
| stabilization_energy_J | 6.46298672971e-15 |
| kinetic_energy_J | 3.25552193412e-17 |
| external_work_J | 2.16972305732e-13 |
| damping_dissipation_J | 1.17191318963e-13 |
| damage_dissipation_estimate_J | 2.1121269959e-17 |
| transport_work_J | 0 |
| numerical_energy_residual_J | -2.74776553105e-17 |


R=U_elastic+U_stabilization+K−E_initial+D_damping+D_damage−W_surface−W_transport. Surface work follows actual trapezoidal force/displacement integration; damping/transport use actual kinetic-energy changes. Final coarse residual is −0.012664%, dt/2 −0.003166%. Small residuals verify bookkeeping, not local physical accuracy. Energy histories are supplied in CSV and PNG/PDF.

## F–G. Fragments and transport

Resolved-fragment diagnostics require at least 8 particles, volume 1.5625e−11 m³, 12 active internal bonds and 80% volume with rank≥2. These are numerical thresholds, not physical clot-size thresholds. Singletons are a separate category; other unsupported components are debris. No particles are deleted/merged. P_singleton>.20 or P_lowrank>.30 triggers a visible inadequacy warning. Intrinsic reduced-rank continuation is diagnosed, not silently accepted as 3-D NOSB material. All new detached/clearance volumes are zero; no new fragment-size distribution exists. Legacy count/volume/equivalent-diameter distributions are shown with different-loading labels.

Optional fragment_drag computes component volume/mass/COM, equivalent sphere diameter/area, same-field velocity sampling and Re. It uses Stokes drag with Schiller–Naumann correction below Re=1000, and Cd=.44 above; see [OpenFOAM Foundation source](https://cpp.openfoam.org/v12/SchillerNaumann_8C_source.html). Frozen-coefficient COM relaxation distributes total impulse conservatively by mass; no torque is implemented. Under-resolved debris retains labeled per-particle relaxation. The main case retains relaxation and has no detached transport to exercise. An independent preclassified component ODE comparison gives 1.13e−8/5.65e−9 m/s errors at dt=10/5 µs, without manufacturing a fragment in the PD simulation. Clearance is classified once at first crossing; plane crossing is not clinical thrombolysis.

## Parameter sensitivity

| Run | Gc demo (J/m²) | m_E | C_E | U (m/s) | Final mean D | Max D | First failure / detach N |
|---|---|---|---|---|---|---|---|
| surface_Gc_0.005_U_0.005 | 0.005 | 1.0 | 0.001 | 0.005 | 5.36094e-05 | 0.0002399158 | Not observed |
| surface_Gc_0.005_U_0.02 | 0.005 | 1.0 | 0.001 | 0.02 | 0.0008426697 | 0.003843757 | Not observed |
| surface_Gc_0.02_U_0.005 | 0.02 | 1.0 | 0.001 | 0.005 | 1.301727e-05 | 5.824638e-05 | Not observed |
| surface_Gc_0.02_U_0.02 | 0.02 | 1.0 | 0.001 | 0.02 | 0.0002044748 | 0.0009310573 | Not observed |
| sweep_C_E_0.0005 | 0.01 | 1.0 | 0.0005 | 0.01 | 8.28509e-05 | 0.0003757651 | Not observed |
| sweep_C_E_0.002 | 0.01 | 1.0 | 0.002 | 0.01 | 0.0001315396 | 0.0005966815 | Not observed |
| sweep_Gc_demo_J_m2_0.005 | 0.005 | 1.0 | 0.001 | 0.01 | 0.0002113513 | 0.0009589593 | Not observed |
| sweep_Gc_demo_J_m2_0.02 | 0.02 | 1.0 | 0.001 | 0.01 | 5.131418e-05 | 0.000232709 | Not observed |
| sweep_m_E_0.8 | 0.01 | 0.8 | 0.001 | 0.01 | 0.001426735 | 0.005333286 | Not observed |
| sweep_m_E_1.2 | 0.01 | 1.2 | 0.001 | 0.01 | 7.726234e-06 | 4.258811e-05 | Not observed |
| sweep_streaming_velocity_scale_m_s_0.005 | 0.01 | 1.0 | 0.001 | 0.005 | 2.648195e-05 | 0.000118501 | Not observed |
| sweep_streaming_velocity_scale_m_s_0.02 | 0.01 | 1.0 | 0.001 | 0.02 | 0.0004160712 | 0.001895644 | Not observed |


Gc/m_E/C_E changes each trigger a new coupon calibration; U changes reuse the same scale. m_E/C_E curves therefore compare jointly recalibrated models at fixed target coupon Gc, not partial derivatives at fixed ψcrit. Four extra corner runs complete the 3×3 Gc–U sample grid. All outcomes are retained, no response fit or visual selection is used. No run fails/detaches; largest particle D≈.00533 at m_E=.8. [All configurations, hashes and metrics](../../verification/regularization/comparisons/study_summary.csv).

## H. Importer and limitations

ResolvedStreamingFieldProvider supports original linear-tetrahedral VTU/legacy VTK P1 interpolation and explicitly declared convex-hull VTP/CSV interpolation. Velocity and pressure are nodal. Missing volume gradients are consistently derived; planar input must supply a full 3-D gradient. Explicit units are mandatory. NaN, out-of-domain queries, missing coverage, excessive divergence and gradient inconsistency are rejected. Observed ranges and optional range thresholds are reported. Outward-solid normal metadata is a declaration, not geometric proof. Static snapshots only; no temporal interpolation or nonconvex scattered-domain reconstruction.

Independent four-format affine tests pass, maximum SI errors u=8.67e−19, p=1.14e−13, grad=7.22e−16. A separate weak uniform synthetic VTU also passed one actual PD macro step through the config-selected importer; [verification](../../verification/regularization/resolved_import_with_runner/IMPORT_VERIFICATION.json), [run audit](../../verification/regularization/resolved_import_with_runner/runner_case/AUDIT.json). These fixtures are not real microbubble CFD.

Smallest next step: obtain one actual resolved-microbubble tetrahedral VTU covering the loaded surface and anticipated motion, with nodal u,p, explicit pressure reference/coordinate registration and preferably grad; provide the metadata below. In a copied config select `streaming.provider="resolved_field"`, `field_path`, `metadata_path`, then inspect import diagnostics and one macro step. The PD solver need not change. The meaning of cyclic forcing from a static mean field still needs definition.

```json
{
  "length_unit": "m", "velocity_unit": "m/s", "pressure_unit": "Pa", "time_unit": "s",
  "normal_convention": "outward_solid",
  "arrays": {"velocity": "velocity", "pressure": "pressure", "velocity_gradient": "velocity_gradient"},
  "validation": {"maximum_relative_divergence": 0.05, "maximum_relative_gradient_inconsistency": 0.25}
}
```


Clot material and fatigue law are not experimentally/sonothrombolysis calibrated. Manufactured streaming is not resolved microbubble flow, and no two-way fluid–PD feedback is solved. Surface quadrature/normals, fixed-slab discretization, stabilization and low-rank continuation remain limitations.

## Visual deliverables and reproduction

[Viewer](../../visualization/streaming_regularized_demo/index.html), [MP4](../../visualization/streaming_regularized_demo/regularized_damage_transport.mp4), [GIF](../../visualization/streaming_regularized_demo/regularized_damage_transport.gif), [QC](../../visualization/streaming_regularized_demo/QC.json). Black/Arial regular, single damage view, no arrows/large blue bubbles/Force panel/N header. Fixed camera, displacement scale one. 1920×1080, 60 fps, 751 frames, 12.5167 s. Diamond position is fixed and does not encode phase. Damage colorbar is explicitly 0–.0005; yellow is not complete damage one. Particle damage derives from neighboring bond integrity; D_break applies to bonds.

Every 30 display frames exactly matches one of 26 original states. Positions alone are linearly interpolated, damage holds its recorded macro value; all particles and exact anchors are preserved. Tiny real displacements (~.261 µm maximum) can produce identical adjacent decoded frames. Constant video cadence does not turn display interpolation into new solver samples. Four-stage panels use N=0/8000/17000/25000 and explicitly state that failure/transport stages were not reached; streamlines are an illustrative representative-peak-phase overlay.

- [cycle_jump_convergence PNG](../../visualization/streaming_regularized_demo/figures/cycle_jump_convergence.png) · [PDF](../../visualization/streaming_regularized_demo/figures/cycle_jump_convergence.pdf)
- [energy_history PNG](../../visualization/streaming_regularized_demo/figures/energy_history.png) · [PDF](../../visualization/streaming_regularized_demo/figures/energy_history.pdf)
- [four_stage_summary PNG](../../visualization/streaming_regularized_demo/figures/four_stage_summary.png) · [PDF](../../visualization/streaming_regularized_demo/figures/four_stage_summary.pdf)
- [fracture_energy_calibration PNG](../../visualization/streaming_regularized_demo/figures/fracture_energy_calibration.png) · [PDF](../../visualization/streaming_regularized_demo/figures/fracture_energy_calibration.pdf)
- [fragment_size_distribution PNG](../../visualization/streaming_regularized_demo/figures/fragment_size_distribution.png) · [PDF](../../visualization/streaming_regularized_demo/figures/fragment_size_distribution.pdf)
- [fragment_size_distribution_detail PNG](../../visualization/streaming_regularized_demo/figures/fragment_size_distribution_detail.png) · [PDF](../../visualization/streaming_regularized_demo/figures/fragment_size_distribution_detail.pdf)
- [localization_metrics PNG](../../visualization/streaming_regularized_demo/figures/localization_metrics.png) · [PDF](../../visualization/streaming_regularized_demo/figures/localization_metrics.pdf)
- [manufactured_field PNG](../../visualization/streaming_regularized_demo/figures/manufactured_field.png) · [PDF](../../visualization/streaming_regularized_demo/figures/manufactured_field.pdf)
- [mesh_convergence PNG](../../visualization/streaming_regularized_demo/figures/mesh_convergence.png) · [PDF](../../visualization/streaming_regularized_demo/figures/mesh_convergence.pdf)
- [parameter_sensitivity PNG](../../visualization/streaming_regularized_demo/figures/parameter_sensitivity.png) · [PDF](../../visualization/streaming_regularized_demo/figures/parameter_sensitivity.pdf)
- [response_surface PNG](../../visualization/streaming_regularized_demo/figures/response_surface.png) · [PDF](../../visualization/streaming_regularized_demo/figures/response_surface.pdf)
- [time_step_sensitivity PNG](../../visualization/streaming_regularized_demo/figures/time_step_sensitivity.png) · [PDF](../../visualization/streaming_regularized_demo/figures/time_step_sensitivity.pdf)
- [topology PNG](../../visualization/streaming_regularized_demo/figures/topology.png) · [PDF](../../visualization/streaming_regularized_demo/figures/topology.pdf)

Use fresh output directories. Current-suite testing explicitly excludes preserved copies:

```bash
cd /home/lzy/projects/clot_clearance
export PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$PWD/build_regularization/numba_cache"
export MPLCONFIGDIR="$PWD/build_fragmentation/matplotlib"
export PYTHONPATH="$PWD"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
CL_PY=/home/lzy/projects/computation_examples/_shared/python_environment/bin/python

# Coupon: use a new output directory; units are m and J/m².
"$CL_PY" -B -m pd_clot.regularization.coupon --config configs/streaming_regularized_demo_coarse.json --spacing 0.000125 --Gc-demo 0.01 --output verification/regularization/replay_coupon
# Demo using that new calibration.
"$CL_PY" -B -m pd_clot.regularization.runner --config configs/streaming_regularized_demo_coarse.json --calibration verification/regularization/replay_coupon/CALIBRATION.json --output results/streaming_regularized_demo_replay/coarse
# Medium/fine coupons are independently recalibrated by the mesh command.
"$CL_PY" -B -m pd_clot.regularization.campaign --study mesh --output-root results/streaming_regularized_demo_replay
# Cycle tests retain the baseline material and coarse calibration.
"$CL_PY" -B -m pd_clot.regularization.campaign --study cycle --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study timestep --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study sweep --output-root results/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.campaign --study sweep_surface --output-root results/streaming_regularized_demo_replay
# Current tests only; preserved snapshot directories are not a test suite.
"$CL_PY" -B -m pytest tests -q -p no:cacheprovider
"$CL_PY" -B scripts/verify_regularized_interfaces.py --mode drag --output verification/regularization/drag_replay
"$CL_PY" -B scripts/verify_regularized_interfaces.py --mode import --output verification/regularization/import_replay
# Presentation replay requires a fresh output directory.
"$CL_PY" -B scripts/render_regularized_damage.py --run results/streaming_regularized_demo/coarse --output visualization/streaming_regularized_demo_replay
"$CL_PY" -B scripts/check_regularized_damage_v2.py --output visualization/streaming_regularized_demo_replay
"$CL_PY" -B -m pd_clot.regularization.figures --output visualization/streaming_regularized_demo_replay/figures
```
