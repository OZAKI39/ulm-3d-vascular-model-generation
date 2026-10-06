# Data conventions

All coordinates and volumes refer to the undeformed configuration. The primary physical cube has 6,859 particles. Node order is the existing Cartesian `block_family` order; no node permutation is introduced.

| Array in `NODEWISE_RESIDUAL_DECOMPOSITION.npz` | Units / meaning |
| --- | --- |
| `X`, `x_target` | m; reference coordinates and manufactured current positions |
| `volume` | m³; physical trapezoid quadrature, including half/quarter/eighth boundary weights |
| `distance_faces` | m; columns x−, x+, y−, y+, z−, z+ |
| `d_min_over_delta`, `d_min_over_h` | Dimensionless; minimum of the six face distances |
| `q_bulk` | N/m³; closed-stencil BAQP internal force density at the same absolute Z |
| `b_local` | N/m³; independently evaluated −Div P |
| `e_model` | N/m³; q_bulk + b_local |
| `f_FD` | N; unchanged finite-domain BAQP internal nodal force |
| `f_bulk_target` | N; q_bulk times the actual physical volume |
| `f_BC_localized` | N; f_bulk_target − f_FD |
| `f_BC_naive` | N; −f_FD, diagnostic only |
| `f_body_local` | N; b_local times physical volume |
| `f_model` | N; e_model times physical volume |
| `f_external_NL` | N; −q_bulk V + f_BC_localized, frozen |
| `f_external_classical` | N; b_local V + f_BC_localized, frozen |
| `Fbar_target` | Dimensionless; actual discrete BAQP reconstructed gradient at x_target |
| `F_local` | Dimensionless; independent analytic diag(1,1,a(Z)) |

`NL_RECOVERED_STATE.npz` and `classical_RECOVERED_STATE.npz` retain current positions, nonzero perturbed initial positions, fixed-node mask, frozen external force, final internal force, reaction, reconstructed Fbar and per-node errors. Reactions are `−(internal+external)` on fixed nodes and zero elsewhere. The whole z− and z+ faces are prescribed; side nodes away from those faces remain free.

`F_recovery_error` is norm(Fbar(solved)−Fbar(target))/0.004. `F_analytic_error` is norm(Fbar(solved)−F_local)/0.004. These differ because a truncated nonaffine shape reconstruction is not analytically exact even at the prescribed nodal solution. `u_error` is norm(x−x_target) in m.

The solver residual is the free-node force L2 norm divided by max(free-node external-force L2 norm, F_ref). Local correction gates use the same fixed F_ref=3.143774999999835e−8 N on both N19 and N27; specimen area is not the normalization for local corrections.

The reused adapter's legacy `surface_traction_mean_Pa` / `surface_traction_max_Pa` fields are zero placeholders, not physical traction measurements. Actual analytic tractions are in the exact reference, and actual applied volume forces are in `external`. The legacy diagnostic label `phase=SCIPY_REFERENCE` names an equilibrium reporting stage; no fracture phase field is enabled. Mechanical/physical time remains null for these static solves.

`N27_DEPTH_CHECK.json` retains signed vector components and norms. Its points approach from a lateral x face at fixed Z=0, at depths 4h, 8h, 12h and 13h. Only the last point has a fully closed force dependency graph. The bulk comparator uses an independently shifted N29 scaffold; no N29 equilibrium is solved.

`FORCE_BALANCE.json` vectors are in N; `MOMENT_BALANCE.json` vectors are in N m. Moments use the reference centroid X=0. `VIRTUAL_WORK_DECOMPOSITION.csv` work values are J. Continuum integrations and nodal volume quadrature are separately labeled.

`BOUNDARY_LOCALIZATION.csv` profiles are weighted by nodal force magnitude. Quantiles measure conversion-load support in this finite specimen; they do not define a physical BC thickness. The bulk-model difference is retained in the classical residual and never moved into the boundary correction.
