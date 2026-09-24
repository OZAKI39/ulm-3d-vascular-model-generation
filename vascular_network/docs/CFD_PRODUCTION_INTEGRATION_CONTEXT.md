# CFD production integration context

Read this file first. The CFD research contract is validated and the current
branch has been reduced to a small integration surface. This phase did **not**
promote that contract into the production runner and made zero Seeder or
Musubi calls.

## Scientific status

- Base: `CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS`.
- Coarse: accepted; Coarse-to-Base two-grid sensitivity is available.
- Fine: mesh/q, controller root cause/fix, and fresh 5,000-step safety are
  accepted; Fine steady was not completed because of the resource budget.
- Formal three-grid Richardson/GCI: not completed; do not claim it passed.
- Conservation owner: `MUSUBI_ONE_STEP_DISCRETE_MASS_IDENTITY_V2`.
- Physical-Q owner: V3 physical planes plus continuous-aperture Gauss/MLS V2.

## Production integration goal

The next phase may promote the validated Tau=1 contract with one low-cost
solver smoke and a new visual acceptance package. Until then, preserve the
current behavior of `cfd_flow.py` and `configs/cfd_flow.yaml`.

The current production path still uses the historical baseline (`mfr_eq`,
the existing production scaling/schema/status). These are explicit
`NEXT_PHASE_REPLACE` targets, not defects to repair in this slimming commit.

## Canonical active modules

- `cfd_flow.py`: production CLI entry point.
- `utils/cfd_flow/config.py`: current production config schema/parser.
- `utils/cfd_flow/apes.py`: APES scaling, Lua rendering, and launch plumbing.
- `utils/cfd_flow/geometry.py`: production surface partition and mesh geometry.
- `utils/cfd_flow/pipeline.py`: current production orchestration; next-phase
  rewrite target, with behavior frozen in this phase.
- `utils/cfd_flow/qc.py`: current production QC.
- `utils/cfd_flow/visualization.py`: production visualization infrastructure.
- `utils/cfd_flow/io.py`: common paths, inputs, JSON, and hashes.
- `utils/cfd_flow/validated_contract.py`: sole validated Tau=1 inputs/formulas.
- `utils/cfd_flow/validated_evidence.py`: compact evidence/artifact locator.
- `utils/cfd_flow/restart_decode.py`: canonical restart/PDF/macro decode.
- `utils/cfd_flow/exact_link_flux.py`: generic D3Q19 boundary/PDF primitives.
- `utils/cfd_flow/musubi_boundary_mass_referee.py`: diagnostic boundary replay.
- `utils/cfd_flow/full_timestep_mass_referee.py`: final Full V2 conservation.
- `utils/cfd_flow/port_grid_sensitivity.py`: continuous physical port geometry.
- `utils/cfd_flow/physical_port_flux.py`: V3 planes and physical-flux V2.
- `utils/cfd_flow/steady_state.py`: generic physical-time steady audit.
- `utils/cfd_flow/steady_export.py`: decoded Cartesian VTK reconstruction.
- `utils/cfd_flow/adaptive_target_population.py`: 375/376 regression oracle.
- `utils/cfd_flow/port_flux_audit.py`: source-proven APES boundary readers.

Optional solver-free research utilities, excluded from the default production
context, are `dimensionless_geometry_kernel.py`, `tau1_grid_convergence.py`,
and `wall_qvalue_oracle.py`.

## Canonical formulas

All are implemented only in `validated_contract.py`:

- `dt = dx^2 / (6 nu)`; Base `dt = 2.038735983690112e-9 s`.
- `P_ref = rho0 * cs^2 * (dx/dt)^2`; Base `3387510.7199999993 Pa`.
- `targetFlux = targetMassFlow/rho0 * dt/dx^3`.
- gauge/absolute pressure conversion.
- physical-to-lattice relaxation (`nu_lattice`, `tau`, `omega`).

The pressure reference is a numerical LBM offset, not physiological absolute
pressure. `23622.320128 Pa` is retained only as a historical regression value.
Accepted Qin, pressure drops, and flow fractions remain in evidence JSON, not
Python constants.

## Canonical evidence

Use `utils/cfd_flow/validated_evidence.py` or
`docs/CFD_VALIDATED_EVIDENCE_INDEX.md`; do not search the full output tree.
The compact set covers geometry freeze, Base final, Coarse final/two-grid
sensitivity, plane V3, adaptive 375/376 root cause/fix, Fine 5,000 safety,
Fine resource termination, and Full V2.

## Binary provenance and accepted artifacts

- Base restart SHA256:
  `ffcd98b2dc684d1569d937d915b603805809c581d5341e71b17afac2ac64c39f`.
- Base mesh `elemlist.lsb`:
  `f7d7b1d55273c78c336ac04e39bc018dd9ebb470a9f29ce833ff01711de8c386`.
- Base mesh `bnd.lsb`:
  `520d7dd1e4a46a45f9b1218a5807cfd89d6f054e0a247872362b130ff6bcfe69`.
- Base mesh `qval.lsb`:
  `35884406b5f0111cd4ab471f7b08ac3df00e478d3458a57636d1bd8921cb0fe6`.
- Coarse restart SHA256:
  `9eda88e685e5eaa6af650757b8c97accba392c88988df2766d2bb4de825fcfbb`.
- Plane contract SHA256:
  `ffaa49bdb6e43fb7208ff29df07a90d4e92ef9bfa4b96ca4f997d4f453a7f005`.

Large binaries are preserved but excluded from default text context. Check
existence, size, SHA256, and headers only unless targeted replay is authorized.

## Retired implementations

One-shot Tau1/Base/Fine runners, long-run monitors, fallback-q and fixed-pressure
harnesses, CELL_CUBE V1, boundary-only final-referee claims, high-tau benchmark
orchestration, duplicated restart/pressure/scaling code, and their redundant
CLI/tests were removed. They remain recoverable from Git history. Consult
`retired_cfd_components.json` for the path-level inventory; do not recreate
them as active alternatives.

## Exact next-phase integration surface

Read in this order:

1. This file.
2. `utils/cfd_flow/validated_contract.py`.
3. `outputs/cfd_flow/code_slimming_for_production_20260902_002647/qc/production_integration_surface.json`.
4. `cfd_flow.py`.
5. `configs/cfd_flow.yaml`.
6. `utils/cfd_flow/config.py`.
7. `utils/cfd_flow/apes.py`.
8. `utils/cfd_flow/pipeline.py`.
9. `utils/cfd_flow/qc.py`.
10. `utils/cfd_flow/visualization.py`.
11. `utils/cfd_flow/validated_evidence.py`.
12. `utils/cfd_flow/restart_decode.py`.
13. `utils/cfd_flow/full_timestep_mass_referee.py`.
14. `utils/cfd_flow/physical_port_flux.py`.
15. `utils/cfd_flow/steady_state.py`.
16. `utils/cfd_flow/adaptive_target_population.py`.

Modify the production YAML/schema/scaling/Lua/orchestration/QC coherently,
retain the 0.001 boundary-window gate, run only the separately authorized
low-cost smoke, and create fresh visual regression evidence. Do not rerun the
historical long Base or resume Fine steady.

NEXT: `PROMOTE VALIDATED TAU1 CFD CONTRACT TO PRODUCTION PIPELINE WITH LOW-COST SOLVER SMOKE AND NEW VISUAL ACCEPTANCE PACKAGE`.
