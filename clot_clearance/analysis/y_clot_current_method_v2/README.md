# Current V3 method on the original irregular Y-vessel clot

Task: `y_clot_current_method_v2`, run: `20261019T131047Z`.

This is one corrected authorized application, not a parameter search. A preceding boundary-adaptation attempt was stopped at 1.03375 ms and invalidated; its states are excluded from scientific conclusions. The 1927 material
points, 1561 original fixed anchors, 157637 original bonds, reference weights,
and unnotched initial geometry come verbatim from the accepted historical
`lps_variational_softening_v3_integration_and_single_pilot_v1` input. It is an
engineered wall-attached three-lobe shape, not patient-segmented anatomy.

Current variational softening, the original weak/mixed/strong thresholds, L1/L2
surface quadrature, beta=40 Pa s/m, and the uniform +X nominal 1000 Pa / 3.6 ms
protocol are retained. The V2 material field is translated by (-6,.05,1) mm and
clipped without scaling or refitting. This leaves 882 strong points (45.7706%)
and eight source-domain labels. The phase fraction differs from the rectangular
specimen; this is not a controlled geometry-only comparison.

The original Y lumen and sidewall geometry are active. Contact uses the old
20 µm effective radius, 50 N/m stiffness and 0.5 damping ratio. The elastic
normal is the closest-distance gradient, checked against the wall potential.
Wall forces and work are separate from the current-surface fluid forces.
No old Gaussian, CFD nodal load, COM drag or numerical damping is added.
There is no new CFD solve, bubble simulation or two-way FSI.

## Adaptations and preflight

- `supported_surface.py`: retain only simplices with all four existing lattice
  vertices and all six original material edges. The original graph rejects
  void-crossing chords. Do not manufacture missing edges or tetra vertices.
- `core_surface.py`: isolated copy of current L1/L2 code, with the support module
  and restart geometry hash pointing to this application's files.
- `current_surface.py`: require current quadrature centroids and 3 µm outward
  probes inside the original Y tetra lumen; evaluate energy-consistent wall
  repulsion and one-sided wall damping. This remains partial approximate
  surface loading, not a validated continuous wetted surface.
- `native.py`: include wall force in the velocity-dependent second half-kick.
- `run_cases.py`: separate wall work, retain the native material, contact,
  damage, identity and energy gates. The full original live graph determines
  components. No healing, particle replacement or displacement rescaling.

An initial preflight failed when a missing original edge was looked up while
building a support tetrahedron. It was corrected by rejecting unsupported
simplices before the one formal launch. Its initialization directory remains
under the preceding V1 run’s `checks/initialization_missing_reference_edge`. No physical parameters
were changed to resolve that implementation error.

`checks.py` validates zero initial material force, no initial wall contact,
native threshold installation, force/torque/power pairing, wall energy gradient,
three microsteps and one restarted continuation. These are not full protocols.

## Execution and outputs

The isolated server environment is `/workspace/y_clot/env.sh`. Native execution
also prepends the frozen V2 library directory to `LD_LIBRARY_PATH`, as recorded
in `run_formal.sh`. `prepare.py`, `checks.py`, then `launch.py` establish the
inputs, checks, source hash freeze and single supervisor-managed formal run.
Existing result directories must never be rerun or overwritten; reproduction
requires a new task/run directory and its own freeze.

`wait_then_post.py` waits for the formal summary, then runs `postprocess.py` and
`finish.py`. These generate real accepted-state red/Z-up video, Y-vessel
overview, initial/final point clouds, surface load field, scalar histories,
fragment size CSV, scientific report and validation records. No synthetic
physics frames or smooth replacement clot shells are created. Displayed short
bonds are a rendering subset; all original live bonds determine fragments.

`status.py` reads progress. `update_records.py` is WSL-only and updates the
scientific notebook and current entry point with a backup, hash guard and diff
after completed outputs have been copied and verified locally.

The scientific energy ledger is
`Wdrive + Wdrag + Wwall = ΔU + ΔK + D + ΔUc − Jc + R`.
Wall potential is diagnostic because the whole wall action is already in
`Wwall`; adding it again to the right-hand side would double count it.
The residual is never absorbed into damage dissipation.

Final conclusions must be read from `APPLICATION_RESULT.json` and
`REVIEW_Y_CLOT_CURRENT_METHOD_V1.md`; an animation alone does not establish
complete surface loading, material calibration or successful fragmentation.

The independent original-boundary ray check found 315 attachment-facing L1 facets incorrectly loaded in V1 (about 19% of its initial area). V2 classifies the original faces from the reconstructed historical exposed/attachment surfaces. An attachment face has zero fluid load while anchor-connected. It may activate only after genuine component release and an additional wall gap greater than 1 µm, smooth to 5 µm, plus current solid/lumen checks. Locally peeled but still globally attached back faces remain conservatively unrepresented. Geometry-only probes verify zero attached-back loading and positive loading after release/opening. No material or load parameter was changed for this correction.
