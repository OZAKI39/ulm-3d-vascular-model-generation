# Isolated straight-pipe clot prototype

## Repository inspection (2026-09-29)

- Working/output root: `/home/lzy/projects/clot_clearance` (initially empty).
- Existing BraVa case: `/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin`.
- Actual local svMultiPhysics source: `/home/lzy/projects/computation_examples/_shared/formal_3D_flow_solver/FEM_SimVascular/external/flow_solver_source`.
- This source has CMake build infrastructure, `Code/Source/solver/mat_models.cpp`, `ustruct.cpp`, `fsi.cpp`, `set_bc.cpp`, and VTK readers (`vtk_xml*`). Searches found no peridynamics implementation in its solver code.
- BraVa's `vendor/flow_solver_support/wss.py` provides the existing P1 velocity-gradient/traction convention. Its original bytes will be copied before any reuse. Existing population mechanics describe rigid bubbles, not a deformable clot; they will not be repurposed as PD mechanics.
- Shared Python already supplies NumPy, SciPy, Numba, VTK/PyVista, Matplotlib and pytest. Use these libraries read-only; put caches/build/results inside this new root.
- No applicable AGENTS.md was found at this root or its ancestors. The shared particle workflow's AGENTS.md protects that separate project; no writes there are planned.
- User-supplied paper: `/home/lzy/projects/1-s2.0-S2666496826000567-main.pdf`; inspect and retain a local reference copy before selecting detailed parameters.

## Integration decision

Implement an independent Python/Numba NOSB-PD tool with a CMake build/CTest wrapper, standard VTK outputs, and JSON inputs. This avoids compiling or changing the existing SimVascular application. Use the existing VTK runtime to import external CFD fields. All new files, dependency additions, compiled caches, reports and visualizations remain under this directory. Original artifacts are checked against a read-only provenance manifest; nothing is written to BraVa or the solver source.

## Ordered phases

0. Inspect the supplied paper, snapshot relevant input/source identities, finish model notes.
1. Generate a rectangular clot attached to a small flat support inside a straight cylindrical pipe. Implement neighbor lists, active shape tensors, compressible Neo-Hookean correspondence mechanics, fixed base, energy-consistent provisional non-affine stabilization, and numerical guards. Build and test affine/rigid motion, force signs and zero load.
2. Implement analytic Poiseuille pipe background and localized synthetic surface traction. Resolve a low-frequency representative envelope with explicit damped dynamics; keep the MHz carrier as a distinct cycle-count clock. Verify small-load deformation.
3. Add irreversible instantaneous-stretch damage and verify thresholds/unloading.
4. Add history-based cycle-jump damage with every parameter in JSON. Verify more cycles and zero-rate controls.
5. Detect components from active bonds, distinguish connectivity from damage, and verify artificial fragmentation. Audit any loss of 3-D shape support explicitly; never use a singular inverse or invent a deformation gradient.
6. Implement VTK/CSV streaming import, normal/tangential stress transfer, coverage rejection, and analytic/manufactured import comparisons.
7. Run requested A–F tests plus mechanics/import/safety checks; execute a small straight-pipe demonstration; export particle/bond VTK, summary CSV/JSON, PNG/PDF figures, interactive HTML and animation; write exact reproduction commands and scientific limitations.

## Acceptance evidence

Store phase checks and final test logs, raw trajectories/bond integrity, source/config hashes and original-file preservation checks. Report actual damage, failed bonds and disconnected components. Synthetic loading, time-scale substitution, provisional stabilization, damage and erosion rules are research/software verification choices, not a validated thrombolysis prediction. Do not start or interfere with the existing remote BraVa study queues.

## Completed evidence (2026-09-30 local time)

- Phase 0: `references/PAPER_READING.md`, original PDF copy/text and 2,640 protected original hashes.
- Phase 1: `logs/build_phase1.log`, `logs/tests_phase1.log`; affine/objective mechanics, energy-gradient/sign checks, base/no-load and shape/J guards passed.
- Phases 2–6: `logs/tests_phases2to6.log`; traction, both damage modes, connectivity and CSV/VTK import checks passed.
- Phase 7: `logs/tests_final.log` / `logs/ctest_results.xml`: 18 automated tests passed. `verification/ACCEPTANCE.json`: A–F passed. `verification/IMPORT_AND_DT_CHECK.json`: importer trajectory parity and limited dt refinement passed.
- Executed main result: `results/straight_pipe_demo/`, 384 particles / 12,804 bonds / 8 macro steps, mild damage and no fragmentation. Artificial cleavage is an explicitly separate verification result.
- Smallest documented 96-particle configuration also executed at `runs/minimal_smoke/`.
- `OPEN_RESULTS.html`, 300 dpi PNG/PDF, offline Plotly view, 108-frame MP4 and GIF delivered. Corrected an animation label overlap after visual inspection; no simulation values were changed.
- `provenance/FINAL_QC.json`: original-file hashes, independent damage/connectivity, VTK equality and complete movie decoding passed. No browser engine was available for an automated interactive UI test.
- No original files changed and no remote queue/solver actions taken.
