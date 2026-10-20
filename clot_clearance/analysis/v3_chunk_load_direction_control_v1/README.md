# V3 load direction control

Task `V3_CHUNK_FRAGMENTATION_LOAD_DIRECTION_CONTROL_V1`, run `20261016T105642Z`.

This isolated experiment changes only the direction operand of the existing Gaussian surface force from the historical reference tangent to the original clot's reference outward exposed-surface normal. Original loaded atom IDs, reference areas, current-position Gaussian, pulse, materials, topology, drag, contact, damping and event-resolved native integration remain unchanged. This is a synthetic mechanical control, not a verified microbubble force.

The actual status and data are under `verification/v3_chunk_fragmentation_load_direction_control_v1/20261016T105642Z/`. Figures and their static index are under the matching `results/` directory. Always inspect the actual status before describing completion.

Server: activate `/workspace/y_clot/env.sh`; project `/workspace/y_clot/project`. The unchanged frozen library is `verification/lps_variational_softening_v3_integration_and_single_pilot_v1/20261014T232612Z/build/liblammps.so`. Prepend that build directory to `LD_LIBRARY_PATH` for native runs. The selected Git snapshot does not contain all required native binaries and raw accepted trajectories.

Call chain:

1. `prepare.py`: verifies both protected historical manifests; copies exact original input files; maps exposed triangles onto the 453 original loaded IDs; checks closed-solid and lumen sidedness and same-input force magnitudes.
2. `freeze.py`: records immutable directions, scalar parameters, source/library hashes and original scientific gates.
3. `run_protocol.py` through `start_managed.py run_protocol`: one native case, original initial state. Imports the frozen Native, Integrator, Mechanics, Journal, Checks and TypedLPS modules. `external.py` subclasses the historical External and only replaces its private reference direction array. No historical file or original reference array is edited.
4. `protocol_support.py`: preserves additional exact component observations; at 0.8 ms writes a permanent checkpoint, records the continuation decision, and actually reloads that checkpoint while retaining accumulated work and numerical defects.
5. `analyze.py`: verifies actual topology against every accepted history row; reuses the exact previous candidate classifier; preserves fragment identity evolution and force/energy semantics. Regular nodal sampling is not claimed to be every-step full state.
6. `figures.py`: actual points and original bonds, common scales, displacement ×1, PNG/SVG and `OPEN_RESULTS.html`.
7. Final review, actual WSL document update with backup/diff/readback, and payload sealing follow numerical completion.

Do not rerun preparation or native dynamics into this frozen directory. A future replication requires its own run ID and a fresh output directory; it is not a continuation of an already completed accepted trajectory. The prior C1 baseline is read only and is never replayed here. The root `ACTIVE_DEVELOPMENT_ENTRYPOINTS.md` is protected by the original manifest; this README supplies the new entrypoint without altering that historical artifact.

Scalar-rule equality is tested at identical coordinates and times. Different trajectories can produce different instantaneous scalar forces and work. No equal-energy efficiency claim, calibrated Gc, clinical fragment-size claim, or new material-law validation follows from this single comparison.
