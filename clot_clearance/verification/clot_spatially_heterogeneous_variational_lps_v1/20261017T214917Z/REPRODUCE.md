# Reproduce CLOT_SPATIALLY_HETEROGENEOUS_VARIATIONAL_LPS_V1

Run ID: 20261017T214917Z. Existing accepted outputs are immutable; create a new RUN before executing another physical case. This delivery contains exactly two formal full native shear cases, not a parameter search.

## Environment

```bash
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH=/workspace/y_clot/project/verification/clot_spatially_heterogeneous_variational_lps_v1/20261017T214917Z/build:${LD_LIBRARY_PATH:-}
cd /workspace/y_clot/project
```

Additional isolated-environment packages: GSTools 1.7.0, pandas 3.0.6. PyVista, VTK, Matplotlib, NumPy, scipy, mpi4py and FFmpeg are available in the existing environment. Exact package inventory is saved in ENVIRONMENT.json. No system GPU driver changes.

The native library is an isolated relink of locked LAMMPS objects. BUILD_AUDIT.json records every object and compile/link command. The modified bond class requires recompiling the original force/style factory and the dynamic coordinator as well as its own source; do not reuse factories compiled with the old class size. The old native library remains unchanged.

## Preparation and focused validation (new output directory only)

1. Change RUN in common.py to a new unique identifier. Keep the frozen original reference and material card. prepare.py copies these exact bytes, freezes the old accepted checkerboard amplitude trace, and records current Git state without modifying unrelated files.
2. Run spatial_heterogeneity_field.py once. It generates GSTools Gaussian z with rescale=1, 4096 RandMeth modes, seed 20261017 and nominal length 300 µm. q=tanh(z−c) has zero material-volume mean. sc is symmetric and normalized using original bulk reference weights. Preserve HETEROGENEITY_FIELD.npz; do not regenerate a favorable realization after outcomes.
3. Run build.py. It creates a new isolated liblammps.so and source diff.
4. In separate Python processes run uniform_probe.py old_uniform with HET_LIB and LD_LIBRARY_PATH pointing to the old library; then run uniform_probe.py new_uniform with the new library. Never load both native libraries into one process.
5. Run checks.py and path_audit.py. The latter is a prescribed opening path using the complete independent discrete energy, not a new native full-time case and not calibrated Gc.
6. Run freeze_sources.py only after all new-path tests pass. Do not modify frozen runtime code or inputs after formal calculations start.
7. Run manage.py uniform and manage.py heterogeneous. Each installs one supervisor program, autostart=false and autorestart=false. Supervisor logs and configs are retained. No RSS/wall-time cap is imposed; numerical energy and damage admission remain active.

## Native constitutive semantics

U(x,d) and its minimized endpoint-star force are unchanged. The sole new material variable is sc_b, frozen by original bond ID. It enters the native irreversible damage law. Local IDs are never tied to current LAMMPS array slots. Restart V4 stores sc; reading old V2/V3 defaults to typed sc. Interface and anchor properties remain unchanged.

Initial uncut kappa=0 and d=0 are identical. Original precut d=1 is preserved. The inherited initializer uses kappa=sc as the already-failed sentinel; thus H precut kappa values differ from U when local sc differs. This has no change in precut force or stored energy because those bonds remain d=1, but it is not bitwise equality of the entire initial kappa arrays. See INITIAL_HISTORY_AUDIT.json.

Predictor rejections happen before a native commit. There is no assertion that a postcommit failure is rolled back; a failure is archived as failed_unaccepted, and last_accepted is reported separately. No unsupported retry or altered material amplitude was used to seek fragmentation.

## Drive and evidence

SHEAR_DRIVE.npz is the original previous accepted amplitude history, linearly interpolated at time, with the unchanged half-cosine envelope applied separately. Both materials use the identical file and function. No work feedback, outcome-driven tuning or force rescaling.

Each accepted step records actual endpoint force-displacement work, signed drive and drag, U/K/D/Uc/Jc/R and actual surface diagnostics. R is not included in D. Original damage changes and exact surviving-graph component memberships are retained every accepted step. Full FP64 x/v/kappa/d snapshots are saved every approximately 25 µs and at named events, rather than at every accepted step.

## Analysis and visual replay

```bash
python analysis/clot_spatially_heterogeneous_variational_lps_v1/literature.py
python analysis/clot_spatially_heterogeneous_variational_lps_v1/analyze.py
python analysis/clot_spatially_heterogeneous_variational_lps_v1/initial_history_audit.py
python analysis/clot_spatially_heterogeneous_variational_lps_v1/figures.py
python analysis/clot_spatially_heterogeneous_variational_lps_v1/render_videos.py
python analysis/clot_spatially_heterogeneous_variational_lps_v1/report.py
```

For presentation-only replay on a local copy set Y_CLOT_PROJECT_ROOT to that project root. Native recomputation requires the frozen server library/dependencies. No network deployment is required: open results/.../OPEN_RESULTS.html locally.

Static images are PNG/PDF. The three videos are 1920×1080 H.264 at 24 FPS, with approximately 145 regular accepted frames plus real onset/failure/secondary event frames and 60 repeats of the actual final state. Events create a shared variable physical-time sampling schedule; both materials use the same playback schedule. Nearest actual saved states are aligned, their distinct actual timestamps are shown, and FRAME_TIMES.csv lists both paths. No trajectory interpolation, artificial displacement scaling or invented fragment shell. Full original surviving bonds determine topology even when only r0≤1.01dx bonds are drawn. Per-video FFmpeg command, cover, source state indices and hash are retained.

The original current-surface module is unchanged and remains a partial supported-geometry proxy. Hydrodynamics of fine components without a supported 3D surface remains unresolved. Reported correlations use one frozen realization and must not be interpreted as independent tissue samples or causality.

## Record preservation

update_records.py runs in WSL and updates /home/lzy/projects/clot_rupture.md plus the local entrypoint after a hash check, retaining a byte-for-byte before backup and exact diff. Sync the entrypoint to the server only if its previous hash matches. seal.py records delivery files; local_readback.py verifies the downloaded manifest and HTML/video links. Historical evidence and all original topology remain intact. There is no GitHub push in this task.
