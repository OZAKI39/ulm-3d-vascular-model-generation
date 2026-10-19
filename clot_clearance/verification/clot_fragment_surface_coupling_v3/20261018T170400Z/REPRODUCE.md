# Reproduce / replay

Task clot_fragment_surface_coupling_v3, run 20261018T170400Z. One formal H attempt only.

Activate `/workspace/y_clot/env.sh`; set OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1. Runtime root defaults to `/workspace/y_clot/project`; `Y_CLOT_PROJECT_ROOT` overrides it.

Read-only replay, with full archived data present:

```bash
python analysis/clot_fragment_surface_coupling_v3/analyze.py
python analysis/clot_fragment_surface_coupling_v3/compare_motion.py
python analysis/clot_fragment_surface_coupling_v3/coverage_audit.py
python analysis/clot_fragment_surface_coupling_v3/figures.py
python analysis/clot_fragment_surface_coupling_v3/render_videos.py --smoke
# Review the five actual-frame smoke previews before final video rendering.
python analysis/clot_fragment_surface_coupling_v3/render_videos.py
python analysis/clot_fragment_surface_coupling_v3/report.py
```

These commands read accepted dynamics and do not restart the physical experiment. The launch script deliberately refuses a second formal attempt in the frozen directory. A justified new experiment requires a new independent RUN and a new frozen request.

V2 dependency root: `verification/clot_biphasic_heterogeneous_shear_fragmentation_v2/20261017T231900Z`. Required original geometry / material / restart / states and native library hashes are in SOURCE_FREEZE.json and V2_DEPENDENCY_FREEZE.json. The raw library, dense trajectories, large native restart/data and paper PDFs are archival dependencies, not incidental Git payload. The native library SHA is `cd9abb11749025d398d17abb69cab68072dd1bf5e4c465e9709e98fa0144b970`. V2 BUILD_AUDIT.json records compiler/object provenance; upstream LAMMPS commit is `8de817dd79bfe4525d5d39246a212d833e6dee07`. This round changes no LAMMPS material/force kernel. The small ray helper is compiled from ray_cells.cpp with g++ -O3 -std=c++17 -fPIC -shared and its source SHA is recorded.

The binary restart is accompanied by surface_restart.npz and fluid_work_restart.npz; reference ID/time/step checks are mandatory. A failed-unaccepted native restart, if present, is explicitly not a checkpoint of last_accepted.npz. Continuation equivalence of the new surface sidecar is tested in NATIVE_SURFACE_HISTORY_RESTART_AUDIT.json.

Coordinates are metres in FP64 archives; figures/videos display micrometres without displacement magnification. Frame-time CSV files point to actual accepted native keys. No synthetic trajectory interpolation.

The physical-scale candidates are offline random fields with fixed physical Fourier coefficients, never formal dynamics. The active material remains the frozen V2 BIPHASIC_FIELD.npz.

After the full sealed archive has been synchronized locally, verify every archived file and relative link with the standard library only:

```bash
python3 analysis/clot_fragment_surface_coupling_v3/verify_local_delivery.py --root /path/to/clot_clearance
```

This checks DELIVERY_MANIFEST.json and writes LOCAL_DELIVERY_READBACK.json. It requires the full archive; a compact Git source export intentionally lacks large native states and is not a substitute for that archive.
