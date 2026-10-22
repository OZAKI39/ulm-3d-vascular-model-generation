# CPU/GPU resources

| Method | Observations | Complete | COM-feasible | Raw hits | First hit | Charged wall (s) | Ask/tell (s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Sobol | 32 | 32 | 12 | 1 | 26 | 98.200 | 0.000 |
| Ax_BoTorch | 32 | 32 | 23 | 9 | 10 | 106.082 | 10.303 |

Scientific workflow wall: 750.990 s. It includes baseline, search, extension, recovery and validation orchestration, but not initial environment downloads, radial audit, later artifact rendering or the additional low-pressure admission. Final cost ledger records those separately. Per-run solver/integration/postprocessing wall, process CPU and RSS are in ALL_RUNS.json. Resource samples are in RESOURCES.jsonl.

Four independent processes; OMP/OpenBLAS/MKL/NumExpr each one thread. RTX 4090 runs Warp acoustics; Ax uses CPU Torch. An unrelated CPU task was present; this is an observed campaign throughput, not an exclusive hardware benchmark. Effective scientific target yield is based on independent confirmation, not raw high scores.

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false


Final measured ledger: core scientific wall 1094.491 s (18.24 min); latest figures and videos 65.919 s. Median / p95 / max complete solver wall: {'median': 6.520963991002645, 'p95': 23.734557223896253, 'max': 39.061073250020854} s. Incomplete fraction 0.0130. Independently confirmed target count 0; confirmed targets/core-hour 0.000. Per-process peak RSS 1023.6 MiB. Resource samples: {"samples": 249, "interval_s": 3, "gpu_utilization_mean_percent": 22.32128514056225, "gpu_utilization_max_percent": 65.0, "gpu_memory_max_MiB": 2186.0, "cpu_count": 32, "worker_processes": 4, "BLAS_threads_each": 1, "scope": "Whole GPU snapshot, not exclusive per-process attribution; CPU RSS recorded per trajectory"}. Exclusions: Core = original radial jobs incl failed refinement + full scientific workflow + extra sub-0.1 Pa admission + full independent audit + incremental audit + single frozen low-domain followup. Excludes installs, developer time, network copies, idle job waits, visual/manual review, and repeated render versions. Artifact rendering reported separately.
