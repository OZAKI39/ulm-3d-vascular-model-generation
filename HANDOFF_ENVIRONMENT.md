# Handoff environment and lightweight reproducibility

2026-10-06. Current lightweight queries are saved in [ENVIRONMENT.json](github_handoff_20261006/ENVIRONMENT.json). No packages were installed or upgraded, and no server environment is copied into this branch.

| Component | Version / role |
|---|---|
| Main Python | 3.13.11, Anaconda; `/home/lzy/projects/computation_examples/_shared/python_environment/bin/python` |
| NumPy / SciPy | 2.5.3 / 1.18.1 |
| matplotlib / pytest / psutil | 3.11.2 / 9.1.1 / 7.2.2 |
| Optional legacy packages | numba 0.67.0; pyvista 0.49.0 |
| PETSc environment | Separate Python 3.12.13 at `/home/lzy/.cache/codex-dolfinx-env/bin/python` |
| PETSc / petsc4py | 3.25.3 / 3.25.3, queried without a solve; CUDA disabled |
| Julia | 1.10.12, recorded in accepted RK cross-check; not rerun for handoff |
| Peridynamics.jl | 0.6.0-DEV, commit `b6e86563845a83b5cb06c805af19a817bde6b106`, current local commit checked |
| Host | WSL2 Linux x86_64; Intel Core i7-13700HX, 24 logical CPUs visible |
| RAM | 7,991,452 KiB total (about 7.62 GiB); solver concurrency must follow available memory, not CPU count |
| Backend | NumPy/SciPy CPU, `GPU_NOT_USED_BY_DESIGN` |

The current PD2 intact suites need NumPy, SciPy and pytest. PETSc is a historical phase/VI dependency and is not needed by these PD2 checks. Historical provenance may record different environments; those records are not rewritten.

All commands below start in `clot_clearance/`. Set the interpreter to an environment with the listed packages; the absolute server path is descriptive, not a portable requirement.

```bash
export PYTHONPATH=.
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

Server check already performed during handoff:

```bash
python -B -m pytest -q -p no:cacheprovider \
  verification/higher_order_baqp_rk_oracle_v1/tests/test_pd2_oracle.py \
  verification/higher_order_baqp_energy_force_candidate_v1/tests/test_variational_candidate.py \
  verification/higher_order_baqp_boundary_reintegration_v1/tests/test_boundary_reintegration.py
```

Result: **36 passed in 4.42 s**, using the existing local artifacts. See [server test log](github_handoff_20261006/SERVER_TARGETED_TESTS.log). These tests read saved scientific results and perform small live mechanics checks; no equilibrium campaign is run.

Portable subset using only selected branch files, with unchanged tests and assertions:

```bash
python -B -m pytest -q -p no:cacheprovider \
  verification/higher_order_baqp_energy_force_candidate_v1/tests/test_variational_candidate.py \
  -k 'not test_01_cached_PD2_kinematics_preserved'
python -B -m pytest -q -p no:cacheprovider \
  verification/higher_order_baqp_boundary_reintegration_v1/tests/test_boundary_reintegration.py \
  -k 'not test_01_pd2_mechanics_cache_and_force_regression_unchanged and not test_04_nonaffine_nodewise_decomposition and not test_12_classical_shear_target_residual'
```

The first command selects 11 tests; the second selects 9. The omitted tests require the excluded N19 cache, original protected local files, or the 3.44 MB nodewise boundary array. They are explicitly deselected for this portable subset, not weakened, skipped by source edits, or claimed verified from incomplete input. Results of an isolated file-copy check are recorded in [PORTABLE_TARGETED_TESTS.json](github_handoff_20261006/PORTABLE_TARGETED_TESTS.json) and its log.

Full scientific reruns and all original plotting entrypoints require the local-only inputs listed in [LOCAL_ARTIFACTS_NOT_SYNCED.md](LOCAL_ARTIFACTS_NOT_SYNCED.md). Reconstructing a reference cache is possible from included source/configuration, but byte identity is not promised across environments; this handoff does not rebuild it. Do not run `run_audit.py` merely to inspect the project: some stages launch long equilibrium solves and write scientific output paths.

External source provenance: [EXTERNAL_DEPENDENCIES.md](github_handoff_20261006/EXTERNAL_DEPENDENCIES.md).
