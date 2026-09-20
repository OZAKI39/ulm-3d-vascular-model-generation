# PETSc 3.25 compatibility adapter — repair_01

WSL isolated source: `external/petsc325/svMultiPhysics-compat`, base commit `c3f0bb892b765b718f61069ecd9726dbc6d177fd`.

Initial PETSc 3.25 CPU compile and official smoke passed without source changes. Source inspection shows the application still omits PetscFinalize. Stage M demonstrated the resulting GPU teardown error with this same fixed svMP source (see frozen `exit_failure_diagnosis.json` and `finalize_ksp_probe.json`). This is an existing lifecycle defect exposed by GPU resource cleanup, not a claim that upgrading PETSc introduced it.

Official contract: https://petsc.org/release/manualpages/Sys/PetscFinalize/ and pinned v3.25.5 `src/sys/objects/pinit.c`. PETSc must clean up before application-owned MPI is finalized. Because svMP calls MPI_Init first, PetscFinalize does not take ownership of MPI_Finalize.

The patch wires the existing finalize interface to checked cleanup, remembers the allocated equation count, and performs process-wide PETSc cleanup exactly once across equation interfaces. Solver context allocation is zero-initialized so unused/non-PETSc equation slots or KSP-only slots have safe null handles. PETSc null-handle destroy and null-pointer free semantics are verified in pinned source. Every newly invoked destruction/finalization operation checks errors with PetscCallAbort; no error is ignored. The existing invalid-LHS abort remains.

No equations, assembly, residual, mesh, physics, boundary conditions, timestep, convergence threshold or solver options change. Original PETSc source and historical solver trees remain unmodified. The adapter was clean-built against new 3.25.5 CUDA. The latest user strategy cancelled further CPU builds and comparisons; old 3.19.6 compatibility is not claimed as tested.

The single official GPU smoke includes read-only GDB tracing: at application MPI_Finalize, PetscInitializeCalled=0, PetscFinalizeCalled=1, actual PetscFinalize call count=1; the inferior exits normally. The same run independently passes solver convergence, finite VTU and fresh reload. A separate tiny CUDA-object lifetime probe also exits0. Scientific comparisons are DEFERRED BY USER DECISION under USER_GPU_PRIORITY_UPDATE.txt. Actual statuses and affected hunks are in `reports/sv1_3n/compatibility_adapter.json`.

## repair_02 — build configuration and observation parser

The new Stage N build script initially used CUDAC_FLAGS as a configure argument. Pinned `config/BuildSystem/config/setCompilers.py:151` declares CUDAFLAGS. The corrected clean source build configured, compiled, installed and passed CPU1/MPI2/CUDA checks. A case-sensitive parser then misread uppercase CUDA; complete_gpu_build_remote.py revalidated the original exit0 log and finished link/source provenance without repeating any passed build/check. No PETSc source patch or scientific setting changed. Original rejected command/log and all later successful logs remain in reports/sv1_3n/remote and logs/sv1_3n/remote.

[build_config_parser_repair.patch](build_config_parser_repair.patch) records the two expression corrections and the added completion script; exact affected-line counts and file hashes are in `compatibility_adapter.json`. Its before-side is explicitly a minimal reconstruction of those expressions from retained failure evidence, not a byte-for-byte historical script snapshot. `scripts/sv13n/record_build_repair.py` reproduces this documentation without running a build or solver. This completes the existing repair_02 record; it introduces no further compatibility change.
