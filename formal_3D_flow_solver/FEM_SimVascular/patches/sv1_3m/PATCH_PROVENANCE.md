# PETSc 3.19.6 + upstream CUDA ghost compatibility backport

Base archive SHA256: `6045e379464e91bb2ef776f22a08a1bc1ff5796ffd6825f15270159cbb2464ae`.

repair_01 upstream: `8ba8e283c7d3ba8c2bd73cc674b79590a073ab5a`, first fixed tag `v3.25.0`.

Reason: derived CUDA Vec types are rejected by concrete ghost dispatch.

Modified file: `src/vec/vec/impls/mpi/commonmpvec.c`. Functions: GetLocalForm, UpdateBegin, UpdateEnd. Five removed lines / five added lines, exact upstream substitutions. No RestoreLocalForm change.

Scientific impact: **NONE — compatibility only**. svMultiPhysics source and all model inputs remain unchanged. Original PETSc tree/archive are read-only; clean extraction is patched in a separate tree.

Revert: in a disposable clean source tree, `patch -R -p1 < petsc319_cuda_ghost_backport.patch`; then discard all build objects and rebuild. Never apply/revert in historical trees.

## Cumulative repair_03 patch

repair_02 is a CUDA MPI-only adaptation of upstream `ca0374c5dd2e044b693bca7c31c716bd3be32b2a` (v3.20.0). It uses the existing 3.19 CUPM initializer to retain host allocation and ghost metadata. repair_03 adds local-reproducer-supported host-view coherence using official GetArray/RestoreArrayWrite semantics; this part is **not an exact upstream cherry-pick**. Full rationale is in `reports/sv1_3m/repair_03.json`.

Cumulative scope: 4 files, 6 functions, +54/-5 lines. Modified files: src/vec/vec/impls/mpi/commonmpvec.c, include/petsc/private/vecimpl.h, src/vec/vec/interface/vecreg.c, src/vec/vec/impls/mpi/cupm/cuda/vecmpicupm.cu. A header declaration is included. Per-iteration patches are preserved under repair_01/02/03; the top-level patch is cumulative. Modified hunk line locations are recorded in each unified diff. Scientific impact remains **NONE — compatibility only**.
