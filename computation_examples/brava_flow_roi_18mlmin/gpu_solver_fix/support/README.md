# Frozen GPU production support

These are the original deployed GPU wrapper, MPI launcher, RA header and policy, plus the recorded PETSc/svMultiPhysics build recipe. SHA and original server locations are in PROVENANCE.json and production_workflow/GPU_MANUFACTURED_RUN_IDENTITY.json. No environment, binary, shared library or full third-party checkout is shipped.

The RA build executable SHA 0e509fe2… is the **pre host-sync fix** baseline. The formal BraVa and frozen V2 executable SHA is **43cd7d7b9ae9bde91275982225c6bf084859141cdaf2ec1030e1b83c3b94691d** after gpu_host_sync.patch. PETSc 3.25.5 library SHA is b79e98d278e105aa8a1848bd08d66ef36a3e2425458610f97213d3fc8e45c5ef. Preserve patch order, PETSC_OPTIONS and recorded original paths; these wrappers are not a portable installer.

The removed third-party checkout is recoverable from GitHub parent 0b808fb6c90e27674d79078c49358ced2e9258ff at brava_flow_roi_18mlmin/gpu_solver_fix/source/. Use the recorded build/configure/link commands only in a separately authorized environment setup; this cleanup ran none. Project-owned before/after petsc_impl.cpp, the patch and backend_equivalence.json remain outside source/. The older build.json BUILT_NOT_YET_VALIDATED label describes its build-time stage; later backend_equivalence and formal GPU identities are the execution evidence.

Gmsh reference generator imports lazily: native Gmsh 4.15.2 and its recorded binary SHA are external runtime dependencies. A future environment setup may install gmsh==4.15.2 into a separate runtime and verify its SHA against candidate_v2/IDENTITY.json; no runtime or mesh was created by this sync.
