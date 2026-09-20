# Official PETSc ghost-vector audit

Official repository: https://gitlab.com/petsc/petsc.git

First base-type dispatch fix: [8ba8e283c7d3ba8c2bd73cc674b79590a073ab5a](https://gitlab.com/petsc/petsc/-/commit/8ba8e283c7d3ba8c2bd73cc674b79590a073ab5a), 2026-03-13, first released tag **v3.25.0**. Commit rationale: “Vec: add device support for ghost vectors”; closes official issue [1096](https://gitlab.com/petsc/petsc/-/issues/1096). The issue identifies both strict type matching and destroyed ghost metadata during backend conversion. Full git log, blame, commit diff and tag diff are in `upstream/`, with commands/hashes in `upstream_audit.json`.

Old 3.19.6:
```c
PetscCall(PetscObjectTypeCompare((PetscObject)g, VECMPI, &ismpi));
PetscCall(PetscObjectTypeCompare((PetscObject)g, VECSEQ, &isseq));
```
First fixed source:
```c
PetscCall(PetscObjectBaseTypeCompare((PetscObject)g, VECMPI, &ismpi));
PetscCall(PetscObjectBaseTypeCompare((PetscObject)g, VECSEQ, &isseq));
```

`seqcuda`/`mpicuda` are derived types. Concrete equality misses both; the existing base-aware helper matches their seq/mpi families. The Stage M pre-patch CPU probes pass; CUDA probes record NULL local form and UpdateBegin error.

repair_01 changes only five calls in GetLocalForm, UpdateBegin, UpdateEnd. RestoreLocalForm has no concrete type dispatch; it only synchronizes state and releases a reference, so needs no change. `VecGhostIsLocalForm` is unused in the solver/probe path and is not expanded speculatively. The upstream control-flow style rewrite is omitted because it is semantically identical.

This fixes recognition, but does **not** yet establish that MPI CUDA conversion retains ghost metadata. That is explicitly checked by the two-rank data probe; a second controlled repair requires independent evidence. Sequential solver runs have zero ghost entries; two-rank tests verify nonzero remote ghosts and reverse contributions element by element.

## repair_02 — conversion retains ghost metadata

Official [ca0374c5dd2e044b693bca7c31c716bd3be32b2a](https://gitlab.com/petsc/petsc/-/commit/ca0374c5dd2e044b693bca7c31c716bd3be32b2a), first tag v3.20.0, adds in-place CUDA conversion. The complete 388-line upstream change is **not** imported. Only the MPI CUDA route is adapted using the existing 3.19 initializer. It preserves the existing host allocation and Vec_MPI metadata, and sets the CUDA operations/type. Forward/reverse data checks remain mandatory. See `repair_02.json` and `upstream/conversion_*.txt`.

## repair_03 — host local-form coherence (local adaptation)

The second repair removes API errors but **fails actual reverse data**: rank 0 sees [1,2] instead of [111,113]. The source shows a host `seq` localrep sharing MPI's host allocation while CUDA owner updates invalidate it. Ghost state-number synchronization alone does not transfer data. The local reproducer proves this independently of svMultiPhysics.

Official PETSc 3.19.6 `VecGetLocalVector` maps data with `VecGetArray`, and its restore counterpart publishes writes. The third repair adapts these semantics inside ghost Get/Restore for `mpicuda` with a non-null localrep. Get synchronizes owned data to the existing host view; Restore uses write access to publish local writes without overwriting them from device. Global Vec remains `mpicuda`; no Vec type conversion to CPU, ghost-communication bypass, or solver change is used. This preserves the native CPU localrep architecture, with potential transfers explicitly subject to profiling. CPU paths are unchanged.

This coherence addition is **not** represented as a verbatim upstream commit. It is an independently reproduced, local compatibility repair under user section 0, grounded in official source (`upstream/official_host_mapping_3.19.6.c`). The two true upstream commits remain separately attributed. Added device-arithmetic/owned-writeback/duplicate checks prevent accepting merely suppressed errors.
