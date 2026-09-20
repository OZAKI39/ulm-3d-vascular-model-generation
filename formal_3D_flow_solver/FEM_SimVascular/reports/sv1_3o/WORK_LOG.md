# Stage SV1.3O completed

Closed 2026-09-19T21:20:22.165657+00:00. Status: PASS — SIMPLE_GPU_TUNING_EXHAUSTED.

Only the authorized A-D fixed windows were attempted. All used the original
Stage N step60 checkpoint and requested steps61-70. A and B completed; B's
0.209572% observed reduction fell below the frozen clear-gain threshold, so A
was retained. C (step61) and D (step62) remain rejected DIVERGED_BREAKDOWN
failures. No confirmation repeat was triggered. One PROFILE_A window served
both baseline and retained winner and was excluded from timing selection.

The unique full REAL_VASCULAR_GPU_PERF run started at t=0, first passed the
unchanged five-interval production gate at step70,
and stopped natively at step71. Process wall 2430.023318998s;
161 KSP solves; 81757 iterations; zero linear/nonlinear failures. Exactly eight
VTUs at steps10,20,30,40,50,60,70,71 were saved. Final native checkpoint and
VTU TimeValue agree, hashes and fresh-process reload pass. No CFD remains
active; native_solver_stopped.json records the process check.

Stage O: 28 passed, zero failures/errors/skips. Full artifact-only pytest:
462 passed, 8 preserved historical failures, 65 skipped. Historical test
identities and counts exactly match frozen Stage N. Local preservation audit
checks 18174 entries and unchanged old FEM/source/production; native audit
confirms frozen drivers, CUDA, MPI, PETSc and historical GPU stacks. All pass.

Five requested figures were visually reviewed. REPORT.md answers the nine
requested sections. TERMINAL_SUMMARY.txt contains section45; delivery_manifest
freezes Stage O files and references adopted source and native binary hashes.

CPU_EARLY_STOP_PRODUCTION is unchanged. CPU/GPU science equivalence, formal
benchmark, CPU comparisons, multi-rank CUDA and advanced GPU preconditioners
remain DEFERRED. No Stage SV1.3P or further CFD was started.
