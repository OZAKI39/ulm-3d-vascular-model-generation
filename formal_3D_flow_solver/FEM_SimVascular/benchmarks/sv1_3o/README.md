# Stage O bounded GPU tuning evidence

These are single-run development measurements, not a formal benchmark.
The frozen decisions and complete PETSc options are in
`configs/sv1_3o/policy.json`, with observed decisions in
`reports/sv1_3o/winner.json`.

All A-D trials restart the same full native GPU step60 checkpoint and request
the same ten timesteps through step70. The measurement is the new native
process monotonic wall time; the solver's printed cumulative elapsed time
contains inherited checkpoint timing and must not be used as window wall.
Incomplete failed windows are excluded from timing comparisons.

The selection rule accepts a healthy reduction of at least 10%, treats less
than 5% as no clear gain, and permits only one confirmation in the 5%-10%
band. A and B completed; B's 0.21% observed reduction did not replace A.
C and D failed with DIVERGED_BREAKDOWN and were rejected. There were no
confirmation repeats and no parameter sweep beyond A-D.

One separate PROFILE_A window serves both the baseline and retained winner.
Its log_view_gpu_time timing is excluded from winner selection. The sole
full-run plan is REAL_VASCULAR_GPU_PERF; it uses the unchanged production
monitor and native STOP_SIM. Final closed evidence is recorded in
`reports/sv1_3o/optimized_steady_candidate.json` and REPORT.md.

Artifact tests run with `.venv/bin/python -B -m pytest -q -p no:cacheprovider
tests/test_sv13o_*.py`. They read existing evidence and never launch CFD.
Execution scripts deliberately reject overwriting existing run evidence;
do not rerun them as part of report review.
