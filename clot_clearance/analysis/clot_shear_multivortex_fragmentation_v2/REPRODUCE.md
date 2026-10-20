# Reproduce clot_shear_multivortex_fragmentation_v2

Accepted run 20261017T175456Z is immutable. Use a NEW task/run directory; never overwrite its results.

Activate /workspace/y_clot/env.sh. Set single-thread BLAS/OpenMP and LD_LIBRARY_PATH to the frozen library directory. Library SHA: 465f8144c6da5c9de3aec4137b9c325f53c26737b199c536fb4a20845511ef77.

Chain: shear_field_library.py → current_surface.py → native.py → run_cases.py → analyze.py.
work_controller.py changes one global amplitude only.

Change isolated RUN before reproducing. Copy exact original inputs; run prepare.py, checks.py, native_checks.py. Freeze and read back every SHA. Geometry fixtures are not fragmentation evidence; native checks take one microstep per mode.

Run the TWO specified modes once each under supervisor, from the identical original state. No additional uniform case or outcome-based calibration trials.

Read-only postprocessing: python analyze.py; python figures.py; python report.py.

Initialization failures and source records are retained in startup_failure_00/ and preflight_failure_01/. Do not hide process attempts.

Full accepted scalar/damage/component records and FP64 key states are retained. Work is actual endpoint force dot displacement; no nominal-area substitution or curve scaling. Unknown true surface area remains unknown.
