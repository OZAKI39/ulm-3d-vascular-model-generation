# Reproduce this diagnostic review

Use the existing server project Python, with this directory as cwd. Read `/etc/vast-agents-guide.md` first. Runtime package installation stays inside `runtime/rt_dependencies`; no system environment changes.

Run `python supervise.py identity`, `python supervise.py terms`, `python supervise.py install_rt`, `python supervise.py rt_reference`, `python supervise.py source_audit`, `python supervise.py plot_results`, `python supervise.py run_tests` in this order. Each invocation exits; there is no automatic next job. Existing attempts remain in monitoring and count toward the 5400 s cumulative budget. These scripts do NOT run CFD, produce meshes, or integrate particles.

`source_audit.py` documents STOP_C for these frozen term results; it is not a general decision engine for other cases. `write_reports.py` packages this specific frozen review, using the upstream identity receipt. Scientific numbers must not be reused for other field identities. The scripts assume the recorded immutable mesh, field and event protocol.

`python write_reports.py` regenerates text reports. Inspect both PNG files, record their SHA in VISUAL_QC.json, then run `python finalize_delivery.py` after all supervised jobs exit. The finalizer verifies inputs, checks no remaining own processes and packages receipts; it never launches a solver. Local files mirror server results byte-for-byte. Output manifest excludes itself, runtime dependencies, monitoring, and temporary TASK.json. Raw supervisor logs remain available with hashes in EXECUTION_RECEIPT.json.

References and formulas: CITATIONS.json and STABILIZATION_AUDIT_ZH.md. RT is a constrained velocity-L2 diagnostic, not a Darcy hydraulic-gradient reproduction or a production particle input. Independent Git commit only, no push.
