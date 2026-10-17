# Current entry index

| Entry | Purpose |
| --- | --- |
| src/distribution.py | Read-only adapter to original OpenCV detector; center, density, covariance; no tracker |
| src/methods.py | Official OpenPIV and pycpd calls, quality masks |
| src/affine.py | Huber affine displacement fit, homogeneous matrix logarithm, rank/coverage/uncertainty |
| src/state.py | M4 messages, independent valid masks, stale/discontinuity gate and local ZMQ mock |
| src/synthetic.py | Real image patches plus explicit analytic visual motion |
| configs/protocol.json | Preregistered seeds, windows, limits and budgets |
| configs/collective_state.schema.json | Machine-readable state and future action-event schema |
| scripts/benchmark.py | Development selection and frozen independent synthetic test |
| scripts/real_observation.py | Offline causal analysis, finite camera check, latest-slot replay |
| scripts/ai_comparison.py | One official SEA-RAFT-S checkpoint, 24 pairs, no training |
| scripts/simulation_bridge.py | Read-only actual F10S HDF5 projection |
| scripts/window_response.py | Causal step-response diagnosis |
| scripts/supplement_checks.py | Local defocus, false candidates, ZMQ fault injection and schema |
| scripts/analyze.py | Truth metrics with explicit accepted/rejected denominators |
| scripts/visualize.py | Five figure groups and two PTS-preserving videos |
| scripts/write_reports.py | Reports assembled from recorded data |
| scripts/finalize.py | Integrity checks, progress-document backup/diff, output hashes |

No temporary historical entry was proven both unused and disposable, so no old code or data was deleted. Old GUI/annotation/physics entry points remain in their original index; this index is the F10U start point.
