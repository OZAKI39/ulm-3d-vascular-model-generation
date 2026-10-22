# F5 delivery

Status and scope: see [FINAL_GATE.json](FINAL_GATE.json). Full development and interpretation: [report](../../../docs/FAST_PZT_F5_FLOW_RELEASE_DEVELOPMENT.md).

The actual main run forms a 32-bubble cluster under continuous flow, then switches **all ultrasound OFF at 300 ms** while flow continues to 600 ms. The release moves the centroid 24.66 µm and increases fixed-ID Rg from 4.92 to 12.64 µm. Eleven persistent graph splits are observed; the last additional split is right-censored. The same-state no-flow control remains essentially stationary. Graph splitting is not physical bubble fragmentation.

- [60 fps MP4](visualization/Video_F5_acoustic_release_flow_60fps.mp4)
- [Flow and geometry](visualization/Fig_F5_1_flow_original_Y.png)
- [Initial → ON → OFF](visualization/Fig_F5_2_attraction_release.png)
- [Flow/no-flow release comparison](visualization/Fig_F5_3_release_comparison.png)
- [Main actual trajectory](science_v1/main/trajectory.h5)
- [Exact switch state](science_v1/SWITCH_STATE.npz)
- [Branch A, exact main suffix](science_v1/branch_A_actual_suffix.npz)
- [Branch B, no-flow continuation](science_v1/branch_B/trajectory.h5)
- [Comparison](ON_OFF_FLOW_COMPARISON.json), [transport and shape](CLUSTER_TRANSPORT_METRICS.json), [performance](PERFORMANCE.json)

All 1,201 frames use saved positions without animation interpolation. Sphere radius is 1.78 µm, with fixed cameras and actual coordinate scales. The complete video decode and uniform presentation timestamps are verified. Exact CFD, basis, response and geometry inputs are included; runtime environments and source hashes are recorded.

Retained attempts: `checks_v1` stopped for missing PyVista; `checks_v2` passed. `science_v1` completed both trajectories before a timing-only exception. `performance_v1` recovered only that instrumentation, preserving trajectory hashes. `checks_v3` passed final diagnostics and archived-trajectory prefix regressions. No scientific rerun or parameter search was performed.

`EXPLORATORY / PRESCRIBED_STEADY_FLOW / GEOMETRIC_WALL_ONLY`. The inherited CFD uses water parameters and has no mesh-independent local validation. See the report for mobility/contact and physical limitations.
