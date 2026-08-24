# v9 strict reconstruction report

- ROI: `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cfd_domain`
- status: **FAIL**
- decision: **KEEP_V7_HARD_MIN_V9_ACCEPTANCE_FAILED**
- root cause: `PIECEWISE_LINEAR_CENTERLINE_SEGMENT_CREASE`
- selected spline theta: `0.5 deg`
- selected spline sagitta: `0.01 r`
- selected k/r: `None`

## Required answers

1. S0 defect-near same-branch segment-switch fraction=92.636% (adjacent=91.702%, nonadjacent=5.336%); root cause=PIECEWISE_LINEAR_CENTERLINE_SEGMENT_CREASE.
2. Maximum resampled-centerline tangent kink=49.9188 deg at branch 4, source neighborhood 49;83;84, segments 974/975.
3. C1 knot tangent discontinuity=0 deg; selected dense-polyline maximum turn=0.49752 deg.
4. Smooth-centerline-only changed same-branch gradient P99 from 18.4846 deg to 14.526 deg, visible count from 313 to 705, and silhouette roughness from 1.44074 to 1.42147.
5. Competition-aware union at measured k/r=0.50 changed cross-branch gradient P99 from 93.8152 deg to 60.7908 deg, visible count from 705 to 304, and silhouette roughness from 1.42147 to 1.43197.
6. Final fairing requirement/status=APPLIED; after fairing Newton/hard-min reprojection=False.
7. Same-camera flat/wireframe/silhouette generated=True; quantitative visible-defect disappearance=False.
8. Source arrays/topology modified=False; maximum Hausdorff=0.172397 um (0.172397r), maximum bidirectional P95=0.143636 um, maximum length change=0.199193%, junction/endpoint/port error=0 um.
9. Final radius P95=0.00678531; hydraulic error=0; acceptance radius/hydraulic checks=True/True.
10. Recommended backend=unified_polyball; decision=KEEP_V7_HARD_MIN_V9_ACCEPTANCE_FAILED; strict multi-criterion pass=False.

## Final acceptance

- [x] same_camera_flat_wireframe_generated
- [ ] visible_sawtooth_disappeared
- [ ] same_branch_segment_switch_defect_significantly_removed
- [ ] cross_branch_switch_smoothed
- [x] radius_p95_below_one_percent
- [x] hydraulic_error_within_tolerance
- [x] junction_volume_within_tolerance
- [x] topology_all_pass
- [x] surface_qc_pass
- [x] source_centerline_fidelity_pass
- [x] junction_position_exact
- [x] endpoint_position_exact
- [x] port_position_exact
