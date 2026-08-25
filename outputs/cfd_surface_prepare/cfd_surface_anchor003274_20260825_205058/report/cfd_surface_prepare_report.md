# CFD surface preparation report

This run derives a CFD-only surface from the validated Ultraliser lumen. The source STL was treated as immutable; only the four saved boundary neighborhoods were cut, extruded, and capped.

- Input preprocess run: `E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_preprocess\global_to_roi_anchor003274_20260825_183628`
- Original STL unchanged: `True`
- Whole-surface reconstruction performed: `False`
- Volume mesh created: `False`
- 3D CFD run: `False`

## Boundary construction and numerical pressure correction

The outlet correction accounts only for the predicted loss in each artificial straight extension. It is not a physiological outlet model, resistance boundary condition, or Windkessel model. Negative gauge pressure is valid where produced.

- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_000`: ASSUMED_INLET, CUT_PORT, L=15.8739 um, r_eq=1.57544 um, P_original=216.456832 Pa, Q_expected=7.69350848e-16 m3/s, P_solver=None
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_001`: ASSUMED_OUTLET, CUT_PORT, L=11.7354 um, r_eq=1.18228 um, P_original=16.9774127 Pa, Q_expected=4.84870237e-17 m3/s, P_solver=14.416492037919472
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_002`: ASSUMED_OUTLET, CUT_PORT, L=11.5286 um, r_eq=1.17065 um, P_original=152.029792 Pa, Q_expected=3.81521158e-16 m3/s, P_solver=131.4357629035375
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__terminal_000`: ASSUMED_OUTLET, TRUE_TERMINAL, L=12.55 um, r_eq=1.27427 um, P_original=0 Pa, Q_expected=3.39342665e-16 m3/s, P_solver=-14.20315070469081

## Quality control

- Surface topology: `PASS`
- Core preservation: `PASS`; max 4.26326e-14 um, P95 2.00972e-14 um
- Extension collision QC: `PASS`; count 0

## Manual review

- STL: `E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_surface_prepare\cfd_surface_anchor003274_20260825_205058\geometry\cfd_surface_extended_um.stl`
- Tagged VTP: `E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_surface_prepare\cfd_surface_anchor003274_20260825_205058\geometry\cfd_surface_extended_um.vtp`
- Status: `CFD_SURFACE_PREPARE_PASS_PENDING_MANUAL_REVIEW`
- Next stage: `MANUAL_CFD_SURFACE_REVIEW`
