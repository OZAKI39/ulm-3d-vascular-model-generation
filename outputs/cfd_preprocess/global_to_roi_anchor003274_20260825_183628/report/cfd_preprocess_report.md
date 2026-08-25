# CFD boundary preprocessing report

Final status: **CFD_PREPROCESS_BASELINE_PASS**

ROI: `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274`

Global source model: `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01`

## Formal baseline and its limits

SWC parent→current direction is used only as the **simulation direction**. It is not an
experimentally measured blood-flow direction and is not physiological ground truth. The single
structural root is temporarily treated as `ASSUMED_GLOBAL_INLET`. Its configured
0.7 mm/s velocity is a literature-derived baseline assumption,
not a measurement from the current fMOST mouse.

All structural leaves are assigned 0 Pa gauge pressure as the global baseline
reference. This
does not mean that physiological blood pressure is zero. The 1D calculation uses original
analysis-SWC source radii; the Ultraliser radius scale of 0.91 is only a surface-reconstruction
compensation and is not used in hydraulic resistance.

ROI inlet flow and outlet pressure are transferred from the global 1D model. The intended outlet
condition is **DIRECT 1D PRESSURE**, not a resistance or Windkessel condition. Version-2
physiological refinements are disabled.

TRUE_TERMINAL boundaries are treated as `ASSUMED_OUTLET` in the current baseline. This is a
simulation convention used to keep the boundary model simple. It does not assert that these
terminals are experimentally verified physiological outlets. Their pressure is read from the
global 1D structural-leaf solution; the current gauge value is a reference pressure, not a claim
that physiological blood pressure equals zero.

## Global solution

- Structural root: 2410
- Structural leaves: 122
- Root flow: 3.79491284792e-14 m³/s
- Global relative mass error: 2.381396e-13
- Maximum internal relative residual: 1.143511e-13

## Complete CFD boundary transfer

- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_000` — CUT_PORT / ASSUMED_INLET; P=216.456832 Pa; Q=7.69350848e-16 m³/s
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_001` — CUT_PORT / ASSUMED_OUTLET; P=16.9774127 Pa; Q=4.84870237e-17 m³/s
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_002` — CUT_PORT / ASSUMED_OUTLET; P=152.029792 Pa; Q=3.81521158e-16 m³/s
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__terminal_000` — TRUE_TERMINAL / ASSUMED_OUTLET; P=0 Pa; Q=3.39342665e-16 m³/s

ROI relative boundary mass error: 4.780737e-14

Readiness failure reasons: none

## Geometry reference

- Model run: `ultraliser_anchor003274_20260825_133350`
- Geometry status: PASS
- Radius scale: 0.91
- Radius P95 absolute relative error: 0.040113538

No surface was modified, no volume mesh was created, and no 3D CFD or microbubble simulation was
run in this stage.

## Next stage

CFD SURFACE / VOLUME MESH PREPARATION.
