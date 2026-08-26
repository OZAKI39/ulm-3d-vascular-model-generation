# Refined CFD extension surface

This run keeps the validated Ultraliser core and every proximal cut ring locked. Only the artificial extension surface is refined with multiple axial rings, a one-diameter shape transition, quality-driven quad diagonals, and constrained Taubin smoothing.

- Input preprocess run: `E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_preprocess\global_to_roi_anchor003274_20260825_183628`
- Original STL unchanged: `True`
- Previous direct-extrusion reference unchanged: `True`
- Whole-surface smoothing/remeshing/reconstruction: `False`
- Volume mesh created: `False`
- 3D CFD run: `False`

## Mesh method and quality

Aspect ratio definition: `(sqrt(3)/2) * longest_edge / altitude_to_longest_edge; equilateral=1`. The refined extension is intended to improve surface triangle quality before volume meshing.

- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_000`: rings=50, target=0.322426 um, triangles=5880, P95 aspect=2.28016, bad fraction=0.00538194, P_solver=None
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_001`: rings=46, target=0.25948 um, triangles=5580, P95 aspect=2.62257, bad fraction=0.00549853, P_solver=14.416492037919468
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_002`: rings=52, target=0.223693 um, triangles=6630, P95 aspect=2.33819, bad fraction=0.00230769, P_solver=131.4357629035375
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__terminal_000`: rings=43, target=0.297669 um, triangles=4620, P95 aspect=2.30156, bad fraction=0.00221729, P_solver=-14.203150704690803

## Quality controls

- Extension mesh quality: `PASS`
- Surface topology: `PASS`
- Core closest-point preservation: `PASS`; max 4.26326e-14 um, P95 2.00972e-14 um
- Direct locked-original-vertex motion: `PASS`; max 0 um
- Extension collision QC: `PASS`; count 0

## Manual review

- Refined STL: `E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_surface_prepare\cfd_surface_refined_anchor003274_20260826_125703\geometry\cfd_surface_refined_um.stl`
- Tagged refined VTP: `E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_surface_prepare\cfd_surface_refined_anchor003274_20260826_125703\geometry\cfd_surface_refined_um.vtp`
- Status: `CFD_EXTENSION_MESH_REFINED_PENDING_MANUAL_REVIEW`
- Next stage: `MANUAL_REFINED_CFD_SURFACE_REVIEW`
