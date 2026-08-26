# VMTK TPS boundary-normal direction experiment

- Status: `VMTK_TPS_BOUNDARY_NORMAL_INTERFACE_WARNING_PENDING_MANUAL_REVIEW`
- VMTK runtime: `1.5.0`
- VTK runtime: `9.2.6`
- Official filter: `vtkvmtkPolyDataFlowExtensionsFilter`
- Interpolation: `thinplatespline`
- Custom TPS implementation: `false`
- Previous extension mode: `centerlinedirection`
- Current extension mode: `boundarynormal`
- Official direction API: `SetExtensionModeToUseNormalToBoundary`
- Centerlines used by VMTK: `false`
- Transition ratio: `0.5`
- Preserve cross-section shape: `false`
- Parameter sweeps/fallbacks: `none`
- Visible fold assessment: `MANUAL_REVIEW_REQUIRED`

## Boundary measurements

- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_000`: axial length=15.7082 um; area error=0.3573%; direction dot=0.999999999; interface P95/P99=9.5172/11.06 deg
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_001`: axial length=11.7974 um; area error=0.864%; direction dot=0.999999998; interface P95/P99=6.7318/7.4201 deg
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__cut_002`: axial length=11.6963 um; area error=0.6234%; direction dot=1.000000000; interface P95/P99=8.0209/8.7107 deg
- `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274__terminal_000`: axial length=12.6924 um; area error=0.3275%; direction dot=0.999999999; interface P95/P99=8.6732/9.1418 deg

## Final checks

- RAW topology: `PASS`
- Remeshed/capped topology: `PASS`
- Radius P95: `0.03999745507192995`
- Core P95/max (um): `0.009902434252990231` / `0.02449571667958121`
- Manual-review STL: `E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_surface_prepare\vmtk_tps_boundarynormal_anchor003274_20260826_153709\geometry\cfd_surface_vmtk_tps_boundarynormal_um.stl`
- Tagged VTP: `E:\ULM\hatimb-particle_flow_simulator\ulm_3D_vascular\outputs\cfd_surface_prepare\vmtk_tps_boundarynormal_anchor003274_20260826_153709\geometry\cfd_surface_vmtk_tps_boundarynormal_um.vtp`
- Volume mesh created: `false`
- CFD run: `false`

NEXT: `MANUALLY REVIEW VMTK BOUNDARY-NORMAL TPS SURFACE`
