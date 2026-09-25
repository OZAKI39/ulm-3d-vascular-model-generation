# BG001 RMCA manufacturing candidate

NOT manufacturer guaranteed limits. Must be calibrated experimentally.

Status: **PRINT_READY_CANDIDATE**. Physical printing/PDMS/ABS dissolution have NOT been validated. No print jobs submitted.

- Layer A: immutable `../refined_roi/semantic_refined/`. Historical `diameter_075/` is retained and not used here.
- Layer B: `BG001_RMCA_print_roi_original_radius.swc`: native manufacturing topology and raw radii.
- Layer C: `BG001_RMCA_print_roi_compensated.swc`: MANUFACTURING_COMPENSATED_GEOMETRY.
- `VascularMD/`: both native fits, QC and original native surfaces.
- `BG001_RMCA_print_candidate.stl`: selected compensated print derivative, ports capped, rotated/translated, scale=1.
- `BG001_RMCA_print_candidate.3mf`: present only when real Bambu slicing succeeded.
- `BambuStudio/`: exact commands, profile inheritance, actual slice metadata and all candidate 3MFs.
- `orientation/`: preserved identity baseline, Top 5 STL, all transforms and scores.
- `bambu_abs_vessel_calibration_coupon.stl`: 25 separate calibration cylinders; CSV maps diameter/angle to row/column.

Run from project root (choose a new output directory):

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tools/manufacture_bg001_rmca.py --output outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/manufacturing_run02
```

Alternative detail profile: `--config config/bambu_abs_0p2_detail.yaml`. This does not change your physical nozzle. Stages `--stage prepare`, `--stage model`, `--stage finish` allow checkpointed operation; `all` is the default. A new full run refuses an existing output directory. Existing semantic/refinement results are always protected.

Optional contexts can be requested by component ID in `roi.optional_context_components`; they remain context and are never relabelled as core.
