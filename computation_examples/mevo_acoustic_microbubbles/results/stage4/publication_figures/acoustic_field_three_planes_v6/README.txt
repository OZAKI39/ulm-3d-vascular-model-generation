B2 XY/YZ/XZ acoustic pressure sections, undistorted tall panels

Fig1b_acoustic_field: PNG/PDF/SVG/TIFF, 180 x 90 mm, 600 dpi.
XY: Z=0 mm; YZ: X=0 mm; XZ: Y=0 mm.
Original magma map; common linear 0-39 kPa color scale.

All three planes have equal horizontal/vertical physical scale, verified from the final drawing transform.
The vessel sections are unchanged from v2 and are not stretched.
YZ/XZ Z-axis windows are expanded to +/-2.3 mm to retain taller panels.
The stored field covers only +/-1.6 mm: light-grey bands denote no data, not zero pressure.
There is no field extrapolation. Panel widths and overall upper/lower boundaries remain aligned with XY.
The gap between the right panels is 12.92065 mm.

CAPTION.txt and PROVENANCE.json describe display choices and source identity.
source_data.h5 contains exact original complex sections and vessel intersections.
EXPORT_VERIFICATION.json records source data, vessel shape, equal-aspect, export and alignment checks.
Earlier versions remain archived; no simulation was started.

Reproduce from project root:
PYTHONPATH=src:. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.stage4.paper_acoustic_planes --output results/stage4/publication_figures/acoustic_field_three_planes_v6
The paper_acoustic_volume companion supplies common style and paths; no volume rendering is invoked.

Revision v5: removed the two no-data text labels from the figure at user request. The pale-grey empty bands, equal spatial aspect, alignment, source data and color scale are unchanged.

Revision v6: pressure amplitude symbol uses proper absolute-value delimiters with explicit thin spaces around P. All other display settings and data are unchanged.
