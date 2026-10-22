B2 XY/YZ/XZ acoustic pressure sections, taller right panels

Fig1b_acoustic_field: PNG/PDF/SVG/TIFF, 180 x 90 mm, 600 dpi.
XY: Z=0 mm; YZ: X=0 mm; XZ: Y=0 mm.
Original magma map, common linear 0-39 kPa color scale.

YZ and XZ are vertically enlarged at user request, by 38.13% and 49.27% versus v2.
Both right panels are 21.8164 mm high; their gap is 13 mm instead of about 26.2 mm.
Their top and bottom align with XY; their left/right edges coincide.
Only the display aspect changes. Coordinates, section samples, color limits and source geometry are unchanged.
The XY plot retains equal spatial aspect. YZ and XZ have vertical exaggeration 1.38135 and 1.49275; axes remain in mm.

CAPTION.txt and PROVENANCE.json record the interpretation and data identity.
source_data.h5 contains exact original complex sections and vessel intersections.
EXPORT_VERIFICATION.json records data, export and alignment checks.
Earlier versions remain archived; no simulation was started.

Reproduce from project root:
PYTHONPATH=src:. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.stage4.paper_acoustic_planes --output results/stage4/publication_figures/acoustic_field_three_planes_v3
The paper_acoustic_volume companion supplies common style and paths; no 3D rendering is invoked.
