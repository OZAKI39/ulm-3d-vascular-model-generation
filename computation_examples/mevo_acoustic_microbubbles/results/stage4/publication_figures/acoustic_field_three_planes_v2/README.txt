B2 orthogonal acoustic pressure sections, aligned revision

Fig1b_acoustic_field: PNG/PDF/SVG/TIFF at 180 x 90 mm, 600 dpi.
XY: Z=0 mm; YZ: X=0 mm; XZ: Y=0 mm.
Original magma (purple-to-yellow) map; common linear 0-39 kPa scale.
White vessel outlines and background; English Arial labels.

Layout: YZ top aligns with XY top; XZ bottom aligns with XY bottom.
The two right panels have identical left/right edges.
All four alignment residuals are zero at 600 dpi.
Each plane preserves equal spatial aspect. YZ magnification is 13.4/12.4 times that of XY/XZ; coordinates are in mm.

Data are exact sections of the existing B2 complex pressure field.
No new simulation, no smoothing, no source-data changes.
CAPTION.txt explains the figure; PROVENANCE.json records data identities.
source_data.h5 contains the original complex sections and vessel outlines.
EXPORT_VERIFICATION.json records the completed checks.
Earlier figure versions remain archived.

Run from the project root:
PYTHONPATH=src:. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.stage4.paper_acoustic_planes --output results/stage4/publication_figures/acoustic_field_three_planes_v2
The companion paper_acoustic_volume module supplies font style, paths and hashing only; no volume render is called.
