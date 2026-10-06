B2 orthogonal acoustic pressure sections

Fig1b_acoustic_field: PNG/PDF/SVG/TIFF at 180 x 90 mm, 600 dpi.
XY: Z=0 mm; YZ: X=0 mm; XZ: Y=0 mm.
Original magma (purple-to-yellow) map; common linear 0-39 kPa scale.
White vessel outlines and background; English Arial labels.
All three panels preserve the same physical scale and equal aspect.
The original 3D-volume figure remains archived in acoustic_field_volume_reference_v1.

Data are exact sections of the existing B2 complex pressure field.
No new simulation, no smoothing, no source-data changes.
CAPTION.txt explains the figure; PROVENANCE.json records data identities.
source_data.h5 contains the original complex sections and vessel outlines.
EXPORT_VERIFICATION.json records the completed checks.

Run from the project root:
PYTHONPATH=src:. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.stage4.paper_acoustic_planes --output results/stage4/publication_figures/acoustic_field_three_planes_v1
The companion paper_acoustic_volume module supplies font style, paths and hashing only; no volume render is called.
