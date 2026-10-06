B2 acoustic field: three-dimensional volume revision

Main figure: Fig1b_acoustic_field.png / .pdf / .svg / .tif
180 x 104 mm, 600 dpi; Arial regular, white background.

Left: native 3D amplitude volume in the original focal ROI, plus 28/35 kPa isosurfaces.
Right: Z=0 section with the same linear 0-39 kPa scale and an explicit ROI locator.
Volume opacity is a documented display choice; use the right plane for direct color-value lookup.
The reference style is PZT_control_reference_v6. No solver or trajectory run was started.

CAPTION.txt: paper caption and rendering interpretation.
PROVENANCE.json: source identities, native-grid and rendering parameters.
source_data.h5: exact native complex data used by the figure.
EXPORT_VERIFICATION.json: source/export/occlusion checks.

Reproduce from the project root:
PYTHONPATH=src:. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.stage4.paper_acoustic_volume --output results/stage4/publication_figures/acoustic_field_volume_reference_v1

The original workflow_B2_B3_20261005_v1 figure remains archived without edits.
