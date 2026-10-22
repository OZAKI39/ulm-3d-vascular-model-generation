B2 PZT control figure, version 4.
White background, regular Arial, 3D coordinate grid, muted red vessel, grey-blue translucent sources.
Outward source faces use reflective shading. Z+ and Y- have reduced transparency.
Use Fig1a_PZT_control.png or .tif as the 600 dpi master; _preview.png is a browser preview.
PDF and SVG embed the 3D rendering and retain editable text. CAPTION.txt documents display-only backing/rim dimensions and all opacity settings.
PROVENANCE.json records the exact source files, materials, camera, reference video hash and geometry/transparency checks. Original source data are unchanged.
Reproduce in the existing local project with its scientific Python environment:
PYTHONPATH=src:. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -m scripts.stage4.paper_pzt_reference --output results/stage4/publication_figures/PZT_control_reproduction
The reference MP4 and Arial font files are external local resources and are not redistributed in this figure directory.
