MEVO: data-grounded PZT-to-microbubble figures
Generated 2026-10-05 from sealed Stage 4I B2/B3 simulation results.

Open OPEN_FIGURES.html for the complete gallery and download links.

CONTENTS
Fig1_B2_workflow: combined B2 workflow, 180 x 218 mm.
Fig1a_PZT_control: actual four-source geometry and controls, 180 x 70 mm.
Fig1b_acoustic_field: 3D pressure slices and z=0 section, 180 x 70 mm.
Fig1c_microbubble_forces: primary vectors and an actual secondary-force pair, 180 x 70 mm.
Fig1d_particle_tracking: actual 3D trajectories and branch counts, 180 x 70 mm.
FigS1_B2_B3_comparison: matched-seed pressure and trajectory comparison, 180 x 200 mm.
Each figure: PDF, editable-text SVG, 600 dpi PNG, 600 dpi RGB LZW TIFF.
Files ending _preview.png are 180 dpi browser previews, not submission masters.

STYLE
White background; English labels; Arial 8-11 pt with STIX mathematical glyphs.
PDFs embed the used TrueType fonts; SVG text is editable and requires those fonts.
Common pressure scale: 0-39 kPa. Common trajectory time scale: 0-6 s.
Geometric axes have explicit units and equal physical scaling in 3D.
These are general manuscript defaults. A specific journal may require different
page dimensions, typography, color mode or export formats. Consult the selected
journal before final submission. Source guidance consulted:
https://journals.plos.org/plosone/s/figures
https://www.nature.com/nature/for-authors/final-submission
The figures are not claimed to satisfy every journal's exact requirements.

SCIENTIFIC IDENTITY
The user selected current B2/B3; the old Stage-1 title in the reference diagram
is a workflow outline, not the identity of these source runs. Both cases use
2.5 MHz, four ideal 6 mm acoustic sources, 256 initial bubbles of radius 1.78 um,
and 6 s of observation. Source amplitudes are pressure, not electrical voltage.
B2: branch A 76; branch B 0; unassigned at observation cutoff 180.
B3: branch A 0; branch B 88; unassigned at observation cutoff 168.
Counts are not independent experimental replicates or all-success outcomes.
The secondary pair example uses all eight production response harmonics.
The free Stokes velocity equation precedes geometric wall/contact constraints.
Trajectory lines connect recorded states at 20 ms spacing; the nominal solver
step is 0.5 ms with adaptive constraints. They do not reconstruct substeps.
No new simulation was run. Read CAPTIONS.txt for complete model/display limits.

DATA AND AUDIT
FIGURE_PROVENANCE.json: exact run paths, protected input SHA256 hashes,
reconstruction checks, pair IDs and geometric values, renderer hash and layout.
EXPORT_VERIFICATION.json: embedded fonts, dimensions, DPI, source-data equality.
source_data.h5: actual trajectory positions/times/IDs/radii/active states;
initial force/velocity and 8-harmonic response arrays; complex control weights;
branch outcomes/crossing times; acoustic grid axes and three complex pressure
sections; cropped vessel wire and z=0 outline; PZT configuration metadata.
Coordinates ending _m are metres; vessel arrays ending _mm are millimetres.
Force/velocity arrays and response_harmonics are at t=0 only. See HDF5 attributes.
pzt_controls.csv and initial_forces.csv provide convenient tabular plot data.
Original basis/field/mesh/response files remain in the project and are not
repeated in this compact figure bundle; their hashes/paths are in provenance.
DELIVERY_SHA256.json hashes every other file in this directory.

REPRODUCE FIGURES
Use the existing project checkout and its scientific Python environment.
The copied paper_workflow_figures.py is the exact script used; the executable
module is scripts/stage4/paper_workflow_figures.py in the project.
From the project root:

PYTHONPATH=src:. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
 /home/lzy/projects/computation_examples/_shared/python_environment/bin/python \
 -m scripts.stage4.paper_workflow_figures --output results/stage4/publication_figures/reproduction

Add --preview for fast 180 dpi layout review. Use a new output directory.
The renderer checks sealed inputs, pressure-basis reconstruction, initial forces
and free velocity, then checks input hashes again after rendering. It requires
existing project dependencies; packaging validation also used Pillow and pypdf.
For identical typography, Arial is loaded from /mnt/c/Windows/Fonts/arial*.ttf;
otherwise the renderer explicitly falls back to DejaVu Sans. Fonts are not
redistributed in this bundle. Exact fonts are recorded in the export audit.
