# clot_clearance — isolated straight-pipe NOSB-PD demo

Open [OPEN_RESULTS.html](OPEN_RESULTS.html) for the offline interactive result, animation and figures. Read [REPORT_ZH.md](REPORT_ZH.md) for the actual results. All development and outputs are in this folder; the existing BraVa case, shared physical code, original PDF and SimVascular source remain untouched.

The workflow is `prescribed/imported streaming traction → fixed-base NOSB-PD clot → irreversible damage → connectivity analysis`. It is a research prototype, **not a validated ultrasound thrombolysis simulation**. The supplied paper concerns aspiration thrombectomy; see [paper mapping](references/PAPER_READING.md).

## Configure and compile

The verified local interpreter already contains NumPy, SciPy, Numba, VTK/PyVista, Matplotlib, imageio, imageio-ffmpeg, Plotly and pytest. No changes to this environment were made. The optional PDF extraction package lives only in `.tools/`.

```bash
cd /home/lzy/projects/clot_clearance
export PYTHONDONTWRITEBYTECODE=1
export NUMBA_CACHE_DIR="$PWD/build/numba_cache"
export MPLCONFIGDIR="$PWD/build/matplotlib"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export CLOT_PY=/home/lzy/projects/computation_examples/_shared/python_environment/bin/python
cmake -S . -B build -DPython3_EXECUTABLE="$CLOT_PY"
cmake --build build
ctest --test-dir build --output-on-failure --output-junit "$PWD/logs/ctest_results.xml"
```

This builds the Python/Numba mechanics kernels; it does **not** rebuild SimVascular or claim a new C++ solver. Versions are in `provenance/ENVIRONMENT.json`. Cached compiled code stays inside `build/`.

## Run the smallest example

```bash
"$CLOT_PY" -B -m pd_clot.runner --config configs/minimal.json --output runs/minimal_new
```

The 96-particle, one-macro-step example takes seconds after compilation. Existing output directories are refused. Use a new name for each run.

## Run the main straight-pipe example

```bash
"$CLOT_PY" -B -m pd_clot.runner --config configs/straight_pipe.json --output runs/straight_pipe_new
```

The delivered run is already at `results/straight_pipe_demo/`: 384 particles, 12,804 bonds, eight damage macro steps. It completed in about 86 seconds in the recorded local environment. Do not relabel it as a BraVa population result. `CONFIG.json`, `IDENTITY.json`, `SUMMARY.json`, `history.json`, `summary.csv` and `states.npz` provide provenance and raw data. The exact source snapshot matching the run is `provenance/run_source_v1/`.

The default JSON gives every scientific/solver parameter. `streaming.traction_scale_Pa` is a **synthetic input**, independent of the metadata bubble radius/carrier frequency. Both damage modes (`instantaneous_stretch`, `cyclic_accumulation`) are available. `damage.enabled=false` disables changes to integrity. See [MODEL_ASSUMPTIONS.md](MODEL_ASSUMPTIONS.md) for formulas, dimensions and safeguards.

## Imported-field workflow

```bash
"$CLOT_PY" -B scripts/prepare_inputs.py
"$CLOT_PY" -B -m pd_clot.runner --config configs/file_traction.json --output runs/imported_traction_new
```

The preparation script generates `inputs/straight_pipe_poiseuille.vtu` (analytical empty-pipe solution), `straight_pipe_wall.vtp`, and `synthetic_traction.csv`. It recreates generated input fixtures/configs **only inside this new project**. It never reads or writes an active BraVa run.

For real CFD, copy your file under `inputs/`, duplicate `configs/file_traction.json`, set `streaming.path`, and use the duplicate. Supported contracts:

- CSV direct traction: `x,y,z,tx,ty,tz` (m, Pa).
- CSV volume velocity/pressure: `x,y,z,ux,uy,uz,p` (m, m/s, Pa).
- Linear tetrahedral `.vtu`/`.vtk`: point arrays `Velocity` and `Pressure`; customize through `streaming.arrays` if needed.
- Point `.vtp` or tetrahedral files with direct point vector `traction_Pa`.

All imported fields are static and SI. Curved surfaces, higher-order elements, point clouds, and volumetric fields have different interpolation limits; read the model notes. Unsupported geometry and outside-domain queries fail rather than extrapolate. `verification/IMPORT_AND_DT_CHECK.json` records the executed import and timestep comparisons.

## View/export results

```bash
# Regenerate the delivered demo's figures and offline HTML from its saved states:
"$CLOT_PY" -B scripts/visualize.py
xdg-open /home/lzy/projects/clot_clearance/OPEN_RESULTS.html
```

Alternatively open the HTML directly in a browser; no web service, network assets or login are required. The rotatable 3-D view uses saved macro states, a common damage scale and true displacement scale. The bubble marker is enlarged and explicitly labeled. MP4/GIF hold each saved state; they do not resolve the mechanical cycle or MHz carrier.

ParaView: open `results/straight_pipe_demo/particles.pvd`, `bonds.pvd`, `traction.pvd`, and `inputs/straight_pipe_wall.vtp`. Particle arrays include displacement, velocity, damage, maximum bond damage, active degree, strain/stabilization energy density, J/validity and fragment IDs. Surface output separates normal/tangential traction; bond output retains active and broken bonds. `PVD timestep` is **represented cycle count**, not seconds. Coordinates remain in meters; the HTML/figures display millimeters.

Example direct VTK file: `results/straight_pipe_demo/vtk/particles_0008.vtp`. Publication figures are 300 dpi PNG and PDF in `visualization/figures/`.

## Verification and preservation

```bash
ctest --test-dir build --output-on-failure --output-junit "$PWD/logs/ctest_results.xml"
# Revalidate and assemble the already completed A–F evidence without recomputation:
"$CLOT_PY" -B scripts/verify_campaign.py --reuse-completed
"$CLOT_PY" -B scripts/final_qc.py
```

On a fresh copy without `verification/`, run `scripts/verify_campaign.py` without `--reuse-completed`, then `scripts/compare_runs.py`. Both refuse to overwrite existing run directories. The initial acceptance-report serialization error and its correction are retained in `logs/acceptance.log` and `logs/acceptance_report_fix.log`; computed controls were reused after exact configuration checks, not rerun or altered.

## File map

| Path | Role |
|---|---|
| `PLAN.md`, `references/PAPER_READING.md` | Inspected local source tree and paper-specific integration choices |
| `pd_clot/geometry.py` | Cell-centered cloud, fixed/exposed sets, surface quadrature, neighbor list |
| `pd_clot/mechanics.py` | Active shape tensor, Neo-Hookean correspondence, stabilization, dt estimate |
| `pd_clot/damage.py` | Two irreversible damage modes and graph fragments |
| `pd_clot/streaming.py` | Analytic/file providers, Poiseuille field, stress and force transfer |
| `pd_clot/runner.py`, `output.py`, `build.py` | Damped explicit integration, cycle jumping, VTK, compilation |
| `configs/` | Main/minimal/import configurations |
| `vendor/brava_wss_reference.py` | Byte-identical local copy of the existing BraVa P1 helper |
| `tests/`, `scripts/` | Mechanics/import/safety tests, A–F controls, input and figure generation, final QC |
| `inputs/`, `results/`, `verification/` | Generated straight-pipe fields, main result and distinct software controls |
| `visualization/`, `OPEN_RESULTS.html` | Offline interactive view, MP4/GIF, PNG/PDF |
| `provenance/`, `logs/`, `references/` | Source hashes, protected-original checks, versions, logs, local paper/prompt copies |

The complete new-file index is `provenance/DELIVERY_FILES.json`. No original source or case files were modified. No remote jobs were started, stopped or resumed.

## Smallest next step

Provide one real microbubble-resolved **velocity+pressure tetrahedral snapshot or clot-surface traction file**, in SI with pressure reference and normal convention documented. Reuse `FileStreaming`, first check field coverage, gradients and integrated load, then replace the synthetic provider in a copied JSON. A steady BraVa background field does not contain acoustic microstreaming; obtaining that input is a separate CFD task. The paper's two-way OpenFOAM coupling and Ogden model are not implemented by this demo.
