# ULM 3-D vascular lumen model generation

This private repository contains the current `model_generate.py` implementation and the
minimal code/data closure needed to reproduce the active v8 ROI test. It reconstructs a
CFD-ready lumen surface from saved representative-ROI centerlines and SWC radii.

## Included scope

- `model_generate.py` and `cfd_lumen_config.yaml`;
- the `utils/cfd_lumen` reconstruction, QC, export, and visualization modules;
- only the sampling/SWC adapters required by that pipeline;
- `tests/test_cfd_lumen.py`;
- one test ROI: `raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274`;
- the single normalized source SWC component and edge manifest required for strict source-geometry verification.

The full `vessel_model` dataset, raw TIFF/Mask images, unrelated sampling candidates, and all
generated `outputs/model_generate` directories are intentionally excluded.

## Environment

Python 3.10 is the validated interpreter. From the repository root:

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Run unit and synthetic tests

```powershell
python -m pytest -q
```

## Reproduce the current ROI test

Run from the repository root so the portable paths in the test manifest resolve correctly:

```powershell
python model_generate.py `
  --sampling-run test_data/sampling `
  --rodent-run test_data/rodent_run `
  --roi-id raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274 `
  --surface-version v8 `
  --workers 1 `
  --headless `
  --run-id roi003274_v8
```

Results are created under `outputs/model_generate/roi003274_v8/` and remain ignored by Git.
The formal v8 reconstruction is computationally and memory intensive; the checked-in test data
is small, but the generated implicit fields and meshes can be much larger.

## Data boundary

The fixture inventory and SHA-256 hashes are documented in `test_data/README.md`. The normalized
SWC fixture is a derived single connected component (7,419 nodes / 7,418 edges); it is not the
original full mouse-brain dataset. TIFF intensity and segmentation-mask data are not used as
surface geometry and are not present here.
