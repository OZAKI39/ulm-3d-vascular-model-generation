# CFD boundary preprocessing

The project has three separate, sequential responsibilities. Each stage consumes saved output
from the preceding stage; the CFD preprocessing stage does not regenerate an ROI or a surface.

## Step 1 — SWC preprocessing and representative ROI data

```powershell
D:\anaconda3\envs\pmp\python.exe swc_roi_generate.py
```

This produces the normalized analysis SWC, representative ROI records, exact `CUT_PORT`
coordinates, and the `CUT_PORT`-to-global-edge mapping.

## Step 2 — calibrated watertight ROI surface

```powershell
D:\anaconda3\envs\pmp\python.exe swc_stl_model_generate.py
```

This reuses a saved ROI and produces the calibrated Ultraliser STL/VTP surface. The current
surface calibration uses `radius_scale=0.91`; that factor is a reconstruction compensation and
is not applied again to the source radii used by the 1D hydraulic model.

## Step 3 — global 1D solution and ROI boundary transfer

```powershell
D:\anaconda3\envs\pmp\python.exe cfd_preprocess.py
```

An alternative strict YAML can be supplied as the only positional argument:

```powershell
D:\anaconda3\envs\pmp\python.exe cfd_preprocess.py configs\my_cfd_preprocess.yaml
```

The command loads the complete source-edge analysis SWC, validates its stable edge identities
against the sampling manifest, solves a Newtonian sparse 1D resistor network, and transfers the
resulting pressure and flow to every CFD boundary. SWC parent→current is a simulation direction
only; it is neither measured nor asserted to be physiological ground truth.

The formal baseline treats the single structural root as `ASSUMED_GLOBAL_INLET`, uses the
configured literature-derived root velocity or flow, and assigns all structural leaves a 0 Pa
gauge reference. A CFD-ready ROI receives one volumetric-flow inlet with parabolic-profile
metadata and direct 1D pressure at each assumed outlet. Resistance and Windkessel outlet models
are not part of this baseline.

ROI CFD boundaries have two traceable origins:

1. `CUT_PORT` boundaries are classified as `ASSUMED_INLET` or `ASSUMED_OUTLET` from their saved
   local SWC parent→current orientation.
2. `TRUE_TERMINAL` boundaries are treated as `ASSUMED_OUTLET` under the current formal baseline.
   Their pressure and expected flow come from the matching global structural leaf and its unique
   incoming edge.

The TRUE_TERMINAL rule is only a simulation convention that keeps the baseline boundary model
simple. It does not identify an experimentally verified physiological outlet. The structural-leaf
pressure is a global gauge reference, not a statement that real blood pressure is zero.

Readiness mass conservation compares all assumed inlet flow with both CUT_PORT and TRUE_TERMINAL
outlet flow. A failed check reports `CFD_ROI_NOT_READY` and does not create a solver boundary
package. The next-stage instruction then requests failure review; surface/volume mesh preparation
is recommended only after `CFD_PREPROCESS_BASELINE_PASS`.

```text
swc_roi_generate.py
        ↓
sampling ROI + global mapping
        ↓
swc_stl_model_generate.py
        ↓
Ultraliser STL/VTP
        ↓
cfd_preprocess.py
        ↓
global 1D flow + ROI boundary conditions
        ↓
PASS only: NEXT: CFD surface / volume mesh preparation
        ↓
3D CFD
```

This stage does not modify STL files, create a volume mesh, solve Navier–Stokes equations, or run
microbubble tracking.
