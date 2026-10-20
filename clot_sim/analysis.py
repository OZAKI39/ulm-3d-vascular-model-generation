"""Reference summary and read-only energy/fragment diagnostics."""

import numpy as np
from .fragmentation import components
from .geometry import load_reference
from .config import INPUTS, RHO, DX, HORIZON
from .io import load_npz


def reference_summary():
    ref = load_reference()
    n = int(ref["clot_count"])
    field = load_npz(INPUTS / "material.npz")
    return dict(
        material_nodes=n,
        anchors=len(ref["X"]) - n,
        reference_bonds=len(ref["bonds"]),
        reference_volume_m3=float(ref["V"][:n].sum()),
        mass_kg=float(RHO * ref["V"][:n].sum()),
        nominal_spacing_m=DX,
        horizon_m=HORIZON,
        displaced_reference_nodes=int(
            np.any(ref["X"] != ref["regular_X"], axis=1).sum()
        ),
        material_thresholds=sorted(np.unique(field["bond_s0"]).tolist()),
    )


def energy_residual(work, state, initial, contact_jump):
    sc = state["scalars"]
    return work - (
        sc[0]
        - initial["scalars"][0]
        + state["K"]
        - initial["K"]
        + sc[1]
        + state["Uc"]
        - initial["Uc"]
        - contact_jump
    )
