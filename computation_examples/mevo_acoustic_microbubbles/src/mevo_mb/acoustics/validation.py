"""Independent direct excitation versus complex linear superposition."""
from pathlib import Path

import numpy as np

from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.acoustics.kwave_backend import KWaveCW, acoustic_identity
from mevo_mb.coupling.geometry import Cylinder
from mevo_mb.io.artifacts import write_json


def validate_superposition(config: dict, basis_path: Path, run: Path) -> dict:
    db = AcousticBasisDatabase.load_basis(basis_path)
    if acoustic_identity(config) != db.metadata["acoustic_identity"]:
        raise ValueError("Basis identity does not match config")
    a = config["acoustics"]
    u = a["basis_source_amplitude_Pa"] * np.asarray(a["superposition_amplitudes"]) * np.exp(1j*np.asarray(a["superposition_phases_rad"]))
    combined = db.combine(u)
    direct, _, metadata = KWaveCW(config, run/"acoustic").run(u, "direct")
    mask = Cylinder(**config["vessel"]).mask(db.axes_m)
    reference = direct[mask]
    difference = combined.pressure_Pa[mask] - reference
    peak = np.max(np.abs(reference))
    significant = np.abs(reference) > peak*a["phase_mask_relative_amplitude"]
    phase = np.angle(combined.pressure_Pa[mask][significant]*np.conj(reference[significant]))
    l2 = float(np.linalg.norm(difference)/np.linalg.norm(reference))
    result = {"status": "PASS" if l2 <= a["validation_relative_L2_max"] else "FAIL", "relative_L2": l2,
              "max_error_over_peak": float(np.max(np.abs(difference))/peak),
              "max_pointwise_relative_error_masked": float(np.max(np.abs(difference[significant])/np.abs(reference[significant]))),
              "phase_RMS_rad_masked": float(np.sqrt(np.mean(phase**2))), "phase_max_abs_rad_masked": float(np.max(np.abs(phase))),
              "phase_mask_relative_amplitude": a["phase_mask_relative_amplitude"], "masked_count": int(significant.sum()),
              "excluded_near_node_count": int((~significant).sum()), "relative_L2_limit": a["validation_relative_L2_max"],
              "basis_file": str(basis_path), "direct_metadata": metadata}
    write_json(run/"acoustic/SUPERPOSITION.json", result)
    np.savez_compressed(run/"acoustic/comparison.npz", basis=combined.pressure_Pa, direct=direct, mask=mask)
    from mevo_mb.visualization.plots import plot_superposition
    plot_superposition(run, db.axes_m, combined.pressure_Pa, direct, result, config)
    write_json(run/"metadata.json", {"status": result["status"], "task": "superposition", "relative_L2": l2})
    if result["status"] != "PASS":
        raise RuntimeError(f"Superposition validation failed: {l2}")
    return result
