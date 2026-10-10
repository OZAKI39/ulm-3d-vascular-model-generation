"""Fixed physical domain refinement with common-coordinate field comparison."""
import copy
from pathlib import Path

import numpy as np
import yaml

from mevo_mb.acoustics.basis import AcousticBasisDatabase, ComplexAcousticField
from mevo_mb.acoustics.kwave_backend import KWaveCW
from mevo_mb.bubbles.forces import primary_bjerknes
from mevo_mb.bubbles.response import LinearResonantBubbleModel
from mevo_mb.control.pzt import PZTArray
from mevo_mb.coupling.geometry import Cylinder
from mevo_mb.io.artifacts import write_json


def run_refined(config: dict, run: Path) -> dict:
    c = copy.deepcopy(config)
    factor = 1.25
    c["grid"]["shape"] = [int(n*factor) for n in c["grid"]["shape"]]
    c["grid"]["spacing_m"] = [d/factor for d in c["grid"]["spacing_m"]]
    c["acoustics"]["pml_size"] = int(c["acoustics"]["pml_size"]*factor)
    (run/"config.yaml").write_text(yaml.safe_dump(c, sort_keys=False))
    adapter = KWaveCW(c, run/"acoustic")
    control = PZTArray.from_config(c).get_control_vector()
    p, g, metadata = adapter.run(control, "refined_direct")
    np.savez_compressed(run/"acoustic/refined_field.npz", pressure_Pa=p, gradient_Pa_m=g,
                        x_m=adapter.roi_axes[0], y_m=adapter.roi_axes[1], z_m=adapter.roi_axes[2], control_Pa=control)
    write_json(run/"metadata.json", {"status": "GRID_FIELD_COMPLETE", "solver_metadata": metadata})
    return {"refined_field": str(run/"acoustic/refined_field.npz"), "refined_config": c, "solver_metadata": metadata}


def load_refined(path: Path, frequency_Hz: float, scale: float = 1.) -> ComplexAcousticField:
    with np.load(path) as f:
        return ComplexAcousticField(tuple(f[f"{a}_m"].copy() for a in "xyz"), f["pressure_Pa"]*scale, f["gradient_Pa_m"]*scale, frequency_Hz)


def compare_grid(config: dict, basis_path: Path, refined_path: Path, run: Path) -> dict:
    db = AcousticBasisDatabase.load_basis(basis_path)
    field = db.combine(PZTArray.from_config(config).get_control_vector())
    fine = load_refined(refined_path, field.frequency_Hz)
    mask = Cylinder(**config["vessel"]).mask(db.axes_m)
    xyz = np.stack(np.meshgrid(*db.axes_m, indexing="ij"), axis=-1)[mask]
    p0 = field.pressure_Pa[mask]
    p1 = fine.interpolate_complex_pressure(xyz)
    g0 = field.gradient_Pa_m[mask]
    g1 = fine.interpolate_complex_gradient(xyz)
    response = config["bubble_response"]
    model = LinearResonantBubbleModel(response["resonance_frequency_Hz"], response["damping_ratio"])
    r = config["particles"]["size_distribution"]["radius_m"]
    f0 = primary_bjerknes(model.response(r, field.frequency_Hz, p0, config["medium"]), g0)
    f1 = primary_bjerknes(model.response(r, field.frequency_Hz, p1, config["medium"]), g1)
    def rel(a, b):
        return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), np.finfo(float).tiny))
    result = {"status": "TWO_GRID_DIAGNOSTIC_NOT_PROOF_OF_GRID_INDEPENDENCE", "frequency_Hz": field.frequency_Hz,
              "coarse_spacing_m": config["grid"]["spacing_m"], "fine_spacing_m": [float(np.diff(a)[0]) for a in fine.axes_m],
              "common_sample_count": len(xyz), "complex_pressure_relative_L2": rel(p0, p1),
              "complex_gradient_relative_L2": rel(g0, g1), "primary_force_relative_L2": rel(f0, f1),
              "coarse_peak_pressure_Pa": float(np.abs(p0).max()), "fine_at_common_points_peak_pressure_Pa": float(np.abs(p1).max()),
              "coarse_center_pressure_amplitude_Pa": float(np.abs(field.interpolate_complex_pressure([[0, 0, 0]])[0])),
              "fine_center_pressure_amplitude_Pa": float(np.abs(fine.interpolate_complex_pressure([[0, 0, 0]])[0])),
              "comparison": "Identical physical coordinates; fine trilinear interpolation; no smoothing; PML physical thickness preserved"}
    write_json(run/"acoustic/GRID_CONVERGENCE.json", result)
    import matplotlib.pyplot as plt
    from mevo_mb.visualization.plots import save
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axs[0].bar(["Pressure", "Gradient", "Primary force"], 100*np.array([rel(p0, p1), rel(g0, g1), rel(f0, f1)]))
    axs[0].set(ylabel="Common-point relative L2 difference (%)", title="dx vs dx/1.25")
    axs[1].bar(["Coarse peak", "Fine peak", "Coarse center", "Fine center"], [result[k]/1000 for k in ["coarse_peak_pressure_Pa", "fine_at_common_points_peak_pressure_Pa", "coarse_center_pressure_amplitude_Pa", "fine_center_pressure_amplitude_Pa"]])
    axs[1].set(ylabel="Pressure amplitude (kPa)")
    axs[1].tick_params(axis="x", rotation=20)
    save(fig, run/"figures/grid_convergence", config["output"]["figure_dpi"])
    write_json(run/"metadata.json", {"status": "DIAGNOSTIC_COMPLETE", "result": result})
    return result
