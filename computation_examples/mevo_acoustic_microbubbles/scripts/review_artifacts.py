"""Postprocess closed raw data: pair CDF, geometry, and all-frame integrity checks."""
from pathlib import Path
import json

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
from scipy.spatial.distance import pdist

from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.coupling.geometry import Cylinder
from mevo_mb.io.artifacts import write_json
from mevo_mb.utils.config import load_config
from mevo_mb.visualization.plots import save


def integrity(run):
    config = load_config(run/"config.yaml")
    summary = json.loads((run/"particles/SUMMARY.json").read_text())
    cylinder = Cylinder(**config["vessel"])
    minimum_gap = np.inf
    with h5py.File(run/"particles/trajectories.h5") as f:
        assert f.attrs["complete"]
        times, radii = f["time_s"][:], f["radii_m"][:]
        assert times[0] == 0 and times[-1] == config["particles"]["duration_s"]
        assert np.all(np.diff(times) > 0)
        assert len(times) == summary["accepted_substeps"]+1
        tri = np.triu_indices(len(radii), 1)
        pair_radii = (radii[:, None]+radii[None, :])[tri]
        for key, ds in f.items():
            if key in {"radii_m", "particle_id", "pair_distance_histogram_edges_m"}:
                continue
            assert len(ds) == len(times), key
            for offset in range(0, len(times), 256):
                block = ds[offset:offset+256]
                assert np.isfinite(block).all(), key
                if key == "positions_m":
                    for x in block:
                        assert cylinder.contains(x, radii-1e-12).all()
                        minimum_gap = min(minimum_gap, float((pdist(x)-pair_radii).min()))
        assert minimum_gap >= -config["particles"]["contact_tolerance_m"]*1.001
        vtps = sorted((run/"vtk").glob("particles_*.vtp"))
        assert np.array_equal(pv.read(vtps[-1]).points, f["positions_m"][-1])
        fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
        for label, index in [("Initial", 0), ("Final", -1)]:
            d = np.sort(pdist(f["positions_m"][index]))*1e6
            ax.step(d, np.arange(1, len(d)+1)/len(d), where="post", label=label)
        ax.set(xlabel="Pair center distance (µm)", ylabel="Empirical cumulative fraction", title="All distinct pairs; no smoothing")
        ax.legend()
        ax.grid(alpha=0.25)
        save(fig, run/"figures/pair_distance_distribution")
    return {"status": "PASS", "states": len(times), "minimum_surface_gap_m": minimum_gap,
            "checked": ["all dynamic arrays finite and equal length", "strict time order and exact final time", "every state in finite cylinder", "all pair gaps within configured tolerance", "last VTP points equal HDF5"]}


def geometry(run):
    c = load_config(run/"config.yaml")
    field = ComplexAcousticField.load(run/"acoustic/combined_field.h5")
    fig = plt.figure(figsize=(10, 6), layout="constrained")
    ax = fig.add_subplot(projection="3d")
    blocks = pv.MultiBlock()
    theta = np.linspace(0, 2*np.pi, 65)
    for e in c["pzt"]:
        n = np.array(e["normal"], float)
        n /= np.linalg.norm(n)
        a = np.cross(n, [1, 0, 0] if abs(n[0]) < .9 else [0, 1, 0])
        a /= np.linalg.norm(a)
        b = np.cross(n, a)
        center = np.array(e["position_m"])
        points = center + e["size_m"][0]/2*(np.cos(theta)[:, None]*a + np.sin(theta)[:, None]*b)
        ax.plot(*points.T*1e3, lw=2, color="tab:orange")
        ax.text(*center*1e3, e["id"], fontsize=10)
        blocks[e["id"]] = pv.Disc(center=center, inner=0, outer=e["size_m"][0]/2, normal=n, c_res=64)
    radius, length = c["vessel"]["radius_m"], c["vessel"]["length_m"]
    xx, tt = np.meshgrid(np.linspace(-length/2, length/2, 12), theta)
    ax.plot_wireframe(xx*1e3, radius*np.cos(tt)*1e3, radius*np.sin(tt)*1e3, color="gray", alpha=.16, rstride=8, cstride=4)
    blocks["vessel_lumen"] = pv.Cylinder(center=(0, 0, 0), direction=(1, 0, 0), radius=radius, height=length)
    blocks.save(run/"vtk/source_and_vessel_geometry.vtm")
    iy, iz = np.meshgrid(np.arange(0, len(field.axes_m[1]), 2), np.arange(0, len(field.axes_m[2]), 2), indexing="ij")
    ix = np.argmin(abs(field.axes_m[0]))
    image = ax.scatter(np.zeros(iy.size), field.axes_m[1][iy].ravel()*1e3, field.axes_m[2][iz].ravel()*1e3,
                       c=np.abs(field.pressure_Pa[ix][iy, iz]).ravel()/1e3, s=7, cmap="viridis")
    fig.colorbar(image, ax=ax, shrink=.6, label="Local pressure amplitude at x = 0 (kPa)")
    ax.set(xlabel="x (mm)", ylabel="y (mm)", zlabel="z (mm)", title="Four acoustic discs and cylindrical flow region")
    ax.set_box_aspect((2, 1, 1))
    save(fig, run/"figures/geometry_and_field")


def main():
    root = Path("results")
    demo = json.loads((root/"DEMO_CAMPAIGN.json").read_text())
    convergence = json.loads((root/"CONVERGENCE_STATE.json").read_text())
    paths = [Path(c["run"]) for c in demo["scenarios"]]
    for c in convergence["results"].values():
        paths.extend(Path(c[k]["fine_run"]) for k in ["timestep", "grid_particle"])
    paths.append(Path(json.loads((root/"SECONDARY_OFF_STATE.json").read_text())["run"]))
    effective = root/"EFFECTIVE_STEP_STATE.json"
    if effective.exists():
        state = json.loads(effective.read_text())
        assert state["status"] == "COMPLETE"
        paths.extend(Path(x["run"]) for x in state["results"].values())
    checks = {}
    for p in paths:
        print("AUDIT", p, flush=True)
        checks[str(p)] = integrity(p)
    geometry(Path(demo["scenarios"][0]["run"]))
    write_json(Path("docs/RESULT_INTEGRITY.json"), checks)


if __name__ == "__main__":
    main()
