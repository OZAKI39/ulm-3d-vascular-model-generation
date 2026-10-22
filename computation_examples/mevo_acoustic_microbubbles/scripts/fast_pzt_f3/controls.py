"""Two preselected controls and two bounded initial-state sensitivities, all frozen physics."""

import argparse
import numpy as np
from mevo_mb.fast_sim.morphodynamics import topology, morphology_series
from mevo_mb.fast_sim.kinematics import rms_vectors
from scripts.fast_pzt_f2.analyze import read, adjacency
from .common import BASE, PROTOCOL, write, js, identity, state
from .experiment import run_case
from .repeatability import pair_comparison


def postprocess(path, reference):
    d = read(path)
    m = morphology_series(
        d["positions_m"][:, 0],
        np.arange(32),
        d["labels"][:, 0],
        PROTOCOL["observation"]["principal_axis_relative_eigengap"],
    )
    np.savez_compressed(path / "MORPHOLOGY.npz", time_s=d["time_s"], **m)
    np.savez_compressed(path / "ADJACENCY.npz", **adjacency(d))
    top, _, _ = topology(
        d["positions_m"][0, 0],
        d["time_s"],
        d["labels"][:, 0],
        d["edge_times"],
        d["edge_ids"],
        d["time_s"][-1],
        PROTOCOL["observation"]["topology_persistence_s"],
    )
    write(path / "TOPOLOGY_PERSISTENCE.json", top)
    c = pair_comparison(reference, d)
    np.savez_compressed(path / "CONTRAST.npz", time_s=d["time_s"], **c)
    return dict(
        endpoint={k: np.asarray(v[-1]).tolist() for k, v in c.items()},
        topology=top,
        final_aspect_ratio=float(m["fixed_aspect_ratio"][-1]),
        membership_changed=bool(np.any(m["dynamic_member_mask"] != m["dynamic_member_mask"][0])),
    )


def main(label, repeats, driving):
    r = js(BASE / repeats / "RESULT.json")
    assert r["status"] == "PASS_SAME_DEVICE_BITWISE"
    selection = js(BASE / driving / "SELECTED.json")
    mode = r["selected_execution_mode"]
    reference_path = BASE / repeats / (mode + "_0")
    reference = read(reference_path)
    cfg = js(reference_path / "CONFIG.json")
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    rows = {}
    analysis = {}
    # Only offline sidecars are added to the unchanged reference; its raw evidence is not rewritten.
    analysis["unchanged"] = postprocess(reference_path, reference)
    for name in ["phase", "amplitude"]:
        bc = dict(cfg, events=[selection["selected"][name]])
        rows[name] = run_case(out / name, bc)
        np.testing.assert_array_equal(read(out / name)["positions_m"][0], reference["positions_m"][0])
        analysis[name] = postprocess(out / name, reference)
        write(out / "PARTIAL.json", dict(completed=list(rows), rows=rows, analysis=analysis))
    x = state()
    rng = np.random.default_rng(PROTOCOL["initial_sensitivity"]["seed"])
    direction = rng.normal(size=(32, 3))
    direction -= direction.mean(0)
    direction /= rms_vectors(direction)
    perturbations = {}
    for i, size in enumerate(PROTOCOL["initial_sensitivity"]["RMS_displacements_m"]):
        name = f"initial_perturb_{i}"
        position = x + size * direction[None]
        dist = np.linalg.norm(position[0, :, None] - position[0, None, :], axis=-1)
        np.fill_diagonal(dist, np.inf)
        overlap = max(0.0, 2 * 1.78e-6 - dist.min())
        perturbations[name] = dict(
            RMS_m=size,
            max_initial_overlap_m=float(overlap),
            seed=PROTOCOL["initial_sensitivity"]["seed"],
            mean_displacement_m=(size * direction).mean(0),
        )
        if overlap > 0.01 * 1.78e-6:
            perturbations[name]["status"] = "INVALID_INITIAL_OVERLAP_NOT_PROJECTED"
            continue
        bc = dict(
            cfg,
            positions_m=position.tolist(),
            initial_condition=f"Predeclared zero-mean {size} m RMS perturbation of actual mature state",
        )
        rows[name] = run_case(out / name, bc)
        analysis[name] = postprocess(out / name, reference)
        perturbations[name]["status"] = rows[name]["numerical_quality"]
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            identity=identity(),
            rows=rows,
            analysis=analysis,
            perturbations=perturbations,
            reference=str(reference_path.relative_to(BASE)),
            selection_source=driving,
            selected_controls=selection,
            execution_mode=mode,
            physics_scope="EXPLORATORY",
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--repeats", required=True)
    ap.add_argument("--driving", required=True)
    a = ap.parse_args()
    main(a.label, a.repeats, a.driving)
