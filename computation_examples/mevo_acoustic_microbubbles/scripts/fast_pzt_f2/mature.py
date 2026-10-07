"""Mature-state matched PZT branches and limited same-state sensitivity controls."""

import argparse
import json
import h5py
import numpy as np
from mevo_mb.fast_sim.morphodynamics import topology, mature_index, morphology_series, contrast
from .common import BASE, PROTOCOL, write
from .experiment import run_case
from .fields import diagnose


def load(path):
    with h5py.File(path / "trajectory.h5") as h:
        d = {k: h[k][:] for k in ["time_s", "positions_m", "control_index"]}
        d["labels"] = h["metrics/labels"][:, 0]
        d["largest_rg_m"] = h["metrics/largest_rg_m"][:, 0]
        d["quality_samples"] = h["quality_samples/raw_F1_diagonal"][:, 0]
        d["edge_times"] = h["edge_events/time_s"][:]
        d["edge_ids"] = h["edge_events/batch_i_j_linked"][:]
    return d


def analyze_topology(path, d):
    result, tt, ll = topology(
        d["positions_m"][0, 0],
        d["time_s"],
        d["labels"],
        d["edge_times"],
        d["edge_ids"],
        d["time_s"][-1],
        PROTOCOL["persistence_s"],
    )
    write(path / "TOPOLOGY_PERSISTENCE.json", result)
    return result, tt, ll


def main(label, adoption):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    chosen = json.loads((BASE / adoption / "SELECTED.json").read_text())
    field = diagnose(out)
    base = field["baseline"]
    cfg = dict(
        count=32,
        seed=20261008,
        cloud_radius_m=30e-6,
        center_m=[0, 80e-6, 40e-6],
        basis=True,
        dt_s=chosen["dt_s"],
        sample_interval_s=0.0002,
        physics=chosen["physics"],
        events=[base],
        duration_s=PROTOCOL["maturity"]["max_physical_s"],
    )
    formation = run_case(out / "formation", cfg)
    data = load(out / "formation")
    top, ht, hl = analyze_topology(out / "formation", data)
    idx, ids = mature_index(data, ht, hl, PROTOCOL["maturity"], 1.78e-6)
    if idx is None or formation["numerical_quality"] != "PASS":
        write(
            out / "RESULT.json",
            dict(
                complete=True,
                status="MATURE_CLUSTER_NOT_FORMED_OR_UNQUALIFIED",
                formation=formation,
                maturity_index=idx,
            ),
        )
        return
    time = float(data["time_s"][idx])
    state = data["positions_m"][idx]
    maturity = dict(
        status="MATURE",
        frame_index=idx,
        absolute_time_s=time,
        particle_ids=ids,
        member_count=len(ids),
        criterion=PROTOCOL["maturity"],
        raw_topology_checked=True,
    )
    write(out / "MATURITY.json", maturity)
    np.savez_compressed(
        out / "MATURE_STATE.npz",
        positions_m=state,
        particle_id=np.arange(32),
        initial_cluster_ids=ids,
        time_s=time,
    )
    branches = {
        "unchanged": dict(event=base, physics=chosen["physics"]),
        "phase": dict(event=field["selected"]["phase"]["event"], physics=chosen["physics"])
        if field["selected"]["phase"]
        else None,
        "amplitude": dict(event=field["selected"]["amplitude"]["event"], physics=chosen["physics"])
        if field["selected"]["amplitude"]
        else None,
        "unchanged_repeat": dict(event=base, physics=chosen["physics"]),
        "F1_near_field": dict(
            event=base, physics=dict(chosen["physics"], force_model="f1_smooth_near_field")
        ),
        "mobility4": dict(event=base, physics=dict(chosen["physics"], mobility_iterations=4)),
    }
    rows = {}
    positions = {}
    topologies = {}
    for name, branch in branches.items():
        if branch is None:
            continue
        bc = dict(
            cfg,
            positions_m=state.tolist(),
            duration_s=PROTOCOL["continuation_s"],
            events=[dict(branch["event"], time_s=0)],
            physics=branch["physics"],
            absolute_time_offset_s=time,
            initial_condition="Exact saved mature state; no new random initialization",
        )
        rows[name] = run_case(out / name, bc)
        d = load(out / name)
        np.testing.assert_array_equal(d["positions_m"][0], state)
        t, _, _ = analyze_topology(out / name, d)
        topologies[name] = t
        series = morphology_series(
            d["positions_m"][:, 0],
            ids,
            d["labels"],
            PROTOCOL["observation"]["principal_axis_relative_eigengap"],
        )
        np.savez_compressed(
            out / name / "MORPHOLOGY.npz", time_s=d["time_s"], fixed_particle_ids=ids, **series
        )
        positions[name] = d["positions_m"][:, 0]
        write(out / "PARTIAL.json", dict(maturity=maturity, completed=list(rows), rows=rows))
    reference = positions["unchanged"]
    differences = {}
    for name, x in positions.items():
        dif = contrast(reference[-1], x[-1], ids)
        metrics = [contrast(a, b, ids) for a, b in zip(reference, x)]
        np.savez_compressed(
            out / name / "CONTRAST.npz",
            time_s=d["time_s"],
            **{k: np.asarray([m[k] for m in metrics]) for k in metrics[0]},
        )
        dif["quality"] = rows[name]["numerical_quality"]
        differences[name] = dif
    limits = PROTOCOL["observation"]
    noise = differences["unchanged_repeat"]
    model = [differences["F1_near_field"], differences["mobility4"]]
    thresholds = {
        "rg_relative": limits["rg_relative"],
        "eigenvalues_relative": limits["eigenvalues_relative"],
        "aligned_internal_rmsd_m": limits["aligned_rmsd_over_radius"] * 1.78e-6,
        "rigid_rotation_deg": limits["rigid_rotation_deg"],
    }
    effects = {}
    for name in ["phase", "amplitude"]:
        if name not in differences:
            continue
        dif = differences[name]
        detected = {
            k: bool(dif[k] > v and dif[k] > limits["numerical_envelope_multiplier"] * noise[k])
            for k, v in thresholds.items()
        }
        robust = {k: bool(detected[k] and dif[k] > max(m[k] for m in model)) for k in thresholds}
        shape = any(detected[k] for k in ["rg_relative", "eigenvalues_relative", "aligned_internal_rmsd_m"])
        effects[name] = dict(
            status="PASS"
            if shape
            and rows[name]["numerical_quality"] == "PASS"
            and rows["unchanged"]["numerical_quality"] == "PASS"
            else "NO_SIGNIFICANT_SHAPE_RESPONSE"
            if rows[name]["numerical_quality"] == "PASS"
            else "NUMERICAL_QUALITY_FAIL",
            final_contrast=dif,
            above_fixed_and_numerical_thresholds=detected,
            exceeds_tested_model_sensitivity=robust,
            persistent_topology=topologies[name]["counts"],
            criterion=limits,
        )
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            status="COMPLETE",
            formation=formation,
            maturity=maturity,
            rows=rows,
            comparisons=differences,
            effects=effects,
            field_selection=field,
            omitted_contact_branch="Candidate rejected by angular-balance audit; same-dt contact sensitivity retained in group B",
            topology=topologies,
            scope="EXPLORATORY; same mature state, physical translation/rotation/internal deformation reported separately",
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--adoption", required=True)
    a = ap.parse_args()
    main(a.label, a.adoption)
