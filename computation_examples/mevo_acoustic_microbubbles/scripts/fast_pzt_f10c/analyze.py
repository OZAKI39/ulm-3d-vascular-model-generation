"""Offline numerical admission; gates are read from the preregistered protocol."""

import h5py
import numpy as np
from mevo_mb.fast_sim.morphodynamics import contrast, rigid_fit
from .common import BASE, ROOT, F10A, js, write, sha


def read(path):
    with h5py.File(path / "trajectory.h5") as h:
        x = h["positions_m"][:, 0]
        result = dict(
            x=x,
            t=h["time_s"][:],
            center=x.mean(1),
            edge_ids=h["edge_events/batch_i_j_linked"][:],
            edge_times=h["edge_events/time_s"][:],
        )
        for key in [
            "global_rg_m",
            "global_shape_tensor_m2",
            "global_eigenvalues_m2",
            "largest_size",
            "minimum_gap_m",
        ]:
            result[key] = h["metrics/" + key][:, 0]
        e = result["global_eigenvalues_m2"]
        result["aspect"] = np.sqrt(e[:, 0] / e[:, -1])
    return result


def errors(a, b):
    complete = a["x"].shape == b["x"].shape
    count = min(len(a["x"]), len(b["x"]))
    last = float(b["t"][count - 1])
    a = {k: v[:count] if k not in ["edge_ids", "edge_times"] else v for k, v in a.items()}
    b = {k: v[:count] if k not in ["edge_ids", "edge_times"] else v for k, v in b.items()}
    rows = [contrast(x, y, np.arange(32)) for x, y in zip(a["x"], b["x"])]
    return dict(
        complete=complete,
        accepted_duration_s=last,
        max_centroid_distance_m=float(np.linalg.norm(a["center"] - b["center"], axis=-1).max()),
        max_centroid_axis_m=np.abs(a["center"] - b["center"]).max(0),
        max_particle_error_m=float(np.linalg.norm(a["x"] - b["x"], axis=-1).max()),
        max_aligned_rmsd_m=max(row["aligned_internal_rmsd_m"] for row in rows),
        max_rg_relative=float(np.max(np.abs(b["global_rg_m"] / a["global_rg_m"] - 1))),
        max_tensor_relative=max(row["tensor_relative"] for row in rows),
        max_aspect_relative=float(np.max(np.abs(b["aspect"] / a["aspect"] - 1))),
        largest_fraction_difference=float(np.max(np.abs(a["largest_size"] - b["largest_size"])) / 32),
        edge_events_exact=bool(
            np.array_equal(a["edge_ids"], b["edge_ids"]) and np.array_equal(a["edge_times"], b["edge_times"])
        ),
        edge_counts=[len(a["edge_times"]), len(b["edge_times"])],
        minimum_gap_m=float(b["minimum_gap_m"].min()),
        final_displacement_m=b["center"][-1] - b["center"][0],
    )


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    proto = js(ROOT / "configs/fast_pzt/f10c/protocol.json")
    entries = js(BASE / "multirate_v1/RESULTS.json")["rows"]
    paths = {(v["mode"], v["control"]): BASE / v["path"] for v in entries}
    controls = ["OFF", "S_1000", "O_1000"]
    baseline = {key: read(paths["Q0", key]) for key in controls}
    rows = []
    curves = {}
    hashes = {}
    for (mode, key), path in paths.items():
        value = read(path)
        base = baseline[key]
        err = errors(base, value)
        ref = read(F10A / "extended_5_v1" / key)
        fine = errors(ref, value)
        meta = js(path / "SUMMARY.json")
        q = meta["per_sample_quality"][0]["quality"]
        group = (
            err["complete"]
            and err["max_centroid_distance_m"] <= proto["multirate_group_gates"]["max_centroid_distance_m"]
        )
        morphology = err["complete"] and all(
            err[k] <= v if not isinstance(v, bool) else err[k] == v
            for k, v in proto["multirate_morphology_gates"].items()
        )
        physical = (
            meta["status"] == "COMPLETE"
            and q["error_bits"] == 0
            and q["max_accepted_penetration_m"] <= 0.01 * 1.78e-6
            and q["max_contact_pass_correction_m"] <= 0.1 * 1.78e-6
        )
        timing = js(path / "TIMING.json")
        row = dict(
            mode=mode,
            control=key,
            path=str(path.relative_to(BASE)),
            errors_vs_Q0=err,
            errors_vs_5us=fine,
            group_pass=bool(group and physical),
            morphology_pass=bool(morphology and physical),
            physical_quality_pass=physical,
            quality=q,
            timing=timing,
            cache_statistics=meta.get("cache_statistics"),
        )
        rows.append(row)
        hashes[str((path / "trajectory.h5").relative_to(ROOT))] = sha(path / "trajectory.h5")
        prefix = mode + "_" + key
        for name in [
            "t",
            "x",
            "center",
            "global_rg_m",
            "global_shape_tensor_m2",
            "aspect",
            "largest_size",
            "minimum_gap_m",
        ]:
            curves[prefix + "_" + name] = value[name]
        curves[prefix + "_aligned_rmsd_m"] = np.array(
            [rigid_fit(x, y)[1] for x, y in zip(base["x"], value["x"])]
        )
    modes = {}
    for mode in proto["modes"]:
        selected = [r for r in rows if r["mode"] == mode]
        paired = {}
        for other in ["OFF", "S_1000"]:
            dis = lambda m, k: curves[m + "_" + k + "_center"][-1, 2] - curves[m + "_" + k + "_center"][0, 2]
            if any(
                not row["errors_vs_Q0"]["complete"] for row in selected if row["control"] in ["O_1000", other]
            ):
                paired[other] = dict(
                    delta_z_m=None, error_vs_Q0_m=None, status="INCOMPLETE_NO_COMMON_10MS_ENDPOINT"
                )
                continue
            delta = dis(mode, "O_1000") - dis(mode, other)
            reference = dis("Q0", "O_1000") - dis("Q0", other)
            paired[other] = dict(delta_z_m=delta, error_vs_Q0_m=delta - reference)
        group = all(r["group_pass"] for r in selected) and all(
            v["error_vs_Q0_m"] is not None
            and abs(v["error_vs_Q0_m"]) <= proto["multirate_group_gates"]["max_paired_z_error_m"]
            for v in paired.values()
        )
        morphology = all(r["morphology_pass"] for r in selected)
        modes[mode] = dict(
            group_pass=group,
            morphology_pass=morphology,
            paired=paired,
            label="FULL_PARTICLE_BASELINE"
            if mode == "Q0"
            else "MULTIRATE_MORPHOLOGY_PASS"
            if group and morphology
            else "MULTIRATE_GROUP_PASS"
            if group
            else "MULTIRATE_NOT_ADMITTED",
            production_admitted=False,
        )
    batches = js(BASE / "batch_v1/RESULTS.json")["rows"] + js(BASE / "batch64_v1/RESULTS.json")["rows"]
    result = dict(
        batch=batches,
        multirate=rows,
        modes=modes,
        decomposition=js(BASE / "multirate_v1/DECOMPOSITION.json"),
        q4="NOT_RUN: no q2 ON-control mode passed full morphology"
        if not any(v["group_pass"] and v["morphology_pass"] for k, v in modes.items() if k != "Q0")
        else "CONDITIONAL",
        training_ready=False,
        physics_scope="EXPLORATORY",
    )
    exclusive = BASE / "exclusive_v1/COMPLETE.json"
    if exclusive.exists():
        completion = js(exclusive)
        audit = js(BASE / "exclusive_v1/RESOURCE_AUDIT.json")
        assert completion["status"] == audit["status"] == "PASS"
        result["exclusive_performance"] = dict(
            completion=completion,
            resource_audit=audit,
            batch=js(BASE / "exclusive_v1_batch/RESULTS.json")["rows"],
            cold=js(BASE / "exclusive_v1_cold/RESULTS.json")["rows"],
            scope="EXCLUSIVE_GPU_SAMPLED; one observation per size, not confidence intervals",
        )
    write(out / "RESULTS.json", result)
    write(out / "SOURCES.json", hashes)
    np.savez_compressed(out / "curves.npz", **curves)
    print(
        {
            k: dict(group=v["group_pass"], morphology=v["morphology_pass"], paired=v["paired"])
            for k, v in modes.items()
        }
    )
    for r in rows:
        print(
            r["mode"],
            r["control"],
            r["errors_vs_Q0"]["max_particle_error_m"],
            r["errors_vs_Q0"]["max_aligned_rmsd_m"],
            r["errors_vs_Q0"]["edge_counts"],
            r["errors_vs_Q0"]["max_centroid_distance_m"],
        )


if __name__ == "__main__":
    main()
