"""Independent finite-feature and candidate-coverage audit on the actual B2 events."""

from time import perf_counter
import numpy as np
from scripts.stage5a_r3.common import BASE, load, write, sha, identity
from scripts.stage5a_r3.run import restore
from mevo_mb.bubbles.motion_sweep_r3 import WallSweepR3
from mevo_mb.bubbles.geometry_exact_r3 import segment_triangle, record


def main():
    begin = perf_counter()
    tracker, old, _, _, _ = load("B2")
    sweep = WallSweepR3(tracker.cylinder, old.path)
    rows, normals = [], []
    for variant in ["nominal", "refined"]:
        folder = BASE / "diagnosis/reproduction_v1" / variant
        state, _, _ = restore(folder / "terminal.npz")
        active = state["active"]
        x, r, ids = state["positions_m"][active], state["radii_m"][active], state["particle_id"][active]
        for path in sorted(folder.glob("candidate_*.npz")):
            with np.load(path) as z:
                data = {k: z[k] for k in z.files}
            # Candidate files bind the actual full coupled QP trial, not a
            # single-particle synthetic trajectory.
            y = data["end_m"] if "end_m" in data else data["y"]
            proof = sweep(x, y, r, 1e-11)
            completeness = []
            exact_relevant = []
            for i in range(len(x)):
                reach = np.nextafter(r[i] + 1e-11, np.inf)
                low = np.nextafter(np.minimum(x[i], y[i]) - reach, -np.inf)
                high = np.nextafter(np.maximum(x[i], y[i]) + reach, np.inf)
                brute = np.flatnonzero(
                    np.all(sweep.triangle_high >= low, axis=1) & np.all(sweep.triangle_low <= high, axis=1)
                )
                chosen = sweep.candidates(x[i], y[i], reach)
                missing = sorted(set(brute) - set(chosen))
                completeness.append(
                    dict(
                        particle_id=int(ids[i]),
                        brute_count=len(brute),
                        locator_count=len(chosen),
                        missing=missing,
                    )
                )
                if ids[i] == 152:
                    for cell in brute:
                        value = record(segment_triangle(x[i], y[i], *sweep.triangles[cell]), r[i], 1e-11)
                        exact_relevant.append(dict(wall_triangle_id=int(cell), **value))
            rows.append(
                dict(
                    variant=variant,
                    candidate=str(path),
                    candidate_sha256=sha(path),
                    new_certificate=proof,
                    exact_responsible_features=exact_relevant,
                    broad_phase=completeness,
                    complete=all(not z["missing"] for z in completeness),
                )
            )
            for label, coords in [("start", x), ("end", y)]:
                cached = tracker.cylinder.wall_gap_and_normal(coords, r)
                tracker.cylinder._distance_cache = None
                cold = tracker.cylinder.wall_gap_and_normal(coords, r)
                i = int(np.flatnonzero(ids == 152)[0])
                candidates = sweep.candidates(coords[i], coords[i], r[i] + 1e-9)
                features = []
                exact_objects = []
                for cell in candidates:
                    obj = segment_triangle(coords[i], coords[i], *sweep.triangles[cell])
                    exact_objects.append((obj["distance2"], int(cell)))
                    rec = record(obj, r[i], 1e-11)
                    features.append(dict(wall_triangle_id=int(cell), **rec))
                normals.append(
                    dict(
                        variant=variant,
                        candidate=path.name,
                        endpoint=label,
                        gap_m=float(cold[0][i]),
                        outward_normal=cold[1][i].tolist(),
                        cold_gap_difference_max_m=float(np.max(abs(cold[0] - cached[0]))),
                        cold_normal_difference_max=float(np.max(abs(cold[1] - cached[1]))),
                        exact_nearest_triangle_id=min(exact_objects)[1],
                        features=features,
                    )
                )
    report = dict(
        rows=rows,
        normal_switch=normals,
        identity=identity(),
        wall_s=perf_counter() - begin,
        candidate_completeness="Full-mesh triangle AABB enumeration versus VTK locator for every particle and candidate; exact predicates on all responsibility features",
        signed_scope="Signed initial-centre membership from original VTK geometry; continuous positive separation from finite WALL preserves side unless leaving an open port. This +5 ms campaign has no port crossings; caps are excluded from WALL candidates.",
        unchanged_tolerance_m=1e-11,
        projection_m=0,
        classification_policy="SAFE is continuous proof; UNSAFE is exact distance below original tolerance; unresolved floating bounds fall back to exact rational arithmetic",
    )
    write(BASE / "final/WALL_CERTIFICATE_AUDIT.json", report)
    print(
        [
            (r["variant"], r["candidate"].split("/")[-1], r["new_certificate"]["pass_gate"], r["complete"])
            for r in rows
        ],
        flush=True,
    )


if __name__ == "__main__":
    main()
