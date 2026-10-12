"""Compare exactly the common accepted horizon, excluding rejected diagnostics."""

import json
import h5py
from mevo_mb.fast_sim.f10d_metrics import process
from scripts.fast_pzt_f10d.validate import arrays_equal, compare_results


def compare(a, b, member=0):
    checks = []
    with h5py.File(a) as old, h5py.File(b) as new:
        count = min(len(old["time_s"]), len(new["time_s"]))
        end = float(old["time_s"][count - 1])
        for group in ["trajectory", "metrics", "morphology"]:
            assert set(old[group]) == set(new[group])
            for key in old[group]:
                aa = old[group + "/" + key][:count, member]
                bb = new[group + "/" + key][:count, 0]
                if group == "trajectory":
                    exact = aa.dtype == bb.dtype and aa.shape == bb.shape and aa.tobytes() == bb.tobytes()
                    checks.append(
                        dict(key=group + "/" + key, status="PASS" if exact else "FAIL", bitwise=exact)
                    )
                else:
                    checks.append(arrays_equal(aa, bb, group + "/" + key))
        for key in ["time_s", "control_index", "particle_id", "radii_m"]:
            aa = old[key][:count] if key in ["time_s", "control_index"] else old[key][:]
            bb = new[key][:count] if key in ["time_s", "control_index"] else new[key][:]
            exact = aa.tobytes() == bb.tobytes()
            checks.append(dict(key=key, status="PASS" if exact else "FAIL", bitwise=exact))
        metadata = [json.loads(h.attrs["metadata_json"]) for h in [old, new]]
        controls = [
            [e for e in m["effective_control_events"] if e["time_s"] <= end + 1e-14] for m in metadata
        ]
        checks.append(
            dict(key="effective_control_events", status="PASS" if controls[0] == controls[1] else "FAIL")
        )
        processed = []
        edges = []
        for h, bi in [(old, member), (new, 0)]:
            ids = h["edge_events/batch_i_j_linked"][:]
            et = h["edge_events/time_s"][:]
            mask = (ids[:, 0] == bi) & (et <= end + 1e-14)
            ids, et = ids[mask].copy(), et[mask]
            ids[:, 0] = 0
            edges.append((ids, et))
            processed.append(
                process(
                    h["positions_m"][:count, bi : bi + 1],
                    float(h["radii_m"][0]),
                    h["time_s"][:count],
                    ids,
                    et,
                )
            )
        exact = all(x.tobytes() == y.tobytes() for x, y in zip(*edges))
        checks.append(dict(key="accepted_edge_events", status="PASS" if exact else "FAIL", bitwise=exact))
        topology = compare_results(*processed)
    return dict(
        status="PASS"
        if all(c["status"] == "PASS" for c in checks) and topology["status"] == "PASS"
        else "FAIL",
        common_accepted_end_s=end,
        common_frames=count,
        checks=checks,
        recomputed_common_horizon=topology,
        quality_scope="Saved accepted diag/wall/near/free_path compared; failed-interval cumulative quality is excluded",
    )
