"""Saved-state pair-force reconstruction and finite runtime guard checks.

No further scientific trajectory integration. Pair values are independently
reconstructed from recorded positions, original fields and response table.
"""

import argparse
import numpy as np
from mevo_mb.fast_sim.consistent_field import query_numpy
from mevo_mb.fast_sim.nearfield import coupled_fundamental, feedback_weight
from mevo_mb.fast_sim.f10b_surrogate import CentroidSurrogate, admit
from .common import ROOT, BASE, case, js, write, identity, objects


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    obj = objects()
    p, b, r, g, f = obj
    ev = case("O_3000", "pilot_v3")["events"][0]
    field, _, _ = b.combine(ev["amplitude_Pa"], ev["phase_rad"], ev["enabled"])
    chi = r.query(np.array([r.knots[0] + 0j]))[0, 0] / r.knots[0]
    pair_records = {}
    reports = []
    for folder, labels in [("contact_v1", ["A0", "A1"]), ("contact_v2", ["A2", "A3"])]:
        for label in labels:
            d = js(BASE / folder / (label + ".json"))
            data = []
            for t in d["traces"]:
                x = np.array(t["initial_positions_m"])
                pp, grad = query_numpy(field.pressure_Pa, b.origin, b.spacing, x)
                v = r.query(pp)
                primary = -0.5 * np.real(v[:, 0, None] * grad.conj())
                net = np.zeros_like(x)
                rows = []
                for i in range(32):
                    for j in range(i + 1, 32):
                        dx = x[j] - x[i]
                        distance = np.linalg.norm(dx)
                        corr = np.sum(np.arange(1, 9) ** 2 * np.real(v[i] * v[j].conj()))
                        ui, uj, _, _, _ = coupled_fundamental(
                            v[i, 0], v[j, 0], pp[i], pp[j], distance, p, chi
                        )
                        corr += feedback_weight(distance, p.radius_m) * (
                            np.real(ui * uj.conjugate()) - np.real(v[i, 0] * v[j, 0].conjugate())
                        )
                        q = min(distance / (p.transition_over_sum * 2 * p.radius_m), 1)
                        core = p.soft_core_over_sum * 2 * p.radius_m * (1 - q) ** 2
                        ff = (
                            p.density_kg_m3
                            * (2 * np.pi * p.frequency_Hz) ** 2
                            / (8 * np.pi)
                            * corr
                            * dx
                            / (distance**2 + core**2) ** 1.5
                        )
                        net[i] += ff
                        net[j] -= ff
                        rows.append([i, j, distance, distance - 2 * p.radius_m, *ff])
                np.testing.assert_allclose(primary, t["primary_force_N"], rtol=2e-10, atol=1e-23)
                np.testing.assert_allclose(net, t["secondary_force_N"], rtol=2e-10, atol=1e-23)
                data.append(rows)
                reports.append(
                    dict(
                        case=label,
                        step=t["step"],
                        primary_max_error_N=float(abs(primary - t["primary_force_N"]).max()),
                        secondary_max_error_N=float(abs(net - t["secondary_force_N"]).max()),
                        maximum_pair_force_N=float(np.linalg.norm(np.array(rows)[:, 4:], axis=1).max()),
                        secondary_rms_N=float(np.sqrt(np.mean(np.sum(net**2, axis=1)))),
                    )
                )
            pair_records[label] = np.array(data)
    np.savez_compressed(out / "reconstructed_pairs.npz", **pair_records)
    limits = js(ROOT / "configs/fast_pzt/f10b/protocol.json")["admission"]
    cfg = case("O_1000")
    x = np.array(cfg["environment"]["positions_m"])
    ev = cfg["events"][0]
    adm = admit(
        x,
        ev,
        obj,
        limits,
        mobility="f5_local",
        duration_s=0.01,
        maturity=True,
        validated_initial=True,
        validated_control=True,
    )
    guards = []
    # Zero-duration software probes: no accepted motion or fabricated trajectory.
    for mode, change in [("G1", dict(maximum_pressure_Pa=1.0)), ("G0", dict(maximum_shear_s_inv=0.0))]:
        sim = CentroidSurrogate(x, p, r, b, f, ev, mode, limits | change, adm)
        sim.step(0.0)
        try:
            sim.read_center()
        except RuntimeError as exc:
            guards.append(dict(mode=mode, probe=change, error=str(exc)))
        else:
            raise AssertionError("Runtime scope guard failed")
        np.testing.assert_array_equal(sim.center.numpy()[0], x.mean(0))
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            pair_reconstruction=reports,
            pair_columns=["i", "j", "distance_m", "gap_m", "Fx_N", "Fy_N", "Fz_N"],
            pair_source="Independent CPU reconstruction at stored diagnostic start states; validated against recorded production net forces",
            runtime_guards=guards,
            no_new_scientific_integration=True,
        ),
    )
    print("Saved-state near-pair reconstruction and runtime guards PASS", flush=True)


if __name__ == "__main__":
    main()
