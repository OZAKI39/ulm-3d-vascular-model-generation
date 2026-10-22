"""Finite known-state probes, historical ledger and explicit rejection cases."""

import argparse
import json
import h5py
import numpy as np
from mevo_mb.fast_sim.f10b_surrogate import admit, instantaneous, CentroidSurrogate
from scripts.fast_pzt_f7.fixed import probe
from .common import ROOT, BASE, F10A, case, js, write, sha, identity, objects


def maturity():
    path = ROOT / "results/fast_pzt/f5/science_v1/main/trajectory.h5"
    with h5py.File(path) as h:
        t = h["time_s"][:]
        sel = (t >= 0.29 - 1e-12) & (t <= 0.3 + 1e-12)
        rg = h["metrics/global_rg_m"][:, 0][sel]
        size = h["metrics/largest_size"][:, 0][sel]
        x = h["positions_m"][int(np.argmin(abs(t - 0.3))), 0]
    return dict(
        source=str(path.relative_to(ROOT)),
        sha256=sha(path),
        window_s=[0.29, 0.3],
        frames=int(sel.sum()),
        minimum_largest_size=int(size.min()),
        rg_relative_range=float(np.ptp(rg) / rg.mean()),
        mature=bool(size.min() == 32 and np.ptp(rg) / rg.mean() < 0.05),
    ), x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    proto = js(ROOT / "configs/fast_pzt/f10b/protocol.json")
    limits = proto["admission"]
    obj = objects()
    p, b, r, g, f = obj
    mat, oldx = maturity()
    cfg = case("O_1000")
    x = np.array(cfg["environment"]["positions_m"])
    np.testing.assert_allclose(x - x.mean(0), oldx - oldx.mean(0), rtol=0, atol=5e-19)
    write(out / "MATURITY.json", dict(mat, translation_preserves_template=True))
    rows = []
    scope = []
    ledgers = []
    gpu_checks = []
    for key in proto["surrogate"]["controls"]:
        cfg = case(key)
        ev = cfg["events"][0]
        adm = admit(
            x,
            ev,
            obj,
            limits,
            mobility="f5_local",
            duration_s=0.01,
            maturity=mat["mature"],
            validated_initial=True,
            validated_control=True,
        )
        assert adm["status"] == "ADMITTED_FOR_SHADOW_TEST", adm
        scope.append(dict(case="F10A_" + key, source="initial CONFIG", **adm))
        for mode in ["G0", "G1"]:
            sim = CentroidSurrogate(x, p, r, b, f, ev, mode, limits, adm)
            sim.sample()
            sim.read_center()
            gpu = (
                sim.flow_query.velocity.numpy()
                + sim.force.numpy() / (6 * np.pi * p.viscosity_Pa_s * p.radius_m)
            ).mean(0)
            cpu = instantaneous(x, ev, obj, mode)
            np.testing.assert_allclose(gpu, cpu, rtol=1e-10, atol=1e-14)
            gpu_checks.append(
                dict(control=key, mode=mode, max_velocity_difference_m_s=float(abs(cpu - gpu).max()))
            )
        path = F10A / "extended_5_v1" / key / "trajectory.h5"
        with h5py.File(path) as h:
            for frame in [0, 100, 200]:
                state = h["positions_m"][frame, 0]
                v = probe(state, ev, obj, "M0", True, True)
                mean = v["velocity_m_s"].mean(0)
                g1 = instantaneous(state, ev, obj, "G1")
                g0 = instantaneous(state, ev, obj, "G0")
                np.testing.assert_allclose(g1, mean, rtol=1e-10, atol=1e-14)
                rows.append(
                    dict(
                        control=key,
                        frame=frame,
                        time_s=float(h["time_s"][frame]),
                        source=str(path.relative_to(ROOT)),
                        source_sha256=sha(path),
                        full_M0_mean_velocity_m_s=mean,
                        G0_m_s=g0,
                        G1_m_s=g1,
                        G0_error_m_s=float(np.linalg.norm(g0 - mean)),
                        G1_error_m_s=float(np.linalg.norm(g1 - mean)),
                        secondary_net_N=v["secondary_N"].sum(0),
                        lubrication_mean_difference_m_s=(v["velocity_m_s"] - v["free_m_s"]).mean(0),
                        primary_N=v["primary_N"],
                        secondary_N=v["secondary_N"],
                    )
                )
            terms = json.loads(h["motion_ledger"].attrs["terms"])
            sums = h["motion_ledger/displacement_m"][:].sum(0)
            ledger = {t: sums[i].tolist() for i, t in enumerate(terms)}
            for term in ["secondary", "mobility", "pair_contact", "wall", "boundary"]:
                assert np.linalg.norm(ledger[term]) < 1e-12
            ledgers.append(
                dict(
                    control=key,
                    source=str(path.relative_to(ROOT)),
                    sums_m=ledger,
                    all_active=bool(h["active"][:].all()),
                    maximum_rg_relative_range=float(
                        np.ptp(h["metrics/global_rg_m"][:, 0]) / h["metrics/global_rg_m"][0, 0]
                    ),
                )
            )
    examples = [
        ("F9 near-wall", "results/fast_pzt/f9/y1_10_v1/run/trajectory.h5", 0, "f5_local"),
        ("F5 OFF dispersed", "results/fast_pzt/f5/science_v1/main/trajectory.h5", -1, "f5_local"),
        ("F7 RPY", "results/fast_pzt/f7/replay_v1/M2/trajectory.h5", 0, "rpy_free_space"),
    ]
    for name, path, frame, mobility in examples:
        path = ROOT / path
        with h5py.File(path) as h:
            xx = h["positions_m"][frame, 0]
            events = (
                json.loads(h["control"].attrs["events_json"])
                if "control" in h
                else json.loads(h.attrs["controls_json"])
            )
            ev = events[-1] if frame == -1 else events[0]
            active = bool(h["active"][frame].all()) if "active" in h else True
            rg = h["metrics/global_rg_m"][:, 0]
            change = float((rg[-1] - rg[0]) / rg[0])
        decision = admit(
            xx,
            ev,
            obj,
            limits,
            mobility=mobility,
            duration_s=0.01,
            maturity=frame == 0,
            validated_initial=False,
            validated_control=False,
            all_active=active,
            morphology_requested=frame == -1,
        )
        assert decision["status"] == "SURROGATE_NOT_ADMITTED"
        scope.append(
            dict(
                case=name,
                source=str(path.relative_to(ROOT)),
                sha256=sha(path),
                frame=frame,
                observed_full_reference_rg_change=change,
                **decision,
            )
        )
    cfg = case("O_3000", "pilot_v3")
    decision = admit(
        x,
        cfg["events"][0],
        obj,
        limits,
        mobility="f5_local",
        duration_s=0.01,
        maturity=True,
        validated_initial=True,
        validated_control=False,
        known_contact_rejection=True,
    )
    assert decision["status"] == "SURROGATE_NOT_ADMITTED"
    scope.append(dict(case="F10A 3 kPa", source="pilot_v3/O_3000/CONFIG.json", **decision))
    # Defensive interface tests: valid geometry does not override unsupported use.
    ev = case("O_1000")["events"][0]
    rejections = []
    for changed in [
        dict(duration_s=0.02),
        dict(all_active=False),
        dict(morphology_requested=True),
        dict(mobility="rpy_free_space"),
        dict(maturity=False),
    ]:
        args = dict(
            mobility="f5_local",
            duration_s=0.01,
            maturity=True,
            validated_initial=True,
            validated_control=True,
        )
        args.update(changed)
        dec = admit(x, ev, obj, limits, **args)
        assert dec["status"] == "SURROGATE_NOT_ADMITTED"
        rejections.append(dict(change=changed, reasons=dec["reasons"]))
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            one_step=rows,
            scope=scope,
            ledger=ledgers,
            gpu_cpu=gpu_checks,
            defensive_rejections=rejections,
            training_ready=False,
        ),
    )
    print("F10B one-step and admission checks PASS", flush=True)


if __name__ == "__main__":
    main()
