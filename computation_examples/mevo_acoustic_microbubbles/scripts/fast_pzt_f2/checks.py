"""Small F2 formula checks and exactly two two-bubble motion comparisons."""

import argparse
from dataclasses import replace
import numpy as np
from mevo_mb.fast_sim.models import Physics, ResponseSection, random_cloud
from mevo_mb.fast_sim.nearfield import F2Physics, pair_reference, coupled_fundamental
from mevo_mb.fast_sim.contact_models import mobility_tensors, contact_pass
from mevo_mb.fast_sim.engine import FastSimulator, complex_vectors
from mevo_mb.fast_sim.f2_engine import F2Simulator
from mevo_mb.fast_sim import f2_kernels as fk
from .common import BASE, write, source_identity
from .experiment import run_case


def main(label):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    p = F2Physics()
    r = ResponseSection(BASE / "inputs/response.h5", p)
    checks = []
    distances = np.r_[np.linspace(2 * p.radius_m, 8 * p.radius_m, 81), 12 * p.radius_m]
    x = np.zeros((len(distances), 2, 3))
    x[:, 1, 0] = distances
    sim = F2Simulator(x, p, r)
    sim.control({"time_s": 0, "pressure_Pa": 1000})
    scans = {}
    for name, phase in [("in_phase", 0.0), ("counterphase", np.pi)]:
        pressure = np.tile([1000 + 0j, 1000 * np.exp(1j * phase)], len(distances))
        v = r.query(pressure)
        sim.fields()
        sim.harmonics.assign(complex_vectors(v))
        sim.pressure.assign(complex_vectors(pressure))
        sim.pairs()
        actual = sim.force.numpy().reshape(-1, 2, 3)
        reference = np.array(
            [
                pair_reference(xx, v[2 * b : 2 * b + 2], pressure[2 * b : 2 * b + 2], p, sim.chi_zero)
                for b, xx in enumerate(x)
            ]
        )
        np.testing.assert_allclose(actual, reference, rtol=2e-12, atol=1e-25)
        sign = 1 if phase == 0 else -1
        assert np.all(sign * actual[:, 0, 0] > 0)
        np.testing.assert_allclose(actual.sum(axis=1), 0, atol=1e-25)
        scans[name + "_f2_N"] = actual[:, 0, 0]
        baseline = replace(p, force_model="f1_smooth_near_field")
        scans[name + "_f1_N"] = np.array(
            [
                pair_reference(xx, v[2 * b : 2 * b + 2], pressure[2 * b : 2 * b + 2], baseline, sim.chi_zero)[
                    0, 0
                ]
                for b, xx in enumerate(x)
            ]
        )
        for distance in [2 * p.radius_m, 5 * p.radius_m, 8 * p.radius_m]:
            v0 = v[:2, 0]
            pp = pressure[:2]
            chi = v0 / pp
            k = p.density_kg_m3 * (2 * np.pi * p.frequency_Hz) ** 2 / (4 * np.pi * distance)
            direct = np.linalg.solve(np.array([[1, k * chi[0]], [k * chi[1], 1]]), v0)
            ans = coupled_fundamental(*v0, *pp, distance, p, sim.chi_zero)
            np.testing.assert_allclose(ans[:2], direct, rtol=1e-13, atol=1e-32)
        for boundary in [2 * p.radius_m, 8 * p.radius_m]:
            values = []
            for d in [boundary - 1e-7 * p.radius_m, boundary + 1e-7 * p.radius_m]:
                values.append(
                    pair_reference(np.array([[0.0, 0, 0], [d, 0, 0]]), v[:2], pressure[:2], p, sim.chi_zero)[
                        0, 0
                    ]
                )
            assert abs(values[1] - values[0]) / max(abs(values[0]), 1e-30) < 1e-5
    checks.append(dict(check="phase_sign_reciprocity_continuity_and_2x2_reference", status="PASS"))
    np.savez_compressed(out / "PAIR_FORCE_CURVES.npz", separation_m=distances, **scans)
    near = np.array([[0, 0, 0], [3.4e-6, 0, 0], [1.7e-6, 2.8e-6, 0], [1.7e-6, 1e-6, 2.5e-6]])
    contact = F2Simulator(near, p, r)
    contact.launch(
        fk.mobility,
        4,
        [contact.x, contact.contact_mobility, 4, p.radius_m, p.lubrication_activation, p.lubrication_floor],
    )
    m = mobility_tensors(near, p)
    np.testing.assert_allclose(contact.contact_mobility.numpy(), m, rtol=1e-12, atol=1e-14)
    contact.launch(
        fk.contact,
        4,
        [
            contact.x,
            contact.tmp1,
            4,
            p.radius_m,
            p.contact_relaxation,
            contact.diag,
            contact.error,
            contact.contact_mobility,
        ],
    )
    actual = contact.tmp1.numpy()
    expected = contact_pass(near, m, p)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-20)
    np.testing.assert_allclose(actual.mean(axis=0), near.mean(axis=0), rtol=0, atol=1e-20)
    assert np.sum((actual - actual.mean(axis=0)) ** 2) >= np.sum((near - near.mean(axis=0)) ** 2)
    checks.append(dict(check="anisotropic_contact_numpy_net_balance_and_expansion", status="PASS"))
    cloud = random_cloud(4, p.radius_m, 8e-6, 4321)
    f1p = Physics(contact_iterations=12)
    original = FastSimulator(cloud, f1p, r)
    delegated = F2Simulator(
        cloud, replace(p, force_model="f1_smooth_near_field", contact_model="f1_pbd_contact"), r
    )
    event = [{"time_s": 0, "pressure_Pa": 1000}]
    a = original.run(0.0001, 20e-6, event, sample_interval_s=20e-6)
    b = delegated.run(0.0001, 20e-6, event, sample_interval_s=20e-6)
    np.testing.assert_allclose(a["positions_m"], b["positions_m"], rtol=0, atol=1e-19)
    checks.append(dict(check="frozen_F1_execution_delegation", status="PASS"))
    zero = F2Simulator(cloud, p, r)
    z = zero.run(0.0001, 20e-6, [{"time_s": 0, "pressure_Pa": 0}], sample_interval_s=20e-6)
    np.testing.assert_array_equal(z["positions_m"], np.broadcast_to(cloud, z["positions_m"].shape))
    checks.append(dict(check="PZT_OFF_stationary", status="PASS"))
    batch = F2Simulator(np.array([cloud, cloud + np.array([1e-3, 0, 0])]), p, r)
    batch.control(event[0])
    batch.fields()
    batch.pairs()
    f = batch.force.numpy().reshape(2, 4, 3)
    np.testing.assert_allclose(f[0], f[1], rtol=1e-11, atol=1e-23)
    checks.append(dict(check="independent_batch_members", status="PASS"))
    toggle = F2Simulator(cloud, p, r)
    t = toggle.run(0.0001, 20e-6, [event[0], {"time_s": 17e-6, "pressure_Pa": 0}], sample_interval_s=20e-6)
    idx = np.flatnonzero(abs(t["time_s"] - 17e-6) < 1e-15)[0]
    np.testing.assert_array_equal(
        t["positions_m"][idx:], np.broadcast_to(t["positions_m"][idx], t["positions_m"][idx:].shape)
    )
    checks.append(dict(check="noncommensurate_control_event_alignment", status="PASS"))
    write(out / "GPU_CHECKS.json", dict(status="PASS", checks=checks, source=source_identity()))
    cfg = dict(
        count=2,
        seed=0,
        positions_m=[[-4e-6, 0, 0], [4e-6, 0, 0]],
        duration_s=0.025,
        dt_s=20e-6,
        sample_interval_s=0.0001,
        events=event,
    )
    rows = {}
    for model in ["f1_smooth_near_field", "f2_near_field_candidate"]:
        rows[model] = run_case(
            out / model, dict(cfg, physics=dict(force_model=model, contact_model="f1_pbd_contact"))
        )
    write(out / "RESULT.json", dict(complete=True, checks=checks, pair_runs=rows))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    args = ap.parse_args()
    main(args.label)
