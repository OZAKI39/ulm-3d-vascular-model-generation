"""Comparison against independently captured original V4 outputs; no old imports."""

from pathlib import Path
import json
import gzip
import shutil
import numpy as np
import pytest
from clot_sim.config import INPUTS, K, G, MIN_H, BASE, RHO, V
from clot_sim.io import load_npz, sha, verify_inputs, dump
from clot_sim.geometry import load_reference, write_native_inputs
from clot_sim.material import SofteningLPS
from clot_sim.damage import predict
from clot_sim.surface import CurrentSurface
from clot_sim.contact import Contact, CONTACT
from clot_sim.fragmentation import Lineage
from clot_sim.analysis import components, energy_residual

FIXTURES = Path(__file__).parent / "fixtures"
ERRORS = {}
# Native scalars[7] is std::chrono evaluation time, not a physical state.
# All other native scalars (energy/history/counters/IDs) are compared.


def compare(label, actual, expected, atol, rtol=1e-12):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    assert actual.shape == expected.shape, label
    if actual.dtype.kind in "biuUS":
        assert np.array_equal(actual, expected), label
        ERRORS[label] = 0.0
    else:
        difference = float(np.max(np.abs(actual - expected), initial=0))
        ERRORS[label] = difference
        np.testing.assert_allclose(
            actual, expected, atol=atol, rtol=rtol, err_msg=label
        )


@pytest.fixture(scope="module")
def model():
    ref = load_reference()
    return (
        ref,
        SofteningLPS(ref, K=K, G=G),
        CurrentSurface(ref, "uniform"),
        Contact(ref),
    )


def test_input_identity():
    verify_inputs()
    capture_hashes = json.loads((FIXTURES / "SERVER_CAPTURE_SHA256.json").read_text())
    for name, want in capture_hashes.items():
        assert sha(FIXTURES / name) == want, name
    ref = load_reference()
    field = load_npz(INPUTS / "material.npz")
    assert (
        len(ref["X"]) == 3488
        and int(ref["clot_count"]) == 1927
        and len(ref["bonds"]) == 153986
    )
    assert (
        sha(INPUTS / "reference.npz")
        == "6a3599d1165be6a1d0bbc00b54a5f79917975b8fca317f6726d93c3402905f75"
    )
    assert (
        sha(INPUTS / "material.npz")
        == "00a890c039dc3c90b2b29d64f37b006135f70feebac2ba6b70b9e81c8d187212"
    )
    assert (
        sha(INPUTS / "boundary.npz")
        == "e974eacce84882932281f1bd40c2aa7b5fdac94531cc0c9310c00a5bd5cc31b0"
    )
    assert np.array_equal(ref["atom_ids"], field["atom_ids"])
    assert np.isclose(ref["V"][:1927].sum(), 2.40875e-10, rtol=0, atol=1e-25)
    oracle = json.loads((FIXTURES / "ORACLE.json").read_text())
    for item in oracle["states"]:
        assert sha(FIXTURES / (item["name"] + ".npz")) == item["sha256"]


@pytest.mark.parametrize("name", ["initial", "first_failure", "final"])
def test_saved_state_equivalence(model, name):
    ref, material, surface, contact = model
    s = load_npz(FIXTURES / (name + ".npz"))
    expected = load_npz(FIXTURES / (name + "_expected.npz"))
    m = material.evaluate(s["x"], s["d"])
    for key, value in m.items():
        compare(
            name + "/material/" + key,
            value,
            expected["material_" + key],
            1e-18 if key.startswith("U") or key == "dual_energy_error_J" else 1e-12,
        )
    z = surface.evaluate(s["x"], s["velocity"], s["d"], float(s["time"]))
    for key in [
        k.removeprefix("surface_") for k in expected if k.startswith("surface_")
    ]:
        # Across server Python 3.13 and WSL Python 3.12, LAPACK inverse/normal
        # roundoff was 1.9e-14 dimensionless; categorical masks remain exact.
        atol = (
            1e-12
            if key == "normal"
            else 1e-18 if key in ["wall_U", "wall_damping_power", "area"] else 1e-14
        )
        compare(name + "/surface/" + key, z[key], expected["surface_" + key], atol)
    uc, fc, pairs = contact.evaluate(s["x"], s["d"])
    compare(name + "/contact/U", uc, expected["contact_U"], 1e-18)
    compare(name + "/contact/force", fc, expected["contact_force"], 1e-12)
    compare(name + "/contact/pairs", pairs, expected["contact_pairs"], 0)
    top = components(ref, s["d"])
    for key in ["labels", "free", "attached"]:
        compare(name + "/topology/" + key, top[key], expected["topology_" + key], 0)
        compare(name + "/surface_topology/" + key, z["topology"][key], top[key], 0)
    field = load_npz(INPUTS / "material.npz")
    kk, dd = predict(
        material, ref, s["x"], s["kappa"], field["bond_s0"], field["bond_sc"]
    )
    compare(name + "/kappa", kk, expected["predicted_kappa"], 1e-12)
    compare(name + "/damage", dd, expected["predicted_d"], 1e-12)
    assert np.all(kk >= s["kappa"])
    lineage = Lineage(ref)
    lineage.update(top, s)
    records = lineage.records(s, top, np.zeros(len(ref["X"])))
    assert sum(r["nodes"] for r in records) == int(top["free"][:1927].sum())
    assert np.isclose(
        sum(r["reference_volume_m3"] for r in records),
        ref["V"][top["free"]].sum(),
        rtol=1e-12,
        atol=1e-30,
    )
    if name == "final":
        assert len(records) == 912 and max(r["nodes"] for r in records) == 461


def test_surface_virtual_work_and_wall_gradient(model):
    ref, material, surface, contact = model
    x = ref["X"].copy()
    d = ref["initial_cut"].astype(float)
    n = int(ref["clot_count"])
    v = np.random.default_rng(20261019).normal(0, 0.003, x.shape)
    v[n:] = 0
    z = surface.evaluate(x, v, d, 0.0024)
    compare("surface/resultant", z["force"].sum(0), z["face_force"].sum(0), 1e-16)
    compare(
        "surface/moment",
        np.cross(x, z["force"]).sum(0),
        np.cross(z["centroid"], z["face_force"]).sum(0),
        1e-18,
    )
    compare(
        "surface/power",
        np.sum(z["force"] * v),
        np.sum(z["face_force"] * z["face_velocity"]),
        1e-16,
    )
    passive = surface.evaluate(x, v, d, 0.0036)
    assert passive["drag_power"] <= 0 and not passive["drive_force"].any()
    q = surface.wall_geometry(x)
    probe = x.copy()
    probe[0] = q["nearest_wall_point_m"][0] + q["energy_normal"][0] * 10e-6
    q = surface.wall_geometry(probe)
    analytic = CONTACT["stiffness_N_m"] * q["overlap"][0] * q["energy_normal"][0]
    numeric = []

    def energy(xx):
        z = surface.wall_geometry(xx)
        return 0.5 * CONTACT["stiffness_N_m"] * np.dot(z["overlap"], z["overlap"])

    for j in range(3):
        plus = probe.copy()
        minus = probe.copy()
        plus[0, j] += 1e-9
        minus[0, j] -= 1e-9
        numeric.append(-(energy(plus) - energy(minus)) / (2e-9))
    compare("wall/potential_gradient", analytic, numeric, 1e-9)


def test_native_steps_and_restart(tmp_path):
    from clot_sim.solver import NativeSolver

    ns = NativeSolver(tmp_path / "native")
    s = ns.dynamics()
    ref = ns.ref
    material = SofteningLPS(ref, K=K, G=G)
    states = []
    field = load_npz(INPUTS / "material.npz")
    compare("native/installed_onset", ns.field(13), field["bond_s0"], 0, 0)
    compare("native/installed_complete", ns.field(14), field["bond_sc"], 0, 0)
    try:
        ns.tool.commit(ns.tool.evaluate(s["x"], s["velocity"], s["d"], s["time"]))
        for i in range(3):
            p = s
            xp = p["x"].copy()
            xp[: ns.n] += MIN_H * (
                p["velocity"][: ns.n]
                + MIN_H * p["force"][: ns.n] / (2 * ns.mass[: ns.n, None])
            )
            kp, dp = predict(material, ref, xp, p["kappa"], ns.field(13), ns.field(14))
            s = ns.step(MIN_H)
            compare(f"step{i+1}/predicted_x", s["x"], xp, 2e-18, 0)
            compare(f"step{i+1}/predicted_d", s["d"], dp, 1e-12)
            expected = load_npz(FIXTURES / f"native_step{i+1}.npz")
            for key in [
                "x",
                "velocity",
                "force",
                "material_force",
                "d",
                "kappa",
                "scalars",
                "external_force",
                "external_wall_force",
                "atom_ids",
                "native_bonds",
            ]:
                compare(
                    f"step{i+1}/{key}",
                    np.delete(s[key], 7) if key == "scalars" else s[key],
                    np.delete(expected[key], 7) if key == "scalars" else expected[key],
                    2e-18 if key == "x" else 1e-12,
                    0,
                )
            compare(
                f"step{i+1}/VV",
                s["velocity"][: ns.n] - p["velocity"][: ns.n],
                MIN_H
                * (p["force"][: ns.n] + s["force"][: ns.n])
                / (2 * ns.mass[: ns.n, None]),
                1e-10,
                0,
            )
            ns.tool.commit(ns.tool.evaluate(s["x"], s["velocity"], s["d"], s["time"]))
            states.append(s)
            if i == 1:
                ns.l.command(f"write_restart {tmp_path}/native.restart")
                ns.tool.save_checkpoint(tmp_path / "surface.npz", s)
        on = ns.field(13)
        off = ns.field(14)
    finally:
        ns.close()
    rs = NativeSolver(tmp_path / "restart", restart=tmp_path / "native.restart")
    try:
        p = rs.dynamics()
        rs.tool.restore_checkpoint(tmp_path / "surface.npz", p)
        compare("restart/onset", rs.field(13), on, 0, 0)
        compare("restart/complete", rs.field(14), off, 0, 0)
        for key in ["x", "velocity", "d", "kappa", "atom_ids"]:
            compare("restart/" + key, p[key], states[1][key], 0, 0)
        continued = rs.step(MIN_H)
        for key in ["x", "velocity", "force", "d", "kappa"]:
            compare(
                "restart/continued/" + key,
                continued[key],
                states[2][key],
                2e-18 if key == "x" else 1e-12,
                0,
            )
    finally:
        rs.close()


def test_damaged_native_checkpoint(tmp_path):
    from clot_sim.solver import NativeSolver

    checkpoint = FIXTURES / "pilot_checkpoint"
    expected = load_npz(FIXTURES / "pilot_restart_expected.npz")
    restart = tmp_path / "native.restart"
    with gzip.open(checkpoint / "native.restart.gz", "rb") as src, restart.open(
        "wb"
    ) as dst:
        shutil.copyfileobj(src, dst)
    assert (
        sha(restart)
        == json.loads((checkpoint / "RESTART_HASH.json").read_text())["raw_sha256"]
    )
    ns = NativeSolver(tmp_path / "damaged", restart=restart)
    try:
        s = ns.dynamics()
        ns.tool.restore_checkpoint(checkpoint / "surface.npz", s)
        assert np.any(s["d"] == 1) and np.any((s["d"] > 0) & (s["d"] < 1))
        for key in [
            "x",
            "velocity",
            "d",
            "kappa",
            "atom_ids",
            "scalars",
            "material_force",
            "external_force",
            "external_wall_force",
            "native_bonds",
            "Uc",
            "K",
        ]:
            compare(
                "damaged_restart/" + key,
                np.delete(s[key], 7) if key == "scalars" else s[key],
                np.delete(expected[key], 7) if key == "scalars" else expected[key],
                2e-18 if key in ["x", "Uc", "K", "scalars"] else 1e-12,
            )
        ledger = json.loads((checkpoint / "ledger.json").read_text())
        initial = load_npz(FIXTURES / "initial.npz")
        residual = energy_residual(ledger["W"], s, initial, ledger["Jc"])
        assert abs(residual) / ledger["scale"] <= 0.01
        ERRORS["damaged_restart/energy_relative_residual"] = float(
            abs(residual) / ledger["scale"]
        )
        # One diagnostic step from a genuinely damaged accepted checkpoint.
        # The original implementation produced five new failures in this step.
        metadata = json.loads((FIXTURES / "DAMAGED_STEP.json").read_text())
        old = s
        model = SofteningLPS(ns.ref, K=K, G=G)
        h = BASE
        while True:
            xp = old["x"].copy()
            xp[: ns.n] += h * (
                old["velocity"][: ns.n]
                + h * old["force"][: ns.n] / (2 * ns.mass[: ns.n, None])
            )
            kp, dp = predict(
                model, ns.ref, xp, old["kappa"], ns.field(13), ns.field(14)
            )
            if np.max(dp - old["d"]) <= 0.1:
                break
            assert h > MIN_H * (1 + 1e-12)
            h /= 2
        assert h == metadata["dt_s"]
        s = ns.step(h)
        expected_step = load_npz(FIXTURES / "damaged_step_expected.npz")
        for key in [
            "x",
            "velocity",
            "force",
            "material_force",
            "d",
            "kappa",
            "scalars",
            "external_force",
            "external_wall_force",
            "native_bonds",
            "Uc",
            "K",
        ]:
            compare(
                "damaged_step/" + key,
                np.delete(s[key], 7) if key == "scalars" else s[key],
                (
                    np.delete(expected_step[key], 7)
                    if key == "scalars"
                    else expected_step[key]
                ),
                2e-18 if key in ["x", "scalars", "Uc", "K"] else 1e-12,
            )
        compare("damaged_step/predicted_damage", s["d"], dp, 1e-10)
        compare("damaged_step/predicted_x", s["x"], xp, 2e-18, 0)
        assert np.all(s["d"] >= old["d"]) and np.all(s["kappa"] >= old["kappa"])
        new_failures = np.flatnonzero(
            (s["d"] == 1) & (old["d"] < 1) & ~ns.ref["initial_cut"]
        )
        assert len(new_failures) == metadata["new_failures"] == 5
        assert int((s["d"] > old["d"]).sum()) == metadata["damage_updates"]
        work = ns.tool.work(old, s)
        compare("damaged_step/work_channels", work, metadata["work_channels_J"], 1e-18)
        jump = Contact(ns.ref).energy_for_pairs(s["x"], ns.ref["bonds"][new_failures])
        compare("damaged_step/contact_jump", jump, metadata["contact_jump_J"], 1e-18)
        residual = energy_residual(
            ledger["W"] + work.sum(), s, initial, ledger["Jc"] + jump
        )
        assert abs(residual) / ledger["scale"] <= 0.01
        ERRORS["damaged_step/energy_relative_residual"] = float(
            abs(residual) / ledger["scale"]
        )
    finally:
        ns.close()


def test_no_historical_runtime_dependencies():
    import ast

    root = Path(__file__).resolve().parents[1]
    for p in (root / "clot_sim").glob("*.py"):
        s = p.read_text()
        tree = ast.parse(s)
        assert (
            "/workspace/y_clot" not in s and "/home/lzy/projects/temp_storage" not in s
        ), p
        assert "importlib" not in s, p
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                assert node.func.id not in ["exec", "eval", "compile"], p
            if isinstance(node, ast.ImportFrom):
                assert not any(n.name == "*" for n in node.names), p
    assert not any(p.is_symlink() for p in (root / "clot_sim").rglob("*"))


@pytest.fixture(scope="session", autouse=True)
def write_report():
    yield
    root = Path(__file__).resolve().parents[1]
    dump(root / "results/validation/max_absolute_errors.json", ERRORS)


def test_runner_ledger_checkpoint_and_resume(tmp_path):
    from clot_sim.main import run

    case = tmp_path / "runner"
    first = run("pilot", case, max_steps=2)
    assert first["status"] == "SHORT_VALIDATION_PASS" and first["formal_cases"] == 0
    checkpoint = case / "checkpoints/short_validation_step_000002"
    before = load_npz(checkpoint / "accepted.npz")
    ledger = json.loads((checkpoint / "ledger.json").read_text())
    nodework = load_npz(checkpoint / "node_work.npz")["node_work"]
    compare("runner/node_work_before", nodework.sum(), ledger["W"], 1e-17)
    resumed = run("continue", case, max_steps=1, checkpoint=checkpoint)
    assert (
        resumed["accepted_steps"] == 3 and resumed["status"] == "SHORT_VALIDATION_PASS"
    )
    final = load_npz(case / "checkpoints/short_validation_step_000003/accepted.npz")
    oracle = load_npz(FIXTURES / "native_step3.npz")
    for key in ["x", "velocity", "force", "d", "kappa"]:
        compare(
            "runner/resume/" + key,
            final[key],
            oracle[key],
            2e-18 if key == "x" else 1e-12,
            0,
        )
    after = json.loads(
        (case / "checkpoints/short_validation_step_000003/ledger.json").read_text()
    )
    assert after["W"] >= ledger["W"] and after["max_relative_R"] <= 0.01
    assert not np.any(before["d"] > final["d"])
    import csv

    rows = list(csv.DictReader((case / "history.csv").open()))
    assert [int(row["step"]) for row in rows] == [1, 2, 3]
