"""One analytic RMS-matched phase-only control after the initial patterns were ineligible."""

import argparse
import json
import numpy as np
from mevo_mb.fast_sim.models import Basis
from mevo_mb.fast_sim.morphodynamics import morphology_series, contrast
from .common import BASE, PROTOCOL, write
from .experiment import run_case
from .mature import load, analyze_topology


def main(label, mature):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    source = BASE / mature
    field = json.loads((source / "FIELD_SELECTION.json").read_text())
    base = field["baseline"]
    basis = Basis(BASE / "inputs/basis.h5", 2.5e6)
    points = np.load(source / "FIELD_DIAGNOSTICS.npz")["points_m"]

    def evaluate(phase, enabled):
        f, _, _ = basis.combine(base["amplitude_Pa"], phase, enabled)
        return f.interpolate_complex_pressure(points), f.interpolate_complex_gradient(points)

    p0, g0 = evaluate([0] * 4, [True] * 4)
    a, _ = evaluate([np.pi / 6, 0, 0, 0], [True, False, True, True])
    b, _ = evaluate([0] * 4, [False, True, False, False])
    cross = np.mean(a * np.conj(b))
    target = np.mean(abs(p0) ** 2)
    offset = np.mean(abs(a) ** 2 + abs(b) ** 2)
    cosine = float((target - offset) / (2 * abs(cross)))
    assert -1 <= cosine <= 1
    roots = [float(np.angle(cross) + s * np.arccos(cosine)) for s in [-1, 1]]
    theta = min(roots, key=abs)
    event = dict(base, time_s=0, phase_rad=[float(np.pi / 6), theta, 0.0, 0.0])
    p, g = evaluate(event["phase_rad"], event["enabled"])
    np.testing.assert_allclose(np.mean(abs(p) ** 2), target, rtol=1e-12)
    assert min(abs(p)) >= 0.1 and max(abs(p)) <= 40000
    record = dict(
        event=event,
        analytic_roots_rad=roots,
        selected_root_rule="Smallest absolute phase change; no trajectory outcome used",
        rms_matching_equation="mean|a+b exp(i theta)|^2=D+2 Re(B exp(-i theta))",
        local_RMS_pressure_Pa=float(np.sqrt(np.mean(abs(p) ** 2))),
        pressure_relative_difference=float(abs(np.sqrt(np.mean(abs(p) ** 2) / target) - 1)),
        normalized_gradient_contrast=float(np.linalg.norm(g - g0) / np.linalg.norm(g0)),
        reason="All three initial phase-only patterns failed the unchanged 10% pressure comparability bound; one algebraic iso-RMS construction, not an outcome-driven control sweep",
        source_amplitudes_and_enabled_unchanged=True,
    )
    write(out / "FIELD_SELECTION.json", record)
    np.savez_compressed(
        out / "FIELD_DIAGNOSTICS.npz",
        points_m=points,
        baseline_pressure_Pa=p0,
        baseline_gradient_Pa_m=g0,
        pressure_Pa=p,
        gradient_Pa_m=g,
    )
    cfg = json.loads((source / "unchanged/CONFIG.json").read_text())
    cfg["events"] = [event]
    row = run_case(out / "phase", cfg)
    d = load(out / "phase")
    reference = load(source / "unchanged")
    np.testing.assert_array_equal(d["positions_m"][0], reference["positions_m"][0])
    np.testing.assert_array_equal(d["time_s"], reference["time_s"])
    ids = np.load(source / "MATURE_STATE.npz")["initial_cluster_ids"]
    top, _, _ = analyze_topology(out / "phase", d)
    series = morphology_series(
        d["positions_m"][:, 0], ids, d["labels"], PROTOCOL["observation"]["principal_axis_relative_eigengap"]
    )
    np.savez_compressed(out / "phase/MORPHOLOGY.npz", time_s=d["time_s"], fixed_particle_ids=ids, **series)
    diffs = [contrast(a, b, ids) for a, b in zip(reference["positions_m"][:, 0], d["positions_m"][:, 0])]
    np.savez_compressed(
        out / "phase/CONTRAST.npz",
        time_s=d["time_s"],
        **{k: np.asarray([m[k] for m in diffs]) for k in diffs[0]},
    )
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            run=row,
            comparison=diffs[-1],
            topology=top,
            field_selection=record,
            absolute_switch_time_s=cfg["absolute_time_offset_s"],
            mature_source=mature,
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--mature", required=True)
    a = ap.parse_args()
    main(a.label, a.mature)
