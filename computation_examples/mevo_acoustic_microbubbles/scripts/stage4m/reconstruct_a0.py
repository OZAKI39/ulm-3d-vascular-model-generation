"""Read-only force diagnostics on every existing A0 accepted state; no stepping."""

import argparse
import json
from pathlib import Path
import platform
import time

import h5py
import numpy as np

from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel
from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
from scripts.stage4m.common import BASE, refs, sha, write


def reconstruct(legacy, stage4l, case):
    old = refs(legacy, stage4l, case)["dense_A0"]
    identity = json.loads((old / "IDENTITY.json").read_text())
    config = identity["config"]

    def resolve(raw):
        p = legacy / str(raw).removeprefix("/workspace/mevo_acoustic_microbubbles/")
        if sha(p) != identity["input_sha256"][str(raw)]:
            raise ValueError("A0 frozen input hash mismatch")
        return p

    acoustic = resolve(config["effective_acoustic_field"]["path"])
    response = resolve(config["bubble_response"]["path"])
    field = ComplexAcousticField.load(acoustic)
    table = AdaptiveResponseTableBubbleModel(response, backend="numpy")
    table.validate_medium(config["medium"])
    pair = ClassicalFarFieldBjerknes(enforce_far_field=False, record_pair_diagnostics=True)
    out = BASE / "comparisons/dense" / case / "A0"
    out.mkdir(parents=True, exist_ok=True)
    file = out / "PAIR_FORCE_DIAGNOSTICS.jsonl"
    start = time.perf_counter()
    with h5py.File(old / "trajectories.h5", "r") as h, file.open("x") as stream:
        if not h.attrs["complete"]:
            raise ValueError("Existing A0 probe is incomplete")
        t, radii = h["time_s"][:], h["radii_m"][:]
        positions, masks = h["positions_m"][:], h["active"][:]
        for k, time_s in enumerate(t):
            active = masks[k]
            x, r = positions[k][active], radii[active]
            pressure = field.interpolate_complex_pressure(x)
            v = table.harmonic_response(r, field.frequency_Hz, pressure)["harmonics_m3"]
            force = pair.forces(x, r, v, field.frequency_Hz, config["medium"]["density_kg_m3"])
            if k == 0:
                with np.load(old / "INITIAL_FORCES_AND_VELOCITIES.npz") as initial:
                    reference = initial["secondary_force_N"]
                comparison = dict(
                    bitwise_equal=bool(np.array_equal(force, reference)),
                    max_abs_difference_N=float(np.max(np.abs(force - reference))),
                    relative_L2=float(np.linalg.norm(force - reference) / np.linalg.norm(reference)),
                )
            stream.write(
                json.dumps(dict(time_s=float(time_s), **pair.pair_diagnostics()), allow_nan=False) + "\n"
            )
            if k % 2000 == 0:
                stream.flush()
                print(
                    json.dumps(
                        dict(
                            case=case, analyzed_states=k, total=len(t), elapsed_s=time.perf_counter() - start
                        )
                    ),
                    flush=True,
                )
    result = dict(
        case=case,
        policy="A0",
        status="READ_ONLY_REANALYSIS_COMPLETE",
        evaluated_existing_states=len(t),
        integration_steps=0,
        translation_solver_calls=0,
        initial_force_regression=comparison,
        raw_trajectory_sha256=sha(old / "trajectories.h5"),
        acoustic_sha256=sha(acoustic),
        response_sha256=sha(response),
        diagnostics_sha256=sha(file),
        python=platform.python_version(),
        numpy=np.__version__,
        processor="NumPy postprocessing of stored coordinates with the frozen original field and response table",
    )
    write(out / "READ_ONLY_PAIR_RECONSTRUCTION.json", result)
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--legacy-root", type=Path, required=True)
    p.add_argument("--stage4l-root", type=Path, required=True)
    p.add_argument("--case", choices=["B2", "B3"], required=True)
    a = p.parse_args()
    reconstruct(a.legacy_root, a.stage4l_root, a.case)
