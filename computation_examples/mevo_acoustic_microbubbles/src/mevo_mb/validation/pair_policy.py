"""Explicit N=2 full-reference route; no superposition with vascular flow."""

from dataclasses import asdict
import json

import h5py
import numpy as np

from mevo_mb.flow.fields import ZeroFlowField
from .doinikov_bouakaz_pair import DoinikovBouakazPairReferenceModel, DoinikovBouakazParameters


def run_pair_reference(config, population, run, field, flow):
    if len(population.radii_m) != 2:
        raise ValueError("pair_reference permits exactly N=2; pairwise N-body superposition is forbidden")
    if not isinstance(flow, ZeroFlowField):
        raise ValueError("Published pair reference requires ZeroFlowField")
    if config["bubble_response"]["model"] != "doinikov2015_uncoated_reference":
        raise ValueError("Doinikov 2015 pair reference is uncoated; explicitly select its radial model")
    x = population.positions_m
    p = DoinikovBouakazParameters(**config["pair_reference"]["parameters"])
    if not np.allclose(population.radii_m, [p.radius1_m, p.radius2_m], rtol=1e-12, atol=0):
        raise ValueError("Population radii disagree with reference parameters")
    distance = np.linalg.norm(x[1] - x[0])
    if not np.isclose(distance, p.distance_m, rtol=1e-12, atol=0):
        raise ValueError("Initial distance disagrees with reference parameters")
    # No PZT field is silently ignored. This route uses the reference Gaussian
    # waveform only; the explicit field sentinel is required by the caller.
    if field is not None:
        raise ValueError(
            "pair_reference uses its prescribed uniform Gaussian pulse; pass acoustic_field=None"
        )
    z = DoinikovBouakazPairReferenceModel(p, config["pair_reference"].get("multipole_order", 8)).solve()
    direction = (x[1] - x[0]) / distance
    positions = x[0] + z["position_m"][..., None] * direction
    with h5py.File(run / "particles/trajectories.h5", "x") as h:
        for key, value in dict(
            time_s=z["time_s"],
            positions_m=positions,
            radii_m=population.radii_m,
            instantaneous_radii_m=z["radius_m"],
        ).items():
            h.create_dataset(key, data=value)
        h.attrs["complete"] = not z["contact_reached"]
        h.attrs["model"] = "DoinikovBouakaz2015 uncoated N=2 reference"
    summary = dict(
        status="CONTACT_STOP" if z["contact_reached"] else "COMPLETE",
        parameters=asdict(p),
        physics_validity=dict(
            physics_status="partially validated",
            near_field_policy="pair_reference",
            limitation="Prescribed pulse, spherical potential flow, no wall/lubrication/contact physics; Gaussian width assumed",
        ),
    )
    (run / "particles/SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary
