"""Explicit scientific disposition; never overwrite the first candidate's outcomes."""

import argparse
import numpy as np
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.contact_models import mobility_tensors, contact_pass
from .common import BASE, write
from .sensitivity import compare


def main(label, sensitivity):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    parent = BASE / sensitivity
    p = F2Physics(contact_model="f2_mobility_aware_contact")
    x = np.array([[0.0, 0, 0], [3.4e-6, 0, 0], [1.7e-6, 2.8e-6, 0], [1.7e-6, 1e-6, 2.5e-6]])
    y = contact_pass(x, mobility_tensors(x, p), p)
    r = x - x.mean(0)
    angular = np.cross(r, y - x).sum(0)
    inertia = np.eye(3) * (r * r).sum() - r.T @ r
    rotation = np.linalg.solve(inertia, angular)
    audit = dict(
        status="REJECTED_AS_DEFAULT",
        reason="Anisotropic equal/opposite contact displacements are not central and can supply spurious rigid rotation; 20/10 us also fail overlap budget",
        configuration="Fixed four-bubble irregular overlap state, independent of trajectory target",
        initial_positions_m=x,
        corrected_positions_m=y,
        angular_moment_m2=angular,
        infinitesimal_rigid_rotation_rad=rotation,
        rotation_norm_rad=float(np.linalg.norm(rotation)),
        linear_centroid_shift_m=np.linalg.norm(y.mean(0) - x.mean(0)),
        contact_model="f2_mobility_aware_contact",
        candidate_still_available=True,
        qualified_for_PZT_control_conclusions=False,
    )
    write(out / "CONTACT_INVARIANT_AUDIT.json", audit)
    near = compare(parent / "f1_10us", parent / "near_10us")
    contact = compare(parent / "near_5us", parent / "contact_5us")
    assert near["numerical_quality"] == ["PASS", "PASS"] and contact["numerical_quality"] == ["PASS", "PASS"]
    chosen = dict(
        dt_s=10e-6,
        physics=dict(force_model="f2_near_field_candidate", contact_model="f1_pbd_contact"),
        reason="First qualified near-field step; preserve F1 central contact because the tensor contact candidate fails angular-balance audit",
        near_field_comparison=dict(baseline=str(parent / "f1_10us"), candidate=str(parent / "near_10us")),
        contact_comparison=dict(baseline=str(parent / "near_5us"), candidate=str(parent / "contact_5us")),
    )
    write(out / "SELECTED.json", chosen)
    write(
        out / "RESULT.json",
        dict(complete=True, selected=chosen, near_field=near, contact=contact, contact_audit=audit),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--sensitivity", required=True)
    a = ap.parse_args()
    main(a.label, a.sensitivity)
