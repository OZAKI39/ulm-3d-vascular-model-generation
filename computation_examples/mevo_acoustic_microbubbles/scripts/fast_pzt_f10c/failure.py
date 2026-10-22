"""Single bounded replay to preserve the FIRST rejected Q2-S proposal.

Stops per physical step; never continues a failed trajectory or relaxes a gate.
Save interval remains 50 us for exact comparison to the previous accepted prefix.
"""

import argparse
import h5py
import numpy as np
from .common import BASE, case, identity, objects, write
from mevo_mb.fast_sim.f10c_multirate import MultirateBatch


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    p, b, r, g, f = objects()
    cfg = case("S_1000")
    x = np.array(cfg["environment"]["positions_m"])
    s = MultirateBatch(x, p, r, b, flow=f, vessel=g, field_stride=1, flow_stride=2, graph_chunk=1)
    s.control(cfg["events"][0])
    frames = [x.copy()]
    times = [0.0]
    for step in range(1000):
        try:
            s.advance(1e-5, 1)
        except RuntimeError as exc:
            proposal = s.tmp1.numpy()
            dist = np.linalg.norm(proposal[:, None] - proposal[None, :], axis=-1)
            np.fill_diagonal(dist, np.inf)
            index = np.unravel_index(dist.argmin(), dist.shape)
            write(
                out / "REJECTION.json",
                dict(
                    status="EXPECTED_REJECTION",
                    sample_id="Q2_S_1000_0000",
                    first_rejected_step=step + 1,
                    attempted_time_s=(step + 1) * 1e-5,
                    last_accepted_step_time_s=step * 1e-5,
                    last_saved_accepted_time_s=times[-1],
                    error_bits=int(s.error.numpy()[0]),
                    reason=str(exc),
                    proposed_max_overlap_m=2 * p.radius_m - dist.min(),
                    contact_budget_m=0.01 * p.radius_m,
                    closest_particle_ids=index,
                    quality=s.quality(),
                    note="First error attribution retained. Original chunked result also carried secondary bit16384 because its invalid flow buffer was checked again after the prior contact error; guard now preserves the first error. No accepted state or physical tolerance changed.",
                ),
            )
            np.savez_compressed(
                out / "rejected_proposal.npz",
                proposal_m=proposal,
                last_accepted_m=s.x.numpy(),
                diag=s.diag.numpy(),
                flow_stats=s.flow_query.stats.numpy(),
                times_s=np.array(times),
                positions_m=np.array(frames),
            )
            assert int(s.error.numpy()[0]) == 2048 and 2 * p.radius_m - dist.min() > 0.01 * p.radius_m
            break
        if (step + 1) % 5 == 0:
            frames.append(s.state()[0].copy())
            times.append((step + 1) * 1e-5)
    else:
        raise AssertionError("Expected recorded Q2 S rejection not reproduced")
    with h5py.File(BASE / "multirate_v1/Q2_S_1000/trajectory.h5") as h:
        np.testing.assert_array_equal(h["positions_m"][:, 0], np.array(frames))
    write(
        out / "COMPLETE.json",
        dict(status="PASS", accepted_prefix_bitwise=True, physical_budget_unchanged=True),
    )
    print("PASS first rejected proposal retained", flush=True)


if __name__ == "__main__":
    main()
