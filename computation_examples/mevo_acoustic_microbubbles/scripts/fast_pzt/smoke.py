"""First end-to-end GPU prototype, not the final acceptance campaign."""

import argparse
import json
from dataclasses import replace
import numpy as np
import warp as wp
from mevo_mb.fast_sim.models import Physics, ResponseSection, pair_forces_numpy, random_cloud
from mevo_mb.fast_sim.engine import FastSimulator
from .common import BASE, inputs, write, source_identity


def main(label):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    p = Physics()
    response = ResponseSection(inputs()[0], p)
    x = np.array([[-5e-6, 0, 0], [5e-6, 0, 0]])
    sim = FastSimulator(x, p, response)
    sim.control({"time_s": 0, "pressure_Pa": 1000})
    sim.fields()
    sim.pairs()
    wp.synchronize()
    gpu = sim.force.numpy()
    v = response.query(np.array([1000, 1000]))
    ref = pair_forces_numpy(x, v, p)
    np.testing.assert_allclose(gpu, ref, rtol=2e-12, atol=1e-26)
    write(
        out / "GPU_FORMULA.json",
        dict(
            status="PASS",
            gpu_N=gpu,
            reference_N=ref,
            max_absolute_N=np.max(np.abs(gpu - ref)),
            compile_s=sim.compile_s,
        ),
    )
    sim.reset()
    rows = []
    for n, secondary, duration in [(2, False, 0.005), (2, True, 0.01), (32, True, 0.05)]:
        pp = replace(p, secondary=secondary)
        xx = x if n == 2 else random_cloud(n, p.radius_m, 20e-6, 20261007)
        run = FastSimulator(xx, pp, response).run(
            duration, 20e-6, [{"time_s": 0, "pressure_Pa": 1000}], sample_interval_s=0.0005
        )
        np.savez_compressed(
            out / f"N{n}_secondary{secondary}.npz", time_s=run["time_s"], positions_m=run["positions_m"]
        )
        from scipy.spatial.distance import pdist

        row = dict(
            run["summary"], final_min_gap_m=float(pdist(run["positions_m"][-1, 0]).min() - 2 * p.radius_m)
        )
        rows.append(row)
        print(json.dumps(row), flush=True)
    write(out / "RESULT.json", dict(source=source_identity(), response=response.metadata, rows=rows))


if __name__ == "__main__":
    a = argparse.ArgumentParser()
    a.add_argument("--label", default="smoke_v1")
    main(a.parse_args().label)
