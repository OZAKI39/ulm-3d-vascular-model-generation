"""Measure an exact repeated-position cache before enabling it in new formal runs."""
import h5py
import numpy as np
from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4.common import ROOT, new_run, require, finish
from scripts.stage4.experiments import setup,run_case


def main():
    baseline = require("performance/profile")
    original = BifurcationGeometry.level_and_gradient

    def cached(self,positions):
        x = np.atleast_2d(np.asarray(positions,float))
        previous = getattr(self,"_profile_cache",None)
        if previous is not None and np.array_equal(x,previous[0]):
            return previous[1],previous[2]
        value,gradient = original(self,x)
        self._profile_cache = (x.copy(),value,gradient)
        return value,gradient

    BifurcationGeometry.level_and_gradient = cached
    provided = setup()
    out = new_run("performance/cache_trial",dict(duration_s=.02,method="Exact equality keyed last-position cache; no rounded coordinates"),["performance/profile"])
    row = run_case(out,"B0",provided[1]["particles"]["seed"],duration=.02,provided=provided)
    results = {}
    with h5py.File(ROOT / baseline["row"]["run"] / "particles/trajectories.h5") as a, h5py.File(ROOT / row["run"] / "particles/trajectories.h5") as b:
        assert set(a) == set(b)
        for name in a:
            results[name] = bool(np.array_equal(a[name][:],b[name][:],equal_nan=True))
    finish(out,dict(status="PASS" if all(results.values()) else "FAIL",bitwise_dataset_equality=results,
                   baseline=baseline["row"],candidate=row,
                   note="Baseline cProfile timing includes instrumentation; candidate wall time is an approximate benefit estimate, not a rigorous speedup benchmark"))


if __name__ == "__main__":
    main()
