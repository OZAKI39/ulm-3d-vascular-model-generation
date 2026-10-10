"""Verify unchanged-particle cache reuse during iterative contact projection."""
import h5py
import numpy as np
from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4.common import ROOT,new_run,require,finish
from scripts.stage4.experiments import setup,run_case


def main():
    baseline = require("performance/profile_B1")
    original = BifurcationGeometry.level_and_gradient

    def cached(self,positions):
        x = np.atleast_2d(np.asarray(positions,float))
        previous = getattr(self,"_profile_row_cache",None)
        if previous is not None and previous[0].shape == x.shape:
            changed = np.any(x != previous[0],axis=1)
            value,gradient = previous[1].copy(),previous[2].copy()
            if changed.any():
                value[changed],gradient[changed] = original(self,x[changed])
        else:
            value,gradient = original(self,x)
        self._profile_row_cache = (x.copy(),value.copy(),gradient.copy())
        return value,gradient

    BifurcationGeometry.level_and_gradient = cached
    provided = setup()
    out = new_run("performance/cache_row_trial",dict(duration_s=.2,case="B1",method="Exact per-row equality across contact projection iterations"),["performance/profile_B1"])
    row = run_case(out,"B1",provided[1]["particles"]["seed"],duration=.2,provided=provided)
    with h5py.File(ROOT / baseline["row"]["run"] / "particles/trajectories.h5") as a,h5py.File(ROOT / row["run"] / "particles/trajectories.h5") as b:
        results = {k:bool(np.array_equal(a[k][:],b[k][:],equal_nan=True)) for k in a}
    finish(out,dict(status="PASS" if all(results.values()) else "FAIL",baseline=baseline["row"],candidate=row,
                   datasets_bitwise_equal=results,scope="Short case includes active acoustic forces and iterative geometric contacts; no scientific values approximated"))


if __name__ == "__main__":
    main()
