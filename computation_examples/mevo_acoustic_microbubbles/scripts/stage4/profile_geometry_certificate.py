"""Diagnostic only: prove interior wall projections are inactive before skipping queries."""
import cProfile
import io
import json
import pstats
import h5py
import numpy as np
from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4.common import ROOT,new_run,require,finish
from scripts.stage4.experiments import setup,run_case


def main():
    baseline = require("performance/profile_B2_secondary")
    original = BifurcationGeometry.project
    counts = dict(calls=0,certified_inactive=0,minimum_certified_margin_m=None)

    def certified(self,positions,radii):
        counts["calls"] += 1
        x = np.asarray(positions)
        previous = self._distance_cache
        if self._mesh_distance is not None and previous is not None and previous[0].shape==x.shape:
            # Euclidean distance to a fixed closed surface is 1-Lipschitz.
            # A ball about a known interior point cannot cross that surface if
            # the displacement plus bubble radius is below the old distance.
            # The extra margin makes this MORE conservative than the old
            # projection tolerance; it is never a relaxed wall constraint.
            margin = -previous[1]-np.linalg.norm(x-previous[0],axis=1)-np.asarray(radii)
            unit_gradient = np.abs(np.linalg.norm(previous[2],axis=1)-1.)<1e-10
            if np.all((previous[1]<0)&unit_gradient&(margin>1e-10)):
                counts["certified_inactive"] += 1
                smallest = float(margin.min())
                counts["minimum_certified_margin_m"] = smallest if counts["minimum_certified_margin_m"] is None else min(smallest,counts["minimum_certified_margin_m"])
                # Do not manufacture a distance/gradient or update the exact
                # cache with a bound. Only this no-op projection is skipped.
                return x.copy(),np.zeros(len(x),bool)
        return original(self,positions,radii)

    BifurcationGeometry.project = certified
    provided = setup()
    duration = json.loads((ROOT/baseline["row"]["run"]/"particles/SUMMARY.json").read_text())["duration_s"]
    out = new_run("performance/certified_wall_trial",dict(duration_s=duration,secondary=True,
        method="Conservative interior-distance certificate; diagnostic patch only",roundoff_margin_m=1e-10),
        ["performance/profile_B2_secondary"])
    profiler = cProfile.Profile()
    profiler.enable()
    row = run_case(out,"B2",provided[1]["particles"]["seed"],secondary=True,duration=duration,provided=provided)
    profiler.disable()
    profiler.dump_stats(str(out/"tracker.prof"))
    stream = io.StringIO()
    pstats.Stats(profiler,stream=stream).sort_stats("cumulative").print_stats(35)
    (out/"PROFILE.txt").write_text(stream.getvalue())
    with h5py.File(ROOT/baseline["row"]["run"]/"particles/trajectories.h5") as a,h5py.File(ROOT/row["run"]/"particles/trajectories.h5") as b:
        equality = {key:bool(np.array_equal(a[key][:],b[key][:],equal_nan=True)) for key in a}
    events_equal = (ROOT/baseline["row"]["run"]/"particles/events.jsonl").read_bytes()==(ROOT/row["run"]/"particles/events.jsonl").read_bytes()
    validity_equal = json.loads((ROOT/baseline["row"]["run"]/"physics_validity.json").read_text())==json.loads((ROOT/row["run"]/"physics_validity.json").read_text())
    finish(out,dict(status="PASS" if all(equality.values()) and events_equal and validity_equal else "FAIL",
        baseline=baseline["row"],candidate=row,datasets_bitwise_equal=equality,
        events_bitwise_equal=events_equal,validity_equal=validity_equal,certificate=counts,
        scope="Instrumented short secondary-ON comparison. Production source is unchanged; no field interpolation or distance approximation."))


if __name__=="__main__":
    main()
