import json
import numpy as np
from particle_3d.rbc_orientation import shape_axis_distance


def test_saved_axis_distance_uses_abs_dot_not_raw_quaternions(p2_data,real_cases):
    summaries=json.loads((p2_data/"10_axis_comparison_summary.json").read_text())
    assert len(summaries)==15
    for s in summaries:
        g,d=s["geometry_index"],s["dt_index"]
        if d==2:
            assert s["max_axis_difference_rad"]==s["median_axis_difference_rad"]==0
            continue
        ref=real_cases[g*3+2][0];rows=real_cases[g*3+d][0]
        rt=np.array([r["time_s"] for r in ref]);t=np.array([r["time_s"] for r in rows]);keep=t<=rt[-1]
        p=np.array([[r[f"p_{a}"] for a in "xyz"] for r in rows])[keep]
        reference=np.column_stack([np.interp(t[keep],rt,[r[f"p_{a}"] for r in ref]) for a in "xyz"])
        gap=shape_axis_distance(p,reference)
        assert s["max_axis_difference_rad"]==float(np.max(gap))
        assert s["median_axis_difference_rad"]==float(np.median(gap))
        np.testing.assert_array_equal(gap,shape_axis_distance(-p,reference))
        assert not s["production_particle_timestep_frozen"]
