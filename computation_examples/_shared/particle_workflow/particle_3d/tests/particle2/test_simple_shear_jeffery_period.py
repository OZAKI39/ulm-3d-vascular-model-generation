import json
import numpy as np
from particle_3d.particle2_cases import simple_shear_case,SHEAR_ANGLE_STEPS,NORM_ATOL


def test_five_distribution_geometries_match_independent_axis_period(geometries,p2_data):
    # Re-execute quaternion integration, not just a saved PASS flag.
    _,metrics=simple_shear_case(geometries,SHEAR_ANGLE_STEPS[0])
    saved=json.loads((p2_data/"06_shear_dt0_metrics.json").read_text())
    for geometry,m,record in zip(geometries,metrics,saved):
        analytic=np.pi/20*(geometry.r+1/geometry.r)
        assert m["analytic_period_s"]==analytic
        assert m["passed"] and m["relative_error"]<=m["period_relative_error_bound"]
        assert m["measured_period_s"]==record["measured_period_s"]
        assert m["max_p_norm_error"]<=NORM_ATOL and m["max_quaternion_norm_error"]<=NORM_ATOL
    assert np.all(np.diff([m["measured_period_s"] for m in metrics])<0)


def test_halving_uses_three_validation_only_steps_and_no_directed_period_confusion(p2_data):
    records=[json.loads((p2_data/f"06_shear_dt{i}_metrics.json").read_text()) for i in range(3)]
    for a,b,c in zip(*records):
        assert a["validation_dt_s"]==2*b["validation_dt_s"]==4*c["validation_dt_s"]
        assert all(m["timestep_role"]=="VALIDATION_ONLY" for m in [a,b,c])
        assert "directed vector period is 2*T_axis" in a["period_semantics"]
        assert abs(c["measured_period_s"]/c["analytic_period_s"]-1)<c["period_relative_error_bound"]
