import numpy as np
from stage02_helpers import read_case


def test_compatible_reference_approaches_exact_solution_under_refinement():
    rows=[read_case(p+"_reference","analytic_comparison") for p in ("pipe_coarse","pipe_medium","pipe_fine")]
    for key in ("velocity_L2_relative_error","velocity_global_L2_relative_error","pressure_profile_error","lambda_relative_error"):
        values=np.array([r[key] for r in rows])
        assert np.all(np.diff(values)<0)
        assert values[-1]<values[0]/3
