import numpy as np
from particle_3d.particle2_cases import NORM_ATOL


def test_all_real_quaternion_and_short_axis_norms(real_cases):
    for rows,summary in real_cases:
        q=np.array([[r[f"q_{a}"] for a in "wxyz"] for r in rows])
        p=np.array([[r[f"p_{a}"] for a in "xyz"] for r in rows])
        assert np.max(abs(np.linalg.norm(q,axis=1)-1))<=NORM_ATOL
        assert np.max(abs(np.linalg.norm(p,axis=1)-1))<=NORM_ATOL
        assert np.min(np.sum(q[:-1]*q[1:],axis=1))>=0
