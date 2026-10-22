import numpy as np
from particle_3d.rbc_orientation import short_axis,advance_orientation,rotation_matrix


def test_real_omega_matches_independent_jeffery_formula_and_updates(real_cases,center_samples):
    for rows,summary in real_cases:
        q=np.array([[r[f"q_{a}"] for a in "wxyz"] for r in rows])
        p=np.array([[r[f"p_{a}"] for a in "xyz"] for r in rows])
        omega=np.array([[r[f"Omega_{a}_s_inv"] for a in "xyz"] for r in rows])
        assert np.isfinite(q).all() and np.isfinite(p).all() and np.isfinite(omega).all()
        np.testing.assert_array_equal(p,short_axis(q))
        sample=center_samples[summary["dt_index"]]
        e=.5*(sample.velocity_gradient_s_inv+sample.velocity_gradient_s_inv.transpose(0,2,1))
        expected=.5*sample.vorticity_s_inv+summary["jeffery_lambda"]*np.cross(p,np.einsum("nij,nj->ni",e,p))
        np.testing.assert_allclose(omega,expected,rtol=1024*np.finfo(float).eps,atol=1024*np.finfo(float).eps*np.max(np.abs(expected)))
        dt=summary["validation_dt_s"]
        # Each q step is recomputed from its own old position and geometry.
        for i in range(len(rows)-1):
            effective_dt=dt if i<len(rows)-2 else dt*summary["event"]["segment_fraction"]
            new=advance_orientation(q[i],omega[i],effective_dt)
            np.testing.assert_allclose(rotation_matrix(q[i+1]),rotation_matrix(new),rtol=0,atol=32*np.finfo(float).eps)
