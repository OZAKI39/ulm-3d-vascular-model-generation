import numpy as np
from sv_validation.mesh_diagnostics import tetra_sicn

def test_quality_metric_and_actual_distribution(root, report):
    ideal=np.array([[0.,0.,0.],[1.,0.,0.],[.5,np.sqrt(3)/2,0.],[.5,np.sqrt(3)/6,np.sqrt(2/3)]])
    assert abs(tetra_sicn(ideal,np.array([[0,1,2,3]]))[0]-1) < 1e-14
    with np.load(root/"outputs/sv1/SV_MESH/mesh_arrays.npz") as d:
        q=tetra_sicn(d["points_m"],d["tetra"])
    a=report("mesh_quality")
    assert np.isfinite(q).all() and q.min() > 0
    assert abs(q.min()-a["q_min"]) < 1e-14
    assert np.allclose(np.quantile(q,[.01,.05,.5,.95]),[a[k] for k in ("P1","P5","median","P95")],rtol=1e-13,atol=0)
    assert int((q<.1).sum()) == a["N_low"]
