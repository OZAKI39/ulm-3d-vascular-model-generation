import pytest
from fem3d.cap_selection import require_volume_permission
from stage015_helpers import OUT,CANDIDATES,read


@pytest.mark.parametrize('candidate',CANDIDATES)
def test_rejected_geometry_cannot_reach_volume_meshing(candidate):
    s=read(OUT/candidate/'qc/surface_invariants.json')
    if s['status']=='FAIL':
        with pytest.raises(ValueError,match='Geometry gate failed'): require_volume_permission(s)
        assert not (OUT/candidate/'mesh/fluid.msh').exists()
    else:
        require_volume_permission(s)
        assert (OUT/candidate/'mesh/volume_mesh.npz').is_file()
