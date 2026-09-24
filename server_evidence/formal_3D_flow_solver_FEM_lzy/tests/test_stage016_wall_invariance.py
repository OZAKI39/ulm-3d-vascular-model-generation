from stage016_helpers import *
import pytest
from fem3d.cap_remesh import check_wall

@pytest.mark.parametrize('name',CANDIDATES)
def test_all_67071_wall_triangles_exact(name):
 r=check_wall(source(),candidate(name))
 assert r['wall_triangle_count']==67071 and r['maximum_wall_displacement_m']==0

def test_moved_wall_vertex_rejected():
 d={k:v.copy() for k,v in candidate('sparse_C').items()};i=d['triangles'][d['facet_tags']==1][0,0];d['points_m'][i,0]+=1e-15
 with pytest.raises(ValueError,match='Wall vertex'):check_wall(source(),d)
