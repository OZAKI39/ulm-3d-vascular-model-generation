import numpy as np
import pytest
from fem3d.cap_remesh import check_wall
from fem3d.mesh_input import surface_topology
from stage015_helpers import ROOT,OUT,CANDIDATES,read


@pytest.mark.parametrize('candidate',CANDIDATES)
def test_frozen_wall_and_closed_connected_surface(candidate):
    source=np.load(ROOT/'inputs/stage01/tagged_surface_si.npz')
    derived=np.load(OUT/candidate/'surface/tagged_surface_si.npz')
    wall=check_wall(source,derived)
    assert wall['wall_triangle_count']==67071 and wall['maximum_wall_displacement_m']==0
    assert surface_topology(derived['points_m'],derived['triangles'])['connected_surface_components']==1


def test_moving_wall_vertex_or_changing_wall_tag_fails():
    source=np.load(ROOT/'inputs/stage01/tagged_surface_si.npz')
    for change in ('coordinate','tag'):
        derived={k:source[k].copy() for k in source.files}
        i=np.flatnonzero(source['facet_tags']==1)[0]
        if change=='coordinate': derived['points_m'][source['triangles'][i,0],0]+=1e-15
        else: derived['facet_tags'][i]=4
        with pytest.raises(ValueError): check_wall(source,derived)
