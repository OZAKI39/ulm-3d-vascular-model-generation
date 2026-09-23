import numpy as np
import pytest
from fem3d.mesh_qc import boundary_partition


def test_all_exterior_facets_have_exactly_one_original_tag(stage01_mesh):
    q=stage01_mesh["qc"]["topology"]
    assert q["exterior_facet_count"]==67262
    assert q["unlabelled_exterior_facets"]==q["multiply_labelled_exterior_facets"]==0
    assert q["tagged_nonexterior_facets"]==0
    assert q["boundary_counts"]=={"WALL":67071,"INLET":56,"OUTLET_01":44,"OUTLET_02":42,"OUTLET_03":49}


@pytest.mark.parametrize("fault",["missing","duplicate","unknown_tag"])
def test_incomplete_or_ambiguous_partition_fails(fault):
    points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
    cells=np.array([[0,1,2,3]])
    faces=np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]])
    tags=np.ones(4,dtype=int)
    if fault=="missing": faces,tags=faces[:-1],tags[:-1]
    elif fault=="duplicate": faces,tags=np.vstack([faces,faces[0]]),np.r_[tags,2]
    else: tags[0]=9
    with pytest.raises(ValueError,match="Boundary partition failed"):
        boundary_partition(points,cells,faces,tags)
