import pytest
from stage015_helpers import REPORT,OUT,read


@pytest.mark.parametrize('ranks',[1,2])
def test_selected_mesh_new_process_roundtrip(ranks):
    decision=read(REPORT/'candidate_selection.json')
    if decision['selected_candidate'] is None:
        pytest.skip('Stage 1.5 geometry gate failed: no volume/winner exists; DOLFINx round-trip forbidden, overall stage FAIL')
    result=read(OUT/'selected/qc'/f'reload_r{ranks}.json')
    assert result['status']=='PASS' and result['mpi_ranks']==ranks
    assert result['tdim']==result['gdim']==3 and result['cell_type']=='tetrahedron'
