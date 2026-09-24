from stage016_helpers import *
import pytest
@pytest.mark.parametrize('ranks',[1,2])
def test_selected_new_process_DOLFINx_roundtrip(ranks):
 decision=read(REPORT/'candidate_selection.json')
 if decision['selected_candidate'] is None:
  assert not (OUT/'selected').exists()
  pytest.skip('Stage 1.6 FAIL at surface density gates; no selected volume exists; roundtrip forbidden')
 r=read(OUT/'selected/qc'/f'reload_r{ranks}.json');save=read(OUT/'selected/metadata/dolfinx_save.json')
 assert r['status']=='PASS' and r['mpi_ranks']==ranks and r['tdim']==r['gdim']==3
 assert r['tetra_geometry_exactly_equal'] and r['positive_geometric_volume']
 assert all(v['status']=='PASS' for v in r['ports'].values())
 assert all(v['pid']!=save['pid'] for v in r['ranks'])
