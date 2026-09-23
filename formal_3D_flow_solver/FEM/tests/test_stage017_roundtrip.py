from stage017_helpers import *
import pytest
from fem3d.audit import sha256

@pytest.mark.parametrize('ranks',[1,2])
def test_new_process_reload_full_tetra_tags_contract_and_recomputed_quality(ranks):
 p=OUT/'selected';r=read(p/'qc'/f'reload_r{ranks}.json');save=read(p/'metadata/dolfinx_save.json')
 assert r['status']=='PASS' and r['mpi_ranks']==ranks and r['tdim']==r['gdim']==3
 assert r['tetra_geometry_exactly_equal'] and r['quality_summary_matches'] and r['positive_geometric_volume']
 assert r['no_mesh_generation_in_reload'] and not r['fem_space_created']
 assert all(item['pid']!=save['pid'] for item in r['ranks'])
 assert r['xdmf_sha256']==save['xdmf_sha256']==sha256(p/'mesh/fluid.xdmf')
 assert r['hdf5_sha256']==save['hdf5_sha256']==sha256(p/'mesh/fluid.h5')
 assert sha256(p/'planar_port_contract_v2.json')==BASE['contract_sha256']
 assert all(v['status']=='PASS' for v in r['ports'].values())
