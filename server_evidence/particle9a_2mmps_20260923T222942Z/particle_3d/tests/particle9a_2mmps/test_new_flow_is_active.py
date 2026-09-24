from particle_3d.particle9a_provenance import require_current_flow,FLOW_SHA,FEM
from particle_3d.audit import read_frozen

def test_new_flow_is_active():
    assert require_current_flow()['flow_sha256']==FLOW_SHA
    p,mesh,flow,boundaries=read_frozen(FEM)
    assert p['case_role']=='MEAN_2P0_MMPS'
    assert mesh.n_cells==371402 and mesh.n_points==70363
    assert set(boundaries)=={'WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03'}
