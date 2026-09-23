from stage03_helpers import *

def test_actual_spaces_equal_topological_p2_p1_cost_and_one_global_real():
    c=config();mesh=np.load(ROOT/c['mesh']['source']/'mesh/volume_mesh.npz')
    tetra=mesh['tetra'];edges=np.unique(np.sort(np.concatenate([tetra[:,[i,j]] for i,j in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]]),axis=1),axis=0)
    nv=len(mesh['points_m']);a=preflight()
    assert a['velocity_dofs']==3*(nv+len(edges)) and a['pressure_dofs']==nv and a['real_dofs']==1
    assert sum(row['real_owned_dofs'] for row in a['ranks'])==1
