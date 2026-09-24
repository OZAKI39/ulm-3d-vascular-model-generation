import numpy as np

def test_shared_wall_edge_stays_active(rim_case):
 base,new,g=rim_case([5e-6,5e-6,1.1e-7])
 assert g.wall_feature=='EDGE'
 assert new.planar_diagnostics[0]['wall_weight']>0
 assert np.linalg.norm(new.planar_matrix.toarray())>0
