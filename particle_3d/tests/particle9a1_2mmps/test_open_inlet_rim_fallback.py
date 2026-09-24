import numpy as np

def test_true_open_edge_and_vertex_fallback(rim_case):
 for x,feature in [([5e-6,-1e-8,1.1e-7],'EDGE'),([-1e-8,-1e-8,1.1e-7],'VERTEX')]:
  base,new,g=rim_case(x)
  assert g.wall_feature==feature
  np.testing.assert_array_equal(base.matrix.toarray(),new.matrix.toarray())
  np.testing.assert_array_equal(base.rhs,new.rhs)
  assert new.planar_diagnostics[0]['open_rim_fallback']['open_boundary_roles']==['INLET']
