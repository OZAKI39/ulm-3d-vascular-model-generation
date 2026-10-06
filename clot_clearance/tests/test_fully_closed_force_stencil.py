"""Seven targeted intact-mechanics checks; no historical campaign dependencies."""
import unittest
import numpy as np
from pd_clot.pfpd_candidate.force_stencil_audit import (
    H, DELTA, FORCE_TOL, COEFFICIENT_TOL, RECONSTRUCTION_TOL, AffineStencil, lattice_offsets,
    dependency_offsets, bulk_oracle, validate_selective)


class ClosedStencilTests(unittest.TestCase):
    def test_1_bulk_lattice_offsets_deterministic(self):
        a=lattice_offsets();b=lattice_offsets.__wrapped__()
        np.testing.assert_array_equal(a,b)
        self.assertEqual(len(a),256)
        self.assertFalse(np.any(np.all(a==0,axis=1)))
        self.assertLessEqual(np.linalg.norm(a*H,axis=1).max(),DELTA)
        self.assertEqual(set(map(tuple,a)),set(map(tuple,-a)))

    def test_2_selective_equals_full_N11(self):
        result=validate_selective(11)
        self.assertTrue(result['PASS'],result)

    def test_3_graph_closure_deterministic(self):
        a=dependency_offsets();b=dependency_offsets.__wrapped__()
        self.assertEqual([len(x) for x in a],[1,256,1797,5821])
        for x,y in zip(a,b):np.testing.assert_array_equal(x,y)
        self.assertIn((12,0,0),set(map(tuple,a[3])))
        self.assertNotIn((13,0,0),set(map(tuple,a[3])))
        self.assertAlmostEqual(np.linalg.norm(a[2]*H,axis=1).max(),8*H)

    def test_4_decompositions_reconstruct_force(self):
        grid=AffineStencil(27)
        for depth in (4,9,13):
            pid=grid.id_at_offset((0,0,depth-grid.n//2))
            result=grid.evaluate_point(pid)
            self.assertLessEqual(max(result['decomposition_errors_local_scaled'].values()),RECONSTRUCTION_TOL)

    def test_5_zero_strain_zero_force(self):
        grid=AffineStencil(11,epsilon=0.)
        for offset in ((0,0,0),(0,0,1),(1,-1,-1)):
            result=grid.evaluate_point(grid.id_at_offset(offset))
            np.testing.assert_array_equal(result['force_N'],np.zeros(3))

    def test_6_N27_full_stencil_geometry_and_volume(self):
        grid=AffineStencil(27);center=grid.closure(grid.id_at_offset((0,0,0)))
        self.assertTrue(center['FULL_FORCE_STENCIL_CLOSED'])
        self.assertTrue(center['every_D2_node_has_complete_family'])
        self.assertTrue(all(x['closed_with_bulk_volumes'] for x in center['levels']))
        # Three-hop nodes at the boundary exist but carry half volume: not bulk.
        near=grid.closure(grid.id_at_offset((0,0,-1)))
        self.assertTrue(near['levels'][3]['all_required_nodes_exist'])
        self.assertFalse(near['levels'][3]['all_required_volumes_are_interior'])
        self.assertFalse(near['FULL_FORCE_STENCIL_CLOSED'])

    def test_7_bulk_oracle_symmetry_zero_force(self):
        result=bulk_oracle()
        self.assertLessEqual(result['center_force_local_normalized'],FORCE_TOL)
        # Matrix comparisons use the predeclared scaled norm, including entries
        # whose exact reference is zero; no arbitrary dimensional zero tolerance.
        self.assertLessEqual(np.linalg.norm(result['Z_inf']-(2/3)*result['P'])/np.linalg.norm(result['Z_inf']),COEFFICIENT_TOL)
        K=result['K_inf']
        self.assertLessEqual(np.linalg.norm(K-np.eye(3)*np.trace(K)/3)/np.linalg.norm(K),COEFFICIENT_TOL)


if __name__=='__main__':unittest.main()
