"""Only the new displacement-component adapter and its physical work pairing."""
import unittest
import numpy as np
from pd_clot.pfpd_candidate.activation_boundary_control import BoundaryProblem,FULL,AXIAL
from scripts.c7_early_branch import problem


class BoundaryControlTests(unittest.TestCase):
    def test_baseline_matches_existing_mapping_and_force(self):
        old=problem(3,.003);p=BoundaryProblem(3,.003,FULL)
        x=p.affine();s=p.lower
        np.testing.assert_array_equal(p.pack(x,s),old.pack(x,s))
        np.testing.assert_array_equal(p.unpack(p.pack(x,s))[0],old.unpack(old.pack(x,s))[0])
        np.testing.assert_array_equal(p.state(x,s).force,old.state(x,s).force)

    def test_axial_gauge_uses_same_geometric_rule_and_removes_six_rigid_modes(self):
        for n in (5,7):
            p=BoundaryProblem(n,.003,AXIAL)
            self.assertEqual(p.gauge_rank(),6)
            self.assertEqual(int(p.fixed_dof[:,:2].sum()),3)
            self.assertEqual(int(p.fixed_dof[:,2].sum()),2*n*n)
            center=p.gauge[0][0];reference=p.gauge[-1][0]
            np.testing.assert_allclose(p.family.X[center],[0,0,-.000125],atol=1e-15)
            np.testing.assert_allclose(p.family.X[reference],[.000125,0,-.000125],atol=1e-15)

    def test_component_roundtrip_preserves_free_lateral_coordinates(self):
        p=BoundaryProblem(3,.003,AXIAL);z=p.pack(p.affine(),p.lower)
        z[:p.nx]+=np.random.default_rng(19).normal(size=p.nx)*1e-3
        x,s=p.unpack(z)
        np.testing.assert_allclose(p.pack(x,s),z,atol=3e-14,rtol=1e-12)
        np.testing.assert_array_equal(x[p.fixed_dof],p.boundary[p.fixed_dof])

    def test_zero_load_force_energy_and_J(self):
        p=BoundaryProblem(3,0,AXIAL);st=p.state(p.affine(),p.lower)
        self.assertLess(np.linalg.norm(st.force),1e-20)
        self.assertLess(abs(st.energy),1e-30)
        np.testing.assert_allclose(np.linalg.det(st.Ftilde),1,atol=1e-14)

    def test_free_component_energy_gradient_and_load_reaction(self):
        p=BoundaryProblem(3,.003,AXIAL);x=p.affine();s=p.lower
        z=p.pack(x,s);rng=np.random.default_rng(82);v=np.r_[rng.normal(size=p.nx),np.zeros(p.ns)];v/=np.linalg.norm(v)
        # At the C3 affine state the free gradient is zero by symmetry. Use a
        # deterministic non-affine state to make a relative derivative test meaningful.
        z=z+.03*v;x,s=p.unpack(z)
        h=3e-6;analytic=p.objective_gradient(z)[1]@v
        fd=(p.objective_gradient(z+h*v)[0]-p.objective_gradient(z-h*v)[0])/(2*h)
        self.assertLess(abs(fd-analytic)/max(abs(analytic),1e-8),5e-6)
        # Hold free components fixed while changing prescribed top uz: dU/deps=Rtop*L.
        de=1e-6;pplus=BoundaryProblem(3,.003+de,AXIAL);pminus=BoundaryProblem(3,.003-de,AXIAL)
        ep=pplus.state(pplus.prescribe(x),s).energy;em=pminus.state(pminus.prescribe(x),s).energy
        fdload=(ep-em)/(2*de);work=-p.state(x,s).force[p.top,2].sum()*.00025
        self.assertLess(abs(fdload/work-1),5e-6)


if __name__=='__main__':
    unittest.main()
