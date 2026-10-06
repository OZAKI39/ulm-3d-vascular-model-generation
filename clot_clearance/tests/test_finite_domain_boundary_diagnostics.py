"""Seven direct scientific checks for intact finite-domain diagnostics only."""
import unittest
import numpy as np
from scipy.spatial import cKDTree
from pd_clot.pfpd_candidate.kernel import Family
from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.boundary_force_diagnostics import (
    make_cube,continuum,force_terms,endpoint_weight,force_stats,required_force,
    self_neighbor_terms,full_lattice_omega0)
from scripts.finite_domain_boundary_audit import material,DELTA,EPS


class BoundaryDiagnosticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.material=material();cls.cubes=[make_cube(n,L,DELTA) for n,L in [(5,.00025),(7,.00025),(13,.00075)]]
        cls.states=[evaluate(c.family,c.family.X@np.diag([1,1,1+EPS]),cls.material) for c in cls.cubes]

    def test_1_direct_Z_reconstruct_original_force(self):
        for c,s in zip(self.cubes,self.states):
            t=force_terms(c.family,s)
            self.assertLess(np.linalg.norm(t['direct']+t['correction']-s.force)/np.linalg.norm(s.force),1e-12)

    def test_2_endpoint_and_coefficient_splits_reconstruct_original(self):
        for c,s in zip(self.cubes,self.states):
            f=c.family;i,j=f.pairs.T;norm=np.linalg.norm(s.force)
            a=force_terms(f,s,.5*f.omega/f.omega0[i])['total'];b=force_terms(f,s,.5*f.omega/f.omega0[j])['total']
            self.assertLess(np.linalg.norm(a+b-s.force)/norm,1e-12)
            a=force_terms(f,s,endpoint_weight(f,c.full_support))['total'];b=force_terms(f,s,endpoint_weight(f,~c.full_support))['total']
            self.assertLess(np.linalg.norm(a+b-s.force)/norm,1e-12)
            self.assertLess(np.linalg.norm(sum(self_neighbor_terms(f,s).values())-s.force)/norm,1e-12)

    def test_3_global_affine_energy_identity(self):
        for c,s in zip(self.cubes,self.states):
            f=c.family;i,j=f.pairs.T
            self.assertLess(abs(np.sum(2*f.wbar*f.volume[i]*f.volume[j])/f.volume.sum()-1),2e-15)
            self.assertLess(abs(s.energy/continuum(self.material,c.L,EPS)['energy_J']-1),2e-11)

    def test_4_large_cube_has_geometric_full_support_core(self):
        c=self.cubes[-1]
        self.assertEqual(int(c.full_support.sum()),27)
        self.assertLess(abs(c.family.volume.sum()/c.L**3-1),2e-15)
        np.testing.assert_allclose(c.face_areas,c.L**2,atol=0,rtol=2e-15)
        self.assertAlmostEqual(c.h,62.5e-6,places=15)

    def test_5_bulk_statistics_are_independent_of_particle_order(self):
        c=self.cubes[-1];f=c.family;ref=continuum(self.material,c.L,EPS);perm=np.arange(len(f.X))[::-1]
        X=f.X[perm];g=Family.from_arrays(X,f.volume[perm],cKDTree(X).query_pairs(DELTA,output_type='ndarray'),DELTA)
        s=evaluate(g,X@ref['F'].T,self.material)
        a=force_stats(self.states[-1].force,c.full_support,ref['reaction_N']);b=force_stats(s.force,c.full_support[perm],ref['reaction_N'])
        np.testing.assert_allclose([a['max_N'],a['rms_N']],[b['max_N'],b['rms_N']],rtol=2e-11,atol=1e-22)

    def test_6_manufactured_forcing_balances_affine_coordinates(self):
        for c,s in [(self.cubes[0],self.states[0]),(self.cubes[-1],self.states[-1])]:
            R=continuum(self.material,c.L,EPS)['reaction_N'];free=(c.ijk[:,2]!=0)&(c.ijk[:,2]!=c.n-1)
            self.assertLess(np.linalg.norm((s.force+required_force(s))[free])/R,1e-10)

    def test_7_continuum_units_and_distributed_virtual_work(self):
        self.assertEqual(self.material.lame,0)
        for c,s in zip(self.cubes,self.states):
            r=continuum(self.material,c.L,EPS)
            # Independent scalar SVK formula: Pa*m^2=N, Pa*m^3=J.
            e=EPS+EPS**2/2;P=(1+EPS)*2*1000*e;W=1000*e*e
            self.assertAlmostEqual(r['reaction_N']/(P*c.L**2),1,places=12)
            self.assertAlmostEqual(r['energy_J']/(W*c.L**3),1,places=12)
            self.assertAlmostEqual(-np.sum(s.force[:,2]*(c.family.X[:,2]+c.L/2))/(P*c.L**3),1,places=11)


if __name__=='__main__':unittest.main()
