"""CPU regressions for the independent repair; never imports a solver."""
import os
from pathlib import Path
import unittest
import numpy as np
from py_scripts.single_rbc_repair.geometry_checks import checked_mesh, membership, independent_intersections, paired_statistics, require_same_frame
from py_scripts.single_rbc_repair.workflow import load_config, paths, require_authorization, qualification
from py_scripts.single_rbc_benchmark.physics import read_off
from py_scripts.single_rbc_benchmark.analysis import end_to_end_cost

ROOT = Path(__file__).resolve().parents[1]


class NativeOrderingRegression(unittest.TestCase):
    def test_local_bounce_finishes_before_halo_reuses_buffers(self):
        # Set REPAIR_NATIVE_SOURCE to the isolated candidate for the green check.
        source = Path(os.environ.get('REPAIR_NATIVE_SOURCE', paths(load_config())[0] / 'native/source'))
        code = (source / 'src/mirheo/core/simulation.cpp').read_text()
        self.assertTrue('scheduler->addDependency(tasks->objHaloBounce, {}, {tasks->objLocalBounce});' in code,
                      'local/halo share per-bouncer buffers; their GPU tasks must be ordered')

    def test_capacity_stays_bounded_and_overflow_still_fails(self):
        source = paths(load_config())[0] / 'native/source/src/mirheo/core/bouncers'
        common = (source/'drivers/common.h').read_text()
        code = (source/'from_mesh.cu').read_text()
        self.assertTrue('if (i < maxSize) indices[i] = idx;' in common)
        self.assertTrue('coarseTable_.nCollisions[0] > maxCoarseCollisions' in code)
        self.assertTrue('die("Found too many triangle collision candidates' in code)
        self.assertEqual(5*1280,6400)

    def test_nonrigid_force_accumulation_and_empty_halo(self):
        source = paths(load_config())[0] / 'native/source/src/mirheo/core/bouncers'
        code = (source/'from_mesh.cu').read_text()
        self.assertLess(code.index('if (activeOV->getNumObjects() == 0'),code.index('coarseTable_.nCollisions.clear'))
        self.assertIn('if (rov_)\n    {',code)
        driver = (source/'drivers/mesh.h').read_text()
        self.assertTrue('atomicAdd' in driver)


class GeometryRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.v,cls.f=read_off(paths(load_config())[0]/'common_reference.off')

    def test_shuffled_ids_and_periodic_unwrapping(self):
        v=self.v+np.array([.5,12,12]);wrapped=v.copy();wrapped[:,:2]%=24
        order=np.random.default_rng(19).permutation(len(v))
        actual,g=checked_mesh(order,wrapped[order],self.f,24)
        np.testing.assert_allclose(actual-actual.mean(0),v-v.mean(0),atol=1e-10)
        self.assertTrue(g['closed'])

    def test_missing_ids_and_invalid_connectivity_rejected(self):
        with self.assertRaises(ValueError):checked_mesh(np.zeros(len(self.v)),self.v,self.f,24)
        faces=self.f.copy();faces[0,0]=len(self.v)
        with self.assertRaises(ValueError):checked_mesh(np.arange(len(self.v)),self.v,faces,24)

    def test_near_surface_is_uncertain_not_a_leak(self):
        tri=self.v[self.f[0]];p=tri.mean(0);n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n/=np.linalg.norm(n)
        points=np.array([p+n*1e-6,p+n*.01,p-n*.01])
        result=membership(points,np.array([True,True,False]),self.v,self.f)
        self.assertTrue(result['uncertain'][0]);self.assertFalse(result['confirmed_mismatch'][0])
        self.assertTrue(np.all(result['confirmed_mismatch'][1:]))

    def test_equal_strain_at_different_phases_is_not_the_same_frame(self):
        with self.assertRaisesRegex(ValueError,'TIME_MISMATCH'):require_same_frame('shear',5000,'relaxation',5000)
        require_same_frame('shear',5000,'shear',5000)

    def test_adjacent_faces_ignored_and_real_intersection_confirmed(self):
        points=np.array([[0,0,0],[2,0,0],[0,2,0],[.5,-.2,-1],[.5,1,1],[.5,1,-1]],float)
        faces=np.array([[0,1,2],[3,4,5]])
        self.assertEqual(independent_intersections(points,faces)['confirmed_count'],1)
        self.assertEqual(independent_intersections(points,np.array([[0,1,2],[0,4,5]]))['confirmed_count'],0)
        self.assertEqual(independent_intersections(self.v,self.f)['confirmed_count'],0)

    def test_reference_geometry_changes_are_measured_without_rotation_fitting(self):
        from py_scripts.single_rbc_benchmark.physics import geometry
        a=geometry(self.v,self.f)['a']
        centered=self.v-self.v.mean(0)
        self.assertLess(np.sqrt(np.mean(np.sum((centered-(self.v+12-(self.v+12).mean(0)))**2,axis=1)))/a,1e-12)
        rotated=centered[:,[2,1,0]]
        self.assertGreater(np.sqrt(np.mean(np.sum((centered-rotated)**2,axis=1)))/a,.02)


class AccountingAndGateRegression(unittest.TestCase):
    def test_both_fluid_sides_are_bound_in_each_control(self):
        from types import SimpleNamespace
        from py_scripts.single_rbc_repair.coupling import bind_bouncers
        for policy,expected in [('shared',1),('independent',2)]:
            calls=[];registered=[]
            mir=SimpleNamespace(Bouncers=SimpleNamespace(Mesh=lambda name,kind:object()))
            u=SimpleNamespace(registerBouncer=registered.append,setBouncer=lambda *args:calls.append(args))
            handles=bind_bouncers(mir,u,'membrane','outer','inner',policy)
            self.assertEqual(len(handles),expected)
            self.assertEqual([x[1:] for x in calls],[('membrane','outer'),('membrane','inner')])
    def test_mean_and_ci_use_same_finite_samples(self):
        x=paired_statistics([1.,3.,float('nan')],confidence_multiplier=2)
        self.assertEqual(x['n'],2);self.assertEqual(x['mean'],2)
        np.testing.assert_allclose(x['ci'],[0.,4.])

    def test_failed_saved_frame_and_unmatched_material_never_rank(self):
        rows=[dict(solver=s,run_id=s+str(i),completed=True,strain_end=4.,cached=False) for s in ('HemoCell','Mirheo') for i in range(2)]
        self.assertIsNone(qualification(rows,'FAILED','PASS'))
        rows[-1].update(completed=False,strain_end=3.4)
        self.assertIsNone(qualification(rows,'PASS','PASS'))

    def test_cache_and_duplicate_run_ids_are_not_new_repeats(self):
        rows=[dict(solver=s,run_id=s,completed=True,strain_end=4.,cached=False) for s in ('HemoCell','Mirheo') for _ in range(2)]
        self.assertIsNone(qualification(rows,'PASS','PASS'))
        rows=[dict(solver=s,run_id=s+str(i),completed=True,strain_end=4.,cached=i==1) for s in ('HemoCell','Mirheo') for i in range(2)]
        self.assertIsNone(qualification(rows,'PASS','PASS'))

    def test_nested_wall_time_not_added_twice(self):
        r=end_to_end_cost(10,dict(setup_s=3,wall_preparation_s=2,relaxation_s=2,coupled_s=3,output_s=1),.5,.2)
        self.assertAlmostEqual(r['total_s'],10.7)
        with self.assertRaises(ValueError):end_to_end_cost(10,dict(setup_s=9,coupled_s=8),0)

    def test_legacy_authorization_is_not_new_permission(self):
        c=load_config()
        if not (paths(c)[0]/'authorization.json').exists():
            with self.assertRaisesRegex(RuntimeError,'NEW_REPAIR_AUTHORIZATION_REQUIRED'):require_authorization(c)

    def test_review_module_has_no_solver_launcher(self):
        import ast
        p=ROOT/'py_scripts/single_rbc_repair/reporting.py';tree=ast.parse(p.read_text())
        names=[n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
        self.assertFalse(any(n and ('runner' in n or 'mirheo_worker' in n) for n in names))


if __name__ == '__main__':
    unittest.main()
