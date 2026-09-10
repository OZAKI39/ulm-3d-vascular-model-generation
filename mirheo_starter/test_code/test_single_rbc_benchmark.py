"""CPU-only tests for geometry, conventions, authority and evidence handling."""
import ast
import copy
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from py_scripts.single_rbc_benchmark.physics import *
from py_scripts.single_rbc_benchmark.analysis import periodic_poiseuille_basis,affine_response,end_to_end_cost
from py_scripts.single_rbc_benchmark.workflow import run,child_env
from py_scripts.single_rbc_benchmark.quality import points_inside,mesh_intersections,strict_comparison

CONFIG=PROJECT_ROOT/'py_scripts/single_rbc_benchmark.yaml'

class CommonCaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=load_config(CONFIG);cls.mesh=paths(cls.c)[0]/'common_reference.off';cls.v,cls.f=read_off(cls.mesh)

    def test_complete_unit_embedding(self):
        u=units();self.assertAlmostEqual(u['M0']*u['L0']**2/u['T0']**2,u['E0'],delta=u['E0']*1e-12)
        self.assertAlmostEqual(8*u['rho0'],1000);self.assertAlmostEqual(u['nu0']*u['rho0'],u['mu0'])

    def test_groups_have_common_time(self):
        g=groups(8,1.5,.02,3,20,8,642,81,129)
        self.assertAlmostEqual(g['Re'],.12);self.assertAlmostEqual(g['Ca'],.036);self.assertAlmostEqual(g['B'],8/180)
        self.assertEqual(g['lambda_viscosity'],1)

    def test_wall_gap_and_periods(self):
        l=lattice(self.c,1.5);self.assertEqual((l['nz']-1)*l['dx'],24);self.assertEqual(l['nx']*l['dx'],24)
        self.assertAlmostEqual(2*(.02*24/2)/24,.02)
        source=(Path(self.c['hemocell']['root'])/'cases/single_rbc_shear_benchmark/benchmark.cpp').read_text()
        self.assertIn('plb::Array<T,3>(-v,0.,0.)',source);self.assertIn('plb::Array<T,3>(v,0.,0.)',source)

    def test_native_single_cell_and_deformable_paths(self):
        cpp=(Path(self.c['hemocell']['root'])/'cases/single_rbc_shear_benchmark/benchmark.cpp').read_text()
        py=(PROJECT_ROOT/'py_scripts/single_rbc_benchmark/mirheo_worker.py').read_text()
        self.assertIn('hc.iterate()',cpp);self.assertIn('RbcHighOrderModel',cpp);self.assertIn('number_of_cells!=1',cpp)
        self.assertIn('MembraneForces',py);self.assertIn('applyObjectBelongingChecker',py);self.assertNotIn('Integrators.Shear',py)
        self.assertIn("u.setInteraction(interaction,inner,fluid)",py);self.assertIn("u.setBouncer(bounce,rbc,inner)",py)

    def test_real_mesh_closure_and_translation(self):
        g=geometry(self.v,self.f);h=geometry(self.v+[234,-13,912],self.f)
        self.assertEqual(g['vertices'],642);self.assertEqual(g['euler_characteristic'],2);self.assertTrue(g['closed'])
        for key in ('area','volume','D'):self.assertAlmostEqual(g[key],h[key],places=9)
        self.assertFalse(geometry(self.v,self.f[:-1])['closed'])

    def test_deformation_definition_is_covariance_not_native_log(self):
        v=np.array([[2,0,0],[-2,0,0],[0,0,1],[0,0,-1]],float)
        f=np.array([[0,1,2],[0,3,1],[0,2,3],[1,3,2]])
        self.assertAlmostEqual(geometry(v,f)['D'],1/3)
        w=v.copy();w[:,0]/=2;g=geometry(w,f);self.assertTrue(g['angle_degenerate']);self.assertIsNone(g['theta_deg'])

    def test_periodic_mesh_unwrapping_uses_edges(self):
        v=self.v+[.25,12,12];wrapped=v.copy();wrapped[:,:2]%=24
        result=unwrap(wrapped,self.f,[24,24,24]);g=geometry(result,self.f);h=geometry(v,self.f)
        self.assertAlmostEqual(g['volume'],h['volume'],places=8);self.assertAlmostEqual(g['area'],h['area'],places=8)

    def test_stable_particle_ids_and_missing_detection(self):
        rng=np.random.default_rng(8);order=rng.permutation(len(self.v));out=order_vertices(order+200,self.v[order],len(self.v));np.testing.assert_array_equal(out,self.v)
        with self.assertRaises(ValueError):order_vertices([0,0],self.v[:2],2)
        with self.assertRaises(ValueError):order_vertices([0,2],self.v[:2],2)

    def test_strain_samples_and_stricter_steps(self):
        for dt in (.001,.0005,lattice(self.c,1.5)['dt'],lattice(self.c,1.5,True)['dt']):
            s=sample_steps(self.c,dt);self.assertAlmostEqual(s['steps']*dt*.02,4);self.assertAlmostEqual(s['sample_steps']*dt*.02,.1);self.assertAlmostEqual(s['prep_steps']*dt,5)

    def test_poiseuille_is_independent_viscosity_test(self):
        y=np.array([0.,2.,4.,6.,8.]);b=periodic_poiseuille_basis(y,8,.1)
        np.testing.assert_allclose(b,[0,-.2,0,.2,0],atol=1e-15)

    def test_quality_gate_never_invents_speedup(self):
        for key in ('workflow','model','screen'):
            args=dict(workflow='PASS',model='PASS',screen='PASS',repeats=2);args[key]='PARTIAL';self.assertIsNone(speedup(20,10,**args))
        self.assertIsNone(speedup(20,10,workflow='PASS',model='PASS',screen='PASS',repeats=1))
        self.assertEqual(speedup(20,10,workflow='PASS',model='PASS',screen='PASS',repeats=2),2)

    def test_total_cost_counts_setup_relaxation_and_analysis_once(self):
        phase=dict(setup_s=2,wall_preparation_s=1,relaxation_s=3,coupled_s=4,output_s=1)
        c=end_to_end_cost(12,phase,.5,.25)
        self.assertEqual(c['total_s'],12.75);self.assertEqual(c['unassigned_process_overhead_s'],2)
        with self.assertRaises(ValueError):end_to_end_cost(8,phase,.5,.25)

    def test_no_gpu_without_new_authorization(self):
        with tempfile.TemporaryDirectory() as tmp:
            with patch('py_scripts.single_rbc_benchmark.workflow.paths',return_value=(Path(tmp),Path(tmp),Path(tmp))),patch('py_scripts.single_rbc_benchmark.workflow.run_attempt') as launch:
                with self.assertRaisesRegex(RuntimeError,'NEW_GPU_AUTHORIZATION_REQUIRED'):run(self.c,'forbidden',{},lambda d:[],10,gpu=True)
                launch.assert_not_called()

    def test_cli_help_does_not_import_gpu(self):
        for module in ('py_scripts.benchmark_single_rbc','py_scripts.single_rbc_benchmark.mirheo_worker','test_code.review_single_rbc_benchmark'):
            r=subprocess.run([str(PROJECT_ROOT/'.venv/bin/python'),'-B','-m',module,'--help'],cwd=PROJECT_ROOT,text=True,capture_output=True,timeout=10)
            self.assertEqual(r.returncode,0,r.stderr);self.assertIn('usage:',r.stdout.lower())

    def test_review_and_preflight_do_not_launch_solver(self):
        from py_scripts.single_rbc_benchmark.workflow import preflight
        with patch('py_scripts.single_rbc_benchmark.workflow.run_attempt') as launch:
            observation=preflight(self.c);launch.assert_not_called()
            if (paths(self.c)[0]/'fluid_measured.json').exists():self.assertEqual(observation['fluid']['status'],'MEASURED')
            if (paths(self.c)[0]/'native_failure_final.json').exists():self.assertEqual(observation['stage'],'COMPUTATION_STOPPED_FINAL_EVIDENCE')
        source=(PROJECT_ROOT/'test_code/review_single_rbc_benchmark.py').read_text();self.assertNotIn('execute(',source);self.assertNotIn('gpu_task',source)

    def test_cache_return_is_not_new_cold_execution(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            with patch('py_scripts.single_rbc_benchmark.workflow.run_attempt',return_value=(d,{'status':'COMPLETED'},True)),patch('py_scripts.single_rbc_benchmark.workflow.write_json') as write,patch('py_scripts.single_rbc_benchmark.workflow.seal') as seal:
                _,_,cached=run(self.c,'cached',{},lambda d:[],10)
                self.assertTrue(cached);write.assert_not_called();seal.assert_not_called()

    def test_old_files_preserved(self):
        r=json.loads((paths(self.c)[0]/'protection_current.json').read_text());self.assertEqual(r['status'],'PASS');self.assertGreater(r['files_checked'],1000)

    def test_threads_explicit(self):
        for name in ('OMP_NUM_THREADS=1','OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1'):self.assertIn(name,child_env(True))

    def test_palabos_archive_identity_is_not_parent_git_head(self):
        from py_scripts.single_rbc_benchmark.workflow import environment
        e=environment(self.c)
        self.assertEqual(e['palabos']['commit'],'05712164d940a42e06afdd705249912fa0c49f14')
        self.assertEqual(e['palabos']['archive_sha256'],'5a9c4f22c169ad16de259a72ef2f128b233afaf81c251e6c3baa8eb34d14838f')
        self.assertNotIn(str(Path(self.c['hemocell']['root'])/'vendor/HemoCell/palabos'),e['git'])

    def test_cpu_probe_is_native_measured(self):
        p=paths(self.c)[1]/'solver/candidate1_native_response/native_response.csv'
        r=affine_response(p,geometry(self.v,self.f)['area']);self.assertGreater(r['Gs_affine_effective'],0);self.assertEqual(len(r['rows']),14)

    def test_no_old_pool_used(self):
        source=(PROJECT_ROOT/'py_scripts/single_rbc_benchmark/workflow.py').read_text();tree=ast.parse(source)
        for node in ast.walk(tree):
            if isinstance(node,ast.Dict):self.assertNotIn('shared_budget_pool',[x.value for x in node.keys if isinstance(x,ast.Constant)])

    def test_solid_angle_handles_concave_closed_rbc(self):
        self.assertTrue(points_inside(np.array([[0.,0.,0.]]),self.v,self.f)[0])
        self.assertFalse(points_inside(np.array([[10.,0.,0.],[0.,0.,5.]]),self.v,self.f).any())

    def test_triangle_intersections_detect_crossing_and_coplanar(self):
        faces=np.array([[0,1,2],[3,4,5]])
        crossing=np.array([[0,0,0],[2,0,0],[0,2,0],[.5,.5,-1],[.5,.5,1],[1,.5,0]],float)
        self.assertEqual(mesh_intersections(crossing,faces)['nonadjacent_intersections'],1)
        crossing[3:]=[[.2,.2,0],[.8,.2,0],[.2,.8,0]]
        self.assertEqual(mesh_intersections(crossing,faces)['nonadjacent_intersections'],1)
        crossing[3:]+=10
        self.assertEqual(mesh_intersections(crossing,faces)['nonadjacent_intersections'],0)
        self.assertEqual(mesh_intersections(self.v,self.f)['nonadjacent_intersections'],0)

    def test_strict_screen_uses_window_means_and_rejects_incomplete_interval(self):
        def fake(offset,end=2.):
            return dict(frames=[dict(phase='shear',strain=float(g),metrics=dict(D=.4+offset,theta_deg=0.,center=[12.,12.,12.])) for g in np.linspace(0,end,21)])
        self.assertEqual(strict_comparison(fake(0),fake(.01),self.c,[.5,2.])['status'],'PASS')
        self.assertEqual(strict_comparison(fake(0),fake(.1),self.c,[.5,2.])['status'],'FAILED')
        self.assertEqual(strict_comparison(fake(0),fake(0,1),self.c,[.5,2.])['status'],'INCOMPLETE_WINDOW')
        self.assertEqual(strict_comparison(fake(0,1),fake(0),self.c,[.5,2.])['status'],'INCOMPLETE_WINDOW')

    def test_phase_handoff_preserves_two_snapshots_at_same_step(self):
        from py_scripts.single_rbc_benchmark.analysis import read_csv
        a=read_csv(paths(self.c)[1]/'solver/coupling_dpd_check/vertices.csv')
        r=a[a['step']==5000]
        self.assertEqual(len(r),1284)
        self.assertEqual(set(r['phase']),{'relaxation','shear'})
        for phase in ('relaxation','shear'):
            self.assertEqual(len(np.unique(r[r['phase']==phase]['vertex'])),642)

    def test_early_failure_with_empty_csv_retains_failure_and_cost(self):
        from py_scripts.single_rbc_benchmark.analysis import _analyze_record
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'execution.json').write_text(json.dumps(dict(status='STOPPED_OR_FAILED',elapsed_monotonic_s=.2)))
            (p/'timings.csv').touch();(p/'local_flow.csv').touch()
            r=_analyze_record(p,self.c)
            self.assertEqual(r['execution']['status'],'STOPPED_OR_FAILED');self.assertEqual(r['frames'],[])
            self.assertIsNone(r['observed_last_saved_strain']);self.assertGreaterEqual(r['costs']['total_s'],.2)

if __name__=='__main__':unittest.main()
