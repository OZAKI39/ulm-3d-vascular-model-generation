"""CPU tests: synthetic fixtures are code checks, never solver measurements."""
import json
import math
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
from py_scripts.solver_benchmark.physics import load_config,definition,lattice_definition,target,profile_metrics,qualified_ratio,valid_blocks,validate_config
from py_scripts.solver_benchmark.workflow import plans,xml_config,cpu_env,gpu_budget,verify_run
from py_scripts.fluid_physics.common import PROJECT_ROOT,sha256_file
from py_scripts.fluid_physics.analysis import periodic_shape


class PhysicsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=load_config('py_scripts/solver_benchmark.yaml');cls.p,cls.u=definition(cls.c)
    def test_identical_SI_materials(self):
        self.assertEqual(self.p['rho_si'],1056);self.assertEqual(self.p['nu_si'],3.27e-6)
        self.assertEqual(self.p['mu_si'],1056*3.27e-6)
        self.assertAlmostEqual(self.u.to_si(self.p['candidate']['viscosity_mu_star'],'dynamic_viscosity'),self.p['mu_si'])
    def test_mapping_not_applied_twice(self):
        self.assertEqual(self.u.t0,self.p['mapping']['t0'])
        self.assertEqual(self.p['box_si'],[4e-6]*3)
        self.assertAlmostEqual(self.u.temperature_K(self.p['mapping']['kBT_star']),298.15)
    def test_force_mass_acceleration_density(self):
        force=self.u.to_si(self.p['particle_force_star'],'force')
        mass=self.u.to_si(self.p['mapping']['m_star'],'particle_mass')
        self.assertAlmostEqual(force/mass,self.p['acceleration_si'])
        self.assertEqual(self.p['force_density_si'],1056*50000)
    def test_lbm_tau_relation(self):
        for dx in self.c['lbm']['dx_candidates_si']:
            l=lattice_definition(self.p,dx)
            self.assertAlmostEqual((l['tau']-.5)*dx*dx/(3*l['dt_si']),self.p['nu_si'])
            self.assertAlmostEqual(l['acceleration_lbm']*dx/l['dt_si']**2,self.p['acceleration_si'])
    def test_time_does_not_mean_equal_steps(self):
        _,_,gpu,tasks=plans(self.c)
        for t in tasks:
            if t['role']=='main':
                self.assertNotEqual(t['steps'],gpu['steps'])
                self.assertLess(abs(t['steps']*t['dt_si']-gpu['steps']*gpu['dt_si']),t['dt_si'])
    def test_unique_periodic_nodes_and_sign(self):
        for N in (16,32):
            y=(np.arange(N)+.5)*self.p['box_si'][1]/N
            f=np.where(y>self.p['box_si'][1]/2,1,-1)
            self.assertEqual(f.sum(),0);self.assertEqual(len(np.unique(y)),N)
            self.assertLess(y[-1],self.p['box_si'][1])
    def test_analytic_bin_integral(self):
        y,v=target(self.p);w=self.p['bin_width_si']
        for yy,value in zip(y,v):
            pts=yy+np.linspace(-w/2,w/2,1001)
            sample=self.p['acceleration_si']/self.p['nu_si']*periodic_shape(pts,self.p['box_si'][1])
            self.assertLess(abs(np.trapz(sample,pts)/w-value),2e-10)
    def test_half_flow_not_full_cancellation(self):
        _,v=target(self.p);r=profile_metrics(self.p,v)
        self.assertAlmostEqual(r['full_signed_flow_si'],0)
        self.assertLess(r['half_flow_si'][0],0);self.assertGreater(r['half_flow_si'][1],0)
        self.assertLess(r['half_flow_relative_error'],1e-12)
    def test_independent_fit_recovers_changed_viscosity(self):
        _,v=target(self.p);r=profile_metrics(self.p,v/1.2+.0005)
        self.assertAlmostEqual(r['apparent_nu_si']/self.p['nu_si'],1.2)
        self.assertAlmostEqual(r['fitted_offset_si'],.0005)
        self.assertGreater(r['profile_relative_l2'],.05)
    def test_no_qualified_ratio_for_failure(self):
        self.assertIsNone(qualified_ratio({'qualified':False,'time_to_qualified_solution_s':300},{'qualified':True,'time_to_qualified_solution_s':1}))
        self.assertIsNone(qualified_ratio({'qualified':True},{'qualified':True}))
    def test_invalid_CI_not_published(self):
        t=np.arange(30)*1e-9
        r=valid_blocks(t,np.arange(30),self.p,self.c)
        self.assertIsNone(r['ci95_halfwidth']);self.assertEqual(r['status'],'WINDOW_INSUFFICIENT')
    def test_thermal_background_is_not_LBM_thermometry(self):
        self.assertIn('NOT_APPLICABLE',self.p['background_temperature_scope'])
    def test_xml_no_cell_loading(self):
        _,_,_,tasks=plans(self.c);xml=xml_config(tasks[0],self.p)
        self.assertNotIn('RBC',xml);self.assertNotIn('PLT',xml)
        self.assertIn('<rhoP>1056.0</rhoP>',xml)
    def test_child_environment_isolated(self):
        self.assertIn('OMP_NUM_THREADS=1',cpu_env());self.assertIn('LD_LIBRARY_PATH',cpu_env())
        self.assertNotIn('oversubscribe',' '.join(cpu_env()))
    def test_old_scope_denied(self):
        _,_,g,_=plans(self.c);b=gpu_budget(self.c,g)
        if b['new_benchmark_currently_authorized_s']==0:self.assertEqual(b['status'],'PENDING_SCOPE_APPROVAL')
    def test_source_uses_bulk_reduction_and_velocity_correction(self):
        cpp=Path(self.c['hemocell_root'])/'cases/pure_fluid_benchmark/benchmark.cpp'
        t=cpp.read_text();self.assertIn('getUniqueBulk',t);self.assertIn('MPI_SUM',t);self.assertIn('MPI_MAX',t)
        self.assertIn('cell.computeVelocity(u)',t);self.assertIn('momentTemplates<T,DESCRIPTOR>::get_j',t)
        self.assertNotIn('hemocell.iterate()',t);self.assertNotIn('loadParticles',t)
    def test_help_no_solver(self):
        for module in ['py_scripts.benchmark_mirheo_hemocell','py_scripts.solver_benchmark.mirheo_worker','test_code.review_solver_benchmark']:
            p=subprocess.run([str(PROJECT_ROOT/'.venv/bin/python'),'-B','-m',module,'--help'],capture_output=True,text=True,timeout=10,cwd=PROJECT_ROOT)
            self.assertEqual(p.returncode,0,p.stderr);self.assertIn('usage:',p.stdout)
    def test_budget_and_campaign_cannot_be_reset(self):
        for key,value in [('cpu_limit_s',1201),('concurrency',2),('automatic_retry',True)]:
            c=json.loads(json.dumps(self.c));c['budget'][key]=value
            with self.assertRaises(ValueError):validate_config(c)
        c=json.loads(json.dumps(self.c));c['campaign_id']='new_free_pool'
        with self.assertRaises(ValueError):validate_config(c)
    def test_statistics_time_keys_are_SI(self):
        out=valid_blocks(np.arange(1000)*1e-8,np.random.default_rng(5).normal(size=1000),self.p,self.c)
        self.assertIn('sample_interval_si',out);self.assertNotIn('sample_interval_star',out)
    def test_cache_tampering_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            d=Path(folder);(d/'data.csv').write_text('original')
            (d/'output_sha256.json').write_text(json.dumps({'data.csv':sha256_file(d/'data.csv')}));verify_run(d)
            (d/'data.csv').write_text('changed')
            with self.assertRaises(ValueError):verify_run(d)
    def test_descendants_in_other_process_groups_are_sampled(self):
        import os,signal
        from py_scripts.solver_benchmark.workflow import host_sample
        code="import subprocess,sys,time; p=subprocess.Popen([sys.executable,'-c','import time; time.sleep(5)'],process_group=0); print(p.pid,flush=True); time.sleep(5)"
        p=subprocess.Popen([str(PROJECT_ROOT/'.venv/bin/python'),'-B','-c',code],stdout=subprocess.PIPE,text=True,start_new_session=True)
        child=int(p.stdout.readline())
        try:
            sample=host_sample(p.pid)
            self.assertIn(child,[x['pid'] for x in sample['processes']])
            self.assertGreater(sample['tree_rss_bytes'],0)
        finally:
            os.kill(child,signal.SIGTERM);p.terminate();p.wait(timeout=5);p.stdout.close()


class ActualEvidenceTests(unittest.TestCase):
    """Readback of sealed real runs, separate from synthetic formula tests."""
    @classmethod
    def setUpClass(cls):
        from py_scripts.solver_benchmark.workflow import freeze
        from py_scripts.solver_benchmark.analysis import analyze_all
        cls.c=load_config('py_scripts/solver_benchmark.yaml');f,_=freeze(cls.c)
        cls.data=analyze_all(cls.c,f);cls.rows=cls.data['results']
    def test_smoke_is_real_zero_cell_nonzero_steps(self):
        row=next(r for r in self.rows if r['role']=='smoke')
        self.assertEqual(row['actual_steps'],64);self.assertEqual(row['completion']['cell_count'],0)
        self.assertFalse(row['completion']['cellfield_created']);self.assertFalse(row['qualified'])
    def test_mpi_unique_bulk_profiles_agree(self):
        rows=[r for r in self.rows if r['candidate_id']=='lbm_N16' and r['role']=='main']
        first=rows[0]['accuracy']['measured_profile_si']
        self.assertEqual({r['rank_count'] for r in rows},{1,2,4})
        for r in rows:np.testing.assert_allclose(first,r['accuracy']['measured_profile_si'],rtol=1e-10,atol=1e-13)
    def test_forced_velocity_halfstep_observed(self):
        for r in self.rows:
            if r['backend']=='HemoCell':self.assertLess(r['velocity_half_force_check_max_si'],1e-12)
    def test_refinement_measured_not_zero_total_error(self):
        self.assertGreater(self.data['refinement_relative_l2'],0)
        for r in self.rows:
            if r['backend']=='HemoCell' and r['role']=='main':
                self.assertGreater(r['accuracy']['profile_relative_l2'],0)
                self.assertIsNone(r['accuracy']['profile_ci95_si'])
    def test_memory_devices_not_confused(self):
        for r in self.rows:
            if r['backend']=='HemoCell':
                self.assertIsNone(r['memory']['device_sampled_peak_used_MiB'])
                host=r['memory']['host']
                if host.get('coverage')=='DESCENDANT_TREE':
                    self.assertGreater(host['sampled_peak_tree_rss_bytes'],0)
                    self.assertEqual(len(host['ranks']),r['rank_count'])
                else:self.assertIsNone(host['sampled_peak_tree_rss_bytes'])
    def test_compile_time_excluded_from_solver(self):
        for r in self.rows:
            if r.get('completion'):
                self.assertAlmostEqual(r['timing']['compute_s'],r['completion']['compute_s'],places=8)
                self.assertLessEqual(r['timing']['compute_s'],r['timing']['total_wall_s'])
    def test_bad_material_never_qualified_speedup(self):
        self.assertTrue(all(r['qualified_speedup'] is None for r in self.data['qualified_speedups']))
    def test_failed_startup_is_zero_steps_and_charged(self):
        r=next(r for r in self.rows if r['task_id']=='sdpd_main')
        self.assertEqual(r['actual_steps'],0);self.assertGreater(r['timing']['total_wall_s'],2)
        self.assertIsNone(r['accuracy']);self.assertFalse(r['qualified'])
    def test_real_native_MPI_transport(self):
        f=PROJECT_ROOT/'data/solver_benchmark/solver_benchmark_20260909/native_mpi_CPU_check.json'
        x=json.loads(f.read_text());self.assertEqual(x['status'],'PASS');self.assertEqual(x['ranks'],2)
        self.assertEqual(x['sum'],3);self.assertGreater(x['maximum_rank_elapsed_s'],.02)
    def test_reused_old_modules_unchanged(self):
        old=json.loads(Path('/home/lzy/projects/hemocell_starter/metadata/protection_before.json').read_text())['hashes']
        for name in ['units.py','analysis.py','runner.py','budget_authorizations.py']:
            p=PROJECT_ROOT/'py_scripts/fluid_physics'/name;self.assertEqual(sha256_file(p),old[str(p)])
    def test_HTML_matches_result_JSON(self):
        latest=json.loads((PROJECT_ROOT/'data/solver_benchmark/LATEST.json').read_text())
        text=Path(latest['html']).read_text();embedded=text.split('<script id="benchmark-audit-data" type="application/json">',1)[1].split('</script>',1)[0]
        self.assertEqual(json.loads(embedded),json.loads(Path(latest['results']).read_text()))
        self.assertNotIn('<script src=',text);self.assertIn('PENDING',text)
    def test_review_is_read_only(self):
        ledger=PROJECT_ROOT/'runs/solver_benchmark/solver_benchmark_20260909/gpu/budget_ledger.json'
        before=sha256_file(ledger)
        p=subprocess.run([str(PROJECT_ROOT/'.venv/bin/python'),'-B','-m','test_code.review_solver_benchmark'],capture_output=True,text=True,timeout=10,cwd=PROJECT_ROOT)
        self.assertEqual(p.returncode,0,p.stderr);self.assertEqual(before,sha256_file(ledger))


if __name__=='__main__':unittest.main()
