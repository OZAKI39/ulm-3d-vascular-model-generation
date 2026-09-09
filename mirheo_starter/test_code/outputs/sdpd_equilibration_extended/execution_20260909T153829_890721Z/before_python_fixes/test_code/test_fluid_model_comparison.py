"""Synthetic CPU tests, never used as experimental measurements or GPU charges."""
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT,write_json,read_json
from py_scripts.fluid_physics.units import Units,thermal_density_units
from py_scripts.fluid_physics.analysis import fit_viscosity,periodic_shape,native_pressure,block_statistics
from py_scripts.fluid_physics.reporting import eos_analysis,sensitivity,evidence_group
from py_scripts.fluid_physics.runner import shared_budget_state,run_attempt,used_budget
from py_scripts.fluid_comparison.models import validate_model,linear_eos,snapshot_density,frozen_targets
from py_scripts.fluid_comparison.experiments import load_config
from py_scripts.fluid_comparison.gpu_worker import thermal_bins
from py_scripts.fluid_comparison.reporting import choose


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.c=load_config('py_scripts/fluid_model_comparison.yaml');self.case,self.mapping,self.u=frozen_targets(self.c)
        self.a=copy.deepcopy(self.c['resolved_candidates'][0])
    def test_exact_frozen_units(self):
        self.assertEqual(self.u.t0,self.mapping['t0']);self.assertAlmostEqual(self.u.temperature_K(1),298.15)
    def test_dynamic_not_kinematic_input(self):
        self.assertEqual(self.a['viscosity_mu_star'],self.mapping['required_mu_star'])
        self.assertAlmostEqual(self.a['viscosity_mu_star']/self.mapping['required_nu_star'],8)
    def test_nonunit_mass_density(self):
        u=thermal_density_units(5e-7,1056,4,2,1,298.15)
        self.assertAlmostEqual(u.to_si(4*2,'mass_density'),1056)
        self.assertNotEqual(u.to_si(4,'mass_density'),1056)
    def test_missing_density_fails(self):
        self.a.pop('density_interaction')
        with self.assertRaisesRegex(ValueError,'REQUIRES_DENSITY'):validate_model(self.a)
    def test_matching_kernel_and_rc(self):
        self.a['density_interaction']['rc_star']=.5
        with self.assertRaisesRegex(ValueError,'MISMATCH'):validate_model(self.a)
    def test_no_double_thermostat(self):
        self.a['also_apply_dpd']=True
        with self.assertRaisesRegex(ValueError,'DUPLICATED'):validate_model(self.a)
    def test_force_time_and_viscosity_conversion(self):
        for quantity in ['force','time','dynamic_viscosity','kinematic_viscosity']:
            self.assertAlmostEqual(self.u.to_star(self.u.to_si(2.5,quantity),quantity),2.5)
    def test_density_self_and_periodicity(self):
        self.assertAlmostEqual(snapshot_density([[0,0,0]],[8,8,8],1)[0],21/(2*np.pi))
        a=snapshot_density([[.1,0,0],[7.9,0,0]],[8,8,8],1)
        b=snapshot_density([[3.1,0,0],[2.9,0,0]],[8,8,8],1)
        np.testing.assert_allclose(a,b);self.assertGreater(a[0],21/(2*np.pi))
    def test_EOS_mass_density_and_negative_gauge(self):
        self.assertAlmostEqual(float(linear_eos(4*2,120,0)),115200)
        p0=100.;negative=self.case['outlet_gauge_pressures_pa']['outlet_03']
        self.assertLess(p0+self.u.to_star(negative,'pressure'),p0)
    def test_measured_viscosity_not_copied_from_input(self):
        y=np.arange(.125,8,.25);truth=2.3;v=.4/(2*truth)*periodic_shape(y,8,.25)+.2
        result=fit_viscosity(y,v,.4,2,8,.25)
        self.assertAlmostEqual(result['nu_star'],truth);self.assertNotEqual(result['nu_star'],self.a['viscosity_mu_star'])
    def test_local_thermal_flow_subtraction(self):
        rng=np.random.default_rng(3);pos=rng.uniform(0,8,(20000,3));v=rng.normal(size=(20000,3));v[:,0]+=2*np.floor(pos[:,1])
        k=thermal_bins(pos,v,8,8,1)[0]
        self.assertLess(abs(k-1),.02);self.assertGreater(np.mean(v*v),2)
    def test_insufficient_statistics(self):
        self.assertNotEqual(block_statistics([0,1,2],[1,1,1],min_duration=1)['status'],'SUFFICIENT')
    def test_help_has_no_mirheo_import(self):
        for module in ['py_scripts.compare_dpd_sdpd','py_scripts.fluid_comparison.gpu_worker','test_code.review_fluid_model_comparison']:
            code="import runpy,sys;sys.modules['mirheo']=None;sys.argv=['x','--help'];runpy.run_module("+repr(module)+",run_name='__main__')"
            r=subprocess.run([sys.executable,'-B','-c',code],cwd=PROJECT_ROOT,capture_output=True,text=True)
            self.assertEqual(r.returncode,0,r.stderr)


class EvidenceTests(unittest.TestCase):
    def rows(self,method='SDPD',candidate='one',group='g'):
        return [{'method':method,'candidate_id':candidate,'comparison_group':group,'task_id':str(n),'kind':'equilibrium',
                 'test_kind':role,'actual_n_star':n,'parameters':{'m_star':2},
                 'pressure':{'mean':20*n,'ci95_halfwidth':.1,'status':'SUFFICIENT'},'temperature_status':'PASS_PROPOSED'}
                for n,role in [(3.8,'eos_low'),(4,'equilibrium'),(4.2,'eos_high')]]
    def test_eos_never_borrows_another_candidate(self):
        with self.assertRaisesRegex(ValueError,'MIXED'):eos_analysis(self.rows()+self.rows(candidate='two'),{},Units(1,1,1),{})
        r=eos_analysis(self.rows()+self.rows(candidate='two'),{'outlet_gauge_pressures_pa':{'x':.1}},Units(1,1,1),{},candidate_id='one')
        self.assertEqual(len(r['points']),3)
    def test_methods_and_groups_are_independent(self):
        for more in [self.rows(method='DPD'),self.rows(group='other')]:
            with self.assertRaisesRegex(ValueError,'MIXED'):evidence_group(self.rows()+more)
    def test_reference_density_not_hardcoded(self):
        r=eos_analysis(self.rows(),{'density_kg_m3':8,'outlet_gauge_pressures_pa':{'p':1}},Units(1,1,1),{})
        self.assertAlmostEqual(r['bulk_modulus_pa'],8*10)
        rows=self.rows()
        for t in rows:t['parameters']['m_star']=1
        r=eos_analysis(rows,{'density_kg_m3':4,'outlet_gauge_pressures_pa':{'p':1}},Units(1,1,1),{})
        self.assertAlmostEqual(r['bulk_modulus_pa'],4*20)
    def test_large_errorbars_do_not_extend_pressure_range(self):
        rows=self.rows()
        for r in rows:r['pressure']['ci95_halfwidth']=1000
        r=eos_analysis(rows,{'outlet_gauge_pressures_pa':{'p':100}},Units(1,1,1),{})
        self.assertEqual(r['target_pressure_coverage_status'],'OUTSIDE_MEASURED_EOS_RANGE')
    def test_nominal_coverage_is_not_confident_coverage(self):
        rows=self.rows()
        for row in rows:row['pressure']['ci95_halfwidth']=2
        r=eos_analysis(rows,{'outlet_gauge_pressures_pa':{'p':3}},Units(1,1,1),{})
        self.assertEqual(r['target_pressure_nominal_coverage_status'],'COVERED')
        self.assertEqual(r['target_pressure_coverage_status'],'NOMINAL_COVERAGE_UNCERTAIN')
    def test_half_test_mixed_candidates_fail(self):
        rows=[{'task_id':'flow','candidate_id':'a'},{'task_id':'flow_half_dt','candidate_id':'b'}]
        with self.assertRaisesRegex(ValueError,'MIXED'):sensitivity(rows,{'sensitivity_relative_difference':.1})
    def test_native_time_offset_must_match_this_run(self):
        s=np.array([(0.001,1,0,0,0,100)],dtype=[(n,float) for n in ['time','kBT','vx','vy','vz','num_particles']])
        v=np.array([(0,30)],dtype=[('time',float),('pressure',float)])
        self.assertAlmostEqual(native_pressure(s,v,1,10,.001)[0][0],13)
        with self.assertRaisesRegex(ValueError,'TIMESTAMP'):native_pressure(s,v,1,10,.002)
    def test_six_digit_rounding_at_half_step(self):
        s=np.array([(5e-7,1,0,0,0,100),(.100001,1,0,0,0,100)],dtype=[(n,float) for n in ['time','kBT','vx','vy','vz','num_particles']])
        v=np.array([(0,30),(.1,30)],dtype=[('time',float),('pressure',float)])
        np.testing.assert_allclose(native_pressure(s,v,1,10,5e-7)[0],[13,13])
        with self.assertRaisesRegex(ValueError,'TIMESTAMP'):native_pressure(s,v,1,10,1e-6)
    def test_failed_candidates_cannot_be_selected(self):
        self.assertIsNone(choose([{'status':'NOT_QUALIFIED','method':'SDPD','candidate_id':'close_to_target','score':100}]))
        self.assertIsNone(choose([{'status':'COST_ONLY','method':'DPD','candidate_id':'fast'}]))
    def test_half_test_wrong_units_cannot_pass(self):
        base={'method':'SDPD','candidate_id':'one','comparison_group':'g','kind':'flow','task_id':'a','test_kind':'flow',
              'parameters':{'locked_units':{'L0':1},'m_star':1,'kBT_star':1,'rc_star':1,'domain_star':[8]*3,'candidate':{'id':'one'},
                            'dt_star':.01,'task':{'n_star':8,'force_star':1}},
              'viscosity':{'nu_star':2,'ci95_star':[1.99,2.01],'sampling_status':'SUFFICIENT'}}
        other=copy.deepcopy(base);other['task_id']='b';other['test_kind']='half_dt';other['parameters']['dt_star']=.005;other['parameters']['locked_units']={'L0':2}
        result=sensitivity([base,other],{'sensitivity_relative_difference':.1})
        self.assertNotEqual(result['status'],'PASS_PROPOSED')
        self.assertFalse(next(x for x in result['comparisons'] if x['task_id']=='half_dt')['locked_same_units'])


class SharedBudgetTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(prefix='synthetic_comparison_',dir=PROJECT_ROOT/'test_code');self.d=Path(self.tmp.name)
        self.old=self.d/'old';self.new=self.d/'new';self.old.mkdir();self.new.mkdir();self.pool=self.d/'pool.json'
        self.b={'task_limit_s':2,'campaign_limit_s':3,'minimum_launch_budget_s':.05,'graceful_margin_s':.1}
        self.c={'campaign_id':'new','budget':self.b,'shared_budget_pool':str(self.pool),'require_shared_budget':True}
        write_json(self.old/'budget_ledger.json',{'campaign_id':'old','campaign_limit_s':3,'attempts':[{'status':'COMPLETED','reserved_s':2,'charged_s':2.5}]})
        write_json(self.pool,{'limit_s':3,'gpu_lock':str(self.d/'gpu.lock'),'members':[
            {'directory':str(self.old),'campaign_id':'old','ledger_required':True},{'directory':str(self.new),'campaign_id':'new','ledger_required':False}]})
    def tearDown(self):self.tmp.cleanup()
    def test_cross_campaign_total_and_cache(self):
        old_bytes=(self.old/'budget_ledger.json').read_bytes()
        state=shared_budget_state(self.new,self.c);self.assertEqual(state['remaining_s'],.5)
        cmd=lambda d:[sys.executable,'-c','pass']
        d,r,hit=run_attempt(self.new,self.c,'test',{'a':1},cmd,2,monitor_gpu=False)
        self.assertLessEqual(r['allocation_s'],.5);self.assertFalse(hit)
        used=shared_budget_state(self.new,self.c)['total_charged_or_reserved_s']
        self.assertGreater(used,2.5)
        _,_,hit=run_attempt(self.new,self.c,'test',{'a':1},cmd,2,monitor_gpu=False)
        self.assertTrue(hit);self.assertEqual(shared_budget_state(self.new,self.c)['total_charged_or_reserved_s'],used)
        self.assertEqual((self.old/'budget_ledger.json').read_bytes(),old_bytes)
    def test_missing_old_ledger_fails_closed(self):
        (self.old/'budget_ledger.json').rename(self.old/'moved.json')
        with self.assertRaisesRegex(ValueError,'LEDGER_MISSING'):shared_budget_state(self.new,self.c)
    def test_unfinished_reservation_cannot_launch(self):
        p=self.old/'budget_ledger.json';l=read_json(p);l['attempts'][0]['status']='RUNNING';p.write_text(json.dumps(l))
        with self.assertRaisesRegex(RuntimeError,'SHARED_RUNNING'):run_attempt(self.new,self.c,'a',{},lambda d:[],1,monitor_gpu=False)
    def test_unregistered_campaign_cannot_reset_allowance(self):
        with self.assertRaisesRegex(ValueError,'NOT_REGISTERED'):shared_budget_state(self.d/'third',self.c)
    def test_timeout_is_charged(self):
        _,r,_=run_attempt(self.new,self.c,'timeout',{},lambda d:[sys.executable,'-c','import time;time.sleep(5)'],.4,monitor_gpu=False)
        self.assertTrue(r['timeout']);self.assertFalse(r['remaining_own_group_processes'])
        self.assertAlmostEqual(shared_budget_state(self.new,self.c)['total_charged_or_reserved_s'],2.5+r['elapsed_monotonic_s'])


if __name__=='__main__':unittest.main()
