"""CPU-only tests. Synthetic data and fake MPI/CUDA are not physical evidence."""
import copy
import ctypes
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT, read_json, write_json, sha256_file
from py_scripts.fluid_physics.runner import shared_budget_state, run_attempt
from py_scripts.fluid_physics.budget_authorizations import account_extensions, authorization_records
from py_scripts.sdpd_diagnostics.equilibration import (load_config, source_bundle, build_plan,
    exact_steps, budget_status, request_budget, execute, read_raw, verify_manifest, export, register_authorization)
from py_scripts.sdpd_diagnostics.equilibration_analysis import analyze_series, completed, empty_result
from py_scripts.sdpd_diagnostics.equilibration_worker import EquilibrationObserver, mechanical_pressure
from py_scripts.sdpd_diagnostics.gpu_worker import continuous_chunks

CONFIG=PROJECT_ROOT/'py_scripts/sdpd_equilibration.yaml'


def successful(plan):
    return ({'status':'COMPLETED','timeout':False,'exit_code':0,'remaining_own_group_processes':[]},
            {'status':'COMPLETED_PLANNED_STEPS','steps':plan['steps'],'time_star':plan['duration_star']})


def synthetic(plan, temperature=1.):
    # Deliberately synthetic alternating signals; do not write into runs/data.
    steps=np.arange(0,plan['steps']+1,plan['parameters']['snapshot_every'])
    t=steps*plan['dt_star'];count=len(t);wave=np.where(np.arange(count)%2,1.,-1.)
    result={'step':steps,'time_star':t,'N':np.full(count,4096.),'total_mass_star':np.full(count,4096.),
            'density_time_star':t-plan['dt_star'],'kBT_COM_star':np.full(count,temperature)+.0002*wave,
            'pressure_star':100000.+.01*wave,'kernel_number_mean':8.+.0001*wave,
            'kernel_number_variance':.05+.00001*wave}
    for k,v in [('kernel_q10',7.8),('kernel_q50',8.),('kernel_q90',8.2)]:result[k]=v+.0001*wave
    for axis in 'xyz':result['mean_v'+axis]=1e-6*wave
    return result


class SourceAndPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=load_config(CONFIG);cls.source=source_bundle(cls.c);cls.plan=build_plan(cls.c,cls.source)

    def test_delivery_not_chosen_by_modification_time(self):
        self.assertEqual(self.c['diagnostic_package'],str(PROJECT_ROOT/'data/sdpd_diagnostics/result_4ad869c67164f241'))
        c=copy.deepcopy(self.c);c['delivery_record_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'DELIVERY_RECORD_CHANGED'):source_bundle(c)

    def test_exact_original_physics_and_units(self):
        before=self.source['spec'];after=self.plan['parameters']
        for key in ['candidate','dt_star','m_star','kBT_star','rc_star','domain_star','locked_units','initial_velocity_seed','precision']:
            self.assertEqual(before[key],after[key],key)
        self.assertEqual(after['candidate']['viscosity_mu_star'],3312.4601269222444)
        self.assertEqual(after['target_temperature_K'],298.15)

    def test_zero_force_native_integrator(self):
        s=self.plan['parameters']
        self.assertEqual(s['task']['force_star'],0)
        self.assertEqual(s['native_integrator'],'VelocityVerlet')
        self.assertFalse(self.plan['temperature_excursion_is_stop_condition'])

    def test_plan_and_nested_task_agree(self):
        p=self.plan;s=p['parameters']
        self.assertEqual(s['steps'],s['desired_steps'])
        self.assertEqual(s['steps'],s['task']['desired_steps'])
        self.assertEqual(s['desired_time_star'],s['task']['desired_time_star'])
        self.assertAlmostEqual(s['steps']*s['dt_star'],.4)
        self.assertEqual(p['formal_window_star'],[.3,.4])
        self.assertEqual(p['planned_formal_samples'],500)
        self.assertIn('historical_template',s)

    def test_main_cannot_change_timestep(self):
        c=copy.deepcopy(self.c);c['plan']['dt_star']/=2
        with self.assertRaisesRegex(ValueError,'FROZEN_PROPOSAL'):build_plan(c,self.source)

    def test_main_cannot_shrink_or_move_window(self):
        for key,value in [('steps',80000),('formal_window_star',[.35,.4]),('duration_star',.08)]:
            c=copy.deepcopy(self.c);c['plan'][key]=value
            with self.assertRaises(ValueError):build_plan(c,self.source)

    def test_different_dt_comparisons_use_same_physical_time(self):
        a,b=exact_steps(.4,1e-6),exact_steps(.4,.5e-6)
        self.assertEqual(b,2*a);self.assertEqual(a*1e-6,b*.5e-6)
        with self.assertRaisesRegex(ValueError,'GRID'):exact_steps(.4,.003)

    def test_cost_uses_actual_unforced_run(self):
        p=self.plan
        self.assertEqual(p['cost_evidence']['actual_steps'],172000)
        self.assertAlmostEqual(p['predicted_total_wall_s'],582.1243803313985)
        self.assertLess(p['predicted_total_wall_s'],588)
        self.assertFalse(p['cost_evidence']['new_observer_cost_measured'])

    def test_all_ledgers_counted_without_new_authorization(self):
        s=budget_status(self.c,self.source);r=request_budget(self.c,self.plan,s)
        self.assertGreaterEqual(len(s['members']),3)
        self.assertEqual(self.c['extra_authorized_gpu_seconds'],0)
        self.assertGreaterEqual(s['total_charged_or_reserved_s'],3492.9193288880488)
        self.assertAlmostEqual(r['additional_required_s'],max(0.,600-s['remaining_s']))
        self.assertEqual(s['original_limit_s'],3600)

    def test_execute_refuses_insufficient_budget_before_runner(self):
        request={'additional_required_s':492.92,'current_scope_remaining_s':107.08}
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT/'test_code/outputs') as tmp:
            d=Path(tmp);write_json(d/'run_plan.json',self.plan);write_json(d/'budget_request.json',request)
            # This rejection case must remain a new attempt after the real run exists.
            c=copy.deepcopy(self.c);c['runs_root']=str(d/'isolated_runs')
            with patch('py_scripts.sdpd_diagnostics.equilibration.export',return_value=d),patch('py_scripts.sdpd_diagnostics.equilibration.cpu_validation',return_value={'status':'PASS'}),patch('py_scripts.sdpd_diagnostics.equilibration.run_attempt') as run:
                with self.assertRaisesRegex(PermissionError,'NO_GPU_STARTED'):execute(c)
                run.assert_not_called()

    def test_yaml_extra_seconds_cannot_authorize(self):
        import yaml
        c=yaml.safe_load(CONFIG.read_text());c['extra_authorized_gpu_seconds']=493
        with tempfile.NamedTemporaryFile('w',suffix='.yaml') as f:
            yaml.safe_dump(c,f);f.flush()
            with self.assertRaisesRegex(ValueError,'CANNOT_GRANT'):load_config(f.name)

    def test_help_and_missing_mode_do_not_import_cuda(self):
        for args,code in [(['--help'],0),([],2)]:
            script="import sys,runpy;sys.argv=['run_sdpd_equilibration']+"+repr(args)+"\ntry:runpy.run_module('py_scripts.run_sdpd_equilibration',run_name='__main__')\nfinally:assert 'mirheo' not in sys.modules and 'libmirheo' not in sys.modules\n"
            p=subprocess.run([sys.executable,'-B','-c',script],cwd=PROJECT_ROOT,capture_output=True,text=True)
            self.assertEqual(p.returncode,code,p.stderr)


class StationarityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        c=load_config(CONFIG);cls.plan=build_plan(c,source_bundle(c))

    def analyze(self,data=None,**kwargs):
        execution,completion=successful(self.plan)
        return analyze_series(synthetic(self.plan) if data is None else data,self.plan,
                              kwargs.get('execution',execution),kwargs.get('completion',completion),
                              kwargs.get('structure',{'measurement_status':'PASS','status':'PASS'}))

    def test_equilibrium_and_thermal_screen_pass_are_not_all_properties(self):
        r=self.analyze();s=r['summary']
        self.assertEqual(s['status'],'EQUILIBRIUM_SCREEN_PASS',s)
        self.assertTrue(s['liquid_stationary']);self.assertTrue(s['stable_temperature_matches_target'])
        self.assertFalse(s['all_liquid_properties_qualified']);self.assertIsNone(s['selection'])
        self.assertIsNone(s['viscosity_new_measurement']);self.assertIsNone(s['EOS_new_measurement'])

    def test_high_stationary_temperature_is_biased(self):
        r=self.analyze(synthetic(self.plan,1.1))
        self.assertEqual(r['summary']['status'],'STATIONARY_BUT_THERMALLY_BIASED')
        self.assertFalse(r['summary']['stable_temperature_matches_target'])

    def test_wide_thermal_ci_is_not_forced_pass_or_bias(self):
        # Borderline stationary plateau, just outside target, with uncertainty overlapping threshold.
        data=synthetic(self.plan,1.02001);rng=np.random.default_rng(812)
        data['kBT_COM_star']+=rng.normal(0,.001,len(data['time_star']))
        r=self.analyze(data)
        self.assertEqual(r['summary']['status'],'STATIONARITY_OR_SAMPLING_INCONCLUSIVE')
        self.assertIsNone(r['summary']['stable_temperature_matches_target'])

    def test_cooling_cannot_pass_and_does_not_claim_stationary_acf(self):
        data=synthetic(self.plan);data['kBT_COM_star']=1.5-data['time_star']
        r=self.analyze(data)
        self.assertEqual(r['summary']['status'],'TRANSIENT_PERSISTS')
        self.assertIsNone(r['summary']['temperature_mean_star'])
        self.assertIsNone(r['sampling']['measured_stationary_ACF_samples'])

    def test_small_coherent_quarter_changes_are_json_serializable(self):
        # A monotonic point trend can remain below both material limits.
        # NumPy's final comparison used to leak np.bool_(False) into JSON.
        data=synthetic(self.plan);t=data['time_star']
        for i,(lo,hi) in enumerate(zip(np.linspace(.3,.4,5)[:-1],np.linspace(.3,.4,5)[1:])):
            mask=(t>lo+1e-12)&(t<=hi+1e-12)
            data['kBT_COM_star'][mask]=1.+.0002*i
        result=self.analyze(data)
        encoded=json.dumps(result,allow_nan=False)
        restored=json.loads(encoded)
        self.assertIs(restored['stationarity']['metrics']['kBT_COM_star']['coherent_drift'],False)
        self.assertIsNone(restored['summary']['stable_temperature_matches_target'])

    def test_pressure_cannot_hide_behind_background(self):
        data=synthetic(self.plan);data['pressure_star']=1e9+1000*data['time_star']
        r=self.analyze(data)
        self.assertEqual(r['summary']['status'],'TRANSIENT_PERSISTS')
        self.assertEqual(r['stationarity']['metrics']['pressure_star']['normalization_scale_star'],8.)

    def test_large_noise_monotonic_quarter_means_are_inconclusive(self):
        data=synthetic(self.plan);t=data['time_star']
        for i,(lo,hi) in enumerate(zip(np.linspace(.3,.4,5)[:-1],np.linspace(.3,.4,5)[1:])):
            mask=(t>lo+1e-12)&(t<=hi+1e-12);count=mask.sum()
            data['pressure_star'][mask]=1e8+100*i+1e4*np.sin(np.arange(count)*2*np.pi/count)
        r=self.analyze(data)
        self.assertEqual(r['summary']['status'],'STATIONARITY_OR_SAMPLING_INCONCLUSIVE')
        self.assertIsNone(r['sampling']['measured_stationary_ACF_samples'])

    def test_kernel_structure_drift_cannot_pass(self):
        data=synthetic(self.plan);data['kernel_number_mean']=8+10*data['time_star']
        self.assertEqual(self.analyze(data)['summary']['status'],'TRANSIENT_PERSISTS')

    def test_momentum_drift_cannot_pass(self):
        data=synthetic(self.plan);data['mean_vx']=.1*data['time_star']
        self.assertEqual(self.analyze(data)['summary']['status'],'TRANSIENT_PERSISTS')

    def test_stationary_looking_insufficient_blocks_have_no_ci(self):
        with patch('py_scripts.sdpd_diagnostics.equilibration_analysis.correlation_time_samples',return_value=82.44):
            r=self.analyze()
        self.assertEqual(r['summary']['status'],'STATIONARITY_OR_SAMPLING_INCONCLUSIVE')
        self.assertIsNone(r['statistics']['kBT_COM_star']['ci95_halfwidth'])
        self.assertGreater(r['sampling']['additional_formal_samples_condition'],0)

    def test_complete_block_mean_and_ci_samples_agree(self):
        data=synthetic(self.plan)
        data['kBT_COM_star'][-8:]+=.001
        r=self.analyze(data);s=r['statistics']['kBT_COM_star']
        self.assertGreater(s['discarded_tail_samples'],0)
        self.assertAlmostEqual(s['mean'],np.mean(s['block_means']),places=13)
        self.assertNotEqual(s['mean'],s['all_sample_mean'])

    def test_old_acf_is_not_written_as_new_measurement(self):
        r=self.analyze()
        values=r['sampling']['measured_stationary_ACF_samples']
        self.assertLess(values['kBT_COM_star']['full'],2)
        self.assertIn('never copied',r['sampling']['historical_ACF_role'])

    def test_fixed_window_does_not_search_for_best_tail(self):
        data=synthetic(self.plan)
        t=data['time_star'];data['kBT_COM_star'][(t>.3)&(t<=.36)]=1.1
        r=self.analyze(data)
        self.assertNotEqual(r['summary']['status'],'EQUILIBRIUM_SCREEN_PASS')
        self.assertEqual(r['sampling']['formal_window_star'],[.3,.4])

    def test_particle_count_error_is_measurement_error(self):
        data=synthetic(self.plan);data['N'][-1]-=1
        self.assertEqual(self.analyze(data)['summary']['status'],'MEASUREMENT_OR_ANALYSIS_ERROR')

    def test_nonfinite_or_wrong_phase_cannot_pass(self):
        for key,value in [('kBT_COM_star',np.nan),('density_time_star',0.)]:
            data=synthetic(self.plan);data[key][-1]=value
            self.assertEqual(self.analyze(data)['summary']['status'],'MEASUREMENT_OR_ANALYSIS_ERROR')

    def test_incomplete_zero_exit_is_partial(self):
        e,k=successful(self.plan);k['steps']=100000;k['time_star']=.1;k['status']='COOPERATIVE_PARTIAL'
        r=self.analyze(execution=e,completion=k)
        self.assertFalse(completed(e,k,self.plan))
        self.assertEqual(r['summary']['execution_status'],'PARTIAL')
        self.assertIsNone(r['summary']['liquid_stationary'])

    def test_timeout_with_full_steps_is_not_completed_experiment(self):
        e,k=successful(self.plan);e['timeout']=True
        r=self.analyze(execution=e,completion=k)
        self.assertEqual(r['summary']['execution_status'],'PARTIAL')

    def test_runtime_failure_is_distinct(self):
        e,k=successful(self.plan);e.update(exit_code=1,status='STOPPED_OR_FAILED');k={}
        self.assertEqual(self.analyze(execution=e,completion=k)['summary']['status'],'NUMERICAL_OR_RUNTIME_FAILURE')

    def test_snapshot_error_not_treated_as_physics(self):
        r=self.analyze(structure={'measurement_status':'FAIL','status':'FAIL'})
        self.assertEqual(r['summary']['status'],'MEASUREMENT_OR_ANALYSIS_ERROR')

    def test_initial_peak_retained_but_not_an_emergency_stop(self):
        data=synthetic(self.plan);data['kBT_COM_star'][1:100]=9.
        r=self.analyze(data)
        self.assertEqual(r['summary']['status'],'EQUILIBRIUM_SCREEN_PASS')
        self.assertEqual(np.max(data['kBT_COM_star']),9.)

    def test_no_run_measurements_are_null(self):
        r=empty_result()
        for k in ['actual_steps','temperature_mean_star','temperature_CI95_star','liquid_stationary','stable_temperature_matches_target']:
            self.assertIsNone(r['summary'][k])
        self.assertIsNone(r['sampling']['measured_stationary_ACF_samples'])


class BudgetAndContinuationTests(unittest.TestCase):
    def test_real_registration_code_only_appends_to_synthetic_pool(self):
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT/'test_code/outputs') as tmp:
            d=Path(tmp);pool=d/'pool.json';old=d/'old';old.mkdir()
            write_json(old/'budget_ledger.json',{'campaign_id':'old','attempts':[{'task_id':'old','status':'COMPLETED','reserved_s':3493,'charged_s':3493}]})
            write_json(pool,{'limit_s':3600,'gpu_lock':str(d/'gpu.lock'),'members':[{'directory':str(old),'campaign_id':'old','ledger_required':True}]})
            c={'runs_root':str(d/'new_runs'),'campaign_id':'synthetic','case_id':'only'}
            scope={'campaign_directory':str(d/'new_runs/synthetic'),'task_id':'only','plan_sha256':'SYNTHETIC_NOT_PHYSICAL'}
            write_json(d/'budget_request.json',{'authorization_scope':scope})
            message=d/'synthetic_user_message.txt';message.write_text('SYNTHETIC CPU TEST ONLY, NOT REAL GPU APPROVAL: 493 seconds for the synthetic scope.')
            confirmation={'authority':'explicit_user_message','user_message_file':str(message),'additional_seconds':493,
                          'scope':scope,'budget_request_sha256':sha256_file(d/'budget_request.json')}
            write_json(d/'confirmation.json',confirmation)
            ledger_before=sha256_file(old/'budget_ledger.json')
            with patch('py_scripts.sdpd_diagnostics.equilibration.export',return_value=d),patch('py_scripts.sdpd_diagnostics.equilibration.source_bundle',return_value={'comparison_config':{'shared_budget_pool':str(pool)}}):
                path=register_authorization(c,d/'confirmation.json')
                with self.assertRaisesRegex(ValueError,'ALREADY_REGISTERED'):register_authorization(c,d/'confirmation.json')
            current=read_json(pool)
            self.assertEqual(current['limit_s'],3600);self.assertEqual(len(current['authorization_records']),1)
            self.assertEqual(authorization_records(current)[0]['additional_seconds'],493)
            self.assertEqual(sha256_file(old/'budget_ledger.json'),ledger_before)
            self.assertTrue(path.with_name('shared_pool_before.json').is_file())
            message.write_text('tampered synthetic evidence')
            with self.assertRaisesRegex(ValueError,'AUTHORIZATION_EVIDENCE'):authorization_records(current)

    def test_registration_cannot_accept_self_approval_flag(self):
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT/'test_code/outputs') as tmp:
            p=Path(tmp)/'confirmation.json';write_json(p,{'approved':True,'additional_seconds':493})
            with self.assertRaisesRegex(ValueError,'EXPLICIT_USER_MESSAGE_REQUIRED'):register_authorization({},p)

    def test_append_accounting_preserves_old_charges_and_scope(self):
        pool={'limit_s':3600};campaign='/synthetic/new'
        record={'scope':{'campaign_directory':campaign,'task_id':'only'},'additional_seconds':493}
        attempts=[('/synthetic/old','old',3492.9193288880488),(campaign,'only',600)]
        with patch('py_scripts.fluid_physics.budget_authorizations.authorization_records',return_value=[record]):
            new=account_extensions(pool,attempts,campaign,'only')
            other=account_extensions(pool,attempts,'/synthetic/other','unapproved')
        self.assertEqual(new['extra_authorized_gpu_seconds'],493)
        self.assertAlmostEqual(new['original_charged_or_reserved_s'],3599.9193288880488)
        self.assertAlmostEqual(other['remaining_s'],.08067111195123)
        self.assertEqual(sum(x[2] for x in attempts),4092.9193288880488)

    def test_unapproved_scope_cannot_spend_extension(self):
        record={'scope':{'campaign_directory':'/synthetic/new','task_id':'only'},'additional_seconds':493}
        with patch('py_scripts.fluid_physics.budget_authorizations.authorization_records',return_value=[record]):
            s=account_extensions({'limit_s':3600},[('/old','old',3492.9)],'/synthetic/new','wrong')
        self.assertAlmostEqual(s['remaining_s'],107.1)

    def test_no_record_means_zero_extra_not_yaml_amount(self):
        s=account_extensions({'limit_s':3600,'extra_authorized_gpu_seconds':999},[('/old','old',3492.9)],'/new','only')
        self.assertEqual(s['extra_authorized_gpu_seconds'],0)

    def test_missing_or_changed_authorization_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'record.json';p.write_text('{}')
            with self.assertRaisesRegex(ValueError,'HASH_MISMATCH'):
                authorization_records({'authorization_records':[{'path':str(p),'sha256':'0'*64}]})

    def test_strict_runner_never_clamps_full_plan(self):
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT/'test_code/outputs') as tmp:
            d=Path(tmp);old=d/'old';new=d/'new';old.mkdir()
            write_json(old/'budget_ledger.json',{'campaign_id':'old','attempts':[{'task_id':'old','status':'COMPLETED','reserved_s':3493,'charged_s':3493}]})
            pool=d/'pool.json';write_json(pool,{'limit_s':3600,'gpu_lock':str(d/'gpu.lock'),'members':[
                {'directory':str(old),'campaign_id':'old','ledger_required':True},
                {'directory':str(new),'campaign_id':'new','ledger_required':False}]})
            c={'campaign_id':'new','require_shared_budget':True,'shared_budget_pool':str(pool),
               'budget':{'campaign_limit_s':3600,'task_limit_s':600,'minimum_launch_budget_s':30}}
            before=sha256_file(pool)
            with patch('py_scripts.fluid_physics.runner.bounded_process') as proc:
                with self.assertRaisesRegex(RuntimeError,'NO_SHORTENED_RUN'):
                    run_attempt(new,c,'full',{},lambda d:[],600,strict_allocation=True,monitor_gpu=False)
                proc.assert_not_called()
            self.assertEqual(sha256_file(pool),before)
            self.assertFalse((new/'budget_ledger.json').exists())

    def test_one_coordinator_rank_matching_and_no_reinitialization(self):
        class Coordinator:
            def __init__(self):self.calls=[];self.time=0
            def run(self,n,dt):self.calls.append((n,dt));self.time+=n*dt
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'control').mkdir();u=Coordinator();other=Coordinator()
            result=list(continuous_chunks(u,d,True,100,30,.001,lambda:None))
            result2=list(continuous_chunks(other,d,False,100,30,.001,lambda:None))
            self.assertEqual(u.calls,other.calls)
            self.assertEqual([x[0] for x in u.calls],[30,30,30,10])
            self.assertEqual(result[-1][1],100);self.assertEqual(result2[-1][1],100)
            self.assertAlmostEqual(u.time,.1)

    def test_cooperative_stop_propagates_to_both_ranks(self):
        class Coordinator:
            def __init__(self):self.calls=[]
            def run(self,n,dt):self.calls.append(n)
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp);(d/'control').mkdir();u=Coordinator()
            it=continuous_chunks(u,d,True,100,20,.001,lambda:None)
            next(it);(d/'STOP_REQUESTED').touch();self.assertEqual(list(it),[])
            other=Coordinator();list(continuous_chunks(other,d,False,100,20,.001,lambda:None))
            self.assertEqual(u.calls,other.calls);self.assertEqual(u.calls,[20])

    def test_pressure_reuses_correct_stress_components_no_extra_eos(self):
        vel=np.array([[1.,0,0],[-1.,0,0]])
        stress=np.array([[3.,100,200,6.,300,9.],[3.,-100,-200,6.,-300,9.]])
        p,v=mechanical_pressure(stress,vel,2.,4.,1e-6,200)
        self.assertAlmostEqual(v,3.);self.assertAlmostEqual(p,3.+1/3)
        with self.assertRaisesRegex(ValueError,'STRESS_LAYOUT'):
            mechanical_pressure(stress[:,:3],vel,2.,4.,1e-6,200)

    def test_native_stress_layout_and_get_ids_exist(self):
        source=PROJECT_ROOT/'vendor/Mirheo/src/mirheo'
        text=(source/'core/datatypes.h').read_text().split('struct Stress',1)[1].split('};',1)[0]
        self.assertLess(text.index('real xx'),text.index('real yy'))
        self.assertLess(text.index('real yy'),text.index('real zz'))
        self.assertIn('MIR_NUMPY_PARTIAL_INFO(Stress, 6', (source/'bindings/cuda_array_interface.cpp').read_text())
        self.assertIn('.def("get_indices"', (source/'bindings/particle_vectors.cpp').read_text())
        self.assertIn('(channel_names::stresses, DataManager::PersistenceMode::None)',(source/'core/interactions/pairwise/sdpd.cu').read_text())


class FakeNativeObserverTests(unittest.TestCase):
    def test_production_worker_lifecycle_with_cpu_fake_native_arrays(self):
        # Runs the real Python worker and observer with a labeled fake device API;
        # no real Mirheo import, CUDA initialization or GPU process.
        from py_scripts.sdpd_diagnostics import gpu_worker
        calls={'coordinator':0,'particle_vector':0,'velocity_set':0,'runs':[],'integrator':[],'savers':[]}
        class Device:
            def __init__(self,a):self.a=a
            @property
            def __cuda_array_interface__(self):
                return {**self.a.__array_interface__,'version':3}
        class PV:
            def __init__(self,name,mass):
                calls['particle_vector']+=1;self.mass=mass
                rng=np.random.default_rng(2);self.x=np.zeros((64,4),np.float32);self.x[:,:3]=rng.uniform(-1,1,(64,3))
                self.v=np.zeros((64,4),np.float32);self.x[:,3]=np.arange(64,dtype=np.int32).view(np.float32)
                self.a={'positions':self.x,'velocities':self.v,'saved_positions':self.x.copy(),
                        'saved_velocities':self.v.copy(),'saved_density':np.full(64,8,np.float32),
                        'saved_forces':np.zeros((64,4),np.float32),'stresses':np.ones((64,6),np.float32),
                        'saved_stresses':np.full((64,6),np.nan,np.float32)}
                self.local=types.SimpleNamespace(per_particle={k:Device(v) for k,v in self.a.items()})
            def getCoordinates(self):return self.x[:,:3].astype(float)+4
            def getVelocities(self):return self.v[:,:3].astype(float)
            def setVelocities(self,v):calls['velocity_set']+=1;self.v[:,:3]=v
            def get_indices(self):return list(range(64))
        class Coordinator:
            def __init__(self,*args,**kwargs):calls['coordinator']+=1
            def registerParticleVector(self,pv,ic):self.pv=pv
            def registerInteraction(self,*args):pass
            def setInteraction(self,*args):pass
            def registerIntegrator(self,*args):pass
            def setIntegrator(self,*args):pass
            def registerPlugins(self,*args):
                if args[0] and args[0][0]=='saver':calls['savers'].append(args[0][1:])
            def run(self,n,dt):
                calls['runs'].append((n,dt))
                for _ in range(n):
                    self.pv.a['stresses'][:]=1
                    self.pv.a['saved_positions'][:]=self.pv.x
                    self.pv.a['saved_velocities'][:]=self.pv.v
                    for source,saved in calls['savers']:
                        if source=='stresses':self.pv.a[saved][:]=self.pv.a[source]
                    self.pv.x[:,:3]+=dt*self.pv.v[:,:3]
                # Nonpersistent stress is not contractually preserved by the final cell-list rebuild.
                self.pv.a['stresses'][:]=np.nan
        class Plugins:
            def createParticleChannelSaver(self,name,pv,source,saved):return ('saver',source,saved)
            def __getattr__(self,name):return lambda *args,**kwargs:None
        def integrator(name):calls['integrator'].append('VelocityVerlet');return None
        mir=types.SimpleNamespace(Mirheo=Coordinator,ParticleVectors=types.SimpleNamespace(ParticleVector=PV),
            InitialConditions=types.SimpleNamespace(Uniform=lambda **kw:None),
            Interactions=types.SimpleNamespace(Pairwise=lambda *a,**kw:None),
            Integrators=types.SimpleNamespace(VelocityVerlet=integrator),Plugins=Plugins())
        def memcpy(dst,src,n,kind):ctypes.memmove(dst,src,n);return 0
        runtime=types.SimpleNamespace(cudaMemcpy=memcpy,cudaDeviceSynchronize=lambda:0)
        spec={'task':{'force_star':0,'kind':'equilibrium','n_star':8},'candidate':{'method':'SDPD','density_interaction':True,
            'viscosity_mu_star':3312.4601269222444,'sound_speed_star':120,'rho_0_star':0},
            'domain_star':[8,8,8],'dt_star':1e-6,'m_star':1.,'kBT_star':1.,'rc_star':1.,'initial_velocity_seed':20260908,
            'snapshot_every':2,'y_bin_star':8.,'steps':16,'snapshot_marks_steps':[2,8,16],'expected_N':64,
            'target_temperature_K':298.15,'locked_units':{'t0':3e-5,'si_per_star':{'mass_density':132.,'pressure':.03}}}
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT/'test_code/outputs') as tmp:
            d=Path(tmp);(d/'control').mkdir();write_json(d/'spec.json',spec);cwd=Path.cwd()
            try:
                os.chdir(d)
                with patch.dict(sys.modules,{'mirheo':mir}),patch.dict(os.environ,{'OMPI_COMM_WORLD_SIZE':'2','OMPI_COMM_WORLD_RANK':'0'}),patch.object(gpu_worker.ctypes,'CDLL',return_value=runtime),patch.object(sys,'argv',['cpu-fake','--spec',str(d/'spec.json')]):
                    gpu_worker.main(observer_factory=EquilibrationObserver)
            finally:os.chdir(cwd)
            self.assertEqual(calls['coordinator'],1);self.assertEqual(calls['particle_vector'],1)
            self.assertEqual(calls['velocity_set'],1);self.assertEqual(calls['integrator'],['VelocityVerlet'])
            self.assertEqual(sum(n for n,_ in calls['runs']),16)
            self.assertEqual(read_json(d/'worker_completion.json')['status'],'COMPLETED_PLANNED_STEPS')
            data=read_raw(d);self.assertEqual(len(data['step']),9);self.assertEqual(data['step'][0],0)
            self.assertTrue(np.isnan(data['pressure_star'][0]));self.assertTrue(np.isfinite(data['pressure_star'][1:]).all())
            self.assertEqual(len(list(d.glob('snapshot_*_ids.npy'))),3)
            self.assertFalse(read_json(d/'initial_state_identity.json')['snapshot_is_checkpoint'])
            self.assertEqual(read_json(d/'equilibration_observer.json')['continuous_coordinator_count'],1)


class ReadOnlyReviewTests(unittest.TestCase):
    def test_repeated_preflight_analysis_and_view_never_launch_gpu_or_fill_measurements(self):
        # Isolate output directories so this remains valid after a real authorized run.
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT/'test_code/outputs') as tmp:
            self.check_repeated_view(Path(tmp))

    def check_repeated_view(self, tmp):
        from test_code.review_sdpd_equilibration import write_page
        c=load_config(CONFIG)
        c.update(data_root=str(tmp/'data'),review_root=str(tmp/'review'),runs_root=str(tmp/'runs'))
        pool=PROJECT_ROOT/'runs/fluid_calibration/shared_budget_pool.json'
        paths=[pool]+[Path(m['directory'])/'budget_ledger.json' for m in read_json(pool)['members']]
        paths += [Path(c['delivery_record']),Path(c['diagnostic_package'])/'report_zh.md',
                  PROJECT_ROOT/'metadata/CMakeCache.txt',PROJECT_ROOT/'vendor/Mirheo/src/mirheo/core/utils/cuda_common.h']
        hashes={str(p):sha256_file(p) for p in paths}
        with patch('py_scripts.fluid_physics.runner.bounded_process') as process,patch('py_scripts.sdpd_diagnostics.gpu_worker.main') as worker,patch('py_scripts.sdpd_diagnostics.equilibration.cpu_validation',return_value=None):
            first=export(c);second=export(c,analyze=True)
            self.assertEqual(first,second)
            page=write_page(c,first);self.assertEqual(write_page(c,second),page)
            process.assert_not_called();worker.assert_not_called()
        self.assertEqual(hashes,{str(p):sha256_file(p) for p in paths})
        self.assertEqual(len((first/'raw_statistics.csv').read_text().splitlines()),1)
        summary=read_json(first/'equilibration_summary.json')
        self.assertEqual(summary['execution_status'],'NOT_RUN_AWAITING_AUTHORIZATION')
        self.assertIsNone(summary['actual_steps']);self.assertIsNone(summary['temperature_mean_star'])
        self.assertFalse((first/'actual_parameters.json').exists())
        self.assertFalse((Path(c['runs_root'])/c['campaign_id']/c['case_id']).exists())

    def test_read_only_sources_reject_writes(self):
        from py_scripts.fluid_physics.common import output
        for path in [PROJECT_ROOT.parent/'ulm_3D_vascular/new.json',PROJECT_ROOT.parent/'bloodflow_starter/new.json']:
            with self.assertRaises(Exception):output(path)


if __name__=='__main__':unittest.main()
