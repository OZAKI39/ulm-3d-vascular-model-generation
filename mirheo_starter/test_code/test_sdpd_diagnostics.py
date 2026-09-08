"""CPU diagnostics; synthetic cases are never material-measurement evidence."""
import unittest
import math
import json
from pathlib import Path
import tempfile
import sys
import subprocess
from unittest.mock import patch
import numpy as np
from py_scripts.fluid_physics.analysis import block_statistics,fit_viscosity,periodic_shape,correlation_time_samples
from py_scripts.fluid_physics.units import thermal_density_units
from py_scripts.fluid_physics.common import PROJECT_ROOT,write_json,read_json,sha256_file
from py_scripts.fluid_physics.calibration import verify_package
from py_scripts.fluid_physics.runner import used_budget,run_attempt,ledger_read
from py_scripts.sdpd_diagnostics.calculations import (cuda_layout,mass_field_semantics,recover_profile_moments,
    profile_arrays,shared_reference_covariance,matched_physical_interval,ou_discrete_temperature_ratio,average_profile,window_audit)
from py_scripts.sdpd_diagnostics.observations import phase_audit,velocity_statistics
from py_scripts.sdpd_diagnostics.probes import completion_state


class StatisticsRepairTests(unittest.TestCase):
    def test_eight_block_budget_respects_actual_sampling_grid(self):
        from py_scripts.sdpd_diagnostics.workflow import load_config,plan_experiments
        from py_scripts.fluid_physics.units import Units
        c=load_config('py_scripts/sdpd_diagnostics.yaml');p=Path(c['comparison_result'])
        rows=read_json(p/'comparison_summary.json')['tasks'];mapping=read_json(p/'unit_mapping.json')
        u=Units(mapping['L0'],mapping['M0'],mapping['t0']);plan=plan_experiments(c,rows,{'remaining_s':0},u,{'block_relaxation_multiples':2})
        flow=next(r for r in rows if r['task_id']=='sdpd_flow');spec=flow['parameters'];spacing=spec['snapshot_every']*spec['dt_star']
        for r in plan['eight_block_cost_scenarios']:
            required=math.ceil(r['block_star']/spacing-1e-9)*8
            available=math.floor((r['duration_star']-r['warmup_star'])/spacing+1e-9)
            self.assertGreaterEqual(available,required)

    def test_mean_and_ci_use_same_complete_blocks(self):
        values=np.r_[np.tile([-1.,1.],50),[50.,50.,50.,50.]]
        result=block_statistics(np.arange(len(values)),values,min_duration=10,min_blocks=8)
        self.assertGreaterEqual(result['block_count'],2)
        self.assertGreater(result['discarded_tail_samples'],0)
        self.assertAlmostEqual(result['mean'],float(np.mean(result['block_means'])),places=12)

    def test_one_block_never_has_CI(self):
        r=block_statistics(np.arange(12.),np.tile([-1.,1.],6),min_duration=10)
        self.assertEqual(r['block_count'],1);self.assertIsNone(r['ci95_halfwidth'])
        self.assertEqual(r['estimator_sample_count'],10);self.assertEqual(r['status'],'WINDOW_INSUFFICIENT')

    def test_no_block_retains_explicit_descriptive_mean(self):
        r=block_statistics(np.arange(6.),np.arange(6.),min_duration=20)
        self.assertEqual(r['block_count'],0);self.assertEqual(r['mean'],2.5)
        self.assertIn('without_complete_block',r['estimator_scope']);self.assertIsNone(r['ci95_halfwidth'])

    def test_drift_and_correlated_noise_are_not_independent(self):
        self.assertGreater(correlation_time_samples(np.arange(1000.)),100)
        rng=np.random.default_rng(812);x=np.zeros(5000)
        for i in range(1,len(x)):x[i]=.98*x[i-1]+rng.normal()
        r=block_statistics(np.arange(len(x)),x,min_duration=1)
        self.assertGreater(r['block_samples'],100);self.assertLess(r['effective_sample_estimate'],500)

    def test_irregular_or_duplicate_time_cannot_be_sufficient(self):
        with self.assertRaisesRegex(ValueError,'regularly'):block_statistics([0,1,2,7],[1,2,3,4],min_duration=1)
        self.assertEqual(block_statistics([0,1,1,2],[1,2,3,4],min_duration=1)['status'],'WINDOW_INSUFFICIENT')

    def test_all_windows_preserved_and_earliest_rule_not_target_search(self):
        tol={'warmup_relaxation_multiples':2,'max_stationarity_drift':.1,'temperature_relative_error':.02}
        t=np.arange(1,1001.)/1000;v=np.full(1000,1.2)
        rows=window_audit(t,v,v,.001,tol,[.2,.3,.4,.5,.6,.7,.8,.9])
        self.assertEqual(len(rows),8);self.assertTrue(rows[0]['selected_by_original_rule'])
        self.assertEqual(sum(x['selected_by_original_rule'] for x in rows),1)
        self.assertTrue(all(r['temperature_mean_all']>1.1 for r in rows))


class MeasurementTests(unittest.TestCase):
    def test_units_mass_mu_and_nu(self):
        u=thermal_density_units(5e-7,1056,4,2,1,298.15)
        nu=u.to_star(3.27e-6,'kinematic_viscosity');mu=u.to_star(1056*3.27e-6,'dynamic_viscosity')
        self.assertAlmostEqual(mu/nu,8);self.assertAlmostEqual(u.temperature_K(1),298.15)
        self.assertAlmostEqual(u.to_si(8,'mass_density'),1056)
        self.assertAlmostEqual(u.scales['force']/2/u.M0,u.scales['acceleration']/2)

    def test_versioned_mass_fields_not_column_guess(self):
        arr=np.array([(4,8)],dtype=[('N',float),('mass_star',float)])
        old="moments.writerow([done,done*dt,N,N*spec['m_star'],raw])"
        self.assertTrue(mass_field_semantics({'m_star':2},old,arr)['matches_executed_formula'])
        new='moments.writerow([done,done*dt,N,mass,raw])'
        self.assertFalse(mass_field_semantics({'m_star':2},new,arr)['matches_executed_formula'])
        arr['mass_star']=2;self.assertTrue(mass_field_semantics({'m_star':2},new,arr)['matches_executed_formula'])
        with self.assertRaisesRegex(ValueError,'UNVERIFIED'):mass_field_semantics({'m_star':2},'unknown',arr)

    def test_equilibrium_and_flow_DOF_from_real_particle_moments(self):
        vel=np.array([[2,1,0],[4,-1,0],[-4,2,1],[-2,-2,-1.]])
        means=np.stack([vel[:2].mean(axis=0),[np.nan]*3,vel[2:].mean(axis=0)])
        sums=[np.sum(vel[:2]**2),0,np.sum(vel[2:]**2)]
        r=recover_profile_moments([2,0,2],means,sums,2)
        self.assertAlmostEqual(r['COM'],2*np.sum((vel-vel.mean(axis=0))**2)/9)
        self.assertAlmostEqual(r['coarse'],2*(np.sum((vel[:2]-means[0])**2)+np.sum((vel[2:]-means[2])**2))/6)
        self.assertEqual(r['occupied_bins'],2)

    def test_low_count_and_nonfinite_live_bins(self):
        with self.assertRaisesRegex(ValueError,'DOF'):recover_profile_moments([1,1],[[0]*3]*2,[0,0],1)
        with self.assertRaisesRegex(ValueError,'NONFINITE'):recover_profile_moments([2,2],[[np.nan,0,0],[0]*3],[0,0],1)

    def test_empty_bins_never_zero_velocity_observations(self):
        eq,weighted=average_profile([[2,np.nan],[4,8]],[[1,0],[3,1]])
        np.testing.assert_allclose(eq,[3,8]);np.testing.assert_allclose(weighted,[3.5,8])
        with self.assertRaisesRegex(ValueError,'WITHOUT_OBSERVATIONS'):average_profile([[1,np.nan]],[[1,0]])

    def test_particle_phase_and_periodic_translation(self):
        x=np.array([[3.99,0,0],[-2,1,0.]]);v=np.array([[2,3,1],[1,-1,0.]]);f=np.ones((2,3));dt=.1;m=2
        postv=v+dt*f/m;postx=(x+dt*postv+4)%8-4
        r=phase_audit(x,v,f,postx,postv,dt,m,[8]*3)
        self.assertLess(r['kick_max_abs'],1e-14);self.assertLess(r['drift_PBC_max_abs'],1e-14)
        wrong=phase_audit(x,v[::-1],f,postx,postv,dt,m,[8]*3)
        self.assertGreater(wrong['kick_max_abs'],1)

    def test_component_temperature_mean_equals_COM(self):
        v=np.random.default_rng(91).normal(size=(100,3))+[2,1,3]
        r=velocity_statistics(v,2);self.assertAlmostEqual(np.mean(r['component_kBT']),r['COM_kBT'])

    def test_force_padding_and_zero_cuda_layout(self):
        interface={'shape':(5,3),'typestr':'<f4','strides':(16,4),'data':(123,False)}
        shape,dtype,strides,size=cuda_layout(interface);self.assertEqual(size,76)
        raw=np.arange(20,dtype=np.float32);a=np.ndarray(shape,dtype,buffer=raw,strides=strides)
        np.testing.assert_array_equal(a[1],[4,5,6]);self.assertEqual(a[-1,-1],18)
        interface.update(shape=(0,3),data=(0,False));self.assertEqual(cuda_layout(interface)[-1],0)
        interface.update(shape=(1,3));
        with self.assertRaisesRegex(ValueError,'NULL'):cuda_layout(interface)
        interface.update(strides=(-16,4))
        with self.assertRaisesRegex(ValueError,'STRIDES'):cuda_layout(interface)

    def test_profile_time_and_bin_order_checked(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp);(d/'moments.csv').write_text('step,time_star\n10,1\n')
            header='step,time_star,bin,count,ux,uy,uz,sum_v2\n'
            (d/'profile_samples.csv').write_text(header+'10,1,0,2,1,0,0,2\n10,1,1,2,-1,0,0,2\n')
            spec={'domain_star':[2,2,2],'y_bin_star':1};profile_arrays(d,spec)
            (d/'profile_samples.csv').write_text(header+'10,1,0,2,1,0,0,2\n11,1,1,2,-1,0,0,2\n')
            with self.assertRaisesRegex(ValueError,'PHASE'):profile_arrays(d,spec)

    def test_partial_spatial_fit_can_recover_but_missing_data_is_not_zero(self):
        y=np.arange(.125,8,.25)[::2];v=.3/(2*4)*periodic_shape(y,8,.25)+.2
        self.assertAlmostEqual(fit_viscosity(y,v,.3,2,8,.25)['nu_star'],4)
        with self.assertRaises(ValueError):fit_viscosity(y[:3],v[:3],.3,2,8,.25)

    def test_shared_reference_covariance_cancels_for_endpoint_width(self):
        cov=shared_reference_covariance([4,9,16]);np.testing.assert_equal(cov,[[13,9],[9,25]])
        self.assertEqual(np.array([-1,1])@cov@np.array([-1,1]),20)

    def test_half_dt_uses_physical_time_not_step_count(self):
        self.assertEqual(matched_physical_interval([1,100],[1,100],.01,.005),(.01,.5))
        self.assertEqual(matched_physical_interval([1,100],[2,200],.01,.005),(.01,1.))

    def test_fixed_OU_reference_bias_not_solver_qualification(self):
        self.assertAlmostEqual(ou_discrete_temperature_ratio(1,.1),1/.95)
        self.assertLess(ou_discrete_temperature_ratio(1,.05),ou_discrete_temperature_ratio(1,.1))
        with self.assertRaises(ValueError):ou_discrete_temperature_ratio(2,1)


class ExecutionProtectionTests(unittest.TestCase):
    def test_diagnostic_force_saver_uses_pinned_native_channel_name(self):
        import re
        native=(PROJECT_ROOT/'vendor/Mirheo/src/mirheo/core/utils/common.cpp').read_text()
        force_name=re.search(r'const std::string forces\s*=\s*"([^"]+)"',native).group(1)
        worker=(PROJECT_ROOT/'py_scripts/sdpd_diagnostics/gpu_worker.py').read_text()
        self.assertIn("('"+force_name+"','saved_forces')",worker)

    def test_cooperative_stop_is_not_complete_comparison(self):
        e={'status':'STOPPED_OR_FAILED','timeout':True};w={'steps':99,'status':'COOPERATIVE_PARTIAL'}
        self.assertIn('PARTIAL',completion_state(e,w,100))
        e.update(status='COMPLETED',timeout=False);w.update(steps=100,status='COMPLETED_PLANNED_STEPS')
        self.assertEqual(completion_state(e,w,100),'COMPLETED_FIXED_PHYSICAL_WINDOW')

    def test_cache_budget_and_immutable_manifest(self):
        with tempfile.TemporaryDirectory(dir=PROJECT_ROOT/'test_code/outputs/sdpd_diagnostics') as temp:
            d=Path(temp);c={'campaign_id':'synthetic_diagnostic','budget':{'campaign_limit_s':2,'task_limit_s':2,'minimum_launch_budget_s':.1,'graceful_margin_s':.1}}
            result,record,hit=run_attempt(d,c,'synthetic',{'v':1},lambda x:[sys.executable,'-c','pass'],.8,monitor_gpu=False)
            charge=used_budget(ledger_read(d,c['campaign_id'],c['budget']))
            with patch('py_scripts.fluid_physics.runner.bounded_process',side_effect=AssertionError('CACHE_MUST_NOT_RUN')):
                self.assertTrue(run_attempt(d,c,'synthetic',{'v':1},lambda x:[],.8,monitor_gpu=False)[2])
            self.assertEqual(charge,used_budget(ledger_read(d,c['campaign_id'],c['budget'])))
            p=d/'old.json';write_json(p,{'original':1});before=sha256_file(p)
            with self.assertRaises(FileExistsError):write_json(p,{'original':2})
            self.assertEqual(before,sha256_file(p))
            write_json(d/'package_sha256.json',{'old.json':before});verify_package(d)
            p.write_text('{}')
            with self.assertRaisesRegex(ValueError,'HASH_MISMATCH'):verify_package(d)

    def test_help_and_no_arguments_do_not_import_GPU(self):
        script="import runpy,sys;sys.modules['mirheo']=None;sys.argv=['diagnose_sdpd']+sys.argv[1:];runpy.run_module('py_scripts.diagnose_sdpd',run_name='__main__')"
        for args,code in [(['--help'],0),([],2)]:
            r=subprocess.run([sys.executable,'-B','-c',script,*args],cwd=PROJECT_ROOT,capture_output=True,text=True)
            self.assertEqual(r.returncode,code);self.assertNotIn('PROBE_RESULT',r.stdout)


class RecordedProbeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=PROJECT_ROOT/'runs/sdpd_diagnostics/thermal_cause_20260909/equilibrium_half_dt_trend_corrected_channel'
        if not (cls.d/'output_sha256.json').exists():raise unittest.SkipTest('Real probe not yet executed; no synthetic substitute.')

    def test_owned_native_snapshots_verify_actual_kick_drift_and_ID_order(self):
        rows=read_json(self.d/'phase_observations.json');spec=read_json(self.d/'actual_parameters.json')
        self.assertEqual(len(rows),5)
        for i,row in enumerate(rows):
            a=[np.load(self.d/(f'snapshot_{i:02d}_'+n+'.npy')) for n in ('positions','pre_velocities','pre_forces','post_positions','post_velocities')]
            measured=phase_audit(*a,spec['dt_star'],spec['m_star'],spec['domain_star'])
            self.assertTrue(row['packed_ID_words_preserved']);self.assertEqual(row['getter_velocity_max_abs'],0)
            self.assertLess(measured['kick_max_abs'],4*np.finfo(np.float32).eps*np.max(abs(a[-1])))
            self.assertLess(measured['drift_PBC_max_abs'],np.finfo(np.float32).eps*max(spec['domain_star']))

    def test_actual_fixed_window_and_failed_job_charges_remain(self):
        e=read_json(self.d/'execution.json');w=read_json(self.d/'worker_completion.json')
        self.assertEqual(completion_state(e,w,160000),'COMPLETED_FIXED_PHYSICAL_WINDOW')
        self.assertFalse(e['remaining_own_group_processes'])
        failed=read_json(self.d.parent/'equilibrium_half_dt_trend/execution.json')
        self.assertGreater(failed['elapsed_monotonic_s'],0);self.assertNotEqual(failed['exit_code'],0)
        ledger=read_json(self.d.parent/'budget_ledger.json');self.assertEqual(len(ledger['attempts']),2)
        self.assertAlmostEqual(used_budget(ledger),e['elapsed_monotonic_s']+failed['elapsed_monotonic_s'])

    def test_initial_state_is_identical_between_observer_attempts(self):
        old=self.d.parent/'equilibrium_half_dt_trend'
        for name in ('initial_positions_global.npy','initial_velocities.npy'):
            self.assertEqual(sha256_file(old/name),sha256_file(self.d/name))

    def test_cache_rejects_changed_config_without_launching(self):
        from py_scripts.sdpd_diagnostics.probes import verify_probe_cache
        with patch('py_scripts.sdpd_diagnostics.probes.run_attempt',side_effect=AssertionError('NO_GPU')):
            with self.assertRaisesRegex(ValueError,'PROBE_CONFIG_CHANGED'):verify_probe_cache(self.d,{'_config_sha256':'changed'})

    def test_cache_rejects_changed_native_environment(self):
        from py_scripts.sdpd_diagnostics.probes import verify_probe_cache
        preflight=read_json(self.d/'preflight.json');prov=read_json(Path(preflight['diagnostic_package'])/'provenance.json')
        env=dict(prov['environment']);env['library_sha256']='changed'
        with patch('py_scripts.sdpd_diagnostics.probes.environment_identity',return_value=env):
            with self.assertRaisesRegex(ValueError,'NATIVE_ENVIRONMENT_CHANGED'):
                verify_probe_cache(self.d,{'_config_sha256':prov['identity']['config_sha256']})


if __name__=='__main__':unittest.main()
