"""CPU regression tests. Synthetic arrays here are never experimental evidence."""
import copy
import json
from pathlib import Path
import tempfile
import io
import contextlib
import unittest
from unittest.mock import patch
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT,read_json,write_json,fingerprint,sha256_file
from py_scripts.fluid_physics.budget_authorizations import account_extensions
from py_scripts.sdpd_diagnostics.extended import load_config,json_safe,require_budget
from py_scripts.sdpd_diagnostics.extended_analysis import observable_analysis,hac_trend,describe
from py_scripts.sdpd_diagnostics.extended_restart import (native_contract,archive_checkpoint,verify_checkpoint,
    compare_by_id,restore_coordinator,merge_segments,require_formal_restart,native_time)
from py_scripts.sdpd_diagnostics.gpu_worker import continuous_chunks


class ArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name);self.src=self.root/'native';self.src.mkdir()
        (self.src/'simulation.state.txt').write_text('0.002\n2000\n1\n')
        (self.src/'pv-1.h5').write_bytes(b'CPU synthetic placeholder, never restored on GPU')
        (self.src/'pv.xmf').write_text('<Xdmf><DataItem Format="HDF">pv-1.h5:/position</DataItem></Xdmf>')
        self.spec={'dt_star':1e-6,'m_star':1};self.cp=self.root/'archive'
        self.make()
    def tearDown(self):self.temp.cleanup()
    def make(self):
        return archive_checkpoint(self.src,self.cp,spec=self.spec,contract={'native_rng_persistence_proven':False},completion={'status':'COMPLETED_PLANNED_STEPS','steps':2001})
    def test_missing_rng_rejected_even_with_particles_and_time(self):
        with self.assertRaisesRegex(ValueError,'MISSING_RANDOM_STATE'):verify_checkpoint(self.cp)
    def test_diagnostic_integrity_is_not_full_restart(self):
        self.assertEqual(verify_checkpoint(self.cp,require_rng=False)['restart_validity'],'RESTART_NOT_VALIDATED')
    def test_missing_file_fails(self):
        (self.cp/'pv-1.h5').unlink()
        with self.assertRaisesRegex(ValueError,'MISSING'):verify_checkpoint(self.cp,require_rng=False)
    def test_changed_hash_fails(self):
        (self.cp/'pv-1.h5').write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError,'HASH'):verify_checkpoint(self.cp,require_rng=False)
    def test_incomplete_marker_fails(self):
        (self.cp/'COMPLETE.json').unlink()
        with self.assertRaisesRegex(ValueError,'INCOMPLETE'):verify_checkpoint(self.cp,require_rng=False)
    def test_manifest_hash_fails(self):
        (self.cp/'checkpoint_manifest.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'MANIFEST_HASH'):verify_checkpoint(self.cp,require_rng=False)
    def test_missing_time_fails_before_copy(self):
        (self.src/'simulation.state.txt').unlink()
        with self.assertRaisesRegex(ValueError,'TIME'):native_time(self.src,1e-6)
    def test_snapshot_only_fails(self):
        directory=self.root/'snapshot';directory.mkdir();np.save(directory/'positions.npy',np.zeros((4,3)))
        with self.assertRaisesRegex(ValueError,'SNAPSHOT'):verify_checkpoint(directory)
    def test_absolute_external_xmf_reference_fails(self):
        (self.src/'pv.xmf').write_text('<Xdmf><DataItem Format="HDF">/etc/passwd:/x</DataItem></Xdmf>')
        with self.assertRaisesRegex(ValueError,'UNSAFE'):archive_checkpoint(self.src,self.root/'bad',spec=self.spec,contract={},completion={'status':'COMPLETED_PLANNED_STEPS','steps':2001})
    def test_symlink_escape_fails(self):
        (self.src/'escape').symlink_to('/etc/passwd')
        with self.assertRaisesRegex(ValueError,'ESCAPES'):archive_checkpoint(self.src,self.root/'bad',spec=self.spec,contract={},completion={'status':'COMPLETED_PLANNED_STEPS','steps':2001})
    def test_existing_version_not_overwritten(self):
        with self.assertRaises(FileExistsError):self.make()
    def test_unfinished_writer_rejected(self):
        with self.assertRaisesRegex(ValueError,'EXIT_NORMALLY'):archive_checkpoint(self.src,self.root/'partial',spec=self.spec,contract={},completion={'status':'COOPERATIVE_PARTIAL','steps':2001})
    def test_restore_does_not_initialize_velocity(self):
        class Coordinator:
            def restart(s,path):s.called=path
            def getState(s):return type('State',(),{'current_step':2000,'current_time':.002})()
            def setVelocities(s,*a):raise AssertionError('velocity reinitialization')
        u=Coordinator();restore_coordinator(u,self.cp,diagnostic_only=True);self.assertEqual(u.called,str(self.cp))
    def test_preadvance_time_mismatch_fails(self):
        class Coordinator:
            def restart(s,p):pass
            def getState(s):return type('State',(),{'current_step':0,'current_time':0.})()
        with self.assertRaisesRegex(ValueError,'TIME_MISMATCH'):restore_coordinator(Coordinator(),self.cp,diagnostic_only=True)


class IdentityAndChainTests(unittest.TestCase):
    def test_id_reordering_compares_correctly(self):
        a={'ids':np.array([7,1,2]),'positions':np.arange(9).reshape(3,3)/10,'velocities':np.arange(9).reshape(3,3)}
        b={k:v[[2,0,1]] for k,v in a.items()};self.assertEqual(compare_by_id(a,b,[8,8,8])['status'],'PASS')
    def test_equal_means_do_not_hide_per_particle_jump(self):
        a={'ids':np.array([1,2]),'positions':np.zeros((2,3)),'velocities':np.zeros((2,3))};b=copy.deepcopy(a);b['velocities'][0,0]=1;b['velocities'][1,0]=-1
        self.assertEqual(compare_by_id(a,b,[8]*3)['status'],'FAIL')
    def test_duplicate_ids_fail(self):
        a={'ids':np.array([1,1]),'positions':np.zeros((2,3)),'velocities':np.zeros((2,3))}
        with self.assertRaisesRegex(ValueError,'IDS'):compare_by_id(a,a,[8]*3)
    def chain(self):
        return [{'status':'COMPLETED','segment_index':i,'parent_checkpoint':None if i==0 else 'cp0',
            'checkpoint_sha256':f'cp{i}','start_step':i*400,'end_step':(i+1)*400,'segment_steps':400,
            'rows':[{'step':n,'time_star':n*1e-6} for n in range(0 if i==0 else 600,(i+1)*400+1,200)]} for i in range(2)]
    def test_absolute_merge_does_not_duplicate_boundary(self):
        s=self.chain();s[1]['rows'].insert(0,{'step':400,'time_star':.0004,'restore_boundary':True})
        self.assertEqual([r['step'] for r in merge_segments(s,1e-6,200,800)],[0,200,400,600,800])
    def test_reset_time_or_duplicate_row_fails(self):
        s=self.chain();s[1]['rows'][0]={'step':0,'time_star':0}
        with self.assertRaisesRegex(ValueError,'SAMPLES'):merge_segments(s,1e-6,200,800)
    def test_broken_parent_chain_fails(self):
        s=self.chain();s[1]['parent_checkpoint']='other'
        with self.assertRaisesRegex(ValueError,'CHAIN'):merge_segments(s,1e-6,200,800)
    def test_partial_segment_cannot_reach_formal_window(self):
        s=self.chain();s[1]['status']='PARTIAL'
        with self.assertRaisesRegex(ValueError,'INCOMPLETE'):merge_segments(s,1e-6,200,800)
    def test_short_chain_not_complete(self):
        with self.assertRaisesRegex(ValueError,'NOT_REACHED'):merge_segments(self.chain()[:1],1e-6,200,800)
    def test_CPU_certificate_cannot_replace_real_restore(self):
        with self.assertRaisesRegex(ValueError,'NO_FORMAL'):require_formal_restart({'native_rng_persistence_proven':True},{'status':'PASS','evidence':'CPU_FAKE'})
    def test_native_rng_gap_blocks_even_close_particle_means(self):
        with self.assertRaisesRegex(ValueError,'NO_FORMAL'):require_formal_restart(native_contract(),{'status':'PASS','evidence':'REAL_GPU_SEPARATE_PROCESSES'})
    def test_formal_segment_entry_never_launches_with_current_native_contract(self):
        from py_scripts.sdpd_diagnostics.extended import run_formal_segments
        with patch('py_scripts.sdpd_diagnostics.extended.launch_task',side_effect=AssertionError('unexpected GPU launch')):
            with self.assertRaisesRegex(ValueError,'NO_FORMAL'):run_formal_segments({}, {}, {}, native_contract(), {'status':'PASS','evidence':'REAL_GPU_SEPARATE_PROCESSES'}, {}, None)
    def test_checkpoint_schedule_is_not_global_step_period(self):
        c=native_contract()
        self.assertIn('nExecutions=0',c['state_contract']['save_phase'])
        self.assertIn('200-step chunk',c['state_contract']['save_phase'])
    def test_absolute_chunk_cuts_keep_sampling_grid_and_same_coordinator(self):
        with tempfile.TemporaryDirectory() as td:
            d=Path(td);(d/'control').mkdir();calls=[]
            class U:
                def run(self,n,dt):calls.append(n)
            chunks=list(continuous_chunks(U(),d,True,400,200,1e-6,lambda:None,start_step=2000,stop_steps=[2001]))
            self.assertEqual(calls,[1,199,200]);self.assertEqual([r[1]+2000 for r in chunks],[2001,2200,2400])


class BudgetTests(unittest.TestCase):
    def account(self,attempts,task='A'):
        records=[{'scope':{'campaign_directory':'/tmp/campaign','task_ids':['A','B']},'additional_seconds':80}]
        with patch('py_scripts.fluid_physics.budget_authorizations.authorization_records',return_value=records):
            return account_extensions({'limit_s':100},attempts,'/tmp/campaign',task)
    def test_segments_share_one_extension_no_reset(self):
        rows=[('/tmp/campaign','A',50),('/tmp/campaign','B',40)]
        a=self.account(rows);self.assertEqual(a['remaining_s'],90);self.assertEqual(a['original_charged_or_reserved_s'],10)
    def test_wrong_task_cannot_spend_scoped_extension(self):
        self.assertEqual(self.account([],task='C')['remaining_s'],100)
        self.assertEqual(self.account([],task='A')['remaining_s'],180)
    def test_cumulative_full_plan_gate(self):
        p={'scope_task_ids':['A','B'],'tasks':[{'id':'A','allocation_s':30},{'id':'B','allocation_s':30}]}
        with self.assertRaises(PermissionError):require_budget(p,{'remaining_s':59,'members':[]})
    def test_single_task_limit_gate(self):
        p={'scope_task_ids':['A'],'tasks':[{'id':'A','allocation_s':601}]}
        with self.assertRaisesRegex(ValueError,'SINGLE_TASK'):require_budget(p,{'remaining_s':10000,'members':[]})
    def test_no_second_concurrent_job(self):
        p={'scope_task_ids':['A'],'tasks':[{'id':'A','allocation_s':30}]}
        with self.assertRaisesRegex(RuntimeError,'CONCURRENCY'):require_budget(p,{'remaining_s':1000,'members':[{'running_reservation_s':1}]})
    def test_scope_plan_does_not_accept_extra_task(self):
        with self.assertRaisesRegex(ValueError,'OUTSIDE'):require_budget({'scope_task_ids':['A'],'tasks':[]},{},['B'])


class PlanMetadataTests(unittest.TestCase):
    def test_new_plan_records_actual_extended_worker_hash(self):
        from py_scripts.sdpd_diagnostics.extended import history,build_plan
        c=load_config(PROJECT_ROOT/'py_scripts/sdpd_equilibration_extended.yaml')
        p=build_plan(c,history(c),native_contract())
        self.assertEqual(p['parameters']['script_sha256'],sha256_file(PROJECT_ROOT/'py_scripts/sdpd_diagnostics/extended_worker.py'))


class StatisticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plan=read_json(PROJECT_ROOT/'data/sdpd_equilibration/result_26300860927d8b87/run_plan.json')
    def data(self,n=2001,temp=None):
        p=self.plan;d={'time_star':np.arange(n)*.0002,'N':np.full(n,4096),'total_mass_star':np.full(n,4096)};d['kBT_COM_star']=np.full(n,1.0) if temp is None else np.asarray(temp)
        d['pressure_star']=np.full(n,100000.);d['kernel_number_mean']=np.full(n,8.)
        d['kernel_number_variance']=np.full(n,.1)
        for k in ['kernel_q10','kernel_q50','kernel_q90']:d[k]=np.full(n,8.)
        for a in 'xyz':d['mean_v'+a]=np.zeros(n)
        return d
    def analyze(self,d):return observable_analysis(d,self.plan,{'status':'PASS','measurement_status':'PASS'})
    def test_stable_high_temperature_is_bias_not_more_waiting(self):
        r=self.analyze(self.data(temp=np.full(2001,1.05)))
        self.assertEqual(r['OVERALL_EQUILIBRIUM']['branch'],'STATIONARY_THERMAL_BIAS')
    def test_coherent_drift_is_transient(self):
        # Within-quarter excursions much smaller than half-window drift.
        t=np.arange(2001)*.0002;y=1+np.floor(np.minimum(t,.39999)/.025)*.02
        r=self.analyze(self.data(temp=y));self.assertEqual(r['OVERALL_EQUILIBRIUM']['branch'],'TRANSIENT_PERSISTS')
    def test_short_data_no_CI_or_formal_window_relocation(self):
        r=self.analyze(self.data(n=1700));item=r['TEMPERATURE_STATIONARITY']
        self.assertEqual(item['status'],'INCONCLUSIVE');self.assertIsNone(item['steady_statistics']);self.assertEqual(item['window_star'],[.3,.4])
    def test_temperature_statistics_survive_pressure_inconclusive(self):
        d=self.data();d['pressure_star']+=np.arange(2001)*10
        r=self.analyze(d);self.assertEqual(r['TEMPERATURE_STATIONARITY']['status'],'PASS');self.assertNotEqual(r['OVERALL_EQUILIBRIUM']['status'],'PASS')
    def test_full_block_mean_and_CI_use_identical_samples(self):
        d=self.data();d['kBT_COM_star']+=np.sin(np.arange(2001))*1e-6
        r=self.analyze(d)['TEMPERATURE_STATIONARITY']['conditional_block_statistics']
        self.assertAlmostEqual(r['mean'],float(np.mean(r['block_means'])),places=12)
        self.assertEqual(r['estimator_sample_count'],r['block_count']*r['block_samples'])
    def test_descriptive_mean_not_steady_estimate(self):
        d=self.data(temp=np.linspace(2,1,2001));r=self.analyze(d)['TEMPERATURE_STATIONARITY']
        self.assertIsNotNone(r['description']['mean']);self.assertIsNone(r['steady_statistics'])
    def test_missing_samples_cannot_pass_even_if_end_time_matches(self):
        d={k:np.delete(v,1510) for k,v in self.data().items()};r=self.analyze(d)
        self.assertEqual(r['TEMPERATURE_STATIONARITY']['status'],'INCONCLUSIVE')
    def test_nonconservation_cannot_pass(self):
        d=self.data();d['N'][-1]=4095;r=self.analyze(d)
        self.assertEqual(r['OVERALL_EQUILIBRIUM']['conservation_check'],'FAIL_OR_MISSING');self.assertNotEqual(r['OVERALL_EQUILIBRIUM']['status'],'PASS')
    def test_HAC_does_not_claim_steady_ACF(self):
        t=np.arange(500)*.0002;r=hac_trend(t,1-.2*t+np.sin(np.arange(500))*.001)
        self.assertEqual(r['direction'],'DECREASING');self.assertIsNone(r['measured_steady_ACF'])
    def test_json_actual_numpy_types_supported(self):
        x=json_safe({'bool':np.bool_(True),'float':np.float64(3),'array':np.array([1,2]),'int':np.int64(7)})
        self.assertEqual(json.loads(json.dumps(x,allow_nan=False)),{'bool':True,'float':3.0,'array':[1,2],'int':7})
    def test_nonfinite_JSON_rejected_not_silently_NULL(self):
        with self.assertRaises(ValueError):json_safe({'value':np.nan})
    def test_fixed_config_uses_original_criteria(self):
        c=load_config(PROJECT_ROOT/'py_scripts/sdpd_equilibration_extended.yaml');self.assertEqual(c['criteria'],self.plan['criteria'])
    def test_imports_do_not_initialize_native_solver(self):
        import sys
        self.assertNotIn('libmirheo',sys.modules);self.assertNotIn('mirheo',sys.modules)
    def test_help_never_invokes_GPU_entrypoints(self):
        from py_scripts import run_sdpd_equilibration as run
        from test_code import review_sdpd_equilibration as review
        from py_scripts.sdpd_diagnostics import extended_worker as worker
        with patch('py_scripts.sdpd_diagnostics.extended.execute',side_effect=AssertionError('GPU')):
            for module in [run,review,worker]:
                with patch('sys.argv',['command','--help']),contextlib.redirect_stdout(io.StringIO()):
                    with self.assertRaises(SystemExit) as cm:module.main()
                self.assertEqual(cm.exception.code,0)


if __name__=='__main__':unittest.main()
