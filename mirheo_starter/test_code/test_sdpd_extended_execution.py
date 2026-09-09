"""CPU reproductions of real 2026-09-09 restore/export failures; no GPU calls."""
import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from py_scripts.fluid_physics.common import PROJECT_ROOT, read_json, write_json
from py_scripts.sdpd_diagnostics import extended as ext
from py_scripts.sdpd_diagnostics.extended_restart import archive_checkpoint, restore_coordinator, checkpoint_channel_forms


class CheckpointReadabilityTests(unittest.TestCase):
    def test_tensor6_stress_accepted_but_one_column_vector_rejected(self):
        with tempfile.TemporaryDirectory() as td:
            xmf = Path(td)/'pv.xmf'
            xmf.write_text('''<Xdmf><Attribute Name="saved_stresses">
              <Information Value="Tensor6" Datatype="Stress"/>
              <DataItem Dimensions="4096 6" Format="HDF">pv.h5:/saved_stresses</DataItem></Attribute></Xdmf>''')
            self.assertEqual(checkpoint_channel_forms(td)['status'],'PASS')
            xmf.write_text(xmf.read_text().replace('Tensor6','Vector').replace('4096 6','4096 1'))
            self.assertEqual(checkpoint_channel_forms(td)['status'],'FAIL')

    def test_unsupported_force_form_rejected_before_native_restart(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td); source = root/'source'; source.mkdir()
            (source/'simulation.state.txt').write_text('0.002\n2000\n0\n')
            (source/'pv.h5').write_bytes(b'CPU fixture, never GPU data')
            (source/'pv.xmf').write_text('''<Xdmf><Domain><Grid>
              <Attribute Name="saved_forces"><Information Value="Other" Datatype="Force"/>
              <DataItem Dimensions="4096 1" Format="HDF">pv.h5:/saved_forces</DataItem></Attribute>
            </Grid></Domain></Xdmf>''')
            archive = archive_checkpoint(source, root/'archive', spec={'dt_star':1e-6,'m_star':1},
                contract={}, completion={'status':'COMPLETED_PLANNED_STEPS','steps':2001})
            class Coordinator:
                called = False
                def restart(self, directory): self.called = True
                def getState(self): return type('State', (), {'current_step':2000,'current_time':.002})()
            coordinator = Coordinator()
            with self.assertRaisesRegex(ValueError, 'UNSUPPORTED_CHECKPOINT_CHANNEL.*saved_forces'):
                restore_coordinator(coordinator, archive, diagnostic_only=True)
            self.assertFalse(coordinator.called)


CAMPAIGN = PROJECT_ROOT/'runs/sdpd_equilibration_extended/fixed_late_restart_20260909'
PREPARATION = PROJECT_ROOT/'data/sdpd_equilibration_extended/preparation_32d12b977de9a64b'


@unittest.skipUnless((CAMPAIGN/'restart_validation.json').is_file(), 'Delivered real GPU evidence not installed')
class ActualExecutionRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = ext.load_config(PROJECT_ROOT/'py_scripts/sdpd_equilibration_extended.yaml')
        cls.history = ext.history(cls.config)
        cls.contract = read_json(PREPARATION/'restart_contract.json')

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.c = copy.deepcopy(self.config); self.c['data_root'] = self.temp.name
        for name, value in [('history', self.history), ('native_contract', self.contract),
                            ('historical_supplement', read_json(PREPARATION/'historical_supplement.json')),
                            ('observable_analysis', read_json(PREPARATION/'historical_observable_status.json'))]:
            p = patch.object(ext, name, return_value=value); p.start(); self.addCleanup(p.stop)
        def output(path):
            path = Path(path)
            if path.is_relative_to(CAMPAIGN/'preparations'):
                return Path(self.temp.name)/'preparations'/path.name
            return path
        p = patch.object(ext, 'output', side_effect=output); p.start(); self.addCleanup(p.stop)
        # Test exports must not write preparation records into the delivered campaign.
        p = patch.object(ext, 'write_preparation_record', create=True); p.start(); self.addCleanup(p.stop)

    def export(self):
        return ext.export(self.c, analyze=True)

    def test_actual_three_GPU_attempts_are_reported(self):
        package = self.export(); summary = read_json(package/'equilibration_summary.json')
        self.assertEqual(summary['execution_status'], 'RESTART_DIAGNOSTIC_FAILED')
        self.assertEqual(summary['GPU_attempt_count'], 3)
        self.assertAlmostEqual(summary['new_GPU_elapsed_s'], 11.711100826971233)
        self.assertIsNone(summary['actual_steps'])
        self.assertEqual(summary['formal_execution_status'], 'NOT_RUN')

    def test_consumed_diagnostics_do_not_create_new_budget_request(self):
        budget = read_json(self.export()/'budget_request.json')
        self.assertEqual(budget['additional_request_whole_seconds'], 0)
        self.assertEqual(budget['extra_authorized_gpu_seconds'], 1379)
        self.assertTrue(budget['user_approval_received'])
        self.assertEqual(budget['status'], 'STOPPED_AFTER_RESTART_FAILURE')

    def test_reanalysis_code_change_preserves_frozen_executed_plan(self):
        from py_scripts.sdpd_diagnostics.equilibration import code_identity
        modified = {**code_identity(), 'CPU_ANALYSIS_REVISION_FOR_TEST':'different'}
        with patch('py_scripts.sdpd_diagnostics.equilibration.code_identity', return_value=modified):
            package = self.export()
        self.assertEqual(read_json(package/'run_plan.json'), read_json(CAMPAIGN/'restart_A/run_plan.json'))

    def test_execute_after_terminal_failure_cannot_retry(self):
        with patch.object(ext, 'launch_task', side_effect=AssertionError('GPU retry forbidden')):
            package = ext.execute(self.c)
        self.assertEqual(read_json(package/'equilibration_summary.json')['execution_status'], 'RESTART_DIAGNOSTIC_FAILED')

    def test_cpu_reanalysis_cannot_change_executed_criteria(self):
        self.c['criteria']['temperature_stationarity_relative_tolerance'] = .1
        with self.assertRaisesRegex(ValueError, 'EXECUTED_DESIGN_CHANGED'):
            self.export()

    def test_own_save_matches_but_A_B_difference_is_not_hidden(self):
        audit = read_json(self.export()/'checkpoint_failure_audit.json')
        self.assertEqual(audit['file_integrity'],'PASS')
        self.assertEqual(audit['native_readability']['status'],'FAIL')
        self.assertEqual(audit['comparisons']['checkpoint_HDF_vs_own_saved_boundary']['max_velocity_error_star'],0.)
        self.assertEqual(audit['comparisons']['independent_A_vs_B_before_any_restore']['status'],'FAIL')
        self.assertEqual(audit['restored_preadvance_state'],'NOT_REACHED')
        self.assertFalse(audit['RNG_persistence_pass'])

    def test_real_bad_archive_rejected_before_mirheo_import(self):
        from py_scripts.sdpd_diagnostics import extended_worker as worker
        import builtins
        original_import = builtins.__import__
        def guarded_import(name,*args,**kwargs):
            if name in ['mirheo','libmirheo']: raise AssertionError('Must reject before Mirheo import')
            return original_import(name,*args,**kwargs)
        with patch.dict('os.environ', {'OMPI_COMM_WORLD_SIZE':'2','OMPI_COMM_WORLD_RANK':'0'}), \
             patch('sys.argv',['worker','--spec',str(CAMPAIGN/'restart_B_restore/actual_parameters.json')]), \
             patch('builtins.__import__',side_effect=guarded_import):
            with self.assertRaisesRegex(ValueError,'UNSUPPORTED_CHECKPOINT_CHANNEL'):
                worker.main()

    def test_viewer_shows_actual_failure_budget_and_independent_curves(self):
        from test_code import review_sdpd_equilibration_extended as viewer
        import json
        package = self.export(); self.c['review_root'] = self.temp.name+'/review'
        with patch.object(viewer,'output',side_effect=lambda path: Path(path)):
            page = viewer.write_page(self.c, package)
        text = page.read_text()
        self.assertTrue('本轮 GPU：NOT_RUN' not in text, 'Executed GPU tasks still labeled NOT_RUN')
        self.assertIn('saved_forces', text)
        data = json.loads(text.split('id="equilibration-audit-data" type="application/json">')[1].split('</script>')[0])
        self.assertEqual(data['GPU_attempt_count'],3)
        self.assertEqual(data['diagnostic_raw_rows'],[21,11,0])
        self.assertEqual(data['requested_additional_s'],0)
        self.assertEqual(data['extra_authorized_gpu_seconds'],1379)


if __name__ == '__main__': unittest.main()
