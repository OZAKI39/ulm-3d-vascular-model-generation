"""Synthetic CPU reader regressions; no native solver or timing claims."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from py_scripts.single_rbc_repair.runtime_analysis import inspect_run
from py_scripts.single_rbc_repair.workflow import load_config,paths
from py_scripts.single_rbc_benchmark.physics import read_off


class RuntimeReaderRegression(unittest.TestCase):
    def setUp(self):
        self.config=load_config();b,_,_=paths(self.config)
        parent=Path(__file__).resolve().parents[1]/'test_code/outputs/single_rbc_repair/cpu_reader_fixtures'
        parent.mkdir(parents=True,exist_ok=True)
        self.temp=tempfile.TemporaryDirectory(dir=parent)
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.run=self.root/'synthetic_cpu_fixture';self.run.mkdir()
        self.mesh=b/'common_reference.off';self.vertices,self.faces=read_off(self.mesh)
        self.positions=self.vertices+12
        self.ids=np.arange(len(self.vertices))
        self.save('execution.json',dict(status='SYNTHETIC_CPU_FIXTURE',elapsed_monotonic_s=0,exit_code=0))
        self.save('spec.json',dict(config=self.config,mesh=str(self.mesh),dt=.001,prep_steps=0,steps=0))
        self.save('completion.json',dict(completed=True,actual_steps=0))
        with (self.run/'vertices.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['step','phase','time_star','strain','vertex','x','y','z','fx','fy','fz'])
            for i,v in enumerate(self.positions):w.writerow([0,'relaxation',0,0,i,*v,0,0,0])
        with (self.run/'moments.csv').open('w') as f:
            w=csv.writer(f);w.writerow(['step','phase','time_star','strain','N','inner_N','temperature','max_speed','wall_crossings']);w.writerow([0,'relaxation',0,0,1,0,1,0,0])
        self.probes()
        np.savez(self.run/'membrane_state_relaxation_00000000.npz',ids=self.ids,positions=self.positions,velocities=np.zeros_like(self.positions),phase='relaxation',global_step=0,strain=0)

    def save(self,name,value):
        (self.run/name).write_text(json.dumps(value))

    def probes(self,phase='relaxation',ids=None):
        np.savez(self.run/'fluid_probes_relaxation_00000000.npz',probes=np.array([[0,99,0,0,0]],float),phase=phase,global_step=0,strain=0,membrane_vertex_ids=self.ids if ids is None else ids)

    def test_consistent_saved_fixture_still_has_no_speedup_or_full_range_claim(self):
        result=inspect_run(self.run,self.root/'review')
        self.assertEqual(result['status'],'SAVED_SHORT_CHECKS_PASS')
        self.assertIsNone(result['qualified_speedup'])
        self.assertIn('not a full-range fix',result['scope'])

    def test_equal_strain_cannot_hide_a_different_phase(self):
        self.probes(phase='shear')
        with self.assertRaisesRegex(ValueError,'TIME_MISMATCH'):
            inspect_run(self.run,self.root/'review')

    def test_sorted_csv_rows_cannot_hide_duplicate_native_ids(self):
        bad=self.ids.copy();bad[-1]=bad[0];self.probes(ids=bad)
        with self.assertRaisesRegex(AssertionError,'ID_CONNECTIVITY'):
            inspect_run(self.run,self.root/'review')

    def test_completion_boolean_cannot_hide_a_wrong_step_count(self):
        self.save('completion.json',dict(completed=True,actual_steps=7))
        result=inspect_run(self.run,self.root/'review')
        self.assertFalse(result['short_screen']['completed_planned_short_protocol'])
        self.assertEqual(result['status'],'SHORT_CONTROL_HAS_UNRESOLVED_QUALITY_FAILURE')


if __name__=='__main__':unittest.main()
