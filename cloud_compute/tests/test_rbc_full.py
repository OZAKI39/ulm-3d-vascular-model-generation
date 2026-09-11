import unittest,tempfile,sys,subprocess,json,os,ast
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common import *
import cloud_run,rbc_plan
from remote_rbc import claim
from rbc_checks import outcome,frame_check,preparation
from rbc_h5 import read_xmf
from rbc_geometry import independent_intersections
from cloud_geometry import read_off,geometry

BASE=Path('/home/lzy/projects/cloud_compute')
M=Path('/home/lzy/projects/mirheo_starter')
S=M/'runs/single_rbc_repair/rbc_repair_20260910T131105Z/gpu/A6_continuous_preparation_half_dt/spec.json'


class FullRBCChecks(unittest.TestCase):
    def spec(self):
        s=read(S);s.update(steps=400000,prep_sample_steps=1000,sample_steps=10000)
        criteria=read(M/'data/single_rbc_repair/rbc_repair_20260910T131105Z/material_matching.json')
        s['preparation_criteria']=criteria['common_preparation_criteria'];v,f=read_off(s['mesh']);s['reference_radius']=geometry(v,f)['a']
        return s

    def test_actual_cloud_library_source_chain(self):
        plan=read(BASE/'build_plan.json');identity=read('/home/lzy/projects/cloud_results/mirheo-sm120-20731713865ae510/compile_identity.json')
        self.assertEqual(identity['runtime_binary_hash'],rbc_plan.EXPECTED_LIBRARY)
        self.assertEqual(identity['original_source_sha256'],plan['original_source_sha256'])
        self.assertEqual(identity['target_architectures'],['120']);self.assertTrue(identity['physics_sources_unchanged'])

    def test_phase_counts_and_old_protocol_rejection(self):
        s=self.spec();r=rbc_plan.validate_protocol(s)
        self.assertEqual((r['prep_steps'],r['shear_steps'],r['wall_hidden_steps']),(60000,400000,4000))
        self.assertEqual(r['shear_steps']*r['dt']*r['shear_rate'],4.)
        s['continuous_observation']=False
        with self.assertRaisesRegex(ValueError,'INCOMPATIBLE'):rbc_plan.validate_protocol(s)

    def test_no_flags_or_new_alone_cannot_start(self):
        with patch('rbc_plan.preflight',side_effect=AssertionError('No remote operation')):
            for a in [(False,False,False),(False,True,False),(False,False,True),(True,True,True)]:
                with self.subTest(a=a),self.assertRaises(ValueError):cloud_run.rbc_full(*a)

    def test_authorization_is_separate_and_single_use(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);plan={'plan_sha256':'fixed'}
            with self.assertRaisesRegex(ValueError,'AUTHORIZATION_REQUIRED'):rbc_plan.authorization(p,plan)
            write(p/'rbc_full/authorization.json',dict(plan_sha256='fixed',purpose=rbc_plan.PURPOSE,process_seconds=1800,attempts=1,user_approval_quote='fixture approval'))
            self.assertEqual(rbc_plan.authorization(p,plan)['attempts'],1)
            write(p/'rbc_full/launch_receipt.json',{'job_id':'already-reserved'})
            with self.assertRaisesRegex(ValueError,'ALREADY_RESERVED'):rbc_plan.authorization(p,plan)
        with tempfile.TemporaryDirectory() as tmp:
            claim(tmp)
            with self.assertRaises(FileExistsError):claim(tmp)

    def test_failed_last_frame_is_only_lower_bound(self):
        s=self.spec();frames=[dict(phase='relaxation',phase_step=60000),dict(phase='shear',phase_step=390000)]
        q=outcome(s,None,1,frames)
        self.assertEqual(q['CLOUD_RBC_RUN_COMPLETE'],'NOT_COMPLETE');self.assertAlmostEqual(q['strain_lower_bound'],3.9)
        self.assertIsNone(q['actual_successful_returned_shear_steps'])
        completed={'completed':True,'shear_steps':400000,'prep_steps':60000,'strain_end':4}
        self.assertEqual(outcome(s,completed,1,frames)['CLOUD_RBC_RUN_COMPLETE'],'NOT_COMPLETE')
        self.assertEqual(outcome(s,completed,0,frames)['CLOUD_RBC_RUN_COMPLETE'],'NOT_COMPLETE')
        frames.append(dict(phase='shear',phase_step=400000,vertices=[[1,2,3]]))
        self.assertEqual(outcome(s,completed,0,frames)['CLOUD_RBC_RUN_COMPLETE'],'PASS')
        self.assertIsNone(q['qualified_speedup'])

    def test_adjacent_faces_and_real_crossing(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float);f=np.array([[0,1,2],[0,1,3]])
        q=independent_intersections(v,f);self.assertEqual(q['confirmed_count'],0);self.assertEqual(q['uncertain_count'],0)
        v=np.array([[-1,-1,0],[1,-1,0],[0,1,0],[0,-.5,-1],[0,-.5,1],[0,.5,.5]],float)
        self.assertEqual(independent_intersections(v,np.array([[0,1,2],[3,4,5]]))['confirmed_count'],1)

    def test_read_existing_native_hdf5_without_solver(self):
        p=Path('/home/lzy/projects/cloud_results/cloud-smoke-20260910T221827Z/C_repair/native_relaxation/rbc00000.xmf')
        r=read_xmf(p);self.assertEqual(r['position'].shape,(642,3));self.assertEqual(r['id'].dtype,np.int64)
        import h5py
        with h5py.File(p.with_suffix('.h5'),'r') as f:
            for key in r:self.assertTrue(np.array_equal(r[key].reshape(-1),np.asarray(f[key]).reshape(-1)),key)
        self.assertNotIn('mirheo',sys.modules);self.assertNotIn('libmirheo',sys.modules)

    def test_original_preparation_window_definition(self):
        s=self.spec();v,f=read_off(s['mesh']);g=geometry(v,f)
        frames=[dict(phase='relaxation',phase_step=step,time_star=step*s['dt'],vertices=(v+[12,12,12]).tolist(),geometry=g,hard_failures=[],soft_failures=[]) for step in range(0,60001,1000)]
        q=preparation(frames,v,s);self.assertTrue(q['passed']);self.assertFalse(q['equilibrium_certified'])
        frames[-1]['soft_failures']=['area_relative_drift'];self.assertFalse(preparation(frames,v,s)['passed'])
        self.assertFalse(preparation(frames[:5],v,s)['passed'])

    def test_periodic_identity_checks_and_real_reference(self):
        s=self.spec();v,f=read_off(s['mesh']);shift=v+[23,12,12];shift[:,0]%=24
        q=frame_check(np.arange(len(v)),shift,np.zeros_like(v),f,v,s)
        self.assertEqual(q['hard_failures'],[]);self.assertLess(q['area_relative_drift'],1e-10)
        bad=np.arange(len(v));bad[-1]=bad[0]
        with self.assertRaisesRegex(ValueError,'DUPLICATE'):frame_check(bad,shift,np.zeros_like(v),f,v,s)

    def test_guard_user_stop_and_closed_failure_archive(self):
        p=subprocess.Popen([sys.executable,'-B','-c','import time;time.sleep(20)'],start_new_session=True)
        q=wait_guarded(p,2,lambda:10*1024**3,stop_reason=lambda:'USER_STOP')
        self.assertEqual(q['stop_reason'],'USER_STOP');self.assertIsNotNone(p.poll())
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'console.log').write_text('retained actual failure')
            seal(root,'test','STOPPED');self.assertEqual(validate_result(root)['run_state'],'STOPPED')

    def test_cpu_view_modules_do_not_import_or_launch_solver(self):
        for name in ['rbc_plan.py','rbc_checks.py','rbc_h5.py','rbc_report.py','rbc_geometry.py']:
            source=(BASE/name).read_text();tree=ast.parse(source)
            imports=[n for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom))]
            for n in imports:self.assertNotIn('libmirheo',ast.unparse(n));self.assertNotEqual(getattr(n,'module',None),'mirheo')
        for args in [['--help'],['rbc-full','--help'],['rbc-full','--new']]:
            q=subprocess.run([sys.executable,'-B',str(BASE/'cloud_run.py'),*args],capture_output=True,text=True)
            self.assertEqual(q.returncode,2 if args[-1]=='--new' else 0)

    def test_extracted_checker_function_bodies_are_unchanged(self):
        got=ast.parse((BASE/'rbc_geometry.py').read_text());functions={n.name:ast.dump(n,include_attributes=False) for n in got.body if isinstance(n,ast.FunctionDef)}
        for relative,record in read(BASE/'rbc_geometry_sources.json').items():
            source=M/'py_scripts'/relative;self.assertEqual(sha(source),record['sha256'])
            for n in ast.parse(source.read_text()).body:
                if isinstance(n,ast.FunctionDef) and n.name in record['functions']:self.assertEqual(functions[n.name],ast.dump(n,include_attributes=False))

    def test_export_repair_reads_existing_archive_without_launch(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'raw';root.mkdir();(root/'simulation').mkdir();s=self.spec()
            s['cloud_full_protocol']='CLOUD_FULL_RBC_V1'
            execution=dict(status='FAILED',exit_code=1,job_id='CPU_EXPORT_REPAIR_FIXTURE',stop_reason='fixture failed before output')
            write(root/'execution.json',execution);write(root/'run_status.json',execution);write(root/'full_spec.json',s)
            write(root/'provenance.json',dict(python_inputs=[dict(remote=s['mesh'],local=s['mesh'],sha256=sha(s['mesh']))]))
            write(root/'disk_after.json',dict(bytes=0,fixture_only=True));seal(root,'fixture','FAILED');before=sha(root/'results_manifest.json')
            from rbc_report import render_review
            with patch.object(cloud_run,'rbc_full',side_effect=AssertionError('No solver')),patch.object(cloud_run,'smoke',side_effect=AssertionError('No smoke')):
                result=render_review(root,Path(tmp)/'review')
            self.assertTrue(Path(result).is_file());self.assertEqual(validate_result(root)['manifest_sha256'],before)
            q=read(Path(tmp)/'review/numerical_checks.json');self.assertEqual(q['CLOUD_RBC_RUN_COMPLETE'],'NOT_COMPLETE')
            self.assertFalse((root/'numerical_checks.json').exists())


if __name__=='__main__':unittest.main()
