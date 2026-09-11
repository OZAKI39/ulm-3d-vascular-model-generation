import unittest,tempfile,sys,subprocess,os,json
from pathlib import Path
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from common import *
import cloud_run
from remote_build import check_upload

class WorkflowChecks(unittest.TestCase):
    def bundle(self,root,state='COMPLETED'):
        root.mkdir();(root/'console.log').write_text('actual task output\n');write(root/'completion.json',dict(status=state,actual_steps=3 if state=='COMPLETED' else 0))
        seal(root,'test-job',state);return validate_result(root)

    def test_archived_source_path_mapping(self):
        u=cloud_run.UPLOAD;n=read(u/'metadata/native_source_locations.json');m=read(u/'metadata/path_mapping.json');manifest=read(u/'transfer_manifest.json')
        _,resolved,key,project=cloud_run.resolve_source(n,m,manifest)
        self.assertEqual(resolved,project+'/data/single_rbc_repair/rbc_repair_20260910T131105Z/native/source')
        self.assertTrue(key.startswith('mirheo_starter/data/'))

    def test_old_vendor_is_rejected(self):
        native={'status':'CONFIRMED_LATEST_NATIVE_SOURCE','latest_source':'/local/vendor/Mirheo','native_trees':[{'path':'/local/vendor/Mirheo'}]}
        with self.assertRaisesRegex(ValueError,'OLD_VENDOR'):cloud_run.resolve_source(native,{}, {})

    def test_upload_modification_is_detected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'physical.cu';p.write_text('original')
            plan={'upload_dir':tmp,'upload_protection':[{'path':'physical.cu','sha256':sha(p),'execute_bits':0}]}
            self.assertEqual(check_upload(plan)['status'],'PASS');p.write_text('changed')
            with self.assertRaisesRegex(RuntimeError,'SNAPSHOT_CHANGED'):check_upload(plan)

    def test_manifest_traversal_and_absolute_paths(self):
        for name in ['/etc/passwd','../secret','x/../secret','x//y','./x','x\\y']:
            with self.subTest(name=name),self.assertRaises(ValueError):safe_relative(name)

    def test_symlink_and_special_file_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'link').symlink_to('/etc/passwd')
            with self.assertRaisesRegex(ValueError,'SYMLINK'):files(root)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);os.mkfifo(root/'pipe')
            with self.assertRaisesRegex(ValueError,'SPECIAL'):files(root)

    def test_unfinished_results_cannot_be_archived(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaisesRegex(ValueError,'NOT_READY'):validate_result(root)
            with self.assertRaisesRegex(ValueError,'CANNOT_SEAL'):seal(root,'job','RUNNING')
            self.assertFalse((root/'RESULTS_READY').exists())

    def test_hash_and_permissions_corruption_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'bundle';self.bundle(root);p=root/'console.log';p.write_text('corruption')
            with self.assertRaisesRegex(ValueError,'CONTENT_MISMATCH'):validate_result(root)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'bundle';self.bundle(root);p=root/'console.log';p.chmod(p.stat().st_mode|0o100)
            with self.assertRaisesRegex(ValueError,'CONTENT_MISMATCH'):validate_result(root)

    def test_repeated_fetch_never_starts_a_process(self):
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'job';check=self.bundle(dest)
            with patch.object(cloud_run,'local_environment'),patch.object(cloud_run,'result_location',return_value=('/remote/job',dest)),patch.object(cloud_run,'result_probe',return_value={'check':check}),patch.object(cloud_run.subprocess,'run',side_effect=AssertionError('No subprocess allowed on completed fetch')):
                self.assertEqual(cloud_run.fetch('job')['status'],'RESULTS_RETURN_VERIFIED')

    def test_cleanup_is_preview_only(self):
        with patch.object(cloud_run,'result_probe',side_effect=AssertionError('Must reject before remote action')):
            with self.assertRaisesRegex(ValueError,'DELETION_NOT_AUTHORIZED'):cloud_run.cleanup('job',False)

    def test_failed_job_still_fetches_closed_logs(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);(base/'records').mkdir();dest=base/'job';partial=base/'job.partial';check=self.bundle(partial,'FAILED')
            manifest=read(partial/'results_manifest.json')
            with patch.object(cloud_run,'BASE',base),patch.object(cloud_run,'local_environment'),patch.object(cloud_run,'result_location',return_value=('/remote/job',dest)),patch.object(cloud_run,'result_probe',return_value={'check':check,'manifest':manifest}),patch.object(cloud_run.subprocess,'run',return_value=subprocess.CompletedProcess([],0,'','')):
                result=cloud_run.fetch('job')
            self.assertEqual(result['run_state'],'FAILED');self.assertTrue((dest/'console.log').is_file());self.assertFalse(partial.exists())

    def test_space_guards_include_growth(self):
        self.assertFalse(disk_allowed(7*1024**3,'build',3*1024**3));self.assertTrue(disk_allowed(9*1024**3,'build',3*1024**3))
        self.assertFalse(disk_allowed(3*1024**3-1,'run'))

    def test_claimed_pass_with_empty_outputs_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'B_liquid').mkdir()
            write(root/'run_status.json',{'status':'CLOUD_SMOKE_PASS'})
            write(root/'B_liquid/completion.json',{'actual_steps':100,'finite':True,'particles':512})
            (root/'B_liquid/particles_final.csv').write_text('id,x,y,z,vx,vy,vz\n')
            seal(root,'job','COMPLETED')
            with self.assertRaisesRegex(ValueError,'EMPTY_PARTICLE_CSV'):validate_result(root)

    def test_timeout_and_disk_limit_stop_only_owned_process(self):
        for seconds,free,expected in [(.3,10*1024**3,'TIMEOUT'),(2,1024,'DISK_LIMIT')]:
            p=subprocess.Popen([sys.executable,'-B','-c','import time;time.sleep(10)'],start_new_session=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
            result=wait_guarded(p,seconds,lambda:free)
            self.assertEqual(result['stop_reason'],expected);self.assertIsNotNone(p.poll());self.assertNotEqual(result['exit_code'],0)

if __name__=='__main__':unittest.main()
