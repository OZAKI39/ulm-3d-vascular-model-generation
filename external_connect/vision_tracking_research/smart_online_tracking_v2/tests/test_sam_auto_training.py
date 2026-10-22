"""Automatic-training boundaries; all generated labels are temporary software fixtures."""
import unittest,tempfile,json,sys,os,hashlib,time,multiprocessing
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sam_auto_training import (snapshot,snapshot_still_valid,config,set_enabled,atomic_json,
                              admission,register_client,unregister_client,has_clients,active_model,
                              attempt_blocks_retraining,rollback,gpu_slot)
from annotation_store import CLASSES,sha

class AutoTraining(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory(prefix='SAM_SOFTWARE_FIXTURE_');self.root=Path(self.temp.name)
  (self.root/'data').mkdir();self.records=[]
 def tearDown(self):self.temp.cleanup()
 def frame(self,name,split,confirmed=True,reviewer='Fixture reviewer'):
  image=self.root/'data'/f'{name}.bin';image.write_bytes(b'SOFTWARE_FIXTURE_NOT_IMAGE')
  r={'id':name,'split':split,'image':str(image.relative_to(self.root)),'image_sha256':sha(image),'resolution_wh':[8,8],'source_sha256':'fixture'}
  self.records.append(r);folder=self.root/'data/annotations'/name;folder.mkdir(parents=True)
  arrays={k:np.zeros((8,8),np.uint32) for k in CLASSES};arrays['SINGLE_BUBBLE'][2:4,2:4]=1
  np.savez(folder/'v1.npz',**arrays)
  m={'candidate_id':name,'split':split,'image_sha256':r['image_sha256'],'revision':'v1','mask_file':'v1.npz',
     'mask_sha256':sha(folder/'v1.npz'),'expert_confirmed':confirmed,'reviewer':reviewer,'notes':'',
     'status':'EXPERT_CONFIRMED' if confirmed else 'DRAFT_UNREVIEWED','review_scope':'FULL_FRAME',
     'instance_counts':{'SINGLE_BUBBLE':1}}
  atomic_json(folder/'latest.json',m);atomic_json(self.root/'data/candidate_manifest.json',{'rows':self.records});return folder
 def test_default_enabled_explicit_opt_out(self):
  self.assertTrue(config(self.root)['enabled']);set_enabled(False,self.root);self.assertFalse(config(self.root)['enabled'])
  set_enabled(True,self.root);self.assertTrue(config(self.root)['enabled'])
 def test_only_human_confirmed_train_validation_and_no_test_reads(self):
  self.frame('train','train');self.frame('val','validation');self.frame('draft','train',False)
  self.frame('software','train',True,'SOFTWARE_TEST')
  p=self.frame('test','test');(p/'v1.npz').unlink()
  self.frame('quarantine','QUARANTINE_UNASSIGNED')
  s=snapshot(self.root);self.assertEqual({r['record']['id'] for r in s['rows']},{'train','val'})
  self.assertEqual(s['excluded'],[])
 def test_draft_revision_invalidates_frozen_snapshot(self):
  p=self.frame('train','train');s=snapshot(self.root);self.assertTrue(snapshot_still_valid(s,self.root))
  m=json.loads((p/'latest.json').read_text());m['expert_confirmed']=False;atomic_json(p/'latest.json',m)
  self.assertFalse(snapshot_still_valid(s,self.root));self.assertEqual(snapshot(self.root)['rows'],[])
 def test_checksum_mismatch_excluded_and_fingerprint_tracks_data(self):
  p=self.frame('train','train');s=snapshot(self.root);self.assertEqual(s['fingerprint'],snapshot(self.root)['fingerprint'])
  (p/'v1.npz').write_bytes(b'corruption');self.assertFalse(snapshot_still_valid(s,self.root))
  new=snapshot(self.root);self.assertNotEqual(s['fingerprint'],new['fingerprint']);self.assertEqual(len(new['excluded']),1)
 def test_reconfirm_same_content_does_not_retrain(self):
  p=self.frame('train','train');s=snapshot(self.root)
  m=json.loads((p/'latest.json').read_text());m['revision']='v2_same_content';atomic_json(p/'latest.json',m)
  self.assertEqual(s['fingerprint'],snapshot(self.root)['fingerprint']);self.assertTrue(snapshot_still_valid(s,self.root))
 def test_no_independent_validation_no_promotion(self):
  self.assertFalse(admission({'validation_frames':0,'validation_objects':0})[0])
  self.assertFalse(admission({'validation_frames':1,'validation_objects':100})[0])
 def test_regression_and_clear_gain_admission(self):
  m={'validation_frames':2,'validation_objects':20,'finite':True,'candidate_mean_iou':.8,
     'original_mean_iou':.6,'active_mean_iou':.7,'worst_frame_delta':.05,
     'candidate_failure_fraction':.1,'active_failure_fraction':.2}
  self.assertTrue(admission(m)[0]);self.assertFalse(admission({**m,'worst_frame_delta':-.1})[0])
  self.assertFalse(admission({**m,'candidate_failure_fraction':.3})[0])
  self.assertFalse(admission({**m,'candidate_mean_iou':.701})[0])
 def test_process_client_identity_and_original_default(self):
  self.assertFalse(has_clients(self.root));register_client(runtime=self.root);self.assertTrue(has_clients(self.root))
  unregister_client(runtime=self.root);self.assertFalse(has_clients(self.root));self.assertEqual(active_model(self.root)['id'],'ORIGINAL_SAM_VIT_B')
 def test_pause_is_resumable_but_completed_data_not_retrained(self):
  self.assertFalse(attempt_blocks_retraining('a',{'fingerprint':'a','status':'CANCELLED'}))
  self.assertTrue(attempt_blocks_retraining('a',{'fingerprint':'a','status':'COMPLETED'}))
  self.assertFalse(attempt_blocks_retraining('b',{'fingerprint':'a','status':'COMPLETED'}))
 def test_manual_rollback(self):
  atomic_json(self.root/'active_model.json',{'id':'SOFTWARE_FIXTURE','admitted':True})
  rollback(self.root);self.assertFalse(active_model(self.root)['admitted'])
 def test_interactive_request_precedes_next_training_step(self):
  ctx=multiprocessing.get_context('fork');queue=ctx.Queue();release=ctx.Event();entered=ctx.Event()
  def first_training():
   with gpu_slot(runtime=self.root):entered.set();release.wait(5)
  def request(interactive):
   with gpu_slot(runtime=self.root,interactive=interactive):queue.put('interactive' if interactive else 'training');time.sleep(.05)
  a=ctx.Process(target=first_training);a.start();self.assertTrue(entered.wait(3))
  b=ctx.Process(target=request,args=(True,));b.start()
  until=time.monotonic()+3
  while not list((self.root/'interactive').glob('*')) and time.monotonic()<until:time.sleep(.01)
  self.assertTrue(list((self.root/'interactive').glob('*')))
  c=ctx.Process(target=request,args=(False,));c.start();release.set()
  try:self.assertEqual(queue.get(timeout=5),'interactive');self.assertEqual(queue.get(timeout=5),'training')
  finally:
   for p in (a,b,c):p.join(5)
   for p in (a,b,c):
    if p.is_alive():p.kill();p.join()

if __name__=='__main__':unittest.main(verbosity=2)
