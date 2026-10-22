"""Replay the actual reported ten points through Qt, in an isolated label sandbox."""
from pathlib import Path
import sys,json,time,argparse,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from annotate import Workbench,QApplication,QTimer,QFont,napari
from qtpy.QtTest import QTest
from qtpy.QtCore import Qt
OUT=ROOT/'evidence/annotation_repair_v4_20261016'
req=json.loads((OUT/'reported_job/request.json').read_text())
records=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows']
idx=next(i for i,r in enumerate(records) if r['id']==req['image_id'])
app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
w=Workbench(argparse.Namespace(index=idx,manual_only=False,smoke=True,annotation_dir=str(OUT/'ui_sandbox')))
report={'research_annotations_modified':False,'checks':{},'PASS':False};beats=[]
timer=QTimer();timer.timeout.connect(lambda:beats.append(time.monotonic()));timer.start(100)
def check(name,value):
 report['checks'][name]=bool(value)
 if not value:raise AssertionError(name)
def run():
 try:
  for layer in w.layers.values():layer.data=np.zeros_like(layer.data)
  w.notes.value='SOFTWARE_UI_FIXTURE; not expert truth'
  w.prompt_history.value=str(ROOT/'results/annotation_jobs'/req['token']/'request.json');w.restore_last_prompts()
  check('reported_ten_prompts_restored',np.array_equal(w.positive.data,np.asarray(req['positive_yx'])))
  check('default_each_point_is_one_bubble',w.sam_semantics.value=='point_instances')
  layer=w.layers['SINGLE_BUBBLE'];w.viewer.layers.remove(layer);QApplication.processEvents()
  check('deleted_layer_restored_same_object',layer in w.viewer.layers and w.layers['SINGLE_BUBBLE'] is layer)
  w.tabs.setCurrentIndex(1);w.zoom_target();beats.clear();w.segment();QTimer.singleShot(100,poll_first)
 except Exception:fail()
def poll_first():
 if w.job is not None:QTimer.singleShot(100,poll_first);return
 try:
  check('real_ten_point_sam_completed',w.last_sam_outcome=='PASS')
  check('ten_independent_targets',len(w.sam_batch)==10)
  labels,items,rejected=w.batch_assembly()
  check('nine_admissible_one_rejected',len(items)==9 and set(rejected)=={6})
  check('separate_preview_ids',len(np.unique(labels))-1==9)
  report['first_wall_s']=w.last_sam_wall_s;report['component_replay_output']=str(OUT/'component_replay/candidate.npz')
  report['event_loop_max_gap_s']=float(max(np.diff(beats)))
  check('ui_responsive_during_real_inference',report['event_loop_max_gap_s']<2)
  w.viewer.screenshot(path=str(OUT/'sam_batch_window.png'),canvas_only=False)
  w.batch_target.value=0;check('per_bubble_candidate_dropdown',w.proposal.value==0)
  w.proposal.value=1;check('candidate_choice_remembered',w.batch_choices[0]==1)
  w.batch_keep.value=False;check('exclude_one_candidate',not w.batch_included[0])
  w.batch_keep.value=True;w.proposal.value=0;w.batch_target.value=-1
  worker=w.proc;w.segment();QTimer.singleShot(100,lambda:poll_warm(worker))
 except Exception:fail()
def poll_warm(worker):
 if w.job is not None:QTimer.singleShot(100,lambda:poll_warm(worker));return
 try:
  check('worker_and_all_ten_embeddings_reused',w.proc is worker and w.sam_details['embedding_cache_hit'])
  report['warm_wall_s']=w.last_sam_wall_s
  layer=w.layers['SINGLE_BUBBLE'];w.accept_sam_batch();accepted=layer.data.copy()
  check('batch_accept_distinct_ids',len(np.unique(accepted))-1==9 and layer.data.max()==9)
  canvas=w.viewer.window._qt_viewer.canvas.native;canvas.setFocus()
  QTest.keyClick(canvas,Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier)
  check('whole_batch_single_keyboard_undo',not np.any(layer.data))
  QTest.keyClick(canvas,Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier|Qt.KeyboardModifier.ShiftModifier)
  check('whole_batch_single_keyboard_redo',np.array_equal(layer.data,accepted))
  w.remove_last_sam();check('remove_whole_batch',not np.any(layer.data));w.undo()
  check('undo_whole_batch_deletion',np.array_equal(layer.data,accepted))
  w.batch_target.value=0
  try:w.accept_sam();raise AssertionError('overlapping instance accepted')
  except ValueError:pass
  check('existing_labels_never_overwritten',np.array_equal(layer.data,accepted))
  w.instance.value=1;w.delete_instance();check('individual_deletion_from_batch',not np.any(layer.data==1))
  w.undo();check('individual_deletion_undo',np.array_equal(layer.data,accepted))
  w.batch_target.value=6
  try:w.accept_sam();raise AssertionError('rejected target accepted')
  except ValueError:pass
  check('failed_target_cannot_accept',np.array_equal(layer.data,accepted))
  w.instance.value=1;w.delete_instance();w.batch_target.value=0;w.accept_sam()
  check('accept_single_target_from_batch',layer.data.max()==10 and len(np.unique(layer.data))-1==9)
  accepted=layer.data.copy();w.remove_last_sam()
  check('remove_single_accept_preserves_other_eight',len(np.unique(layer.data))-1==8)
  w.undo();check('single_target_delete_undo',np.array_equal(layer.data,accepted))
  w.positive.data=np.asarray(req['positive_yx'])+1
  try:w.accept_sam_batch();raise AssertionError('stale accepted')
  except ValueError:pass
  check('changed_points_reject_stale_batch',np.array_equal(layer.data,accepted))
  w.clear_prompts();w.sam_semantics.value='one_object';w.positive.data=np.asarray(req['positive_yx'])
  try:w.segment();raise AssertionError('multiobject singlebubble accepted')
  except ValueError:pass
  check('broad_single_object_prompt_rejected',w.job is None)
  w.clear_prompts();w.save(False);w.go(idx+1);w.go(idx)
  check('batch_save_reopen_exact',np.array_equal(w.layers['SINGLE_BUBBLE'].data,accepted))
  check('software_labels_unconfirmed',not w.meta['expert_confirmed'])
  w.clear_prompts();w.sam_semantics.value='point_instances';w.positive.data=np.array(req['positive_yx'][:1])
  w.segment();w.clear_prompts();check('clear_cancels_pending_batch',w.job is None and w.sam_batch is None and w.proc is None)
  w.viewer.window._qt_window.close();report['PASS']=True;finish()
 except Exception:fail()
def fail():
 report['error']=traceback.format_exc();print(report['error']);w.cancel_sam();w.closing=True;w.viewer.close();finish()
def finish():
 timer.stop();(OUT/'ui_regression.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True);app.quit()
QTimer.singleShot(1500,run);napari.run()
if not report['PASS']:raise SystemExit(1)
