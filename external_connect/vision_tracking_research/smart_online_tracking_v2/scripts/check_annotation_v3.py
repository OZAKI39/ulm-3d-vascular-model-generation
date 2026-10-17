"""Actual Qt regression test; writes only an isolated synthetic-label sandbox."""
from pathlib import Path
import sys, json, time, argparse, traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from annotate import Workbench, QApplication, QTimer, QFont, napari
from annotation_store import load_annotation
from qtpy.QtTest import QTest
from qtpy.QtCore import Qt,QPoint
OUT=ROOT/'evidence/annotation_repair_v3_20261016'
app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
w=Workbench(argparse.Namespace(index=0,manual_only=False,smoke=True,annotation_dir=str(OUT/'ui_sandbox')))
report={'research_annotations_modified':False,'checks':{},'PASS':False};beats=[]
timer=QTimer();timer.timeout.connect(lambda:beats.append(time.monotonic()));timer.start(100)

def check(name,condition):
    report['checks'][name]=bool(condition)
    if not condition:raise AssertionError(name)

def proposal(mask):
    w.sam_masks=mask[None].copy();w.sam_details={'proposals':[{'allowed':True,'area_px':int(mask.sum()),'reasons':[]}]}
    w.preview_signature=w.prompt_signature();w.last_sam_wall_s=0
    w.proposal.choices=[('SOFTWARE FIXTURE',0)];w.proposal.value=0;w.show_proposal()

def run():
 try:
    for layer in w.layers.values():layer.data=np.zeros_like(layer.data)
    layer=w.layers['SINGLE_BUBBLE'];layer.data_setitem((np.array([100,100,120]),np.array([100,101,120])),1)
    before=layer.data.copy();history=len(layer._undo_history)
    w.viewer.layers.remove(layer);QApplication.processEvents()
    check('deleted_class_layer_reinserted_same_object',w.layers['SINGLE_BUBBLE'] is layer and layer in w.viewer.layers)
    check('deleted_layer_pixels_and_undo_retained',np.array_equal(before,layer.data) and len(layer._undo_history)==history)
    w.viewer.layers.remove(w.image);w.viewer.layers.remove(w.preview);w.viewer.layers.remove(w.positive);QApplication.processEvents()
    check('all_required_layers_restored_in_correct_order',list(w.viewer.layers)==w.required_layers())
    w.category.value='SINGLE_BUBBLE';w.instance.value=1;w.select_instance();w.track_ids['SINGLE_BUBBLE:1']='FIXTURE_A'
    w.delete_instance();check('delete_disconnected_instance_with_preserve_labels',not np.any(layer.data))
    w.undo();check('delete_undo_pixels_and_track_id',np.array_equal(layer.data,before) and w.track_ids['SINGLE_BUBBLE:1']=='FIXTURE_A')
    w.redo();check('delete_redo',not np.any(layer.data) and not w.track_ids)
    mask=np.zeros_like(layer.data,bool);mask[200:205,220:226]=True
    proposal(mask);w.accept_sam();check('sam_accept_selects_accepted_id',layer.selected_label==w.instance.value==1)
    canvas=w.viewer.window._qt_viewer.canvas.native;canvas.setFocus()
    QTest.keyClick(canvas,Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier)
    check('sam_accept_keyboard_undo',not np.any(layer.data))
    QTest.keyClick(canvas,Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier|Qt.KeyboardModifier.ShiftModifier)
    check('sam_accept_keyboard_redo',np.array_equal(layer.data>0,mask))
    w.remove_last_sam();check('one_click_remove_last_sam',not np.any(layer.data));check('deleted_sam_pointer_cleared',w.last_sam_instance is None);w.undo()
    w.pick_mode()
    class E:position=(202,223)
    list(w.pick_callback(w.viewer,E()))
    check('click_pick_existing_instance',w.instance.value==1 and not w.picking)
    w.delete_instance();check('picked_instance_deleted',not np.any(layer.data))
    proposal(mask);w.positive.data=np.array([[202.,223.]])
    try:w.accept_sam();raise AssertionError('stale accepted')
    except ValueError:pass
    check('changed_prompt_rejects_stale_candidate',not np.any(layer.data))
    w.clear_prompts();proposal(mask);w.sam_details['proposals'][0]['allowed']=False;w.sam_details['proposals'][0]['reasons']=['TEST channel']
    try:w.accept_sam();raise AssertionError('oversize accepted')
    except ValueError:pass
    check('rejected_proposal_cannot_be_accepted',not np.any(layer.data))
    w.clear_prompts();w.notes.value='SOFTWARE_FIXTURE; never human truth';w.save(False)
    w.go(1);w.go(0);check('saved_edits_roundtrip',not np.any(w.layers['SINGLE_BUBBLE'].data))
    w.box_mode();canvas=w.viewer.window._qt_viewer.canvas.native
    QTest.mousePress(canvas,Qt.MouseButton.LeftButton,pos=QPoint(400,300))
    QTest.mouseMove(canvas,QPoint(460,355),delay=30)
    QTest.mouseRelease(canvas,Qt.MouseButton.LeftButton,pos=QPoint(460,355))
    QApplication.processEvents()
    check('native_rectangle_draw_by_mouse',len(w.box_layer.data)==1 and np.ptp(w.box_layer.data[0],axis=0).min()>2)
    w.clear_prompts()
    w.box_layer.data=[np.array([[382,759],[382,785],[408,785],[408,759]],float)]
    w.zoom_target();check('target_zoom',w.viewer.scene.camera.zoom>2)
    w.positive.data=np.array([[395.,772.]])
    w.tabs.setCurrentIndex(1);beats.clear();w.segment();QTimer.singleShot(200,poll_first)
 except Exception:fail()

def poll_first():
 if w.job is not None:QTimer.singleShot(100,poll_first);return
 try:
    check('real_local_sam_completed',w.last_sam_outcome=='PASS')
    check('real_candidate_has_admitted_option',any(v['allowed'] for v in w.sam_details['proposals']))
    report['first_request_wall_s']=w.last_sam_wall_s;report['first_proposals']=w.sam_details
    report['sam_event_loop_max_gap_s']=float(max(np.diff(beats)));check('ui_responsive_during_sam',report['sam_event_loop_max_gap_s']<2)
    w.viewer.screenshot(path=str(OUT/'sam_local_window.png'),canvas_only=False)
    worker=w.proc;report['worker_pid']=int(worker.processId());w.segment()
    QTimer.singleShot(100,lambda:poll_warm(worker))
 except Exception:fail()

def poll_warm(worker):
 if w.job is not None:QTimer.singleShot(100,lambda:poll_warm(worker));return
 try:
    check('worker_and_embedding_reused',w.proc is worker and w.sam_details['embedding_cache_hit'])
    report['warm_request_wall_s']=w.last_sam_wall_s
    w.accept_sam();check('actual_sam_accept_then_remove',np.any(w.layers['SINGLE_BUBBLE'].data))
    w.remove_last_sam();check('actual_sam_removed',not np.any(w.layers['SINGLE_BUBBLE'].data));w.undo()
    check('actual_sam_delete_undo',np.any(w.layers['SINGLE_BUBBLE'].data))
    old_worker=w.proc;w.box_mode();check('new_box_keeps_idle_worker',w.proc is old_worker)
    w.box_layer.data=[np.array([[430,737],[430,755],[448,755],[448,737]],float)]
    w.segment();QTimer.singleShot(100,poll_new_target)
 except Exception:fail()

def poll_new_target():
 if w.job is not None:QTimer.singleShot(100,poll_new_target);return
 try:
    check('box_only_new_target_uses_same_model',w.last_sam_outcome=='PASS' and not w.sam_details['cold_model'] and not w.sam_details['embedding_cache_hit'])
    report['new_target_request_wall_s']=w.last_sam_wall_s
    w.clear_prompts();w.positive.data=np.array([[395.,772.]]);w.segment();w.clear_prompts()
    check('clear_cancels_and_invalidates_job',w.job is None and w.sam_masks is None and w.proc is None)
    report['whole_test_event_loop_max_gap_s_including_synthetic_batches']=float(max(np.diff(beats)))
    w.save(False);before={k:v.data.copy() for k,v in w.layers.items()};w.go(1);w.go(0)
    check('actual_sam_edit_roundtrip',all(np.array_equal(v,w.layers[k].data) for k,v in before.items()))
    w.viewer.window._qt_window.close();report['PASS']=True;finish()
 except Exception:fail()

def fail():
 report['error']=traceback.format_exc();print(report['error']);w.cancel_sam();w.closing=True;w.viewer.close();finish()
def finish():
 timer.stop();(OUT/'ui_regression.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True);app.quit()
QTimer.singleShot(1500,run);napari.run()
if not report['PASS']:raise SystemExit(1)
