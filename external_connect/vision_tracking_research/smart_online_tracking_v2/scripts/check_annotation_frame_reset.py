"""Real Qt reset/recovery checks in an isolated, non-training annotation folder."""
from pathlib import Path
import argparse,json,sys,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from annotate import Workbench,QApplication,QTimer,QFont,napari
from qtpy.QtWidgets import QMessageBox
from qtpy.QtCore import QProcess
from qtpy.QtTest import QTest
from annotation_store import CLASSES,sha,load_annotation

OUT=ROOT/'evidence/annotation_frame_reset_20261016';OUT.mkdir(exist_ok=True)
RUN=OUT/f'gui_{time.time_ns()}';RUN.mkdir()
app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
args=argparse.Namespace(index=0,smoke=True,manual_only=True,no_auto_train=True,annotation_dir=str(RUN/'annotations'))
w=Workbench(args);report={'software_fixture':True,'checks':{},'PASS':False};others=[]
def hashes(folder):return {str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}
formal_before=hashes(ROOT/'data/annotations')

def check(name,value):
 report['checks'][name]=bool(value)
 if not value:raise AssertionError(name)

def click_dialog(button,accept):
 def answer():
  dialog=QApplication.activeModalWidget()
  if not isinstance(dialog,QMessageBox):QTimer.singleShot(25,answer);return
  role=QMessageBox.ButtonRole.AcceptRole if accept else QMessageBox.ButtonRole.RejectRole
  next(b for b in dialog.buttons() if dialog.buttonRole(b)==role).click()
 QTimer.singleShot(25,answer);w.buttons[button].native.click();QApplication.processEvents()

def run():
 global w
 try:
  check('default_alpha_18_percent',w.label_opacity.value==18 and all(v.opacity==.18 for v in [*w.layers.values(),w.preview,w.bulk_preview]))
  check('no_recovery_without_backup',not w.buttons['restore_frame'].enabled)
  # Preserve another frame as an independent fixture.
  w.go(1);w.layers['SINGLE_BUBBLE'].data_setitem((np.array([30]),np.array([30])),17)
  w.notes.value='SOFTWARE_TEST_OTHER_FRAME';w.save(False)
  other_dir=w.directory/w.record['id'];other_before=hashes(other_dir);w.go(0)
  for i,layer in enumerate(w.layers.values()):
   ys,xs=np.mgrid[50+i*20:60+i*20,80:90]
   layer.data_setitem((ys.ravel(),xs.ravel()),i+1)
  w.notes.value='SOFTWARE_TEST_SAVED';w.save(False)
  folder=w.directory/w.record['id'];history_before={p.name:sha(p) for p in folder.glob('*.*') if p.name!='latest.json'}
  # Include unsaved labels, drafts, prompts and histories in the reset scenario.
  w.layers['SINGLE_BUBBLE'].data_setitem((np.array([300,301]),np.array([300,301])),6)
  w.notes.value='SOFTWARE_TEST_UNSAVED';w.track_ids={'SINGLE_BUBBLE:1':'fixture_track'}
  w.proposal_sources=[{'status':'SOFTWARE_FIXTURE_SOURCE'}]
  w.roi_layer.data=[np.array([[40,40],[40,350],[350,350],[350,40]])]
  w.positive.data=np.array([[60.,60.]]);w.negative.data=np.array([[90.,90.]])
  w.box_layer.data=[np.array([[50,50],[50,100],[100,100],[100,50]])]
  w.split_points.data=np.array([[80.,80.],[85.,85.]])
  w.preview.data_setitem((np.array([60]),np.array([60])),1)
  w.bulk_preview.data_setitem((np.array([80]),np.array([80])),1)
  w.layers['AGGREGATE']._overlays['polygon'].points=[np.array([200.,200.]),np.array([200.,210.])]
  original={k:v.data.copy() for k,v in w.layers.items()};dirty=w.dirty
  w.label_opacity.value=12
  check('all_annotation_and_preview_alpha_adjusted',all(v.opacity==.12 for v in [*w.layers.values(),w.preview,w.bulk_preview]))
  check('opacity_does_not_change_pixels_or_dirty',all(np.array_equal(original[k],v.data) for k,v in w.layers.items()) and w.dirty==dirty)
  click_dialog('reset_frame',False)
  check('cancel_reset_changes_nothing',all(np.array_equal(original[k],v.data) for k,v in w.layers.items())
        and len(w.positive.data)==1 and not (w.reset_backup_dir()/'latest.json').exists())
  # A backup failure must leave all labels in place.
  snapshot=w.reset_snapshot
  def fail_snapshot(*_,**__):raise OSError('SOFTWARE_TEST simulated backup write failure')
  w.reset_snapshot=fail_snapshot;click_dialog('reset_frame',True);w.reset_snapshot=snapshot
  check('backup_failure_prevents_reset',all(np.array_equal(original[k],v.data) for k,v in w.layers.items())
        and 'backup write failure' in w.status.value)
  # Simulated in-flight processes verify cancellation, without invoking a model.
  job_dir=RUN/'pending_jobs';job_dir.mkdir()
  sam_proc=QProcess();sam_proc.start(sys.executable,['-c','import time; time.sleep(30)']);sam_proc.waitForStarted(1000)
  bulk_proc=QProcess();bulk_proc.start(sys.executable,['-c','import time; time.sleep(30)']);bulk_proc.waitForStarted(1000)
  w.proc=sam_proc;w.job={'output':str(job_dir/'sam.npz'),'image_id':w.record['id']}
  w.bulk_proc=bulk_proc;w.bulk_job={'folder':job_dir,'started':time.monotonic()}
  click_dialog('reset_frame',True)
  check('inflight_jobs_cancelled',w.proc is None and w.job is None and w.bulk_proc is None and w.bulk_job is None)
  check('all_five_classes_empty',all(not np.any(v.data) for v in w.layers.values()))
  check('all_previews_empty',not np.any(w.preview.data) and not np.any(w.bulk_preview.data) and w.sam_masks is None and w.bulk_result is None)
  check('all_prompt_ROI_and_split_geometry_empty',all(len(v.data)==0 for v in [w.positive,w.negative,w.box_layer,w.roi_layer,w.split_points]))
  check('pending_polygons_cleared',all(not v._overlays['polygon'].points for v in w.layers.values()))
  check('notes_tracks_review_and_sources_cleared',w.notes.value=='' and w.track_ids=={} and w.review_regions==[] and w.proposal_sources==[] and not w.reviewed.value)
  check('fresh_single_bubble_id_one',w.category.value=='SINGLE_BUBBLE' and w.layers['SINGLE_BUBBLE'].selected_label==1 and w.layers['SINGLE_BUBBLE'].mode=='paint')
  meta,masks=load_annotation(w.directory,w.record)
  check('reset_saved_immediately_as_empty_unreviewed',meta['status']=='DRAFT_UNREVIEWED' and not meta['expert_confirmed'] and all(not np.any(v) for v in masks.values()))
  check('previous_revisions_retained',all(sha(folder/name)==digest for name,digest in history_before.items()))
  state,backup=w.read_reset_snapshot()
  check('backup_includes_unsaved_pixels_and_notes',all(np.array_equal(backup[k],original[k]) for k in CLASSES) and state['notes']=='SOFTWARE_TEST_UNSAVED')
  check('backup_captures_pending_polygon',state['pending_polygons']['AGGREGATE']==[[200.,200.],[200.,210.]])
  for layer in [*w.layers.values(),w.preview,w.bulk_preview]:layer.undo();layer.redo()
  check('old_undo_histories_cannot_resurrect_traces',all(not np.any(v.data) for v in [*w.layers.values(),w.preview,w.bulk_preview]))
  pointer_hash=sha(w.reset_backup_dir()/'latest.json');click_dialog('reset_frame',True)
  check('repeated_empty_reset_preserves_recovery',pointer_hash==sha(w.reset_backup_dir()/'latest.json'))
  w.go(1)
  check('other_frame_has_no_wrong_recovery',not w.buttons['restore_frame'].enabled)
  check('other_frame_unchanged',hashes(other_dir)==other_before)
  w.go(0);check('reset_survives_frame_switch',all(not np.any(v.data) for v in w.layers.values()))
  old=w;old.viewer.window._qt_window.close()
  args.index=None;w=Workbench(args);others.append(old)
  check('reset_and_recovery_survive_restart',all(not np.any(v.data) for v in w.layers.values()) and w.buttons['restore_frame'].enabled)
  check('opacity_preference_restored',w.label_opacity.value==12 and all(v.opacity==.12 for v in [*w.layers.values(),w.preview,w.bulk_preview]))
  click_dialog('restore_frame',False)
  check('cancel_recovery_preserves_empty_frame',all(not np.any(v.data) for v in w.layers.values()))
  # New work after reset is separately backed up before recovery replaces it.
  w.layers['SINGLE_BUBBLE'].data_setitem((np.array([500]),np.array([500])),77)
  backup_count=len(list(w.reset_backup_dir().glob('*.npz')))
  click_dialog('restore_frame',True)
  check('recovery_exact_pixels_and_metadata',all(np.array_equal(v.data,original[k]) for k,v in w.layers.items())
        and w.notes.value=='SOFTWARE_TEST_UNSAVED' and w.track_ids=={'SINGLE_BUBBLE:1':'fixture_track'}
        and w.current_roi()==[40,40,350,350] and w.proposal_sources==[{'status':'SOFTWARE_FIXTURE_SOURCE'}])
  check('new_post_reset_work_backed_up',len(list(w.reset_backup_dir().glob('*.npz')))==backup_count+1)
  check('recovery_pointer_still_original',pointer_hash==sha(w.reset_backup_dir()/'latest.json'))
  meta,masks=load_annotation(w.directory,w.record)
  check('recovery_saved_unreviewed_not_training_truth',meta['status']=='DRAFT_UNREVIEWED' and meta['review_regions']==[] and not meta['expert_confirmed'])
  check('recovery_does_not_restore_candidates_or_pending_polygon',not np.any(w.preview.data)
        and not np.any(w.bulk_preview.data) and not len(w.positive.data)
        and all(not v._overlays['polygon'].points for v in w.layers.values()))
  for width in (440,500):
   w.tabs.setFixedWidth(width);QApplication.processEvents();QTest.qWait(50)
   a=w.buttons['reset_frame'].native.geometry();b=w.buttons['restore_frame'].native.geometry()
   check(f'reset_buttons_aligned_{width}',a.y()==b.y() and abs(a.width()-b.width())<=1)
   for key in ('reset_frame','restore_frame'):
    button=w.buttons[key];check(f'button_text_fits_{key}_{width}',button.native.fontMetrics().horizontalAdvance(button.text)+14<=button.native.width())
  w.label_opacity.value=18
  w.main.native.parentWidget().parentWidget().ensureWidgetVisible(w.buttons['reset_frame'].native)
  QApplication.processEvents();w.viewer.screenshot(path=str(RUN/'reset_and_opacity.png'),canvas_only=False)
  check('formal_annotations_unchanged',hashes(ROOT/'data/annotations')==formal_before)
  check('test_training_disabled',not w.auto_training_allowed)
  report['PASS']=True
 except Exception:
  report['error']=traceback.format_exc();print(report['error'],flush=True)
 finally:
  w.cancel_sam();w.cancel_bulk();w.viewer.window._qt_window.close()
  report['formal_annotations_before']=formal_before;report['formal_annotations_after']=hashes(ROOT/'data/annotations')
  (RUN/'report.json').write_text(json.dumps(report,indent=2));(OUT/'latest_gui_run.json').write_text(json.dumps({'run':str(RUN),'PASS':report['PASS']},indent=2))
  print(json.dumps({'PASS':report['PASS'],'checks':len(report['checks']),'run':str(RUN),'error':report.get('error')}),flush=True)
  app.quit()

QTimer.singleShot(1500,run);napari.run()
if not report['PASS']:raise SystemExit(1)
