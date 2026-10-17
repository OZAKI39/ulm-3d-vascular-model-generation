"""Actual Qt + GPU workflow smoke checks, isolated from formal human labels."""
from pathlib import Path
import argparse,json,sys,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from annotate import Workbench,QApplication,QTimer,QFont,napari
from qtpy.QtCore import Qt,QPoint
from qtpy.QtTest import QTest
from annotation_store import load_annotation
OUT=ROOT/'evidence/annotation_batch_roi_20261016';OUT.mkdir(exist_ok=True)
RUN=OUT/f'gui_{time.time_ns()}';RUN.mkdir()
records=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows']
index=next(i for i,r in enumerate(records) if r['id']=='video_000054')
app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
w=Workbench(argparse.Namespace(index=index,smoke=True,manual_only=False,no_auto_train=True,annotation_dir=str(RUN/'annotations')))
report={'software_fixture':True,'formal_labels_modified':False,'checks':{},'PASS':False}
ticks=[];heartbeat=QTimer();heartbeat.timeout.connect(lambda:ticks.append(time.monotonic()));heartbeat.start(50)
start=time.monotonic()

def check(name,value):
 report['checks'][name]=bool(value)
 if not value:raise AssertionError(name)

def click(key):
 w.buttons[key].native.click();QApplication.processEvents()

def rectangle(box):
 y0,x0,y1,x1=box;w.roi_layer.data=[np.array([[y0,x0],[y0,x1],[y1,x1],[y1,x0]])]

def layout_checks():
 rows=[('roi','zoom_roi'),('generate_bulk','cancel_bulk'),('drop_candidate','undo_candidate'),('incorporate','clear_bulk'),('brush','erase'),('pick','delete_main'),('start_split','finish_split'),('undo','redo'),('save','confirm'),('prev','next')]
 geom={}
 for width in (440,500):
  w.tabs.setFixedWidth(width);QApplication.processEvents();QTest.qWait(100)
  pairs=[]
  for a,b in rows:
   ga=w.buttons[a].native.geometry();gb=w.buttons[b].native.geometry()
   check(f'aligned_{a}_{width}',abs(ga.width()-gb.width())<=1 and ga.y()==gb.y() and ga.height()==gb.height()==32)
   pairs.append([ga.x(),gb.x(),ga.width(),gb.width()])
  check(f'common_columns_{width}',len(set(tuple(x) for x in pairs))==1)
  for name,widget in w.buttons.items():
   # Keep all controls readable, including the less commonly used SAM tab.
   if widget.native.isVisible():
    check(f'text_fits_{name}_{width}',widget.native.fontMetrics().horizontalAdvance(widget.text)+14<=widget.native.width())
  for tab in (1,2):
   w.tabs.setCurrentIndex(tab);QApplication.processEvents();QTest.qWait(60)
   for name,widget in w.buttons.items():
    if widget.native.isVisible():
     check(f'text_fits_{name}_{width}',widget.native.fontMetrics().horizontalAdvance(widget.text)+14<=widget.native.width())
  w.tabs.setCurrentIndex(0);QApplication.processEvents()
  geom[str(width)]=pairs
 report['button_geometry']=geom
 w.tabs.setMinimumWidth(440);w.tabs.setMaximumWidth(500)

def first():
 try:
  check('default_region_and_paint',w.review_scope.value=='REGIONS' and w.layers['SINGLE_BUBBLE'].mode=='paint')
  check('checkbox_matches_scope','框内' in w.reviewed.native.text())
  layout_checks()
  canvas=w.viewer.window._qt_viewer.canvas.native
  w.viewer.window._qt_window.activateWindow();w.viewer.window._qt_window.raise_();QTest.qWait(100)
  click('roi');canvas.setFocus()
  center=QPoint(canvas.width()//2,canvas.height()//2)
  QTest.mousePress(canvas,Qt.MouseButton.LeftButton,Qt.KeyboardModifier.NoModifier,center-QPoint(70,70))
  QTest.mouseMove(canvas,center+QPoint(70,70),100)
  QTest.mouseRelease(canvas,Qt.MouseButton.LeftButton,Qt.KeyboardModifier.NoModifier,center+QPoint(70,70))
  QApplication.processEvents()
  check('mouse_drag_creates_review_rectangle',len(w.roi_layer.data)==1 and len(w.current_roi())==4)
  rectangle([150,900,406,1156]);click('zoom_roi')
  click('generate_bulk');check('async_job_started',w.bulk_job is not None)
  QTimer.singleShot(150,poll)
 except Exception:fail()

def poll():
 if w.bulk_job is not None:
  if time.monotonic()-start>180:
   report['error']='GUI inference timeout';finish();return
  QTimer.singleShot(100,poll);return
 try:
  check('real_gpu_proposals_loaded',w.bulk_result is not None and len(w.bulk_result['objects'])>0)
  report['real_proposal_count']=len(w.bulk_result['objects']);report['real_inference']=w.bulk_result
  check('preview_not_saved_truth',not np.any(w.layers['SINGLE_BUBBLE'].data))
  w.viewer.screenshot(path=str(RUN/'bulk_candidates.png'),canvas_only=False)
  ident=int(np.max(w.bulk_preview.data));old=w.bulk_preview.data.copy()
  click('drop_candidate')
  yy,xx=np.argwhere(old==ident)[0]
  event=argparse.Namespace(position=(float(yy),float(xx)))
  for _ in w.bulk_pick_callback(w.viewer,event):pass
  check('delete_candidate_only',not np.any(w.bulk_preview.data==ident) and not np.any(w.layers['SINGLE_BUBBLE'].data))
  click('undo_candidate');check('candidate_delete_undo_exact',np.array_equal(old,w.bulk_preview.data))
  # Existing human-like fixture pixels survive; overlapping candidate is skipped whole.
  layer=w.layers['SINGLE_BUBBLE'];layer.data_setitem((np.array([yy]),np.array([xx])),900)
  before=layer.data.copy();click('incorporate');accepted=layer.data.copy()
  check('incorporate_preserves_existing_pixels',accepted[yy,xx]==900 and np.count_nonzero(accepted)>1)
  check('incorporate_separate_ids',len(np.unique(accepted))>3 and not np.any(w.bulk_preview.data))
  click('undo');check('bulk_one_step_undo',np.array_equal(before,layer.data))
  click('redo');check('bulk_one_step_redo',np.array_equal(accepted,layer.data))
  check('provenance_recorded',bool(w.proposal_sources and w.proposal_sources[0]['model_sha256']))
  # A joined software mask for the two-seed exact-union test.
  layer.data_setitem((np.repeat(np.arange(420,432),24),np.tile(np.arange(1000,1024),12)),950)
  w.instance.value=950;click('start_split');w.split_points.data=np.array([[425,1005],[425,1019]])
  before=layer.data.copy();click('finish_split');after=layer.data.copy()
  check('split_two_ids_exact_union',np.array_equal(before>0,after>0) and set(np.unique(after[420:432,1000:1024]))=={950,951})
  click('undo');check('split_undo_exact',np.array_equal(before,layer.data));click('redo')
  check('split_redo_exact',np.array_equal(after,layer.data))
  w.notes.value='SOFTWARE_UI_FIXTURE: regional workflow, never human truth';click('save')
  meta,masks=load_annotation(w.directory,w.record)
  check('saved_region_workspace_not_confirmed',meta['workspace_roi_yxyx']==[150,900,406,1156] and not meta['expert_confirmed'] and meta['review_regions']==[])
  w.reviewed.value=True;check('software_confirmation_blocked',w.save(True) is False)
  w.go(index+1);w.go(index)
  check('reload_labels_and_roi_exact',np.array_equal(w.layers['SINGLE_BUBBLE'].data,after) and w.current_roi()==[150,900,406,1156])
  check('restart_retains_batch_provenance',len(w.proposal_sources)==1)
  click('generate_bulk');check('second_job_started',w.bulk_job is not None)
  # Switching frame invalidates and kills in-flight work before it can apply.
  w.go(index+1)
  check('frame_change_cancels_and_clears',w.bulk_job is None and w.bulk_result is None and not np.any(w.bulk_preview.data))
  w.go(index);w.zoom_roi();w.new_bubble()
  click('generate_bulk');w.bulk_diameter.value=17
  check('parameter_change_cancels_and_clears',w.bulk_job is None and w.bulk_result is None)
  w.bulk_diameter.value=16
  rectangle([10,10,900,1000]);click('generate_bulk')
  check('oversized_roi_refused',w.bulk_job is None and '32–768' in w.status.value)
  rectangle([150,900,406,1156]);w.zoom_roi();w.new_bubble()
  check('manual_brush_after_batch',w.layers['SINGLE_BUBBLE'].mode=='paint' and not w.bulk_picking)
  check('software_test_no_training',not w.auto_training_allowed)
  report['heartbeat_max_gap_s']=float(np.max(np.diff(ticks)))
  check('responsive_during_inference',report['heartbeat_max_gap_s']<2)
  w.status.value='软件测试：正式人工标签未改动。';w.bulk_status.value='批量纳入后，可按 N 补漏；点选误标后删除。'
  w.viewer.screenshot(path=str(RUN/'workbench.png'),canvas_only=False)
  report['PASS']=True;finish()
 except Exception:fail()

def fail():
 report['error']=traceback.format_exc();print(report['error'],flush=True)
 try:w.viewer.screenshot(path=str(RUN/'failure.png'),canvas_only=False)
 except Exception:pass
 finish()

def finish():
 heartbeat.stop();w.cancel_bulk();w.viewer.window._qt_window.close()
 report['wall_s']=time.monotonic()-start
 (RUN/'report.json').write_text(json.dumps(report,indent=2));(OUT/'latest_gui_run.json').write_text(json.dumps({'run':str(RUN),'PASS':report['PASS']},indent=2))
 print(json.dumps({'PASS':report['PASS'],'run':str(RUN),'error':report.get('error')}),flush=True)
 app.quit()

QTimer.singleShot(1800,first);napari.run()
if not report['PASS']:raise SystemExit(1)
