"""Actual Qt + GPU workflow smoke checks, isolated from formal human labels."""
from pathlib import Path
import argparse,json,sys,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from annotate import Workbench,QApplication,QTimer,QFont,napari
from qtpy.QtCore import Qt,QPoint
from qtpy.QtTest import QTest
from annotation_store import load_annotation,sha
from annotation_sizes import size_plan,collect_candidates
parser=argparse.ArgumentParser();parser.add_argument('--output-dir',type=Path);args=parser.parse_args()
OUT=args.output_dir or ROOT/'evidence/annotation_diameter_reference_20261016';OUT.mkdir(parents=True,exist_ok=True)
RUN=OUT/f'gui_{time.time_ns()}';RUN.mkdir()
records=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows']
index=next(i for i,r in enumerate(records) if r['id']=='video_000054')
app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
w=Workbench(argparse.Namespace(index=index,smoke=True,manual_only=False,no_auto_train=True,annotation_dir=str(RUN/'annotations')))
report={'software_fixture':True,'formal_labels_modified':False,'checks':{},'PASS':False}
def formal_hashes():
 return {str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'data/annotations').rglob('*')) if p.is_file()}
report['formal_annotations_before']=formal_hashes()
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
 rows=[('roi','zoom_roi'),('generate_bulk','cancel_bulk'),('erase_candidate','drop_candidate'),('undo_candidate','redo_candidate'),('incorporate','clear_bulk'),('brush','erase'),('pick','delete_main'),('start_split','finish_split'),('undo','redo'),('save','confirm'),('prev','next'),('reset_frame','restore_frame')]
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
  ga=w.bulk_diameter_min.native.geometry();gb=w.bulk_diameter_max.native.geometry()
  check(f'diameter_inputs_aligned_{width}',ga.y()==gb.y() and abs(ga.width()-gb.width())<=1 and ga.x()==pairs[0][0] and gb.x()==pairs[0][1])
  geom[str(width)]=pairs
 report['button_geometry']=geom
 w.tabs.setMinimumWidth(440);w.tabs.setMaximumWidth(500)

def first():
 try:
  check('default_region_and_paint',w.review_scope.value=='REGIONS' and w.layers['SINGLE_BUBBLE'].mode=='paint')
  check('default_diameter_range',w.current_diameter_range()==[8,24])
  check('default_reference_only',not w.bulk_strict_size.value and w.current_diameter_filter()=='REFERENCE_ONLY')
  check('reference_help_visible','仅作参考' in w.diameter_help.value)
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
  w.bulk_diameter_min.value=30;click('generate_bulk')
  check('reversed_range_refused',w.bulk_job is None and '最小直径' in w.status.value)
  w.bulk_diameter_min.value=24;click('generate_bulk')
  check('equal_range_refused',w.bulk_job is None and '最小直径' in w.status.value)
  w.bulk_diameter_min.value=8
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
  check('range_recorded_in_result',w.bulk_result['diameter_range_px']==[8,24])
  check('completed_job_status_cleared',w.status.value.startswith('候选已生成'))
  check('three_inference_scales',np.allclose(w.bulk_result['size_plan']['inference_diameters_px'],[8,np.sqrt(8*24),24]))
  check('reference_policy_recorded',w.bulk_result['size_plan']['diameter_filter']=='REFERENCE_ONLY')
  check('outside_reference_candidates_retained',any(o['outside_reference_range'] for o in w.bulk_result['objects']))
  check('retained_sizes_disclosed','超出参考直径，已保留' in w.bulk_status.value)
  q=json.loads(Path(w.bulk_result['request_path']).read_text());z=np.load(q['output']);cy0,cx0,_,_=w.bulk_result['context_yxyx'];y0,x0,y1,x1=w.bulk_result['roi_yxyx']
  strict=size_plan(dict(q,diameter_filter='STRICT_RANGE'))
  raw=[z[f'raw_labels_s{i:02d}'] for i in range(3)]
  strict_labels,strict_objects,_=collect_candidates(raw,[y0-cy0,x0-cx0,y1-cy0,x1-cx0],q['roi_yxyx'],strict)
  old=ROOT/'results/preannotation_diameter_range_20261016/range_8_24_video_000054/proposals.npz'
  with np.load(old) as prior:check('strict_mode_preserves_previous_labels',np.array_equal(strict_labels,prior['labels']))
  check('strict_mode_filters_outside_range',all(8<=o['equivalent_diameter_px']<=24 for o in strict_objects))
  report['strict_comparison_count']=len(strict_objects)
  for obj in w.bulk_result['objects']:
   raw=z[f'raw_labels_s{obj["scale_index"]:02d}'][y0-cy0:y1-cy0,x0-cx0:x1-cx0]
   check(f'whole_mask_preserved_{obj["id"]}',np.array_equal(z['labels']==obj['id'],raw==obj['raw_id']))
  z.close()
  # A mismatched response must be rejected, even when the image and ROI agree.
  w.bulk_job={'request':q,'signature':w.bulk_prompt_signature()}
  try:
   try:w.receive_bulk(dict(w.bulk_result,size_plan=strict))
   except ValueError as e:check('wrong_policy_response_rejected','筛选方式不一致' in str(e))
   else:check('wrong_policy_response_rejected',False)
  finally:w.bulk_job=None
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
  check('reference_policy_in_draft_provenance',w.proposal_sources[0]['size_plan']['diameter_filter']=='REFERENCE_ONLY')
  # A joined software mask for the two-seed exact-union test.
  joined=int(layer.data.max())+50
  layer.data_setitem((np.repeat(np.arange(420,432),24),np.tile(np.arange(1000,1024),12)),joined)
  w.instance.value=joined;click('start_split');w.split_points.data=np.array([[425,1005],[425,1019]])
  before=layer.data.copy();click('finish_split');after=layer.data.copy()
  check('split_two_ids_exact_union',np.array_equal(before>0,after>0) and set(np.unique(after[420:432,1000:1024]))=={joined,joined+1})
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
  click('generate_bulk');w.bulk_diameter_max.value=30
  check('parameter_change_cancels_and_clears',w.bulk_job is None and w.bulk_result is None)
  w.bulk_diameter_max.value=24
  click('generate_bulk');w.bulk_diameter_min.value=6
  check('minimum_change_cancels_and_clears',w.bulk_job is None and w.bulk_result is None)
  w.bulk_diameter_min.value=8
  before=layer.data.copy();click('generate_bulk');w.bulk_strict_size.value=True
  check('policy_change_cancels_and_clears',w.bulk_job is None and w.bulk_result is None and not np.any(w.bulk_preview.data))
  check('policy_change_preserves_manual_labels',np.array_equal(before,layer.data))
  check('strict_help_visible','严格模式' in w.diameter_help.value)
  rectangle([10,10,900,1000]);click('generate_bulk')
  check('oversized_roi_refused',w.bulk_job is None and '32–768' in w.status.value)
  rectangle([150,900,406,1156]);w.zoom_roi();w.new_bubble()
  check('manual_brush_after_batch',w.layers['SINGLE_BUBBLE'].mode=='paint' and not w.bulk_picking)
  check('software_test_no_training',not w.auto_training_allowed)
  report['heartbeat_max_gap_s']=float(np.max(np.diff(ticks)))
  check('responsive_during_inference',report['heartbeat_max_gap_s']<2)
  w.status.value='软件测试：正式人工标签未改动。';w.bulk_status.value='批量纳入后，可按 N 补漏；点选误标后删除。'
  w.viewer.screenshot(path=str(RUN/'workbench.png'),canvas_only=False)
  w.bulk_diameter_min.value=6;w.bulk_diameter_max.value=30;w.save_session()
  check('range_saved_in_session',json.loads(w.sessionfile.read_text())['bulk_diameter_range_px']==[6,30])
  check('policy_saved_in_session',json.loads(w.sessionfile.read_text())['bulk_diameter_filter']=='STRICT_RANGE')
  again=Workbench(argparse.Namespace(index=None,smoke=True,manual_only=True,no_auto_train=True,annotation_dir=str(w.directory)))
  check('range_restored_on_reopen',again.current_diameter_range()==[6,30])
  check('strict_policy_restored_on_reopen',again.current_diameter_filter()=='STRICT_RANGE')
  again.viewer.window._qt_window.close()
  report['PASS']=True;finish()
 except Exception:fail()

def fail():
 report['error']=traceback.format_exc();print(report['error'],flush=True)
 try:w.viewer.screenshot(path=str(RUN/'failure.png'),canvas_only=False)
 except Exception:pass
 finish()

def finish():
 heartbeat.stop();w.cancel_bulk();w.viewer.window._qt_window.close()
 report['formal_annotations_after']=formal_hashes()
 if report['formal_annotations_after']!=report['formal_annotations_before']:
  report['PASS']=False;report['error']='Formal annotations changed during isolated check'
 report['checks']['formal_annotations_unchanged']=report['formal_annotations_after']==report['formal_annotations_before']
 report['wall_s']=time.monotonic()-start
 (RUN/'report.json').write_text(json.dumps(report,indent=2));(OUT/'latest_gui_run.json').write_text(json.dumps({'run':str(RUN),'PASS':report['PASS']},indent=2))
 print(json.dumps({'PASS':report['PASS'],'run':str(RUN),'error':report.get('error')}),flush=True)
 app.quit()

QTimer.singleShot(1800,first);napari.run()
if not report['PASS']:raise SystemExit(1)
