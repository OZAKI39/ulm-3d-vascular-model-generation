"""Actual Qt candidate erasing, circular cursor and sidebar layout; isolated fixtures."""
from pathlib import Path
import argparse,json,sys,time,traceback
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from annotate import Workbench,QApplication,QTimer,QFont,napari
from qtpy.QtCore import Qt,QPoint,QPointF
from qtpy.QtGui import QWheelEvent
from qtpy.QtTest import QTest
from annotation_store import sha,load_annotation
from annotation_sizes import size_plan

OUT=ROOT/'evidence/annotation_candidate_tools_20261016';OUT.mkdir(exist_ok=True)
RUN=OUT/f'gui_{time.time_ns()}';RUN.mkdir()
app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
w=Workbench(argparse.Namespace(index=0,smoke=True,manual_only=True,no_auto_train=True,annotation_dir=str(RUN/'annotations')))
report={'software_fixture':True,'checks':{},'PASS':False}
def formal_hashes():return {str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'data/annotations').rglob('*')) if p.is_file()}
report['formal_before']=formal_hashes()
def check(name,value):
 report['checks'][name]=bool(value)
 if not value:raise AssertionError(name)

def run():
 try:
  canvas=w.viewer.window._qt_viewer.canvas.native
  w.viewer.window._qt_window.activateWindow();w.viewer.window._qt_window.raise_();QTest.qWait(100)
  def key(k,mods=Qt.KeyboardModifier.NoModifier):
   canvas.setFocus();QTest.keyClick(canvas,k,mods);QApplication.processEvents()
  def mouse(dx=0):
   p=QPoint(canvas.width()//2+dx,canvas.height()//2)
   QTest.mouseMove(canvas,p);QTest.mouseClick(canvas,Qt.MouseButton.LeftButton,Qt.KeyboardModifier.NoModifier,p)
   QApplication.processEvents()
  def wheel(delta):
   p=QPoint(canvas.width()//2,canvas.height()//2)
   event=QWheelEvent(QPointF(p),QPointF(canvas.mapToGlobal(p)),QPoint(),QPoint(0,delta),
                     Qt.MouseButton.NoButton,Qt.KeyboardModifier.ControlModifier,Qt.ScrollPhase.NoScrollPhase,False)
   QApplication.sendEvent(canvas,event);QApplication.processEvents()
  check('no_candidates_refuses_eraser',w.safe(w.erase_candidate_mode) is False and '先生成' in w.status.value)
  # Synthetic masks only, intentionally not model predictions or human truth.
  box=[340,440,460,620];y0,x0,y1,x1=box
  w.roi_layer.data=[np.array([[y0,x0],[y0,x1],[y1,x1],[y1,x0]])]
  yy,xx=np.mgrid[y0:y1,x0:x1];local=np.zeros(yy.shape,np.uint32)
  local[(yy-400)**2+(xx-500)**2<=20**2]=1
  local[(yy-400)**2+(xx-560)**2<=16**2]=2
  path=RUN/'original_proposals.npz';np.savez_compressed(path,labels=local)
  req={'model':'cyto3','output':str(path),'diameter_range_px':[8,24],'diameter_filter':'REFERENCE_ONLY','software_test':True}
  request=RUN/'request.json';request.write_text(json.dumps(req))
  data=np.zeros(w.image.data.shape[:2],np.uint32);data[y0:y1,x0:x1]=local
  w.bulk_preview.data=data;w.bulk_preview._reset_history()
  w.bulk_result={'request_path':str(request),'model':'SOFTWARE_FIXTURE','model_sha256':{},
                 'output_sha256':sha(path),'diameter_range_px':[8,24],'size_plan':size_plan(req),'merge_policy':'SOFTWARE_FIXTURE'}
  w.bulk_signature=w.bulk_prompt_signature()
  manual=w.layers['SINGLE_BUBBLE'];manual.data_setitem((np.array([400]),np.array([500])),10)
  original_manual=manual.data.copy();w.notes.value='SOFTWARE_UI_FIXTURE';w.save(False)
  dirty=w.dirty;original_sha=sha(path)
  w.viewer.camera.center=(0,400,500);w.viewer.camera.zoom=4
  key(Qt.Key.Key_R);QTest.qWait(80)
  check('R_selects_candidate_eraser',w.viewer.layers.selection.active is w.bulk_preview and w.bulk_preview.mode=='erase' and not w.bulk_preview.preserve_labels)
  before=w.bulk_preview.data.copy();mouse()
  erased=w.bulk_preview.data.copy()
  check('real_mouse_erases_part_of_candidate',0<np.count_nonzero(erased==1)<np.count_nonzero(before==1))
  check('other_candidate_keeps_exact_pixels',np.array_equal(erased==2,before==2))
  check('manual_labels_and_dirty_unchanged',np.array_equal(manual.data,original_manual) and w.dirty==dirty)
  check('erasure_only_removes_pixels',np.all((erased==before)|(erased==0)))
  key(Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier)
  check('candidate_keyboard_undo_exact',np.array_equal(w.bulk_preview.data,before))
  key(Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier|Qt.KeyboardModifier.ShiftModifier)
  check('candidate_keyboard_redo_exact',np.array_equal(w.bulk_preview.data,erased))
  w.buttons['undo_candidate'].native.click();check('candidate_button_undo_exact',np.array_equal(w.bulk_preview.data,before))
  w.buttons['redo_candidate'].native.click();check('candidate_button_redo_exact',np.array_equal(w.bulk_preview.data,erased))
  size=w.bulk_preview.brush_size;zoom=w.viewer.camera.zoom;wheel(120)
  check('candidate_Ctrl_wheel_resizes_without_zoom',w.bulk_preview.brush_size==size+1 and w.viewer.camera.zoom==zoom)
  key(Qt.Key.Key_BracketLeft);check('candidate_bracket_resizes',w.bulk_preview.brush_size==size)
  check('candidate_size_help_is_local','候选橡皮擦直径' in w.candidate_help.value)
  key(Qt.Key.Key_P);check('P_returns_to_manual_labels',w.viewer.layers.selection.active is manual and manual.mode=='paint')
  key(Qt.Key.Key_R);key(Qt.Key.Key_E);check('E_in_candidates_stays_candidate_eraser',w.viewer.layers.selection.active is w.bulk_preview)
  w.drop_candidate_mode()
  for _ in w.bulk_pick_callback(w.viewer,argparse.Namespace(position=(400.,560.))):pass
  check('whole_candidate_delete_still_works',not np.any(w.bulk_preview.data==2))
  w.bulk_preview.undo();check('whole_candidate_delete_undo_exact',np.array_equal(w.bulk_preview.data,erased))
  w.erase_candidate_mode();w.step_jump.value=3;QApplication.processEvents()
  QTest.mouseMove(canvas,QPoint(canvas.width()//2,canvas.height()//2));QTest.qWait(100)
  w.viewer.screenshot(path=str(RUN/'candidate_eraser.png'),canvas_only=False)
  w.incorporate_bulk()
  edits=w.proposal_sources[-1]['manual_candidate_edits']
  check('edited_candidates_added_as_draft',np.count_nonzero(manual.data)>np.count_nonzero(original_manual) and manual.data[400,500]==10)
  for item in w.proposal_sources[-1]['accepted']:
   check(f'accepted_exact_edited_mask_{item["candidate_id"]}',np.array_equal(manual.data==item['instance_id'],erased==item['candidate_id']))
  check('original_model_file_unchanged',sha(path)==original_sha)
  check('manual_edit_provenance_and_hash',edits['erased_pixels']==np.count_nonzero(before!=erased) and sha(edits['mask_path'])==edits['mask_sha256'])
  with np.load(edits['mask_path']) as z:check('edited_mask_saved_exactly',np.array_equal(z['labels'],erased[y0:y1,x0:x1]))
  w.save(False);meta,_=load_annotation(w.directory,w.record)
  check('draft_provenance_saved_not_reviewed',not meta['expert_confirmed'] and meta['proposal_sources'][-1]['manual_candidate_edits']==edits)
  w.bulk_preview.undo();w.bulk_preview.redo()
  check('old_candidates_do_not_return_after_incorporation',not np.any(w.bulk_preview.data))
  # The circle must remain visible even when the data-space brush is tiny on screen.
  for mode in ('paint','erase'):
   w.edit_mode(mode);manual.brush_size=3
   for zoom in (.25,1,4):
    w.viewer.camera.zoom=zoom;QApplication.processEvents()
    QTest.mouseMove(canvas,QPoint(canvas.width()//2+8,canvas.height()//2));QTest.qWait(50)
    overlay=w.viewer.canvas.overlays._brush_circle
    check(f'circle_visible_{mode}_{zoom}',overlay.visible and canvas.cursor().shape()==Qt.CursorShape.BlankCursor)
    check(f'circle_size_{mode}_{zoom}',overlay.size==max(8,round(3*zoom)) and manual.brush_size==3)
  w.viewer.camera.zoom=.25;key(Qt.Key.Key_N)
  check('N_keeps_circle_in_overview',w.viewer.canvas.overlays._brush_circle.visible and canvas.cursor().shape()==Qt.CursorShape.BlankCursor)
  check('minimum_display_ring_disclosed','小圆环仅定位' in w.brush_help.value)
  w.pan();check('pan_restores_normal_cursor',not w.viewer.canvas.overlays._brush_circle.visible and canvas.cursor().shape()!=Qt.CursorShape.BlankCursor)
  w.prompt_mode(True);check('SAM_point_cursor_not_forced_circle',not w.viewer.canvas.overlays._brush_circle.visible)
  key(Qt.Key.Key_N);w.viewer.camera.center=(0,400,500);w.viewer.camera.zoom=4
  for width in (440,500):
   w.tabs.setFixedWidth(width);QApplication.processEvents();QTest.qWait(60)
   columns=[]
   for a,b in [('roi','zoom_roi'),('generate_bulk','cancel_bulk'),('erase_candidate','drop_candidate'),('undo_candidate','redo_candidate'),('incorporate','clear_bulk'),('brush','erase'),('save','confirm'),('reset_frame','restore_frame')]:
    ga=w.buttons[a].native.geometry();gb=w.buttons[b].native.geometry()
    check(f'aligned_{a}_{width}',ga.y()==gb.y() and abs(ga.width()-gb.width())<=1 and ga.height()==gb.height()==32)
    columns.append((ga.x(),gb.x(),ga.width(),gb.width()))
    for name in (a,b):
     btt=w.buttons[name];check(f'text_fits_{name}_{width}',btt.native.fontMetrics().horizontalAdvance(btt.text)+14<=btt.native.width())
   check(f'common_columns_{width}',len(set(columns))==1)
   check(f'five_step_headers_ordered_{width}',all(a.native.y()<b.native.y() for a,b in zip(w.step_headers,w.step_headers[1:])))
   for step in range(1,6):
    w.step_jump.value=step;QApplication.processEvents()
    header=w.step_headers[step-1].native;scroll=w.main_scroll.verticalScrollBar().value()
    check(f'jump_step_{step}_{width}',0<=header.y()-scroll<w.main_scroll.viewport().height())
   check(f'no_horizontal_scroll_{width}',w.main_scroll.horizontalScrollBar().maximum()==0)
  w.viewer.window._qt_window.resize(1280,720);QApplication.processEvents();QTest.qWait(80)
  w.step_jump.value=3;QApplication.processEvents()
  check('short_window_step_navigation_visible',w.step_jump.native.isVisible() and w.main_scroll.verticalScrollBar().maximum()>0)
  w.viewer.screenshot(path=str(RUN/'compact_sidebar.png'),canvas_only=False)
  w.viewer.window._qt_window.resize(1680,1000);w.step_jump.value=4;QApplication.processEvents()
  QTest.mouseMove(canvas,QPoint(canvas.width()//2+40,canvas.height()//2));QTest.qWait(100)
  w.viewer.screenshot(path=str(RUN/'round_brush.png'),canvas_only=False)
  check('formal_annotations_unchanged',formal_hashes()==report['formal_before'])
  check('no_training_started',not w.auto_training_allowed)
  report['PASS']=True
 except Exception:
  report['error']=traceback.format_exc();print(report['error'],flush=True)
  w.viewer.screenshot(path=str(RUN/'failure.png'),canvas_only=False)
 finally:
  w.viewer.window._qt_window.close();report['formal_after']=formal_hashes()
  (RUN/'report.json').write_text(json.dumps(report,indent=2));(OUT/'latest_gui_run.json').write_text(json.dumps({'run':str(RUN),'PASS':report['PASS']},indent=2))
  print(json.dumps({'PASS':report['PASS'],'checks':len(report['checks']),'run':str(RUN),'error':report.get('error')}),flush=True);app.quit()

QTimer.singleShot(1500,run);napari.run()
if not report['PASS']:raise SystemExit(1)
