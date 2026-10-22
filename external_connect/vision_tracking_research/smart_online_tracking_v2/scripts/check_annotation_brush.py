"""Exercise real Qt brush clicks and shortcuts; never write human annotations."""
from pathlib import Path
import argparse,hashlib,json,sys,time,traceback
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from annotate import Workbench,QApplication,QTimer,QFont,napari
from qtpy.QtCore import Qt,QPoint,QPointF
from qtpy.QtGui import QWheelEvent
from qtpy.QtTest import QTest

parser=argparse.ArgumentParser()
parser.add_argument('--output-dir',type=Path)
args=parser.parse_args()
OUT=args.output_dir or ROOT/f'evidence/annotation_brush_shortcuts_{time.time_ns()}'
OUT.mkdir(parents=True,exist_ok=True)
def formal_hashes():
    return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted((ROOT/'data/annotations').rglob('*')) if p.is_file()}
formal_before=formal_hashes()
app=QApplication.instance() or QApplication([])
app.setFont(QFont('Microsoft YaHei',10))
w=Workbench(argparse.Namespace(index=0,smoke=True,manual_only=True,
    no_auto_train=True,annotation_dir=str(OUT/f'sandbox_{time.time_ns()}')))
report={'software_fixture':True,'formal_labels_modified':False,'checks':{},'PASS':False}

def check(name,value):
    report['checks'][name]=bool(value)
    if not value:raise AssertionError(name)

def run():
    try:
        layer=w.layers['SINGLE_BUBBLE']
        canvas=w.viewer.window._qt_viewer.canvas.native
        w.viewer.window._qt_window.activateWindow();w.viewer.window._qt_window.raise_();QTest.qWait(100)
        def key(k,mods=Qt.KeyboardModifier.NoModifier):
            canvas.setFocus();QTest.keyClick(canvas,k,mods);QApplication.processEvents()
        def paint_at(dx,dy=0):
            QTest.mouseClick(canvas,Qt.MouseButton.LeftButton,
                Qt.KeyboardModifier.NoModifier,QPoint(canvas.width()//2+dx,canvas.height()//2+dy))
            QApplication.processEvents()
        def wheel(delta,mods=Qt.KeyboardModifier.ControlModifier,target=None,pixels=False):
            target=canvas if target is None else target
            pos=QPoint(target.width()//2,target.height()//2)
            event=QWheelEvent(QPointF(pos),QPointF(target.mapToGlobal(pos)),
                QPoint(0,delta) if pixels else QPoint(),QPoint() if pixels else QPoint(0,delta),
                Qt.MouseButton.NoButton,mods,Qt.ScrollPhase.NoScrollPhase,False)
            QApplication.sendEvent(target,event);QApplication.processEvents()
        report['startup']={'category':w.category.value,'mode':str(layer.mode),'selected_label':int(layer.selected_label),
            'active':w.viewer.layers.selection.active.name if w.viewer.layers.selection.active else None,'tab':w.tabs.currentIndex()}
        check('startup_single_bubble_brush',w.category.value=='SINGLE_BUBBLE'
            and layer.mode=='paint' and layer.selected_label==1
            and w.viewer.layers.selection.active is layer and w.tabs.currentIndex()==0)
        w.notes.value='SOFTWARE_UI_FIXTURE; brush shortcuts, not human ground truth'
        canvas.setFocus();paint_at(-50)
        check('real_mouse_brush_writes_pixels',np.any(layer.data==1))
        before=layer.data.copy();paint_at(-35)
        check('multiple_strokes_one_instance',np.count_nonzero(layer.data)>np.count_nonzero(before)
            and np.array_equal(np.unique(layer.data),[0,1]))
        key(Qt.Key.Key_N)
        check('N_allocates_next_single_bubble',layer.selected_label==2 and w.instance.value==2 and layer.mode=='paint')
        paint_at(40)
        check('second_bubble_has_separate_ID',np.array_equal(np.unique(layer.data),[0,1,2]))
        before=layer.data.copy()
        key(Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier)
        check('keyboard_undo_removes_last_stroke_only',np.any(layer.data==1) and not np.any(layer.data==2))
        key(Qt.Key.Key_Z,Qt.KeyboardModifier.ControlModifier|Qt.KeyboardModifier.ShiftModifier)
        check('keyboard_redo_exact',np.array_equal(layer.data,before))
        key(Qt.Key.Key_N);key(Qt.Key.Key_N)
        check('unused_next_ID_not_duplicated_or_saved',layer.selected_label==3 and np.array_equal(layer.data,before))
        paint_at(-50)
        check('painting_does_not_overwrite_existing_ID',np.array_equal(layer.data,before))
        key(Qt.Key.Key_E)
        check('E_eraser_same_ID',layer.mode=='erase' and layer.selected_label==3)
        key(Qt.Key.Key_P)
        check('P_brush_same_ID',layer.mode=='paint' and layer.selected_label==3)
        size=layer.brush_size;key(Qt.Key.Key_BracketRight)
        check('increase_brush_size_key',layer.brush_size==size+1)
        key(Qt.Key.Key_BracketLeft)
        check('decrease_brush_size_key',layer.brush_size==size)
        w.save(False)
        before=layer.data.copy();selected=layer.selected_label
        dirty=w.dirty;zoom=w.viewer.camera.zoom
        # Send the real key press as well, to catch conflicts with native modes.
        QTest.keyPress(canvas,Qt.Key.Key_Control)
        wheel(120)
        check('Ctrl_wheel_up_increases_one_pixel',layer.brush_size==size+1)
        check('size_label_updates_live',f'{size+1:g} px' in w.brush_help.value)
        wheel(-120)
        check('Ctrl_wheel_down_decreases_one_pixel',layer.brush_size==size)
        wheel(60);check('partial_wheel_waits_for_full_step',layer.brush_size==size)
        wheel(60);check('partial_wheel_accumulates',layer.brush_size==size+1)
        QTest.keyRelease(canvas,Qt.Key.Key_Control)
        key(Qt.Key.Key_E);wheel(120)
        check('Ctrl_wheel_works_in_eraser',layer.mode=='erase' and layer.brush_size==size+2)
        check('resize_preserves_labels_ID_and_dirty_state',np.array_equal(layer.data,before)
            and selected==layer.selected_label and dirty==w.dirty)
        check('Ctrl_wheel_does_not_zoom',w.viewer.camera.zoom==zoom)
        layer.brush_size=1;wheel(-120)
        check('minimum_brush_size_one_pixel',layer.brush_size==1)
        wheel(120);check('resize_immediate_after_lower_limit',layer.brush_size==2)
        wheel(40,pixels=True)
        check('pixel_scroll_supported',layer.brush_size==3)
        wheel(120,Qt.KeyboardModifier.NoModifier)
        check('normal_wheel_still_zooms',w.viewer.camera.zoom!=zoom and layer.brush_size==3)
        w.notes.native.setFocus();wheel(120,target=w.notes.native)
        check('wheel_over_text_does_not_resize',layer.brush_size==3)
        w.pan();wheel(120)
        check('pan_mode_does_not_resize',layer.brush_size==3)
        w.prompt_mode(True);wheel(120)
        check('SAM_points_do_not_resize_brush',layer.brush_size==3)
        w.edit_mode('paint');w.viewer.reset_view()
        layer.brush_size=7;key(Qt.Key.Key_BracketRight)
        check('native_keys_update_size_label','8 px' in w.brush_help.value)
        layer.brush_size=9
        w.buttons['new_bubble'].native.click()
        QApplication.processEvents()
        report['focus']={'widget':str(QApplication.focusWidget()),'canvas':str(canvas),
            'active_window':str(QApplication.activeWindow()),'mode':str(layer.mode)}
        check('button_focus_returns_to_canvas',canvas.hasFocus() and layer.mode=='paint')
        w.category.value='AGGREGATE'
        w.notes.native.setFocus();QTest.keyClicks(w.notes.native,'n');QApplication.processEvents()
        check('typing_n_is_not_shortcut',w.notes.value.endswith('n') and w.category.value=='AGGREGATE')
        w.tabs.setCurrentIndex(1);w.prompt_mode(True)
        key(Qt.Key.Key_N)
        check('N_from_SAM_switches_to_single_brush',w.category.value=='SINGLE_BUBBLE'
            and w.tabs.currentIndex()==0 and w.viewer.layers.selection.active is layer and layer.mode=='paint')
        check('brush_size_preserved_for_new_bubble',layer.brush_size==9)
        w.pick_mode();key(Qt.Key.Key_P)
        check('brush_exits_object_picker',not w.picking and layer.mode=='paint')
        w.new_object();points=[[100,100],[100,110],[110,110]]
        layer._overlays['polygon'].points=points
        key(Qt.Key.Key_N)
        check('unfinished_polygon_preserved',layer._overlays['polygon'].points==points
            and '未完成轮廓' in w.status.value and layer.mode=='polygon')
        w.finish_polygon()
        check('optional_polygon_still_works',np.any(layer.data==3))
        before=layer.data.copy();w.save(False)
        w.go(1)
        check('next_frame_brush_ready',layer.mode=='paint' and layer.selected_label==1)
        w.go(0)
        check('saved_labels_reloaded_exact',np.array_equal(layer.data,before)
            and layer.selected_label==4 and layer.mode=='paint' and not w.meta['expert_confirmed'])
        w.viewer.layers.remove(layer);QApplication.processEvents()
        key(Qt.Key.Key_N)
        check('restored_layer_keeps_shortcut_and_labels',layer in w.viewer.layers
            and layer.selected_label==4 and layer.mode=='paint' and np.array_equal(layer.data,before))
        size=layer.brush_size;wheel(120)
        check('restored_layer_keeps_wheel_shortcut',layer.brush_size==size+1 and np.array_equal(layer.data,before))
        w.category.value='BACKGROUND_ARTIFACT';w.buttons['new_brush'].native.click()
        check('new_current_category_brush',w.category.value=='BACKGROUND_ARTIFACT'
            and w.layers['BACKGROUND_ARTIFACT'].mode=='paint')
        other=w.layers['BACKGROUND_ARTIFACT'];size=other.brush_size;single_size=layer.brush_size;wheel(120)
        check('only_active_class_is_resized',other.brush_size==size+1 and layer.brush_size==single_size)
        key(Qt.Key.Key_N)
        check('N_returns_to_single_bubble_without_reclassification',w.category.value=='SINGLE_BUBBLE'
            and np.array_equal(layer.data,before))
        w.main.native.parentWidget().parentWidget().ensureWidgetVisible(w.brush_help.native)
        QApplication.processEvents()
        w.viewer.screenshot(path=str(OUT/'brush_workbench.png'),canvas_only=False)
        check('test_never_starts_auto_training',not w.auto_training_allowed)
        w.viewer.window._qt_window.close()
        report['formal_annotations_before']=formal_before
        report['formal_annotations_after']=formal_hashes()
        check('formal_annotations_unchanged',report['formal_annotations_after']==formal_before)
        report['PASS']=True
    except Exception:
        report['error']=traceback.format_exc()
        w.viewer.screenshot(path=str(OUT/'failed_workbench.png'),canvas_only=False)
        w.closing=True;w.cancel_sam();w.viewer.close()
    finally:
        (OUT/'ui_checks.json').write_text(json.dumps(report,indent=2))
        print(json.dumps(report),flush=True)
        app.quit()

QTimer.singleShot(1800,run)
napari.run()
if not report['PASS']:raise SystemExit(1)
