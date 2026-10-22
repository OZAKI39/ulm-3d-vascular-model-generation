"""Real Qt saved-library regression; labels are isolated software fixtures."""
from pathlib import Path
import argparse, json, sys, time, traceback, subprocess
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'scripts'))
from annotate import Workbench, QApplication, QTimer, QFont, napari
from qtpy.QtWidgets import QMessageBox
from qtpy.QtTest import QTest
from annotation_store import CLASSES, sha, save_annotation, load_annotation
from annotation_library_store import saved_inventory

OUT = ROOT / 'evidence/annotation_saved_library_20261016'; OUT.mkdir(exist_ok=True)
RUN = OUT / f'gui_{time.time_ns()}'; RUN.mkdir()
app = QApplication.instance() or QApplication([]); app.setFont(QFont('Microsoft YaHei', 10))
args = argparse.Namespace(index=0, smoke=True, manual_only=True, no_auto_train=True, annotation_dir=str(RUN / 'annotations'))
w = Workbench(args); report = {'software_fixture': True, 'checks': {}, 'PASS': False}
def hashes(folder): return {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob('*')) if p.is_file()}
formal_before = hashes(ROOT / 'data/annotations'); started = time.monotonic()

def check(name, value):
    report['checks'][name] = bool(value)
    if not value: raise AssertionError(name)

def select(index):
    row = next(i for i, row in enumerate(w.library_visible) if row['index'] == index)
    w.library_table.selectRow(row); QApplication.processEvents()

def action(key, accept=None):
    if accept is not None:
        def answer():
            dialog = QApplication.activeModalWidget()
            if not isinstance(dialog, QMessageBox): QTimer.singleShot(20, answer); return
            role = QMessageBox.ButtonRole.AcceptRole if accept else QMessageBox.ButtonRole.RejectRole
            next(b for b in dialog.buttons() if dialog.buttonRole(b) == role).click()
        QTimer.singleShot(20, answer)
    w.library_buttons[key].click(); QApplication.processEvents()

def run():
    global w
    try:
        check('append_tab_keeps_existing_tab_indices', w.library_tab == 3 and w.tabs.tabText(1) == 'SAM 辅助')
        check('empty_library_controls_disabled', not w.library_visible and not w.library_buttons['edit'].isEnabled()
              and not w.library_buttons['delete'].isEnabled() and not w.library_buttons['restore'].isEnabled())
        for i, cls in enumerate(CLASSES):
            ys, xs = np.mgrid[80+i*20:88+i*20, 100:108]
            w.layers[cls].data_setitem((ys.ravel(), xs.ravel()), i+1)
        w.notes.value = 'Fixture alpha'; w.save(False)
        check('save_immediately_updates_inventory', len(w.library_visible) == 1)
        check('all_class_counts_visible', all(name in w.library_detail.toPlainText() for name in ['单泡 1', '微泡团 1', '背景 / 杂质 1', '不确定 1', '失焦 1']))
        w.go(1); w.layers['SINGLE_BUBBLE'].data_setitem((np.array([30]), np.array([30])), 10)
        w.notes.value = 'Fixture beta'; w.save(False)
        w.go(2); w.save(False)
        check('zero_label_saved_draft_included', len(w.library_visible) == 3 and w.library_visible[2]['meta']['instance_counts']['SINGLE_BUBBLE'] == 0)
        w.tabs.setCurrentIndex(w.library_tab); w.library_search.setText('alpha')
        check('search_saved_notes', len(w.library_visible) == 1 and w.library_visible[0]['index'] == 0)
        w.library_search.setText('NO SUCH FRAME')
        check('no_match_disables_destructive_actions', not w.library_visible and not w.library_buttons['delete'].isEnabled())
        w.library_search.clear(); w.library_scope.setCurrentIndex(2)
        check('regional_filter_excludes_drafts', not w.library_visible)
        w.library_scope.setCurrentIndex(1); check('draft_filter_lists_saved_drafts', len(w.library_visible) == 3)
        w.library_scope.setCurrentIndex(0)
        # Opening a saved record first preserves unsaved edits in the current frame.
        w.layers['SINGLE_BUBBLE'].data_setitem((np.array([40]), np.array([40])), 12); w.notes.value = 'Unsaved gamma'
        select(0); action('edit')
        check('open_edit_loads_selected_and_returns_manual_tab', w.index == 0 and w.tabs.currentIndex() == 0 and w.notes.value == 'Fixture alpha')
        meta, masks = load_annotation(w.directory, w.records[2])
        check('open_edit_saves_current_unsaved_work', masks['SINGLE_BUBBLE'][40, 40] == 12 and meta['notes'] == 'Unsaved gamma')
        w.tabs.setCurrentIndex(w.library_tab)
        w.layers['SINGLE_BUBBLE'].data_setitem((np.array([50]), np.array([50])), 14)
        unsaved = {k: v.data.copy() for k, v in w.layers.items()}
        select(1); pointer = w.directory / w.records[1]['id'] / 'latest.json'; before = sha(pointer)
        action('delete', False)
        check('cancel_delete_preserves_saved_data', sha(pointer) == before)
        action('delete', True)
        check('delete_other_frame_removes_only_active_pointer', not pointer.exists() and w.index == 0 and w.changed()
              and all(np.array_equal(unsaved[k], v.data) for k, v in w.layers.items()))
        check('deleted_not_read_by_annotation_loader', load_annotation(w.directory, w.records[1]) == (None, None))
        w.library_scope.setCurrentIndex(4); select(1)
        check('trash_lists_record_and_enables_restore_only', len(w.library_visible) == 1 and w.library_buttons['restore'].isEnabled()
              and not w.library_buttons['edit'].isEnabled() and not w.library_buttons['delete'].isEnabled())
        action('restore', False); check('cancel_restore_keeps_pointer_absent', not pointer.exists())
        action('restore', True)
        meta, masks = load_annotation(w.directory, w.records[1])
        check('restore_other_frame_preserves_unsaved_current_pixels', all(np.array_equal(unsaved[k], v.data) for k, v in w.layers.items()) and w.changed())
        check('restore_exact_saved_pixels_as_draft', masks['SINGLE_BUBBLE'][30, 30] == 10 and meta['status'] == 'DRAFT_UNREVIEWED')
        # Current-frame deletion captures unsaved labels, then clears all volatile state.
        select(0); w.notes.value = 'Unsaved alpha before deletion'
        w.preview.data_setitem((np.array([60]), np.array([60])), 1)
        w.bulk_preview.data_setitem((np.array([60]), np.array([60])), 1)
        w.positive.data = np.array([[60., 60.]])
        action('delete', True)
        check('delete_current_clears_labels_previews_and_dirty', all(not np.any(v.data) for v in [*w.layers.values(), w.preview, w.bulk_preview])
              and not w.changed() and not w.dirty and w.meta is None and not len(w.positive.data))
        for layer in [*w.layers.values(), w.preview, w.bulk_preview]: layer.undo(); layer.redo()
        w.autosave()
        check('undo_autosave_do_not_resurrect_deleted_record', not (w.directory / w.record['id'] / 'latest.json').exists()
              and all(not np.any(v.data) for v in w.layers.values()))
        old = w; old.viewer.window._qt_window.close(); args.index = None; w = Workbench(args)
        check('deletion_survives_close_and_restart', w.index == 0 and w.meta is None and not w.changed())
        w.tabs.setCurrentIndex(w.library_tab); w.library_scope.setCurrentIndex(4); select(0)
        action('restore', True)
        check('restore_recovers_unsaved_work_captured_at_delete', all(np.array_equal(unsaved[k], v.data) for k, v in w.layers.items())
              and w.notes.value == 'Unsaved alpha before deletion' and not w.meta['expert_confirmed'])
        # Incomplete manual polygons abort deletion instead of discarding work.
        select(0); pointer = w.directory / w.record['id'] / 'latest.json'; before = sha(pointer)
        w.layers['SINGLE_BUBBLE']._overlays['polygon'].points = [np.array([10., 10.]), np.array([10., 15.])]
        action('delete', True)
        check('pending_polygon_blocks_delete', sha(pointer) == before and '未完成轮廓' in w.library_status.text())
        w.layers['SINGLE_BUBBLE']._overlays['polygon'].points = []
        # An external/new save makes a selected row stale; never silently delete it.
        fresh = save_annotation(w.directory, w.records[1], masks, notes='new revision fixture')
        select(1); w.library_action(w.delete_saved_annotation)
        check('stale_selection_rejected', '变化' in w.library_status.text()
              and load_annotation(w.directory, w.records[1])[0]['revision'] == fresh['revision'])
        w.refresh_saved_library(); select(0)
        # Store corrupt pointers visibly rather than silently hiding them.
        corrupt = w.directory / w.records[3]['id']; corrupt.mkdir(); (corrupt / 'latest.json').write_text('{broken')
        w.refresh_saved_library(); select(3)
        check('bad_metadata_visible_and_safe', '读取异常' in w.library_detail.toPlainText() and not w.library_buttons['edit'].isEnabled())
        (corrupt / 'latest.json').unlink(); w.refresh_saved_library(); select(0)
        for width in (440, 500):
            w.tabs.setFixedWidth(width); QApplication.processEvents(); QTest.qWait(40)
            for first, second in [('edit', 'delete'), ('refresh', 'restore')]:
                a = w.library_buttons[first].geometry(); b = w.library_buttons[second].geometry()
                check(f'buttons_aligned_{first}_{width}', a.y() == b.y() and abs(a.width()-b.width()) <= 1)
            check(f'button_text_fits_{width}', all(b.fontMetrics().horizontalAdvance(b.text()) + 16 <= b.width() for b in w.library_buttons.values()))
            check(f'table_fits_{width}', w.library_table.horizontalScrollBar().maximum() == 0)
        w.library_status.setText('选择一条记录查看详情；双击或点“打开编辑”继续。')
        w.tabs.setFixedWidth(440); QApplication.processEvents()
        w.viewer.screenshot(path=str(RUN / 'saved_library.png'), canvas_only=False)
        w.viewer.window._qt_window.resize(1280, 720); QApplication.processEvents(); QTest.qWait(100)
        check('compact_window_library_scrolls_without_horizontal_overflow', w.library_scroll.horizontalScrollBar().maximum() == 0)
        w.viewer.screenshot(path=str(RUN / 'compact_library.png'), canvas_only=False)
        w.library_scope.setCurrentIndex(4)
        check('restored_archives_not_offered_again', not w.library_visible)
        check('formal_annotations_unchanged', hashes(ROOT / 'data/annotations') == formal_before)
        check('test_training_disabled', not w.auto_training_allowed)
        report['PASS'] = True
    except Exception:
        report['error'] = traceback.format_exc(); print(report['error'], flush=True)
    finally:
        w.cancel_sam(); w.cancel_bulk(); w.viewer.window._qt_window.close()
        report['formal_annotations_before'] = formal_before; report['formal_annotations_after'] = hashes(ROOT / 'data/annotations')
        report['wall_s'] = time.monotonic() - started
        report['head'] = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
        report['source_sha256'] = {name: sha(ROOT / name) for name in ['scripts/annotate.py', 'src/annotation_library.py', 'src/annotation_library_store.py', 'scripts/check_annotation_library.py']}
        (RUN / 'report.json').write_text(json.dumps(report, indent=2))
        (OUT / 'latest_gui_run.json').write_text(json.dumps({'run': str(RUN), 'PASS': report['PASS']}, indent=2))
        print(json.dumps({'PASS': report['PASS'], 'checks': len(report['checks']), 'run': str(RUN), 'error': report.get('error')}), flush=True)
        app.quit()

QTimer.singleShot(1000, run); napari.run()
if not report['PASS']: raise SystemExit(1)
