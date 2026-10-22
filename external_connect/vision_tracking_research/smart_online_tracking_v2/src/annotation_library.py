"""Saved-label browser for the existing workbench, without loading any model."""
from datetime import datetime
from zoneinfo import ZoneInfo
from qtpy.QtCore import Qt
from qtpy.QtWidgets import (QWidget, QGridLayout, QLabel, QComboBox, QLineEdit,
                           QPushButton, QTableWidget, QTableWidgetItem,
                           QAbstractItemView, QHeaderView, QPlainTextEdit, QSizePolicy, QScrollArea)
from annotation_library_store import saved_inventory, read_revision, delete_saved, restore_saved

STATUS = {'DRAFT_UNREVIEWED': '草稿', 'REGIONS_CONFIRMED': '区域审核', 'EXPERT_CONFIRMED': '整帧审核'}
CLASS_NAMES = {'SINGLE_BUBBLE': '单泡', 'AGGREGATE': '微泡团', 'BACKGROUND_ARTIFACT': '背景 / 杂质',
               'UNCERTAIN': '不确定', 'OUT_OF_FOCUS': '失焦'}


class SavedAnnotationWorkflow:
    def init_saved_library(self):
        page = QWidget(); grid = QGridLayout(page)
        grid.setContentsMargins(8, 8, 8, 8); grid.setHorizontalSpacing(8); grid.setVerticalSpacing(6)
        grid.setColumnStretch(0, 1); grid.setColumnStretch(1, 1)
        intro = QLabel('查看已保存的标签；双击一行即可继续编辑。\n草稿与审核完成分开显示，未保存笔画不计入。')
        intro.setWordWrap(True); grid.addWidget(intro, 0, 0, 1, 2)
        self.library_scope = QComboBox()
        for text, key in [('全部已保存', 'active'), ('仅草稿', 'DRAFT_UNREVIEWED'),
                          ('仅区域审核', 'REGIONS_CONFIRMED'), ('仅整帧审核', 'EXPERT_CONFIRMED'),
                          ('已删除（可恢复）', 'deleted')]:
            self.library_scope.addItem(text, key)
        grid.addWidget(self.library_scope, 1, 0, 1, 2)
        self.library_search = QLineEdit(); self.library_search.setPlaceholderText('搜索图片序号、ID、审核人或备注')
        grid.addWidget(self.library_search, 2, 0, 1, 2)
        self.library_summary = QLabel(); self.library_summary.setWordWrap(True)
        grid.addWidget(self.library_summary, 3, 0, 1, 2)
        self.library_table = QTableWidget(0, 4)
        self.library_table.setHorizontalHeaderLabels(['图片', '单泡', '团', '状态'])
        self.library_table.verticalHeader().hide()
        self.library_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.library_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.library_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.library_table.setWordWrap(False); self.library_table.setMinimumHeight(150)
        header = self.library_table.horizontalHeader(); header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        grid.addWidget(self.library_table, 4, 0, 1, 2); grid.setRowStretch(4, 1)
        self.library_detail = QPlainTextEdit(); self.library_detail.setReadOnly(True)
        self.library_detail.setFixedHeight(148); grid.addWidget(self.library_detail, 5, 0, 1, 2)
        self.library_buttons = {}
        for key, text, action, row, col in [
                ('edit', '打开编辑', self.edit_saved_annotation, 6, 0),
                ('delete', '删除保存记录…', self.delete_saved_annotation, 6, 1),
                ('refresh', '刷新列表', self.refresh_saved_library, 7, 0),
                ('restore', '恢复为草稿…', self.restore_saved_annotation, 7, 1)]:
            button = QPushButton(text); button.setFixedHeight(32); button.setMinimumWidth(0)
            button.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)
            button.clicked.connect(lambda _=False, fn=action: self.library_action(fn))
            self.library_buttons[key] = button; grid.addWidget(button, row, col)
        self.library_status = QLabel('选择一条记录查看详情。'); self.library_status.setWordWrap(True)
        grid.addWidget(self.library_status, 8, 0, 1, 2)
        self.library_rows = []; self.library_visible = []
        self.library_scope.currentIndexChanged.connect(lambda *_: self.refresh_saved_library())
        self.library_search.textChanged.connect(lambda *_: self.filter_saved_library())
        self.library_table.itemSelectionChanged.connect(self.show_saved_detail)
        self.library_table.itemDoubleClicked.connect(lambda *_: self.library_action(self.edit_saved_annotation))
        self.library_scroll = QScrollArea(); self.library_scroll.setWidgetResizable(True)
        self.library_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.library_scroll.setWidget(page)
        self.library_tab = self.tabs.addTab(self.library_scroll, '已保存标注')
        self.tabs.currentChanged.connect(lambda index: self.refresh_saved_library() if index == self.library_tab else None)

    def library_action(self, fn):
        try:
            return fn()
        except Exception as exc:
            self.library_status.setText('操作未完成：' + str(exc))
            return False

    def selected_saved_row(self):
        index = self.library_table.currentRow()
        return self.library_visible[index] if 0 <= index < len(self.library_visible) else None

    def refresh_saved_library(self):
        self.library_rows = saved_inventory(self.directory, self.records, self.library_scope.currentData() == 'deleted')
        self.filter_saved_library()

    def filter_saved_library(self):
        selected = self.selected_saved_row(); key = selected['key'] if selected else None
        query = self.library_search.text().strip().casefold(); scope = self.library_scope.currentData()
        visible = []
        for row in self.library_rows:
            meta = row['meta']; record = row['record']
            if scope not in ('active', 'deleted') and meta.get('status') != scope:
                continue
            haystack = f'{row["index"]+1} {record["id"]} {meta.get("reviewer", "")} {meta.get("notes", "")}'
            if query and query not in haystack.casefold():
                continue
            visible.append(row)
        table = self.library_table; table.blockSignals(True); table.clearContents(); table.setRowCount(len(visible))
        self.library_visible = visible
        for index, row in enumerate(visible):
            meta = row['meta']; counts = meta.get('instance_counts', {})
            status = '读取异常' if row['error'] else '已删除' if row['deleted'] else STATUS.get(meta.get('status'), '未知')
            values = [f'{row["index"]+1:03d}  {row["record"]["id"]}', str(counts.get('SINGLE_BUBBLE', 0)),
                      str(counts.get('AGGREGATE', 0)), status]
            for col, value in enumerate(values):
                item = QTableWidgetItem(value); item.setToolTip(row['record']['id'] + '\n' + status)
                table.setItem(index, col, item)
        chosen = next((i for i, row in enumerate(visible) if row['key'] == key), 0)
        if visible: table.selectRow(chosen)
        else: table.setCurrentCell(-1, -1)
        table.blockSignals(False)
        self.library_summary.setText(f'显示 {len(visible)} / {len(self.library_rows)} 条' + ('已删除记录' if scope == 'deleted' else '保存记录（每图最新版本）'))
        self.show_saved_detail()

    def show_saved_detail(self):
        row = self.selected_saved_row(); valid = bool(row and not row['error'])
        active = valid and not row['deleted']; deleted = valid and row['deleted']
        self.library_buttons['edit'].setEnabled(active); self.library_buttons['delete'].setEnabled(active)
        self.library_buttons['restore'].setEnabled(deleted and not (self.directory / row['record']['id'] / 'latest.json').exists())
        if not row:
            self.library_detail.setPlainText('暂无符合条件的记录。保存草稿后会自动出现在这里。'); return
        meta = row['meta']; record = row['record']; counts = meta.get('instance_counts', {})
        try: when = datetime.fromisoformat(meta['utc']).astimezone(ZoneInfo('Europe/Zurich')).strftime('%Y-%m-%d %H:%M:%S %Z')
        except (ValueError, KeyError): when = meta.get('utc', '未知')
        details = [f'图片 {row["index"]+1}：{record["id"]}', f'保存时间：{when}',
                   f'状态：{STATUS.get(meta.get("status"), "未知")}；分组：{record["split"]}',
                   '；'.join(f'{name} {counts.get(key, 0)}' for key, name in CLASS_NAMES.items()),
                   f'已确认区域：{len(meta.get("review_regions", []))}；审核人：{meta.get("reviewer") or "未填写"}',
                   f'备注：{meta.get("notes") or "无"}']
        if row['deleted']:
            details.insert(2, '删除编号（UTC）：' + row['pointer'].stem)
            details.append('恢复后为未审核草稿；同图已有记录时不覆盖。')
        if row['error']: details.append('读取异常：' + row['error'])
        self.library_detail.setPlainText('\n'.join(details))

    def edit_saved_annotation(self):
        row = self.selected_saved_row()
        if not row or row['deleted'] or row['error']: return False
        # Preflight before navigation changes the current frame or cancels jobs.
        read_revision(self.directory / row['record']['id'], row['pointer'], row['record'])
        if not self.save_if_changed():
            self.library_status.setText(self.status.value); return False
        self.load_frame(row['index']); self.tabs.setCurrentIndex(0); self.jump_to_step(4)
        self.library_status.setText('已打开；修改后请保存草稿，并重新审核受影响的范围。')
        return True

    def delete_saved_annotation(self):
        row = self.selected_saved_row()
        if not row or row['deleted'] or row['error']: return False
        record = row['record']; current = record['id'] == self.record['id']
        meta, _ = read_revision(self.directory / record['id'], row['pointer'], record)
        if meta['revision'] != row['meta']['revision']:
            raise ValueError('保存记录已变化，请刷新后重试')
        if not self.confirm_frame_action('删除保存记录',
                f'图片 {row["index"]+1}：{record["id"]}\n删除这张图的当前保存记录，移入“已删除（可恢复）”。\n'
                + ('当前未保存标签会先保存，再一起移入已删除记录；候选预览将清空。\n' if current else '')
                + '原图和历史修订保留，其他帧不受影响。后续训练不再读取此记录；已训练模型不回溯撤销。'):
            return False
        if current and not self.save_if_changed():
            self.library_status.setText(self.status.value); return False
        revision = self.meta['revision'] if current else meta['revision']
        delete_saved(self.directory, record, revision)
        if current:
            self.load_frame(self.index)  # Clear dirty data and undo stacks: no autosave resurrection.
            self.preview._reset_history()
        self.refresh_saved_library()
        self.library_status.setText('已删除保存记录。可在上方“已删除（可恢复）”中恢复为草稿；原图保留。')
        return True

    def restore_saved_annotation(self):
        row = self.selected_saved_row()
        if not row or not row['deleted'] or row['error']: return False
        current = row['record']['id'] == self.record['id']
        if current and (self.changed() or any(v._overlays['polygon'].points for v in self.layers.values())):
            raise ValueError('当前图片有新笔画或未完成轮廓，请先保存；恢复不会覆盖现有内容')
        if not self.confirm_frame_action('恢复保存记录',
                f'图片 {row["index"]+1}：{row["record"]["id"]}\n恢复五类标签、备注、短轨迹 ID 和审核选区。\n'
                '恢复为未审核草稿，需要重新检查确认；不会覆盖同图已有保存记录。'):
            return False
        restore_saved(self.directory, row['record'], row['pointer'], self.person.value)
        if current: self.load_frame(self.index)
        self.library_scope.setCurrentIndex(0); self.library_search.clear(); self.refresh_saved_library()
        for i, candidate in enumerate(self.library_visible):
            if candidate['record']['id'] == row['record']['id']: self.library_table.selectRow(i); break
        self.library_status.setText('已恢复为草稿。点“打开编辑”继续；请重新检查并确认审核范围。')
        return True
