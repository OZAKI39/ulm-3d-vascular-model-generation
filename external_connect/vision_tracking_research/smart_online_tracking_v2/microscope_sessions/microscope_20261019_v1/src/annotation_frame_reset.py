"""Frame-local reset with durable recovery of saved and unsaved label pixels."""
from pathlib import Path
from datetime import datetime,timezone
import json,uuid
import numpy as np
from qtpy.QtWidgets import QMessageBox
from annotation_store import CLASSES,sha,validate,save_annotation


class FrameResetWorkflow:
    def reset_backup_dir(self):
        return self.directory/self.record['id']/'reset_backups'

    def reset_snapshot(self,make_current=True):
        """Backup only; never change latest annotation or feed this into training."""
        folder=self.reset_backup_dir();folder.mkdir(parents=True,exist_ok=True)
        name=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'_'+uuid.uuid4().hex[:8]
        masks={k:v.data.copy() for k,v in self.layers.items()}
        path=folder/(name+'.npz')
        with path.with_suffix('.tmp').open('wb') as f:np.savez_compressed(f,**masks)
        path.with_suffix('.tmp').replace(path)
        try:roi=self.current_roi()
        except ValueError:roi=None
        state={'schema':1,'candidate_id':self.record['id'],'image_sha256':self.record['image_sha256'],
               'mask_file':path.name,'mask_sha256':sha(path),'notes':self.notes.value,
               'track_ids':self.track_ids,'review_scope':self.review_scope.value,'workspace_roi_yxyx':roi,
               'proposal_sources':self.proposal_sources,'annotation_metadata_at_reset':self.meta,
               'review_regions_at_reset':self.review_regions,'status':'RESET_BACKUP_NOT_TRAINING',
               'pending_polygons':{k:np.asarray(v._overlays['polygon'].points).tolist() for k,v in self.layers.items()},
               'sam_positive_yx':self.positive.data.tolist(),'sam_negative_yx':self.negative.data.tolist(),
               'sam_boxes':[v.tolist() for v in self.box_layer.data]}
        # Serialize before updating the recovery pointer; a failed backup aborts reset.
        text=json.dumps(state,ensure_ascii=False,indent=2)
        (folder/(name+'.json')).write_text(text)
        if make_current:
            pointer=folder/'latest.tmp';pointer.write_text(text);pointer.replace(folder/'latest.json')
        return state

    def read_reset_snapshot(self):
        folder=self.reset_backup_dir();state=json.loads((folder/'latest.json').read_text())
        if state['candidate_id']!=self.record['id'] or state['image_sha256']!=self.record['image_sha256']:
            raise ValueError('重置备份属于另一张图；当前标签未改变')
        name=state['mask_file']
        if Path(name).name!=name:raise ValueError('重置备份路径无效')
        path=folder/name
        if sha(path)!=state['mask_sha256']:raise ValueError('重置备份校验失败；当前标签未改变')
        with np.load(path,allow_pickle=False) as z:masks={k:z[k].copy() for k in CLASSES}
        validate(masks,self.image.data.shape[:2])
        return state,masks

    def update_reset_recovery(self):
        self.buttons['restore_frame'].enabled=(self.reset_backup_dir()/'latest.json').exists()

    def confirm_frame_action(self,title,text):
        box=QMessageBox(self.viewer.window._qt_window);box.setWindowTitle(title)
        box.setIcon(QMessageBox.Icon.Warning);box.setText(text)
        apply=box.addButton('确认'+title,QMessageBox.ButtonRole.AcceptRole)
        cancel=box.addButton('取消',QMessageBox.ButtonRole.RejectRole)
        box.setDefaultButton(cancel);box.setEscapeButton(cancel);box.exec()
        return box.clickedButton() is apply

    def frame_has_label_content(self):
        return (any(np.any(v.data) for v in self.layers.values()) or bool(self.notes.value.strip())
                or bool(self.track_ids) or bool(self.review_regions)
                or bool((self.meta or {}).get('expert_confirmed')))

    def reload_after_frame_action(self):
        self.load_frame(self.index)
        # Old per-layer undo must not resurrect a discarded candidate preview.
        self.preview._reset_history()
        self.track.value='';self.instance.value=1
        self.category.value='SINGLE_BUBBLE';self.select_class();self.tabs.setCurrentIndex(0)
        self.sam_status.value='本帧候选、提示点和选框已清空；需要时可重新生成。'
        self.update_reset_recovery();self.focus_canvas()

    def reset_current_frame(self):
        if not self.confirm_frame_action('重置当前帧',
                f'只重置当前图片：{self.record["id"]}\n'
                '清空全部五类标签、候选、提示点、选框、未完成轮廓、备注、短轨迹 ID 和审核状态。\n'
                '重置会立即保存为空白草稿。原图、其他帧及历史修订保留；已画标签会先备份，可用旁边按钮恢复。'):
            return False
        # Repeated reset of an already empty frame must not replace its useful backup.
        if self.frame_has_label_content() or not (self.reset_backup_dir()/'latest.json').exists():
            self.reset_snapshot()
        self.cancel_sam();self.cancel_bulk()
        masks={k:np.zeros_like(v.data) for k,v in self.layers.items()}
        save_annotation(self.directory,self.record,masks,self.person.value,False,{},'',review_scope='REGIONS')
        self.reload_after_frame_action()
        self.status.value=f'{self.index+1}/{len(self.records)}  {self.record["id"]}\n当前帧已完全清空并保存为空白草稿；从单泡 #1 开始。'
        self.bulk_status.value='当前帧已重置。可重新框选，或直接用画笔标注；误操作可恢复重置前标签。'
        return True

    def restore_reset_frame(self):
        state,masks=self.read_reset_snapshot()
        if not self.confirm_frame_action('恢复重置前标签',
                f'恢复图片 {self.record["id"]} 最近一次重置前的五类标签、备注、短轨迹 ID 和审核选区。\n'
                '会替换当前草稿；当前内容也会另存备份。恢复内容作为未审核草稿，需要重新确认。\n'
                '未接受候选、SAM 提示点和未完成轮廓不重新加载。'):
            return False
        if self.frame_has_label_content():self.reset_snapshot(make_current=False)
        self.cancel_sam();self.cancel_bulk()
        save_annotation(self.directory,self.record,masks,self.person.value,False,state['track_ids'],state['notes'],
                        review_scope=state['review_scope'],workspace_roi=state['workspace_roi_yxyx'],
                        proposal_sources=state['proposal_sources'])
        self.reload_after_frame_action()
        self.status.value='已恢复重置前标签并保存为草稿；请重新检查和确认审核范围。'
        return True
