"""Local-region review and asynchronous bulk candidates for the napari workbench."""
from pathlib import Path
import json,sys,time,uuid
import numpy as np
from qtpy.QtCore import QProcess
from annotation_store import sha
from annotation_regions import normalize_box,valid_regions
from annotation_batch import merge_proposals,split_instance
from annotation_sizes import diameter_range,size_plan

ROOT=Path(__file__).resolve().parents[1]


class RegionalWorkflow:
    def init_regional_state(self):
        self.bulk_proc=None;self.bulk_job=None;self.bulk_buffer=b'';self.bulk_log=None
        self.bulk_result=None;self.bulk_signature=None;self.bulk_picking=False
        self.review_regions=[];self.proposal_sources=[];self.review_metadata_dirty=False
        self.split_target=None

    def init_regional_layers(self,shape):
        self.bulk_preview=self.viewer.add_labels(np.zeros(shape,np.uint32),name='批量候选（未审核）',opacity=.45)
        self.bulk_preview.brush_size=5;self.bulk_preview.preserve_labels=False
        self.bulk_preview.bind_key('Control-Z',lambda _:self.bulk_preview.undo(),overwrite=True)
        self.bulk_preview.bind_key('Control-Shift-Z',lambda _:self.bulk_preview.redo(),overwrite=True)
        self.bulk_preview.bind_key('E',lambda _:self.safe(self.erase_candidate_mode),overwrite=True)
        self.bulk_preview.bind_key('P',lambda _:self.safe(lambda:self.edit_mode('paint')),overwrite=True)
        self.roi_layer=self.viewer.add_shapes(name='本次审核区域',ndim=2,shape_type='rectangle',
                                             edge_color='#54d9fb',face_color='transparent',edge_width=1)
        self.split_points=self.viewer.add_points(np.empty((0,2)),name='拆分：两个泡各一点',size=5,face_color='#ffcc55')
        self.roi_layer.events.data.connect(self.roi_changed)
        self.viewer.mouse_drag_callbacks.append(self.bulk_pick_callback)

    def current_roi(self):
        if len(self.roi_layer.data)!=1:raise ValueError('请先拖出一个小矩形区域；一次只审核一个区域')
        pts=np.asarray(self.roi_layer.data[0]);shape=self.image.data.shape[:2]
        lo=np.maximum(0,np.floor(pts.min(axis=0)).astype(int))
        hi=np.minimum(shape,np.ceil(pts.max(axis=0)).astype(int))
        return normalize_box([*lo,*hi],shape)

    def roi_changed(self,*_):
        if self.loading:return
        self.cancel_bulk();self.clear_bulk()
        self.reviewed.value=False
        try:
            y0,x0,y1,x1=self.current_roi()
            self.bulk_status.value=f'区域 {x1-x0} × {y1-y0} px。先生成候选，再逐一检查漏标 / 误标。'
        except ValueError:self.bulk_status.value='拖出约 128–256 px 的小区域；建议从清晰、较稀疏处开始。'

    def roi_mode(self):
        self.check_pending_polygon();self.picking=False;self.bulk_picking=False
        self.review_scope.value='REGIONS';self.roi_layer.data=[]
        self.viewer.layers.selection.active=self.roi_layer;self.roi_layer.mode='add_rectangle'
        if not self.roi_layer.visible:self.roi_layer.visible=True
        self.status.value='拖出一个小矩形；只需完整检查框内对象，框外可以暂不标注。';self.focus_canvas()

    def zoom_roi(self):
        y0,x0,y1,x1=self.current_roi()
        self.viewer.scene.camera.center=(0,(y0+y1)/2,(x0+x1)/2)
        self.viewer.scene.camera.zoom=min(self.viewer.window._qt_viewer.canvas.size)/max(y1-y0,x1-x0)*.85

    def scope_changed(self,*_):
        if self.loading:return
        self.reviewed.value=False
        region=self.review_scope.value=='REGIONS'
        self.reviewed.label='我已检查框内全部对象、漏标和背景' if region else '我已检查整幅图，包括漏标和背景'
        self.buttons['confirm'].text='确认当前区域' if region else '确认整幅图'
        self.reviewed.native.setText(self.reviewed.label)

    def revoke_review(self):
        # Deliberate metadata-only change, preserving masks and previous immutable revisions.
        if self.review_scope.value=='REGIONS':
            box=self.current_roi();before=len(self.review_regions)
            self.review_regions=[r for r in self.review_regions if r['box_yxyx']!=box]
            if before==len(self.review_regions):raise ValueError('此框没有已确认记录；可从下拉列表恢复原审核区域')
        else:self.review_regions=[]
        self.review_metadata_dirty=True;self.mark_dirty()
        self.status.value='已撤回所选审核状态，像素标签保留；请保存草稿。'

    def update_region_choices(self):
        entries=valid_regions({'review_regions':self.review_regions},self.frame_arrays())
        self.region_history.choices=[('选择已确认区域','')]+[(f'区域 {i+1}：{r["box_yxyx"]}',r['id']) for i,r in enumerate(entries)]

    def restore_review_region(self):
        record=next((r for r in self.review_regions if r['id']==self.region_history.value),None)
        if record is None:raise ValueError('请先选择一个已确认区域')
        y0,x0,y1,x1=record['box_yxyx']
        self.roi_layer.data=[np.array([[y0,x0],[y0,x1],[y1,x1],[y1,x0]])]
        self.review_scope.value='REGIONS';self.zoom_roi()

    def current_diameter_range(self):
        return diameter_range([self.bulk_diameter_min.value,self.bulk_diameter_max.value])

    def current_diameter_filter(self):
        return 'STRICT_RANGE' if self.bulk_strict_size.value else 'REFERENCE_ONLY'

    def update_diameter_help(self):
        strict=self.bulk_strict_size.value
        model='cyto3 按小、中、大三个尺度检测' if self.bulk_model.value=='cyto3' else '此 SAM 模型使用原生检测尺度'
        policy='严格模式：筛除范围外候选。' if strict else '默认仅作参考：范围外候选保留，需人工核查。'
        self.diameter_help.value='原图像素；直径 = 2 × 半径。'+model+'。\n'+policy
        self.bulk_strict_size.native.setToolTip('关闭后不按所填范围删除候选；仍排除不足 5 像素、触及选区边缘或占选区过半的对象。')

    def bulk_prompt_signature(self):
        return json.dumps([self.record['id'],self.record['image_sha256'],self.current_roi(),
                           self.bulk_model.value,self.current_diameter_range(),self.current_diameter_filter()])

    def clear_bulk(self):
        self.bulk_result=None;self.bulk_signature=None;self.bulk_picking=False
        if hasattr(self,'bulk_preview'):
            self.bulk_preview.mode='pan_zoom'
            self.bulk_preview.data=np.zeros(self.image.data.shape[:2],np.uint32);self.bulk_preview._reset_history()

    def bulk_parameters_changed(self,*_):
        self.update_diameter_help()
        if self.loading:return
        self.cancel_bulk();self.clear_bulk()
        try:
            lo,hi=self.current_diameter_range()
            policy='严格筛选' if self.bulk_strict_size.value else '仅作参考，保留范围外候选'
            self.bulk_status.value=f'直径 {lo:g}–{hi:g} px：{policy}。已有标签保留，请重新生成。'
        except ValueError as e:self.bulk_status.value=str(e)

    def generate_bulk(self):
        self.check_pending_polygon()
        if self.args.manual_only:raise ValueError('当前为纯手工模式，模型未启用')
        box=self.current_roi();y0,x0,y1,x1=box
        if min(y1-y0,x1-x0)<32 or max(y1-y0,x1-x0)>768:raise ValueError('批量候选要求区域边长为 32–768 px，推荐 128–256 px')
        limits=self.current_diameter_range()
        self.cancel_sam();self.clear_bulk();self.cancel_bulk(keep_worker=True)
        token=uuid.uuid4().hex
        folder=ROOT/('evidence/annotation_batch_jobs' if self.args.smoke else 'results/annotation_batch_jobs')/token
        folder.mkdir(parents=True)
        req={'token':token,'image':str(ROOT/self.record['image']),'image_id':self.record['id'],
             'image_sha256':self.record['image_sha256'],'roi_yxyx':box,'model':self.bulk_model.value,
             'diameter_range_px':limits,'diameter_filter':self.current_diameter_filter(),
             'output':str(folder/'proposals.npz'),'software_test':self.args.smoke}
        (folder/'request.json').write_text(json.dumps(req,indent=2))
        signature=self.bulk_prompt_signature()
        self.bulk_job={'request':req,'signature':signature,'started':time.monotonic(),'folder':folder}
        env='.venv_preannotation' if req['model']=='cellpose4' else '.venv' if req['model']=='cyto3' else '.venv_annotation'
        if self.bulk_proc is not None and (getattr(self,'bulk_worker_model',None)!=req['model'] or self.bulk_proc.state()==QProcess.ProcessState.NotRunning):
            proc=self.bulk_proc;self.bulk_proc=None
            proc.readyReadStandardOutput.disconnect();proc.finished.disconnect();proc.errorOccurred.disconnect()
            proc.kill();proc.deleteLater()
        if self.bulk_proc is None:
            self.bulk_proc=QProcess();self.bulk_proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
            self.bulk_proc.setWorkingDirectory(str(ROOT));self.bulk_proc.readyReadStandardOutput.connect(self.read_bulk)
            self.bulk_proc.finished.connect(self.bulk_finished);self.bulk_proc.errorOccurred.connect(self.bulk_process_error)
            self.bulk_proc.start(str(ROOT/env/'bin/python'),['-u','-B',str(ROOT/'scripts/preannotation_worker.py'),'--server'])
            self.bulk_worker_model=req['model']
        self.bulk_log=(folder/'worker.log').open('ab');self.bulk_buffer=b''
        self.bulk_proc.write((json.dumps(req)+'\n').encode())
        self.bulk_status.value='后台生成中；可继续画笔标注。首次加载较慢，可随时取消。'
        self.status.value=f'正在为直径 {limits[0]:g}–{limits[1]:g} px 生成待审候选；已有标签保留。'
        self.buttons['generate_bulk'].enabled=False

    def read_bulk(self):
        proc=self.bulk_proc
        if proc is None:return
        chunk=bytes(proc.readAllStandardOutput())
        if self.bulk_log:self.bulk_log.write(chunk);self.bulk_log.flush()
        self.bulk_buffer+=chunk
        while b'\n' in self.bulk_buffer:
            line,self.bulk_buffer=self.bulk_buffer.split(b'\n',1)
            try:message=json.loads(line)
            except (ValueError,UnicodeDecodeError):continue
            if not self.bulk_job or message.get('token')!=self.bulk_job['request']['token']:continue
            if message.get('event')=='progress':self.bulk_status.value=message['status']
            elif message.get('event')=='error':
                self.bulk_status.value='生成失败，已有标注保留。详见 '+str(self.bulk_job['folder']/'worker.log')
                self.finish_bulk_job('FAILED')
            elif message.get('event')=='result':
                try:self.receive_bulk(message)
                except Exception as e:
                    self.bulk_status.value='候选未加载：'+str(e);self.finish_bulk_job('REJECTED_RESULT')

    def receive_bulk(self,result):
        job=self.bulk_job;req=job['request']
        if self.bulk_prompt_signature()!=job['signature']:raise ValueError('图片或选区已改变，请重新生成')
        output=Path(req['output'])
        if sha(output)!=result['output_sha256'] or result['image_sha256']!=self.record['image_sha256']:raise ValueError('结果校验失败')
        if result['roi_yxyx']!=req['roi_yxyx']:raise ValueError('结果区域不一致')
        if result.get('diameter_range_px')!=req['diameter_range_px']:raise ValueError('结果直径范围不一致，请重新生成')
        if result.get('size_plan')!=size_plan(req):raise ValueError('结果直径筛选方式不一致，请重新生成')
        with np.load(output,allow_pickle=False) as z:local=z['labels']
        y0,x0,y1,x1=req['roi_yxyx']
        if local.shape!=(y1-y0,x1-x0):raise ValueError('结果尺寸不一致')
        data=np.zeros(self.image.data.shape[:2],np.uint32);data[y0:y1,x0:x1]=local
        self.bulk_preview.data=data;self.bulk_preview._reset_history()
        if not self.bulk_preview.visible:self.bulk_preview.visible=True
        self.bulk_signature=job['signature'];self.bulk_result=dict(result,request_path=str(job['folder']/'request.json'))
        duplicates=sum('duplicate_across_scales' in r['reasons'] for r in result['rejected'])
        outside=sum(o.get('outside_reference_range',False) for o in result['objects'])
        extra=f'其中 {outside} 个超出参考直径，已保留待核查。' if not self.bulk_strict_size.value else '已严格按直径筛选。'
        self.bulk_status.value=f'{len(result["objects"])} 个待审候选；{extra}\n去重 {duplicates} 条，另排除 {len(result["rejected"])-duplicates} 条质量 / 边界 / 重叠记录。请检查并补漏。'
        if self.status.value.startswith('正在为直径 '):
            self.status.value='候选已生成；请检查后纳入草稿，再补漏、修正并审核。'
        self.finish_bulk_job('PROPOSALS_READY')

    def finish_bulk_job(self,outcome):
        if self.bulk_job:
            (self.bulk_job['folder']/'gui_outcome.json').write_text(json.dumps({'outcome':outcome,'wall_s':time.monotonic()-self.bulk_job['started']},indent=2))
        if self.bulk_log:self.bulk_log.close();self.bulk_log=None
        self.bulk_job=None;self.buttons['generate_bulk'].enabled=True

    def cancel_bulk(self,keep_worker=False):
        active=self.bulk_job is not None
        if self.bulk_proc and (active or not keep_worker):
            proc=self.bulk_proc;self.bulk_proc=None
            try:proc.readyReadStandardOutput.disconnect();proc.finished.disconnect();proc.errorOccurred.disconnect()
            except (TypeError,RuntimeError):pass
            proc.kill();proc.deleteLater()
        if active:self.finish_bulk_job('CANCELLED');self.bulk_status.value='批量计算已取消，已有标签保留。'

    def bulk_finished(self,*_):
        self.read_bulk()
        if self.bulk_job:
            self.bulk_status.value='批量进程退出，已有标签保留。请查看本次 worker.log。';self.finish_bulk_job('PROCESS_EXIT')

    def bulk_process_error(self,*_):
        if self.bulk_job:
            self.bulk_status.value='模型进程未能运行，已有标签保留。';self.finish_bulk_job('PROCESS_ERROR')

    def bulk_watchdog(self):
        if self.bulk_job and time.monotonic()-self.bulk_job['started']>180:
            self.cancel_bulk();self.bulk_status.value='已超过 180 秒并停止；缩小区域后重试。'

    def drop_candidate_mode(self):
        if self.bulk_result is None:raise ValueError('请先生成候选')
        self.picking=False;self.bulk_picking=True;self.bulk_preview.mode='pan_zoom'
        self.viewer.layers.selection.active=self.bulk_preview
        self.bulk_status.value='单击删除整只候选；“撤销候选编辑”可恢复。';self.focus_canvas()

    def erase_candidate_mode(self):
        self.check_pending_polygon()
        if self.bulk_result is None or self.bulk_signature!=self.bulk_prompt_signature():raise ValueError('请先生成当前区域的批量候选')
        self.picking=False;self.bulk_picking=False;layer=self.bulk_preview
        layer.visible=True;layer.preserve_labels=False;layer.mode='erase'
        self.viewer.layers.selection.active=layer;self.tabs.setCurrentIndex(0)
        self.bulk_status.value='正在擦除批量候选的局部像素；已纳入草稿的标签不变。Ctrl+Z 撤销；修好后点“纳入草稿”。'
        self.update_brush_help();self.focus_canvas()

    def candidate_edit_provenance(self,box):
        """Keep model output immutable and record the exact manually edited mask."""
        request_path=Path(self.bulk_result['request_path']);request=json.loads(request_path.read_text())
        if sha(request['output'])!=self.bulk_result['output_sha256']:raise ValueError('原始候选校验失败；请重新生成')
        with np.load(request['output'],allow_pickle=False) as z:original=z['labels']
        y0,x0,y1,x1=box;edited=self.bulk_preview.data[y0:y1,x0:x1]
        if np.array_equal(original,edited):return None
        path=request_path.parent/f'manual_candidates_{uuid.uuid4().hex}.npz'
        np.savez_compressed(path,labels=edited)
        return {'status':'MANUALLY_EDITED_UNREVIEWED_CANDIDATES','mask_path':str(path),'mask_sha256':sha(path),
                'original_output_sha256':self.bulk_result['output_sha256'],
                'changed_pixels':int(np.count_nonzero(original!=edited)),
                'erased_pixels':int(np.count_nonzero((original>0)&(edited==0)))}

    def bulk_pick_callback(self,viewer,event):
        if not self.bulk_picking:return
        y,x=np.floor(np.asarray(self.image.world_to_data(event.position)[-2:])+.5).astype(int)
        if not (0<=y<self.bulk_preview.data.shape[0] and 0<=x<self.bulk_preview.data.shape[1]):return
        ident=int(self.bulk_preview.data[y,x])
        if ident:self.bulk_preview.data_setitem(np.where(self.bulk_preview.data==ident),0)
        yield

    def incorporate_bulk(self):
        if self.bulk_result is None or self.bulk_signature!=self.bulk_prompt_signature():raise ValueError('没有当前区域的有效候选，请重新生成')
        box=self.current_roi();y0,x0,y1,x1=box
        result,accepted,skipped=merge_proposals(self.frame_arrays(),self.bulk_preview.data[y0:y1,x0:x1],box)
        if not accepted:raise ValueError('没有可加入的对象；已有标签重叠的候选整只跳过，不覆盖手工标注')
        edits=self.candidate_edit_provenance(box)
        self.category.value='SINGLE_BUBBLE';layer=self.layers['SINGLE_BUBBLE']
        changed=result!=layer.data;layer.data_setitem(np.where(changed),result[changed])
        self.proposal_sources.append({'request_path':self.bulk_result['request_path'],'model':self.bulk_result['model'],
                                     'model_sha256':self.bulk_result['model_sha256'],'output_sha256':self.bulk_result['output_sha256'],
                                     'diameter_range_px':self.bulk_result['diameter_range_px'],
                                     'size_plan':self.bulk_result['size_plan'],'merge_policy':self.bulk_result['merge_policy'],
                                     'manual_candidate_edits':edits,
                                     'accepted':accepted,'skipped_overlapping':skipped,'status':'ADDED_AS_DRAFT'})
        self.clear_bulk();self.mark_dirty();self.new_bubble()
        self.bulk_status.value=f'已纳入 {len(accepted)} 个草稿对象，跳过 {len(skipped)} 个重叠候选。N 补画漏泡；点选 / 删除误标；可一次撤销本批。'

    def start_split(self):
        self.select_instance();self.split_target=(self.category.value,self.instance.value)
        self.split_points.data=np.empty((0,2));self.viewer.layers.selection.active=self.split_points
        self.split_points.mode='add';self.picking=False;self.bulk_picking=False
        self.status.value='在粘连对象的两个泡内部各点一次，然后点“完成拆分”。';self.focus_canvas()

    def finish_split(self):
        if self.split_target is None:raise ValueError('先点选粘连对象，再点“两点拆分”')
        category,ident=self.split_target;layer=self.layers[category]
        result,new_id=split_instance(layer.data,ident,self.split_points.data)
        self.category.value=category;changed=result!=layer.data
        layer.data_setitem(np.where(changed),result[changed]);self.split_points.data=np.empty((0,2));self.split_target=None
        self.mark_dirty();self.instance.value=new_id;self.select_instance()
        self.status.value=f'已拆成 #{ident} 和 #{new_id}；原有像素范围不变，可撤销。请检查分界是否正确。'
