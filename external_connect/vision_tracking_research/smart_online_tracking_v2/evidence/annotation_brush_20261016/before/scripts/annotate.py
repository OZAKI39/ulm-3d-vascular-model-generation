"""Persistent napari Labels workbench; optional original SAM in a cancellable process.

The former synchronous micro-sam window is retained in evidence/annotation_repair_20261016.
No custom drawing/segmentation algorithm; napari provides polygon, brush, undo and zoom.
"""
from pathlib import Path
import argparse,sys,json,os,time,uuid,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
os.environ.setdefault('QT_API','pyqt6');os.environ.setdefault('QT_QPA_PLATFORM','xcb')
os.environ.setdefault('PYOPENGL_PLATFORM','glx');os.environ.setdefault('MICROSAM_CACHEDIR',str(ROOT/'models/microsam_cache'))
import numpy as np,imageio.v3 as iio
from annotation_store import CLASSES,load_annotation,save_annotation,sha
from sam_prompt_geometry import prompt_region,point_instance_requests
from sam_batch_review import assemble_candidates
import sam_auto_training as auto_sam
import napari
from napari.layers.labels._labels_key_bindings import complete_polygon,reset_polygon
from magicgui.widgets import Container,Label,ComboBox,LineEdit,CheckBox,PushButton,SpinBox
from qtpy.QtCore import QTimer,QProcess,QObject,QEvent,QLockFile,Qt
from qtpy.QtWidgets import QTabWidget,QApplication,QScrollArea,QSizePolicy
from qtpy.QtGui import QFont

NAMES={'SINGLE_BUBBLE':'单泡','AGGREGATE':'微泡团','BACKGROUND_ARTIFACT':'背景 / 杂质',
       'UNCERTAIN':'不确定','OUT_OF_FOCUS':'失焦'}
COLORS={'SINGLE_BUBBLE':'#3bc4ef','AGGREGATE':'#ffab4a','BACKGROUND_ARTIFACT':'#be9bf3',
        'UNCERTAIN':'#ffe65b','OUT_OF_FOCUS':'#aab3be'}

def label(text):
    w=Label(value=text);w.native.setWordWrap(True);return w

def row(*widgets):return Container(widgets=list(widgets),layout='horizontal',labels=False)

class CloseGuard(QObject):
    def __init__(self,workbench):super().__init__();self.w=workbench
    def eventFilter(self,obj,event):
        if event.type()==QEvent.Type.Close:
            if not self.w.save_if_changed():event.ignore();return True
            self.w.closing=True
            self.w.cancel_sam();self.w.timer.stop();self.w.save_session()
            self.w.watch.stop();self.w.training_timer.stop()
            if self.w.auto_training_allowed:auto_sam.unregister_client()
        return False

class Workbench:
    def __init__(self,args):
        self.args=args;self.records=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows']
        self.directory=Path(args.annotation_dir) if getattr(args,'annotation_dir',None) else ROOT/('evidence/annotation_repair_gui_sandbox' if args.smoke else 'data/annotations')
        if args.smoke and not self.directory.resolve().is_relative_to((ROOT/'evidence').resolve()):
            raise ValueError('软件测试只能写入独立 evidence 沙盒，禁止指向正式标注目录')
        self.directory.mkdir(parents=True,exist_ok=True)
        self.sessionfile=self.directory/'session.json';previous=json.loads(self.sessionfile.read_text()) if self.sessionfile.exists() else {}
        self.index=args.index if args.index is not None else previous.get('index',0)
        if not 0<=self.index<len(self.records):raise ValueError('Frame index out of range')
        self.viewer=napari.Viewer(title='微泡人工标注 — 正在打开图片')
        self.viewer.window._qt_window.resize(1680,1000)
        self.layers={};self.loading=True;self.dirty=False;self.meta=None;self.saved_masks=None
        self.closing=False;self.restoring=False;self.picking=False;self.last_sam_instance=None
        self.sam_masks=None;self.sam_details=None;self.preview_signature=None;self.stdout_buffer=b''
        self.sam_batch=None;self.batch_choices=[];self.batch_included=[];self.updating_sam_widgets=False
        self.last_sam_batch=None
        self.track_ids={};self.proc=None;self.job=None;self.logfile=None;self.job_started=None
        self.last_sam_outcome=None;self.last_sam_wall_s=None;self.autosave_count=0
        self.status=label('正在打开原始图片…');self.instructions=label('选类别 → 新建轮廓 → 沿边界点击 → 完成轮廓 → 保存。\n每个新对象都先点“新建轮廓”。')
        self.category=ComboBox(label='① 对象类别',choices=[(f'{NAMES[k]}  |  {k}',k) for k in CLASSES])
        self.person=LineEdit(label='审核人',value=previous.get('reviewer',''))
        self.notes=LineEdit(label='本帧备注',value='')
        self.reviewed=CheckBox(label='我已检查整幅图，包括漏标和背景',value=False)
        self.jump=SpinBox(label='图片序号',min=1,max=len(self.records),value=self.index+1)
        self.instance=SpinBox(label='要修改的实例编号',min=1,max=100000,value=1)
        self.track=LineEdit(label='短轨迹 ID（选填）',value='')
        self.sam_status=label('SAM 是可选辅助。手工标注不需要加载模型。')
        self.auto_training_allowed=not args.smoke and not args.manual_only and not getattr(args,'no_auto_train',False)
        self.auto_training=CheckBox(label='自动微调（人工确认后后台训练）',value=self.auto_training_allowed and auto_sam.config()['enabled'])
        self.auto_training.enabled=self.auto_training_allowed
        self.training_status=label('自动训练正在检查已确认标注。' if self.auto_training_allowed else '自动训练已关闭（手工模式或软件测试）。')
        self.selection_status=label('当前未选中已有对象；可点选标签后删除 / 修边。')
        self.target_width=SpinBox(label='每个单泡的局部框宽 (px)',min=6,max=256,value=16)
        self.proposal=ComboBox(label='候选轮廓',choices=[('尚未生成',-1)])
        self.sam_semantics=ComboBox(label='提示方式',choices=[('每点一个单泡（可连续点多个泡）','point_instances'),('一个对象：框选 / 多点修边','one_object')])
        self.sam_help=label('单泡模式：每个绿点对应一个独立泡，可连续点多个泡。大框仅限定选区。')
        self.prompt_count=label('0 个绿点 → 0 个独立单泡')
        self.prompt_history=ComboBox(label='历史提示',choices=[('本图暂无历史','')])
        self.batch_target=ComboBox(label='查看哪个泡',choices=[('尚未生成',-1)])
        self.batch_keep=CheckBox(label='保留当前泡，纳入批量接受',value=False)
        self.buttons={}
        def button(key,text,fn):
            b=PushButton(text=text);b.clicked.connect(lambda *_:self.safe(fn));self.buttons[key]=b;return b
        main=Container(widgets=[self.status,self.instructions,self.category,
            row(button('new','② 新建轮廓',self.new_object),button('finish','③ 完成轮廓',self.finish_polygon)),
            row(button('brush','画笔修边',lambda:self.edit_mode('paint')),button('erase','橡皮擦',lambda:self.edit_mode('erase'))),
            row(button('undo','撤销（含 SAM / 删除）',self.undo),button('redo','重做',self.redo)),
            self.selection_status,row(button('pick','点选已有对象',self.pick_mode),button('delete_main','删除选中对象',self.delete_instance)),
            row(button('restore','恢复 / 显示所有图层',self.restore_layers),button('pan','移动 / 缩放',self.pan)),
            row(button('fit','查看整幅图',self.viewer.reset_view),button('clear_polygon','取消未完成轮廓',self.clear_polygon)),
            self.person,self.notes,button('save','④ 保存草稿（Ctrl+S）',lambda:self.save(False)),
            self.reviewed,button('confirm','⑤ 确认本帧人工审核完成',lambda:self.save(True)),
            self.auto_training,self.training_status,button('original_sam','恢复原版 SAM',self.restore_original_sam),
            row(button('prev','上一张',lambda:self.go(self.index-1)),button('next','⑥ 保存并下一张',lambda:self.go(self.index+1))),
            row(self.jump,button('jump','跳转',lambda:self.go(self.jump.value-1))),
            label('草稿每 30 秒保存；换帧 / 关闭时保存修改。\n左右键不换帧；Alt+← / → 换帧。\n空格+拖动平移，滚轮缩放。')],labels=True)
        sam=Container(widgets=[self.sam_status,self.sam_semantics,self.sam_help,self.prompt_count,
            row(button('box','框选一个目标（拖动）',self.box_mode),button('clear_box','清除目标框',self.clear_box)),
            self.prompt_history,button('restore_prompts','恢复所选历史提示',self.restore_last_prompts),
            button('zoom_target','放大查看当前目标',self.zoom_target),
            self.target_width,
            button('positive','1. 放正点（对象内部）',lambda:self.prompt_mode(True)),
            button('negative','2. 放负点（排除区域，可选）',lambda:self.prompt_mode(False)),
            button('segment','3. 生成 / 更新候选',self.segment),button('cancel_sam','取消计算',self.cancel_sam),
            self.batch_target,self.proposal,self.batch_keep,button('accept','4. 接受当前泡 / 对象',self.accept_sam),
            button('accept_batch','接受全部有效单泡（分别编号）',self.accept_sam_batch),
            row(button('undo_sam','删除上次接受的一批 / 对象',self.remove_last_sam),button('sam_undo','撤销删除 / 编辑',self.undo)),
            row(button('sam_pick','点选已有对象',self.pick_mode),button('sam_delete','删除选中对象',self.delete_instance)),
            button('clear_prompts','清空候选和所有点',self.clear_prompts),
            label('多泡预览分别着色，接受后各自编号。\n尺寸异常或互相重叠的候选不纳入批量接受。\n预测轮廓不强行裁小；模型空闲 120 秒释放。')],labels=True)
        extra=Container(widgets=[label('可选：仅给连续 8 帧中能确实对应的对象设置短轨迹 ID。普通标注无需此步骤。'),
            self.instance,button('select','选中此编号进行修改',self.select_instance),
            self.track,button('track','保存短轨迹 ID 到本帧',self.set_track),
            button('delete','删除选中实例（可撤销）',self.delete_instance),
            label('同一类别中的编号对应一个对象。修改旧对象先选编号。\n颜色仅用于辅助辨认，不表示已完成专家审核。')])
        tabs=QTabWidget()
        for widget,title in [(main,'逐张标注'),(sam,'SAM 辅助'),(extra,'实例 / 短轨迹')]:
            widget.native.layout().setAlignment(Qt.AlignmentFlag.AlignTop)
            for child in widget:
                child.native.setSizePolicy(QSizePolicy.Policy.Preferred,QSizePolicy.Policy.Maximum)
            scroll=QScrollArea();scroll.setWidgetResizable(True);scroll.setWidget(widget.native);tabs.addTab(scroll,title)
        tabs.setMinimumWidth(385);tabs.setMaximumWidth(480)
        self.tabs=tabs;self.main=main;self.sam=sam;self.extra=extra
        self.sam_status.native.setMinimumHeight(76)
        self.training_status.native.setMinimumHeight(100)
        self.status.native.setMinimumHeight(68)
        self.selection_status.native.setMinimumHeight(38)
        self.viewer.window.add_dock_widget(tabs,name='标注步骤',area='right')
        self.category.changed.connect(lambda *_:self.category_changed())
        self.auto_training.changed.connect(lambda *_:self.toggle_auto_training())
        self.sam_semantics.changed.connect(lambda *_:self.update_prompt_help())
        self.batch_target.changed.connect(lambda *_:self.show_batch_target())
        self.batch_keep.changed.connect(lambda *_:self.change_batch_keep())
        self.notes.changed.connect(lambda *_:self.mark_dirty())
        self.proposal.changed.connect(lambda *_:self.show_proposal())
        self.viewer.bind_key('Control-S',lambda _:self.save(False),overwrite=True)
        self.viewer.bind_key('Alt-Right',lambda _:self.go(self.index+1),overwrite=True)
        self.viewer.bind_key('Alt-Left',lambda _:self.go(self.index-1),overwrite=True)
        self.guard=CloseGuard(self);self.viewer.window._qt_window.installEventFilter(self.guard)
        self.timer=QTimer();self.timer.timeout.connect(self.autosave);self.timer.start(30000)
        self.watch=QTimer();self.watch.timeout.connect(self.sam_watchdog);self.watch.start(1000)
        self.load_frame(self.index)
        self.training_timer=QTimer();self.training_timer.timeout.connect(self.poll_training)
        if self.auto_training_allowed:
            auto_sam.register_client();auto_sam.ensure_service();self.training_timer.start(2000)
        self.positive.events.data.connect(lambda *_:self.update_prompt_help())
        self.update_prompt_help()
        self.viewer.layers.events.removed.connect(self.layer_removed)
        self.viewer.mouse_drag_callbacks.append(self.pick_callback)
        for layer in self.layers.values():
            layer.bind_key('Control-Z',lambda _:self.undo(),overwrite=True)
            layer.bind_key('Control-Shift-Z',lambda _:self.redo(),overwrite=True)
        self.viewer.bind_key('Control-Z',lambda _:self.undo(),overwrite=True)
        self.viewer.bind_key('Control-Shift-Z',lambda _:self.redo(),overwrite=True)
    def safe(self,fn):
        try:return fn()
        except Exception as e:
            self.status.value='操作未完成：'+str(e)
            self.sam_status.value=self.status.value
            print(traceback.format_exc(),flush=True);return False
    def toggle_auto_training(self):
        if not self.auto_training_allowed:return
        auto_sam.set_enabled(self.auto_training.value);auto_sam.ensure_service();self.poll_training()
    def restore_original_sam(self):
        if not self.auto_training_allowed:return
        auto_sam.rollback();self.cancel_sam();self.invalidate_preview();self.poll_training()
        self.sam_status.value='已恢复原版 SAM；标签未改变，下次生成使用原版。'
    def poll_training(self):
        if self.closing:return
        state=auto_sam.read_json(auto_sam.RUNTIME/'status.json',{})
        names={'WAITING_FOR_CONFIRMED_LABELS':'等待训练集人工确认标注（至少 8 个单泡）',
               'QUEUED':'已排队，稍后自动训练','WAITING_FOR_GPU':'等待 GPU 空闲；可继续手工标注',
               'RUNNING':'后台微调中；可继续标注','PAUSING':'正在停止本轮训练','PAUSED':'自动训练已暂停',
               'UP_TO_DATE':'新模型已通过验证，下一次 SAM 请求自动使用',
               'CANDIDATE_HELD':'候选模型已保存，暂未替换当前模型',
               'LAST_RUN_FAILED':'本轮训练失败，现有模型保留','SERVICE_ERROR':'训练服务异常，现有标签保留',
               'STOPPED':'训练服务已停止'}
        text=names.get(state.get('state'),'正在检查训练状态')
        p=state.get('progress',{})
        stages={'LOADING_MODEL':'加载模型','ENCODING':'准备图像特征','TRAINING_FIT_BEFORE':'检查训练前拟合',
                'TRAINING':'训练','TRAINING_FIT_AFTER':'检查训练后拟合','VALIDATING':'独立验证','COMPLETED':'已完成'}
        if p:text+=f'；{stages.get(p.get("stage"),"处理中")} {p.get("step",0)}/{p.get("max_steps",80)}'
        last=state.get('last_attempt',{})
        if state.get('state')=='CANDIDATE_HELD' and last.get('reason')=='WAITING_FOR_INDEPENDENT_VALIDATION':text+='：需要至少 2 张验证帧 / 10 个单泡'
        current='原版' if state.get('active_model','ORIGINAL_SAM_VIT_B')=='ORIGINAL_SAM_VIT_B' else '已验证微调版'
        self.training_status.value=text+f'。\n当前 SAM：{current}。草稿不训练；未通过独立验证不自动换模型。'
    def required_layers(self):
        return [self.image,*self.layers.values(),self.preview,self.box_layer,self.positive,self.negative]
    def layer_removed(self,event):
        if not self.closing and not self.restoring and event.value in self.required_layers():
            QTimer.singleShot(0,self.restore_layers)
    def restore_layers(self):
        if self.closing or self.restoring:return
        self.restoring=True
        try:
            # Reinsert the SAME objects: unsaved pixels and napari undo history survive.
            for layer in self.required_layers():
                if layer not in self.viewer.layers:self.viewer.add_layer(layer)
                layer.visible=True
            for index,layer in enumerate(self.required_layers()):
                current=self.viewer.layers.index(layer)
                if current!=index:self.viewer.layers.move(current,index)
            self.status.value='必需图层已恢复 / 显示，标签和撤销记录保留。删除对象请用“删除选中对象”。'
        finally:self.restoring=False
    def pick_mode(self):
        self.restore_layers();self.picking=True
        for layer in self.layers.values():layer.mode='pan_zoom'
        self.positive.mode='pan_zoom';self.negative.mode='pan_zoom';self.box_layer.mode='pan_zoom'
        self.viewer.layers.selection.active=self.layers[self.category.value]
        self.status.value='请在已有彩色标签内部单击；选中后可删除或修边。'
        self.sam_status.value=self.status.value
    def pick_callback(self,viewer,event):
        if not self.picking:return
        pos=np.asarray(self.image.world_to_data(event.position)[-2:])
        y,x=np.floor(pos+.5).astype(int)
        if not (0<=y<self.image.data.shape[0] and 0<=x<self.image.data.shape[1]):return
        found=[(k,int(v.data[y,x])) for k,v in self.layers.items() if v.data[y,x]>0]
        if len(found)!=1:
            self.status.value='此处无标签，或有类别重叠；请点在目标标签内部。';return
        k,n=found[0];self.category.value=k;self.instance.value=n;self.select_instance();self.picking=False
        yield
    def mark_dirty(self,*_):
        if not self.loading:
            self.dirty=True;self.reviewed.value=False
            self.status.value=f'{self.index+1}/{len(self.records)}  {self.record["id"]}\n有修改，尚未保存 / 确认；可保存草稿。'
    def frame_arrays(self):return {k:v.data for k,v in self.layers.items()}
    def changed(self):
        if self.saved_masks is None:return any(np.any(v.data) for v in self.layers.values()) or bool(self.notes.value.strip()) or bool(self.track_ids)
        return any(not np.array_equal(self.saved_masks[k],v.data) for k,v in self.layers.items()) or self.notes.value!=(self.meta or {}).get('notes','') or self.track_ids!=(self.meta or {}).get('track_ids',{})
    def save_session(self):
        p=self.sessionfile.with_suffix('.tmp');p.write_text(json.dumps({'index':self.index,'reviewer':self.person.value,'candidate_id':self.record['id']},ensure_ascii=False,indent=2));p.replace(self.sessionfile)
    def save(self,confirm=False):
        try:
            if any(v._overlays['polygon'].points for v in self.layers.values()):
                raise ValueError('仍有未完成轮廓，请先点“完成轮廓”或“取消未完成轮廓”')
            if confirm and self.args.smoke:raise ValueError('软件测试禁止生成正式专家标签')
            if confirm and not self.reviewed.value:raise ValueError('请先检查整幅图，再勾选完整审核；不确定时只保存草稿')
            self.clean_track_ids()
            if not confirm and self.meta is not None and not self.changed():
                self.status.value='内容未改动，已保存的确认状态保留。';self.save_session();return True
            self.meta=save_annotation(self.directory,self.record,self.frame_arrays(),self.person.value,confirm,self.track_ids,self.notes.value)
            self.saved_masks={k:v.data.copy() for k,v in self.layers.items()};self.dirty=False;self.save_session()
            state='人工审核完成' if confirm else '草稿已保存（尚未确认）'
            counts=' / '.join(f'{NAMES[k]} {self.meta["instance_counts"][k]}' for k in CLASSES)
            self.status.value=f'{self.index+1}/{len(self.records)}  {self.record["id"]}\n{state}\n{counts}'
            return True
        except Exception as e:
            self.status.value='未保存：'+str(e)+'\n当前画面保留，请修正后重试。';return False
    def save_if_changed(self):
        pending=any(v._overlays['polygon'].points for v in self.layers.values())
        return self.save(False) if self.changed() or pending else True
    def autosave(self):
        if self.dirty and self.changed():
            if self.save(False):self.autosave_count+=1
    def load_frame(self,index):
        started=time.monotonic();self.loading=True;self.cancel_sam();self.index=index;self.record=self.records[index]
        image=iio.imread(ROOT/self.record['image']);meta,masks=load_annotation(self.directory,self.record)
        if not self.layers:
            self.image=self.viewer.add_image(image,name='原始图像',rgb=True)
            for k in CLASSES:
                layer=self.viewer.add_labels(np.zeros(image.shape[:2],np.uint32),name=k,opacity=.35)
                layer.brush_size=5;layer.preserve_labels=True
                layer.events.data.connect(self.mark_dirty);self.layers[k]=layer
            self.preview=self.viewer.add_labels(np.zeros(image.shape[:2],np.uint8),name='SAM 候选（未接受）',opacity=.45)
            self.box_layer=self.viewer.add_shapes(name='SAM 目标框',ndim=2,shape_type='rectangle',edge_color='#ffe46b',face_color='transparent',edge_width=.7)
            self.positive=self.viewer.add_points(np.empty((0,2)),name='SAM 正点',size=4,face_color='#36ff86')
            self.negative=self.viewer.add_points(np.empty((0,2)),name='SAM 负点',size=4,face_color='#ff5e7e')
        else:self.image.data=image
        for k,layer in self.layers.items():
            layer.data=masks[k].copy() if masks else np.zeros(image.shape[:2],np.uint32)
            layer._reset_history();layer._overlays['polygon'].points=[]
        self.clear_prompts();self.last_sam_instance=None;self.last_sam_batch=None;self.removed_tracks={};self.picking=False
        self.meta=meta;self.saved_masks={k:v.copy() for k,v in masks.items()} if masks else None
        self.track_ids=dict(meta.get('track_ids',{})) if meta else {};self.notes.value=meta.get('notes','') if meta else ''
        self.reviewed.value=False;self.jump.value=index+1;self.dirty=False;self.loading=False
        state='已确认；编辑后需重新确认' if meta and meta['expert_confirmed'] else '未确认 / 草稿'
        self.status.value=f'{index+1}/{len(self.records)}  {self.record["id"]}\n分组：{self.record["split"]}；{state}'
        self.viewer.title=f'微泡人工标注 | {index+1}/{len(self.records)} | {self.record["id"]}'
        self.refresh_prompt_history()
        self.select_class();self.viewer.reset_view();self.save_session()
        self.last_load_seconds=time.monotonic()-started
    def go(self,index):
        if index<0 or index>=len(self.records):self.status.value='已经到达第一张 / 最后一张。';return
        if self.save_if_changed():self.load_frame(index)
    def select_class(self):
        if not self.layers:return
        layer=self.layers[self.category.value];self.viewer.layers.selection.active=layer
        layer.selected_label=int(layer.data.max())+1;self.instance.value=layer.selected_label;layer.mode='pan_zoom'
    def category_changed(self):
        self.select_class()
        self.sam_semantics.value='point_instances' if self.category.value=='SINGLE_BUBBLE' else 'one_object'
        self.sam_semantics.enabled=self.category.value=='SINGLE_BUBBLE'
        self.update_prompt_help()
    def update_prompt_help(self):
        many=self.sam_semantics.value=='point_instances';n=len(self.positive.data) if hasattr(self,'positive') else 0
        self.prompt_count.value=f'{n} 个绿点 → '+(f'{n} 个独立单泡' if many else '同一个对象的提示')
        self.sam_help.value=('单泡：每个绿点对应一个泡，可连续点多个泡。黄色大框仅作选区，分割逐泡进行。'
                             if many else '一个对象：框和所有绿点共同描述同一个对象。多个分散泡请切换“每点一个单泡”。')
        if 'box' in self.buttons:self.buttons['box'].text='画选区框（可选）' if many else '框选一个目标（拖动）'
    def refresh_prompt_history(self):
        candidates=sorted((ROOT/'results/annotation_jobs').glob('*/request.json'),key=lambda p:p.stat().st_mtime,reverse=True)
        choices=[]
        for path in candidates:
            req=json.loads(path.read_text())
            if req.get('image_id')!=self.record['id'] or req.get('software_test',False):continue
            stamp=time.strftime('%H:%M:%S',time.localtime(path.stat().st_mtime))
            choices.append((f'{stamp} | {len(req["positive_yx"])} 个绿点 | {path.parent.name[:6]}',str(path)))
        self.prompt_history.choices=choices or [('本图暂无历史','')]
    def restore_last_prompts(self):
        path=Path(self.prompt_history.value) if self.prompt_history.value else None
        if path is not None:
            req=json.loads(path.read_text())
            if req.get('image_id')!=self.record['id']:raise ValueError('历史提示属于其他图片')
            self.clear_prompts();self.positive.data=np.asarray(req['positive_yx'],float).reshape(-1,2)
            self.negative.data=np.asarray(req['negative_yx'],float).reshape(-1,2)
            if req.get('box_yxyx'):
                y0,x0,y1,x1=req['box_yxyx'];self.box_layer.data=[np.array([[y0,x0],[y0,x1],[y1,x1],[y1,x0]])]
            self.target_width.value=req.get('target_width_px',16);self.update_prompt_help()
            self.sam_status.value=f'已恢复本图 {len(self.positive.data)} 个绿点及原选区；请重新生成，未接受任何标签。';return
        raise ValueError('本图没有可恢复的历史 SAM 请求')
    def new_object(self):
        self.restore_layers();self.select_class();layer=self.layers[self.category.value];layer.mode='polygon';self.picking=False
        if self.last_sam_instance==(self.category.value,layer.selected_label):self.last_sam_instance=None
        if self.last_sam_batch and self.last_sam_batch[0]==self.category.value and layer.selected_label in self.last_sam_batch[1]:self.last_sam_batch=None
        self.removed_tracks.pop(f'{self.category.value}:{layer.selected_label}',None)
        self.status.value=f'{NAMES[self.category.value]}，新实例 {layer.selected_label}：沿边界点击至少 3 点，然后点“完成轮廓”。'
    def finish_polygon(self):
        layer=self.layers[self.category.value]
        if len(layer._overlays['polygon'].points)<3:self.status.value='尚未形成轮廓：先点“新建轮廓”，再沿边界点击至少 3 点。';return
        complete_polygon(layer);self.mark_dirty();self.status.value='轮廓已写入当前类别。下一个对象先点“新建轮廓”；可保存草稿。'
    def clear_polygon(self):reset_polygon(self.layers[self.category.value])
    def pan(self):
        self.viewer.layers.selection.active=self.layers[self.category.value];self.layers[self.category.value].mode='pan_zoom'
    def edit_mode(self,mode):
        layer=self.layers[self.category.value];self.viewer.layers.selection.active=layer;layer.mode=mode
        self.status.value=f'正在编辑 {NAMES[self.category.value]} 编号 {layer.selected_label}。修改旧对象先在“实例 / 短轨迹”选择编号。'
    def undo(self):
        layer=self.layers[self.category.value];layer.undo()
        for key,value in list(self.removed_tracks.items()):
            k,n=key.split(':')
            if np.any(self.layers[k].data==int(n)):self.track_ids[key]=value;del self.removed_tracks[key]
        self.mark_dirty();self.status.value='已撤销当前类别的最近一次编辑（含 SAM 接受 / 删除）。'
    def redo(self):
        self.layers[self.category.value].redo();self.clean_track_ids();self.mark_dirty()
        self.status.value='已重做当前类别的最近一次编辑。'
    def clean_track_ids(self):
        for key,value in list(self.track_ids.items()):
            k,n=key.split(':')
            if not np.any(self.layers[k].data==int(n)):
                self.removed_tracks[key]=value;del self.track_ids[key]
    def select_instance(self):
        layer=self.layers[self.category.value]
        if not np.any(layer.data==self.instance.value):raise ValueError('此类别中不存在这个编号')
        layer.selected_label=self.instance.value;self.viewer.layers.selection.active=layer;layer.mode='pan_zoom';layer.visible=True
        self.selection_status.value=f'已选：{NAMES[self.category.value]} #{self.instance.value}；可删除 / 修边。'
        self.status.value=self.selection_status.value;self.sam_status.value=self.selection_status.value
    def delete_instance(self):
        self.select_instance();layer=self.layers[self.category.value]
        # Direct history-aware edit bypasses preserve_labels, including disconnected pixels.
        n=layer.selected_label;layer.data_setitem(np.where(layer.data==n),0)
        if self.last_sam_instance==(self.category.value,n):self.last_sam_instance=None
        if self.last_sam_batch and self.last_sam_batch[0]==self.category.value:
            ids=[i for i in self.last_sam_batch[1] if i!=n]
            self.last_sam_batch=(self.category.value,ids) if ids else None
        self.clean_track_ids();self.mark_dirty()
        self.selection_status.value=f'已删除 {NAMES[self.category.value]} #{n}；可点“撤销”恢复。'
        self.status.value=self.selection_status.value;self.sam_status.value=self.selection_status.value
    def remove_last_sam(self):
        if self.last_sam_batch is not None:
            k,ids=self.last_sam_batch;layer=self.layers[k];mask=np.isin(layer.data,ids)
            if not np.any(mask):raise ValueError('上次接受的标签已撤销或删除；可用“重做”恢复')
            self.category.value=k;layer.data_setitem(np.where(mask),0)
            self.last_sam_batch=None;self.last_sam_instance=None;self.clean_track_ids();self.mark_dirty()
            self.sam_status.value=f'已删除上次接受的 {len(ids)} 个 SAM 单泡；可一次撤销恢复。';return
        if self.last_sam_instance is None:raise ValueError('本帧还没有接受过 SAM 对象；已有标签可用“点选已有对象”删除')
        k,n=self.last_sam_instance;self.category.value=k;self.instance.value=n;self.delete_instance()
    def set_track(self):
        self.select_instance();key=f'{self.category.value}:{self.instance.value}'
        if self.track.value.strip():self.track_ids[key]=self.track.value.strip()
        else:self.track_ids.pop(key,None)
        self.mark_dirty();self.status.value='短轨迹 ID 已写入草稿；仍需保存 / 确认。'
    def prompt_mode(self,positive):
        self.picking=False;self.restore_layers()
        layer=self.positive if positive else self.negative;self.viewer.layers.selection.active=layer;layer.mode='add'
        self.sam_status.value=('每个泡内部放一个绿点；可连续点多个泡。红点排除误选背景。' if self.sam_semantics.value=='point_instances' else '绿点描述同一个对象内部，红点排除紧邻误选区域。')+'点完后重新生成。'
    def box_mode(self):
        if self.sam_semantics.value=='point_instances':self.clear_box()
        else:self.clear_prompts()
        self.picking=False;self.restore_layers()
        self.viewer.layers.selection.active=self.box_layer;self.box_layer.mode='add_rectangle'
        self.sam_status.value=('拖出选区框，框内每个绿点仍单独分割。无需画框也可逐泡点选。' if self.sam_semantics.value=='point_instances' else '拖出紧贴一个对象的矩形框，稍留边缘。然后生成候选，可补正点。')
    def clear_box(self):
        self.box_layer.data=[];self.invalidate_preview()
    def current_box(self):
        if len(self.box_layer.data)>1:raise ValueError('一次只能框一个对象；点“框选一个目标”清空后重画')
        if not len(self.box_layer.data):return None
        points=np.asarray(self.box_layer.data[0]);return np.r_[points.min(axis=0),points.max(axis=0)].tolist()
    def zoom_target(self):
        crop,box=prompt_region(self.image.data.shape,self.positive.data,self.negative.data,self.current_box(),self.target_width.value)
        cy,cx=(box[:2]+box[2:])/2
        self.viewer.scene.camera.center=(0,float(cy),float(cx))
        size=max(64,float(np.max(box[2:]-box[:2]))*3)
        self.viewer.scene.camera.zoom=min(self.viewer.window._qt_viewer.canvas.size)/size*.85
    def prompt_signature(self):
        return json.dumps({'image_id':self.record['id'],'positive':self.positive.data.tolist(),
                           'negative':self.negative.data.tolist(),'box':self.current_box(),
                           'width':self.target_width.value,'semantics':self.sam_semantics.value,'class':self.category.value},sort_keys=True)
    def invalidate_preview(self):
        self.sam_masks=None;self.sam_details=None;self.preview_signature=None
        self.sam_batch=None;self.batch_choices=[];self.batch_included=[]
        if hasattr(self,'preview'):self.preview.data=np.zeros(self.image.data.shape[:2],np.uint8)
        self.proposal.choices=[('尚未生成',-1)]
        self.batch_target.choices=[('尚未生成',-1)];self.batch_keep.enabled=False
        self.buttons['accept_batch'].enabled=False
    def clear_prompts(self):
        if not hasattr(self,'preview'):return
        # Invalidate outstanding token: a late result cannot resurrect a cleared preview.
        if self.job is not None:self.cancel_sam()
        self.invalidate_preview();self.positive.data=np.empty((0,2));self.negative.data=np.empty((0,2))
        self.box_layer.data=[]
    def segment(self):
        if self.args.manual_only:raise ValueError('当前 --manual-only 模式未启用 SAM')
        if self.job is not None:raise ValueError('正在计算；请等待或点“取消计算”')
        box=self.current_box()
        many=self.sam_semantics.value=='point_instances'
        if many:
            if self.category.value!='SINGLE_BUBBLE':raise ValueError('每点一个泡仅适用于 SINGLE_BUBBLE')
            point_instance_requests(self.image.data.shape,self.positive.data,self.negative.data,box,self.target_width.value)
        else:
            crop,target=prompt_region(self.image.data.shape,self.positive.data,self.negative.data,box,self.target_width.value)
            if self.category.value=='SINGLE_BUBBLE' and np.max(target[2:]-target[:2])>2*self.target_width.value:
                raise ValueError('框或多个绿点覆盖范围过大；多个泡请使用“每点一个单泡”模式')
        self.invalidate_preview();token=uuid.uuid4().hex;folder=ROOT/'results/annotation_jobs'/token;folder.mkdir(parents=True)
        req={'image':str(ROOT/self.record['image']),'image_id':self.record['id'],'token':token,
             'positive_yx':self.positive.data.tolist(),'negative_yx':self.negative.data.tolist(),
             'mode':'point_instances' if many else 'local_box','object_class':self.category.value,'software_test':bool(self.args.smoke),
             'use_validated_adaptation':not self.args.smoke,
             'box_yxyx':box,'target_width_px':self.target_width.value,'context_width_px':128,
             'output':str(folder/'candidate.npz'),'prompt_signature':self.prompt_signature()}
        (folder/'request.json').write_text(json.dumps(req,indent=2));self.job=req;self.job_started=time.monotonic()
        self.refresh_prompt_history()
        self.logfile=(folder/'worker.log').open('ab')
        if self.proc is None or self.proc.state()==QProcess.ProcessState.NotRunning:
            proc=QProcess(self.viewer.window._qt_window);self.proc=proc;self.stdout_buffer=b''
            proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
            proc.readyReadStandardOutput.connect(lambda:self.read_sam(proc))
            proc.started.connect(lambda:self.send_sam_request(proc))
            proc.finished.connect(lambda code,status:self.worker_finished(proc,code))
            proc.start(sys.executable,['-B',str(ROOT/'scripts/sam_worker.py'),'--server'])
        else:self.send_sam_request(self.proc)
        self.sam_status.value='局部分割中，界面可操作；模型会复用。首次加载稍慢，可取消。'
    def send_sam_request(self,proc):
        if proc is self.proc and self.job is not None:
            proc.write((json.dumps(self.job)+'\n').encode())
    def read_sam(self,proc):
        data=bytes(proc.readAllStandardOutput())
        if proc is not self.proc:return
        if self.logfile:self.logfile.write(data);self.logfile.flush()
        self.stdout_buffer+=data
        while b'\n' in self.stdout_buffer:
            line,self.stdout_buffer=self.stdout_buffer.split(b'\n',1)
            try:obj=json.loads(line.decode(errors='replace'))
            except ValueError:continue
            if self.job is None or obj.get('token')!=self.job['token']:continue
            if 'status' in obj:self.sam_status.value=obj['status']
            if obj.get('event')=='result':self.sam_result(obj)
            elif obj.get('event')=='error':
                self.last_sam_outcome='FAILED';self.sam_status.value='SAM 失败，已有标签不受影响。请查看任务日志或手工标注。'
                self.finish_job()
    def finish_job(self):
        self.job=None
        if self.logfile:self.logfile.close();self.logfile=None
    def sam_result(self,obj):
        job=self.job;self.last_sam_wall_s=time.monotonic()-self.job_started
        try:
            with np.load(job['output'],allow_pickle=False) as z:
                if str(z['request_token'])!=job['token'] or str(z['image_id'])!=self.record['id']:raise ValueError('旧图结果已丢弃')
                if job['prompt_signature']!=self.prompt_signature():raise ValueError('提示已改变；请重新生成，旧候选未接受')
                self.sam_details=json.loads(str(z['details']))
                if job['mode']=='point_instances':
                    self.sam_batch=[z[t['array_key']].copy() for t in self.sam_details['targets']]
                else:self.sam_masks=z['masks'].copy()
            self.preview_signature=job['prompt_signature'];self.last_sam_outcome='PASS'
            if self.sam_batch is not None:
                self.batch_choices=[next((j for j,v in enumerate(t['proposals']) if v['allowed']),0) for t in self.sam_details['targets']]
                self.batch_included=[True]*len(self.sam_batch)
                self.batch_target.choices=[('全部单泡预览',-1)]+[(f'泡 {i+1}'+('（需修正）' if not any(v['allowed'] for v in t['proposals']) else ''),i) for i,t in enumerate(self.sam_details['targets'])]
                self.batch_target.value=-1;self.show_batch_target();self.buttons['accept_batch'].enabled=True
            else:
                self.proposal.choices=[(f'{i+1}: {v["area_px"]} px | '+('可检查' if v['allowed'] else '禁止接受'),i)
                                       for i,v in enumerate(self.sam_details['proposals'])]
                self.proposal.value=0;self.show_proposal()
        except Exception as e:
            self.invalidate_preview();self.last_sam_outcome='REJECTED';self.sam_status.value=str(e)
        finally:self.finish_job()
    def show_proposal(self):
        if self.updating_sam_widgets:return
        if self.sam_batch is not None:
            i=self.batch_target.value
            if i is not None and i>=0 and self.proposal.value is not None and self.proposal.value>=0:
                self.batch_choices[i]=self.proposal.value;self.show_batch_target()
            return
        if self.sam_masks is None or self.proposal.value is None or self.proposal.value<0:return
        n=self.proposal.value;self.preview.data=self.sam_masks[n].astype(np.uint8);v=self.sam_details['proposals'][n]
        message='请检查后接受；模型分数不是微泡置信度。' if v['allowed'] else '禁止接受：'+'；'.join(v['reasons'])
        self.sam_status.value=f'候选 {n+1}，{v["area_px"]} 像素；耗时 {self.last_sam_wall_s:.2f} 秒。\n{message}'
    def batch_assembly(self,included=None):
        occupied=np.logical_or.reduce([layer.data>0 for layer in self.layers.values()])
        return assemble_candidates(self.sam_batch,self.sam_details['targets'],self.batch_choices,
                                   self.batch_included if included is None else included,self.image.data.shape[:2],occupied)
    def show_batch_target(self):
        if self.updating_sam_widgets or self.sam_batch is None:return
        i=self.batch_target.value
        if i is None:return
        self.updating_sam_widgets=True
        try:
            self.batch_keep.enabled=i>=0
            if i<0:
                self.proposal.choices=[('选择一个泡后切换候选',-1)]
                self.batch_keep.value=False
                labels,items,rejected=self.batch_assembly();self.preview.data=labels
                failures=', '.join(str(j+1) for j in rejected) or '无'
                self.sam_status.value=f'{len(items)}/{len(self.sam_batch)} 个独立候选可检查；耗时 {self.last_sam_wall_s:.2f} 秒。\n需修正：{failures}。请逐泡检查；通过尺寸检查不代表正确。'
            else:
                t=self.sam_details['targets'][i];n=self.batch_choices[i]
                self.proposal.choices=[(f'{j+1}: {v["area_px"]} px | '+('可检查' if v['allowed'] else '禁止接受'),j) for j,v in enumerate(t['proposals'])]
                self.proposal.value=n;self.batch_keep.value=self.batch_included[i]
                mask=np.zeros(self.image.data.shape[:2],np.uint32);y0,x0,y1,x1=t['crop_yxyx']
                mask[y0:y1,x0:x1]=self.sam_batch[i][n].astype(np.uint32)*(i+1);self.preview.data=mask
                v=t['proposals'][n];message='请检查边界；可改候选或取消保留。' if v['allowed'] else '禁止接受：'+'；'.join(v['reasons'])
                self.sam_status.value=f'泡 {i+1}，候选 {n+1}，{v["area_px"]} 像素。\n{message}'
        finally:self.updating_sam_widgets=False
    def change_batch_keep(self):
        if self.updating_sam_widgets or self.sam_batch is None:return
        i=self.batch_target.value
        if i is not None and i>=0:self.batch_included[i]=self.batch_keep.value;self.show_batch_target()
    def accept_sam_batch(self,only=None):
        if self.sam_batch is None:raise ValueError('请先生成逐泡候选')
        if self.preview_signature!=self.prompt_signature():raise ValueError('提示已改变，请重新生成后再接受')
        included=self.batch_included if only is None else [i==only for i in range(len(self.sam_batch))]
        _,items,rejected=self.batch_assembly(included)
        if not items:raise ValueError('没有可接受的候选：'+'；'.join(f'泡 {i+1}: {reason}' for i,reason in rejected.items()))
        k=self.category.value;layer=self.layers[k];start=int(layer.data.max())+1
        indices=tuple(np.concatenate([ix[d] for ix in items.values()]) for d in (0,1))
        values=np.concatenate([np.full(len(ix[0]),start+j,np.uint32) for j,ix in enumerate(items.values())])
        # One history-aware edit: distinct IDs, atomic undo, no trimming or overwrite.
        layer.data_setitem(indices,values);ids=list(range(start,start+len(items)))
        for n in ids:self.removed_tracks.pop(f'{k}:{n}',None)
        self.last_sam_batch=(k,ids);self.last_sam_instance=None;self.instance.value=ids[-1];self.select_instance()
        for i in items:self.batch_included[i]=False
        self.mark_dirty();self.show_batch_target()
        failed=', '.join(str(i+1) for i in rejected) or '无'
        self.sam_status.value=f'已接受 {len(items)} 个独立单泡（#{start}–#{ids[-1]}），尚未人工确认。\n未接受：{failed}；可整批撤销或点选单独删除。'
    def worker_finished(self,proc,code):
        if proc is not self.proc:return
        self.read_sam(proc);self.proc=None
        if self.job is not None:
            self.last_sam_outcome='FAILED';self.sam_status.value='SAM 进程退出；标签保留，下次生成会自动重启。'
            self.finish_job()
    def cancel_sam(self):
        proc=self.proc;job=self.job;self.proc=None
        if proc is not None:
            proc.kill();proc.waitForFinished(1500)
        self.finish_job();self.stdout_buffer=b''
        if job:
            self.last_sam_outcome='CANCELLED'
            Path(job['output']).parent.joinpath('CANCELLED.json').write_text(json.dumps({'status':'USER_OR_NAVIGATION_CANCELLED','image_id':job['image_id']}))
            self.sam_status.value='已取消计算，已接受标签与草稿未改动。'
    def sam_watchdog(self):
        if self.job is not None and time.monotonic()-self.job_started>90:
            self.cancel_sam();self.sam_status.value='超过 90 秒，已取消。可继续手工标注。'
    def accept_sam(self):
        if self.sam_batch is not None:
            i=self.batch_target.value
            if i is None or i<0:raise ValueError('先在“查看哪个泡”选择单泡；整批请用“接受全部有效单泡”')
            return self.accept_sam_batch(only=i)
        if self.sam_masks is None:raise ValueError('请先生成候选')
        if self.preview_signature!=self.prompt_signature():raise ValueError('提示已改变，请重新生成候选后再接受')
        n=self.proposal.value;v=self.sam_details['proposals'][n]
        if not v['allowed']:raise ValueError('候选未通过范围检查：'+'；'.join(v['reasons'])+'。请重新框选或改用手工。')
        mask=self.sam_masks[n]
        if not np.any(mask):raise ValueError('没有可接受的轮廓')
        if np.any(np.logical_or.reduce([v.data>0 for v in self.layers.values()])&mask):
            raise ValueError('候选与已有实例重叠；请点选原对象修改，或调整框和负点，不会覆盖标签')
        k=self.category.value;layer=self.layers[k];new_id=int(layer.data.max())+1
        layer.data_setitem(np.where(mask),new_id)
        self.removed_tracks.pop(f'{k}:{new_id}',None)
        self.last_sam_instance=(k,new_id);self.last_sam_batch=None;self.clear_prompts();self.instance.value=new_id;self.select_instance()
        self.mark_dirty();self.sam_status.value=f'已接受 {NAMES[k]} #{new_id}，尚未保存。错误可直接“删除上次接受的 SAM 对象”或 Ctrl+Z。'

def smoke(w):
    """Exercise actual Qt buttons and napari polygon API in a separate sandbox."""
    from qtpy.QtTest import QTest
    from qtpy.QtCore import Qt
    report={'annotation_directory':str(w.directory),'research_labels_modified':False,'checks':{},'heartbeats':0};ticks=[]
    timer=QTimer();timer.timeout.connect(lambda:ticks.append(time.monotonic()));timer.start(100)
    def click(key):QTest.mouseClick(w.buttons[key].native,Qt.MouseButton.LeftButton);QApplication.processEvents()
    def run():
        try:
            viewer_id=id(w.viewer);idx=w.index
            # Clear sandbox only, then draw via napari's real Labels polygon tool.
            for layer in w.layers.values():layer.data=np.zeros_like(layer.data)
            w.category.value='SINGLE_BUBBLE';click('new')
            layer=w.layers['SINGLE_BUBBLE'];layer._overlays['polygon'].points=[[480,680],[480,695],[500,695],[500,680]]
            click('finish');assert layer.data.max()==1
            click('new');layer._overlays['polygon'].points=[[520,710],[520,725],[540,725],[540,710]];click('finish');assert layer.data.max()==2
            w.notes.value='SOFTWARE_UI_FIXTURE; not expert labels';click('save')
            before=layer.data.copy();click('next');click('prev')
            assert id(w.viewer)==viewer_id and w.index==idx;np.testing.assert_array_equal(before,w.layers['SINGLE_BUBBLE'].data)
            report['checks'].update(polygon_draw_two_instances=True,save_reopen_equal=True,persistent_window_next_previous=True,
                                   frame_switch_s=w.last_load_seconds,no_gpu_required_for_manual=True)
            w.go(len(w.records)-1);assert w.image.data.shape[:2]==(1024,1280)
            w.go(idx);np.testing.assert_array_equal(before,w.layers['SINGLE_BUBBLE'].data)
            report['checks']['switch_different_resolution_and_restore']=True
            layer=w.layers['SINGLE_BUBBLE'];layer._overlays['polygon'].points=[[100,100],[100,110]]
            click('next');assert w.index==idx
            report['checks']['unfinished_polygon_blocks_navigation']=True
            w.clear_polygon()
            w.notes.value='SOFTWARE_UI_FIXTURE_AUTOSAVE';w.autosave();assert w.autosave_count>=1
            report['checks']['autosave']=True
            w.layers['SINGLE_BUBBLE'].data=np.zeros_like(before);w.save(False)
            w.tabs.setCurrentIndex(1);w.sam_semantics.value='one_object';w.positive.data=np.array([[500.,700.]])
            if w.args.smoke_sam:
                w.segment();QTimer.singleShot(500,check_sam)
            else:finish()
        except Exception:
            report['error']=traceback.format_exc();finish()
    def check_sam():
        if w.job is not None:QTimer.singleShot(250,check_sam);return
        try:
            assert w.last_sam_outcome=='PASS',w.last_sam_outcome
            assert np.any(w.preview.data);report['checks']['sam_async_inference']=True
            report['checks']['sam_wall_s']=w.last_sam_wall_s
            report['checks']['event_loop_max_gap_s']=max(np.diff(ticks))
            assert max(np.diff(ticks))<2,'UI heartbeat blocked'
            w.accept_sam();assert not len(w.positive.data) and not np.any(w.preview.data)
            report['checks']['sam_accept_clears_prompts']=True
            w.positive.data=np.array([[500.,700.]]);w.segment();w.cancel_sam()
            report['checks']['sam_cancel']=w.last_sam_outcome=='CANCELLED'
        except Exception:report['error']=traceback.format_exc()
        finish()
    def finish():
        timer.stop();report['heartbeats']=len(ticks);report['visible']=w.viewer.window._qt_window.isVisible()
        for layer in w.layers.values():layer.data=np.zeros_like(layer.data)
        w.clear_prompts();w.notes.value='软件界面检查；不计入人工标签';w.save(False)
        w.tabs.setCurrentIndex(0);w.viewer.reset_view()
        w.viewer.screenshot(path=str(ROOT/'evidence/annotation_repaired_window.png'),canvas_only=False)
        # Simulate closing with unsaved pixels, through the actual Qt close event.
        try:
            w.layers['SINGLE_BUBBLE'].data[10:13,10:13]=1
            w.viewer.window._qt_window.close()
            meta,masks=load_annotation(w.directory,w.record)
            assert masks['SINGLE_BUBBLE'][11,11]==1 and not meta['expert_confirmed']
            report['checks']['close_autosaves_unsaved_pixels']=True
        except Exception:report['error']=traceback.format_exc()
        report['PASS']='error' not in report
        (ROOT/'evidence/annotation_repair_gui_smoke.json').write_text(json.dumps(report,indent=2))
        kind='sam' if w.args.smoke_sam else 'manual'
        (ROOT/f'evidence/annotation_repair_{kind}_smoke.json').write_text(json.dumps(report,indent=2))
    QTimer.singleShot(2000,run)

def main():
    p=argparse.ArgumentParser();p.add_argument('--index',type=int);p.add_argument('--manual-only',action='store_true')
    p.add_argument('--no-auto-train',action='store_true',help='Disable automatic training for this workbench session')
    p.add_argument('--annotation-dir',help=argparse.SUPPRESS);p.add_argument('--smoke',action='store_true');p.add_argument('--smoke-sam',action='store_true');a=p.parse_args()
    if a.smoke_sam:a.smoke=True
    if a.annotation_dir and not a.smoke:raise ValueError('--annotation-dir is only available for isolated software tests')
    lock=QLockFile(str(ROOT/'data'/('annotation_smoke.lock' if a.smoke else 'annotation_workbench.lock')))
    lock.setStaleLockTime(0)
    if not lock.tryLock(100):raise RuntimeError('标注工作台已在运行。请回到已有窗口；不要同时打开两个写入进程。')
    app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
    w=Workbench(a)
    if a.smoke:smoke(w)
    napari.run();lock.unlock()
if __name__=='__main__':main()
