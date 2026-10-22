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
import napari
from napari.layers.labels._labels_key_bindings import complete_polygon,reset_polygon
from magicgui.widgets import Container,Label,ComboBox,LineEdit,CheckBox,PushButton,SpinBox
from qtpy.QtCore import QTimer,QProcess,QObject,QEvent,QLockFile
from qtpy.QtWidgets import QTabWidget,QApplication
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
            self.w.cancel_sam();self.w.timer.stop();self.w.save_session()
        return False

class Workbench:
    def __init__(self,args):
        self.args=args;self.records=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows']
        self.directory=ROOT/('evidence/annotation_repair_gui_sandbox' if args.smoke else 'data/annotations')
        self.directory.mkdir(parents=True,exist_ok=True)
        self.sessionfile=self.directory/'session.json';previous=json.loads(self.sessionfile.read_text()) if self.sessionfile.exists() else {}
        self.index=args.index if args.index is not None else previous.get('index',0)
        if not 0<=self.index<len(self.records):raise ValueError('Frame index out of range')
        self.viewer=napari.Viewer(title='微泡人工标注 — 正在打开图片')
        self.viewer.window._qt_window.resize(1680,1000)
        self.layers={};self.loading=True;self.dirty=False;self.meta=None;self.saved_masks=None
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
        self.buttons={}
        def button(key,text,fn):
            b=PushButton(text=text);b.clicked.connect(lambda *_:self.safe(fn));self.buttons[key]=b;return b
        main=Container(widgets=[self.status,self.instructions,self.category,
            row(button('new','② 新建轮廓',self.new_object),button('finish','③ 完成轮廓',self.finish_polygon)),
            row(button('brush','画笔修边',lambda:self.edit_mode('paint')),button('erase','橡皮擦',lambda:self.edit_mode('erase'))),
            row(button('undo','撤销本次笔画',self.undo),button('pan','移动 / 缩放',self.pan)),
            row(button('fit','查看整幅图',self.viewer.reset_view),button('clear_polygon','取消未完成轮廓',self.clear_polygon)),
            self.person,self.notes,button('save','④ 保存草稿（Ctrl+S）',lambda:self.save(False)),
            self.reviewed,button('confirm','⑤ 确认本帧人工审核完成',lambda:self.save(True)),
            row(button('prev','上一张',lambda:self.go(self.index-1)),button('next','⑥ 保存并下一张',lambda:self.go(self.index+1))),
            row(self.jump,button('jump','跳转',lambda:self.go(self.jump.value-1))),
            label('草稿每 30 秒保存；换帧 / 关闭时保存修改。\n左右键不换帧；Alt+← / → 换帧。\n空格+拖动平移，滚轮缩放。')],labels=True)
        sam=Container(widgets=[self.sam_status,label('只在需要时使用原版 SAM。先放正点，再分割；结果必须人工检查。'),
            button('positive','1. 放正点（对象内部）',lambda:self.prompt_mode(True)),
            button('negative','2. 放负点（排除区域，可选）',lambda:self.prompt_mode(False)),
            button('segment','3. 生成 / 更新候选',self.segment),button('cancel_sam','取消计算',self.cancel_sam),
            button('accept','4. 接受候选到当前类别',self.accept_sam),
            button('clear_prompts','清空候选和所有点',self.clear_prompts),
            label('接受后自动清空正负点，避免影响下一对象。\n候选过大或选错：清空后改用手工轮廓。\n换帧会取消当前 SAM 任务；旧帧结果不会混入新帧。')])
        extra=Container(widgets=[label('可选：仅给连续 8 帧中能确实对应的对象设置短轨迹 ID。普通标注无需此步骤。'),
            self.instance,button('select','选中此编号进行修改',self.select_instance),
            self.track,button('track','保存短轨迹 ID 到本帧',self.set_track),
            button('delete','删除选中实例（可撤销）',self.delete_instance),
            label('同一类别中的编号对应一个对象。修改旧对象先选编号。\n颜色仅用于辅助辨认，不表示已完成专家审核。')])
        tabs=QTabWidget();tabs.addTab(main.native,'逐张标注');tabs.addTab(sam.native,'SAM 辅助');tabs.addTab(extra.native,'实例 / 短轨迹')
        tabs.setMinimumWidth(385);tabs.setMaximumWidth(480)
        self.tabs=tabs;self.main=main;self.sam=sam;self.extra=extra
        self.viewer.window.add_dock_widget(tabs,name='标注步骤',area='right')
        self.category.changed.connect(lambda *_:self.select_class())
        self.notes.changed.connect(lambda *_:self.mark_dirty())
        self.viewer.bind_key('Control-S',lambda _:self.save(False),overwrite=True)
        self.viewer.bind_key('Alt-Right',lambda _:self.go(self.index+1),overwrite=True)
        self.viewer.bind_key('Alt-Left',lambda _:self.go(self.index-1),overwrite=True)
        self.guard=CloseGuard(self);self.viewer.window._qt_window.installEventFilter(self.guard)
        self.timer=QTimer();self.timer.timeout.connect(self.autosave);self.timer.start(30000)
        self.watch=QTimer();self.watch.timeout.connect(self.sam_watchdog);self.watch.start(1000)
        self.load_frame(self.index)
    def safe(self,fn):
        try:return fn()
        except Exception as e:
            self.status.value='操作未完成：'+str(e)
            print(traceback.format_exc(),flush=True);return False
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
            self.positive=self.viewer.add_points(np.empty((0,2)),name='SAM 正点',size=8,face_color='#36ff86')
            self.negative=self.viewer.add_points(np.empty((0,2)),name='SAM 负点',size=8,face_color='#ff5e7e')
        else:self.image.data=image
        for k,layer in self.layers.items():
            layer.data=masks[k].copy() if masks else np.zeros(image.shape[:2],np.uint32)
            layer._reset_history();layer._overlays['polygon'].points=[]
        self.clear_prompts();self.meta=meta;self.saved_masks={k:v.copy() for k,v in masks.items()} if masks else None
        self.track_ids=dict(meta.get('track_ids',{})) if meta else {};self.notes.value=meta.get('notes','') if meta else ''
        self.reviewed.value=False;self.jump.value=index+1;self.dirty=False;self.loading=False
        state='已确认；编辑后需重新确认' if meta and meta['expert_confirmed'] else '未确认 / 草稿'
        self.status.value=f'{index+1}/{len(self.records)}  {self.record["id"]}\n分组：{self.record["split"]}；{state}'
        self.viewer.title=f'微泡人工标注 | {index+1}/{len(self.records)} | {self.record["id"]}'
        self.select_class();self.viewer.reset_view();self.save_session()
        self.last_load_seconds=time.monotonic()-started
    def go(self,index):
        if index<0 or index>=len(self.records):self.status.value='已经到达第一张 / 最后一张。';return
        if self.save_if_changed():self.load_frame(index)
    def select_class(self):
        if not self.layers:return
        layer=self.layers[self.category.value];self.viewer.layers.selection.active=layer
        layer.selected_label=int(layer.data.max())+1;self.instance.value=layer.selected_label;layer.mode='pan_zoom'
    def new_object(self):
        self.select_class();layer=self.layers[self.category.value];layer.mode='polygon'
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
    def undo(self):self.layers[self.category.value].undo();self.mark_dirty()
    def select_instance(self):
        layer=self.layers[self.category.value]
        if not np.any(layer.data==self.instance.value):raise ValueError('此类别中不存在这个编号')
        layer.selected_label=self.instance.value;self.viewer.layers.selection.active=layer;layer.mode='paint'
    def delete_instance(self):
        self.select_instance();layer=self.layers[self.category.value]
        # Use napari's history-aware public fill method for connected regions.
        ys,xs=np.where(layer.data==layer.selected_label)
        old=layer.contiguous;layer.contiguous=False;layer.fill((int(ys[0]),int(xs[0])),0);layer.contiguous=old
        self.track_ids.pop(f'{self.category.value}:{self.instance.value}',None);self.mark_dirty()
    def set_track(self):
        self.select_instance();key=f'{self.category.value}:{self.instance.value}'
        if self.track.value.strip():self.track_ids[key]=self.track.value.strip()
        else:self.track_ids.pop(key,None)
        self.mark_dirty();self.status.value='短轨迹 ID 已写入草稿；仍需保存 / 确认。'
    def prompt_mode(self,positive):
        layer=self.positive if positive else self.negative;self.viewer.layers.selection.active=layer;layer.mode='add'
        self.sam_status.value='在图上单击放点。绿色为对象内部，红色为排除区域。然后点“生成 / 更新候选”。'
    def clear_prompts(self):
        if not hasattr(self,'preview'):return
        self.preview.data=np.zeros(self.image.data.shape[:2],np.uint8)
        self.positive.data=np.empty((0,2));self.negative.data=np.empty((0,2))
    def segment(self):
        if self.args.manual_only:raise ValueError('当前 --manual-only 模式未启用 SAM；手工描画仍可用')
        if self.proc is not None:raise ValueError('正在计算；请等待或点“取消计算”')
        if not len(self.positive.data):raise ValueError('请先在对象内部放至少一个正点')
        token=uuid.uuid4().hex;folder=ROOT/'results/annotation_jobs'/token;folder.mkdir(parents=True)
        req={'image':str(ROOT/self.record['image']),'image_id':self.record['id'],'token':token,
             'positive_yx':self.positive.data.tolist(),'negative_yx':self.negative.data.tolist(),
             'embedding':str(ROOT/'models/embeddings'/f'{self.record["id"]}.zarr'),'output':str(folder/'candidate.npz')}
        (folder/'request.json').write_text(json.dumps(req,indent=2));self.job=req;self.job_started=time.monotonic()
        self.logfile=(folder/'worker.log').open('ab');proc=QProcess(self.viewer.window._qt_window)
        self.proc=proc;proc.setProcessChannelMode(QProcess.ProcessChannelMode.MergedChannels)
        proc.readyReadStandardOutput.connect(lambda:self.read_sam(proc))
        proc.finished.connect(lambda code,status:self.sam_finished(proc,code))
        proc.start(sys.executable,['-B',str(ROOT/'scripts/sam_worker.py'),str(folder/'request.json')])
        self.sam_status.value='后台计算中，界面可以缩放、保存；可取消。首次加载会较慢。'
    def read_sam(self,proc):
        data=bytes(proc.readAllStandardOutput())
        if self.logfile:self.logfile.write(data);self.logfile.flush()
        if proc is not self.proc:return
        for line in data.decode(errors='replace').splitlines():
            try:
                obj=json.loads(line)
                if 'status' in obj:self.sam_status.value=obj['status']
            except ValueError:pass
    def sam_finished(self,proc,code):
        if proc is not self.proc:return
        self.read_sam(proc);job=self.job;self.last_sam_wall_s=time.monotonic()-self.job_started
        self.proc=None;self.job=None
        if self.logfile:self.logfile.close();self.logfile=None
        if code!=0:
            self.last_sam_outcome='FAILED';self.sam_status.value='SAM 未完成；草稿不受影响。请使用手工轮廓，详见 results/annotation_jobs。';return
        try:
            with np.load(job['output'],allow_pickle=False) as z:
                if str(z['request_token'])!=job['token'] or self.record['id']!=job['image_id']:raise ValueError('旧帧结果已丢弃')
                self.preview.data=z['mask'].astype(np.uint8)
            self.last_sam_outcome='PASS';self.sam_status.value=f'候选已显示（{self.last_sam_wall_s:.1f} 秒）。检查边界后再接受；它还不是人工标签。'
        except Exception as e:self.last_sam_outcome='FAILED';self.sam_status.value=str(e)
    def cancel_sam(self):
        proc=self.proc
        if proc is None:return
        job=self.job;self.proc=None;self.job=None;self.last_sam_outcome='CANCELLED'
        proc.kill();proc.waitForFinished(1500)
        if self.logfile:self.logfile.close();self.logfile=None
        if job:Path(job['output']).parent.joinpath('CANCELLED.json').write_text(json.dumps({'status':'USER_OR_NAVIGATION_CANCELLED','image_id':job['image_id']}))
        self.sam_status.value='已取消当前计算，标签与草稿未改动。'
    def sam_watchdog(self):
        if self.proc is not None and time.monotonic()-self.job_started>90:
            self.cancel_sam();self.sam_status.value='计算超过 90 秒，已自动取消。请先手工标注；故障日志已保留。'
    def accept_sam(self):
        mask=self.preview.data>0
        if not np.any(mask):raise ValueError('还没有可接受的候选，请先生成候选')
        conflict=np.logical_or.reduce([v.data>0 for k,v in self.layers.items() if k!=self.category.value])&mask
        if np.any(conflict):raise ValueError('候选与其他类别重叠，请用负点收缩候选或先修正原标签；未覆盖已有标签')
        layer=self.layers[self.category.value];occupied=(layer.data>0)&mask
        if np.any(occupied):raise ValueError('候选与已有实例重叠；请修改旧实例或清空候选，不会覆盖已有标签')
        arr=layer.data.copy();arr[mask]=int(arr.max())+1;layer.data=arr
        self.mark_dirty();self.clear_prompts();self.select_class();self.sam_status.value='已加入草稿并清空正负点。请保存；专家确认仍需逐帧检查。'

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
            w.tabs.setCurrentIndex(1);w.positive.data=np.array([[500.,700.]])
            if w.args.smoke_sam:
                w.segment();QTimer.singleShot(500,check_sam)
            else:finish()
        except Exception:
            report['error']=traceback.format_exc();finish()
    def check_sam():
        if w.proc is not None:QTimer.singleShot(250,check_sam);return
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
    p.add_argument('--smoke',action='store_true');p.add_argument('--smoke-sam',action='store_true');a=p.parse_args()
    if a.smoke_sam:a.smoke=True
    lock=QLockFile(str(ROOT/'data'/('annotation_smoke.lock' if a.smoke else 'annotation_workbench.lock')))
    lock.setStaleLockTime(0)
    if not lock.tryLock(100):raise RuntimeError('标注工作台已在运行。请回到已有窗口；不要同时打开两个写入进程。')
    app=QApplication.instance() or QApplication([]);app.setFont(QFont('Microsoft YaHei',10))
    w=Workbench(a)
    if a.smoke:smoke(w)
    napari.run();lock.unlock()
if __name__=='__main__':main()
