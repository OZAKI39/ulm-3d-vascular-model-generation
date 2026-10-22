"""One Qt UI shared by replay and LIVE. Screen recordings capture these widgets."""
import time,queue,threading,json
from pathlib import Path
import numpy as np,cv2,av
from PySide6 import QtWidgets as W,QtCore as C,QtGui as G
from ycommon import ROOT,dump,sha,stats
from collections import deque
from pipeline import Pipeline
from live_geometry import new_geometry,validate_geometry
from desktop import open_directory
from pzt_panel import PZTPanel
from piv_display import piv_arrows,arrow_shape,ARROW_RGB,OUTLINE_RGB,REFERENCE_SPEED,GAIN_SECONDS,MAX_LENGTH_PX,STROKE_PX,HALO_PX,display_vectors,display_is_smoothed
from overlay_style import DEFAULT_LAYERS,overlay_opacity,overlay_caption

def observation_text(snapshot):
 state=snapshot.get('state') or {};valid=snapshot.get('valid',False)
 status=snapshot.get('status','WAITING')
 labels={'FRESH':'Active','STOPPED':'Stopped','PAUSED':'Paused','EOF':'Playback finished',
  'WAITING':'Waiting for analysis','WARMUP':'Preparing analysis','STALE':'Waiting for fresh analysis',
  'MASK_NOT_VERIFIED':'Draw or load a channel mask','LOW_PIV_SUPPORT_OR_REGISTRATION':'Motion unavailable',
  'DISCONNECTED':'Camera disconnected','CONNECTING':'Connecting','RECONNECTING':'Reconnecting','ANALYSIS_ERROR':'Analysis error; see records'}
 status=labels.get(status,'Connection error' if status.startswith('ERROR:') else status.replace('_',' ').capitalize())
 aggregation=state.get('aggregation') or {};regions=aggregation.get('regions') or [{}];region=regions[0]
 count=region.get('count',len(state.get('candidates',[]))) if snapshot.get('fresh') else '-'
 velocity=region.get('mean_velocity') if valid else None
 motion=f'x={velocity[0]:.1f}, y={velocity[1]:.1f} px/s' if velocity is not None else 'Unavailable'
 trend=aggregation.get('status') if valid else None
 trend={'VISIBLE_CONVERGENCE':'Visible convergence','VISIBLE_DISPERSION':'Visible dispersion',
  'NO_CONSISTENT_CONVERGENCE':'No consistent convergence'}.get(trend,'Not observable')
 return f'Status: {status}\nCandidates: {count}\nVelocity: {motion}\nTrend: {trend}'

class Canvas(W.QWidget):
 mask_finished=C.Signal()
 mask_cancelled=C.Signal()
 def __init__(self):super().__init__();self.setMinimumSize(700,550);self.snap=None;self.layers={};self.editing=False;self.points=[];self.image_rect=None;self.paint_times=deque(maxlen=6000);self.diagnostics=False;self.active_geometry=None;self.setFocusPolicy(C.Qt.StrongFocus);self._draw_state=None;self._draw_paths=None;self.setAttribute(C.Qt.WA_OpaquePaintEvent,True)
 def image_point(self,e):
  if not self.image_rect:return None
  x,y,w,h,iw,ih=self.image_rect;p=e.position()
  if x<=p.x()<x+w and y<=p.y()<y+h:return [(p.x()-x)*iw/w,(p.y()-y)*ih/h]
 def append_point(self,point):
  if point is not None and (not self.points or np.linalg.norm(np.array(point)-self.points[-1])>1e-6):self.points.append(point)
  self.update()
 def mousePressEvent(self,e):
  if not self.editing:return super().mousePressEvent(e)
  self.setFocus();point=self.image_point(e)
  if e.button()==C.Qt.LeftButton and point is not None:
   scale=self.image_rect[2]/self.image_rect[4]
   if len(self.points)>=3 and np.linalg.norm(np.array(point)-self.points[0])*scale<=8:self.mask_finished.emit()
   else:self.append_point(point)
 def mouseDoubleClickEvent(self,e):
  if self.editing and e.button()==C.Qt.LeftButton and self.image_point(e) is not None:
   self.append_point(self.image_point(e));self.mask_finished.emit()
 def keyPressEvent(self,e):
  if not self.editing:return super().keyPressEvent(e)
  if e.key() in [C.Qt.Key_Return,C.Qt.Key_Enter]:self.mask_finished.emit()
  elif e.key() in [C.Qt.Key_Backspace,C.Qt.Key_Delete]:
   if self.points:self.points.pop();self.update()
  elif e.key()==C.Qt.Key_Escape:self.editing=False;self.points=[];self.unsetCursor();self.update();self.mask_cancelled.emit()
 def paintEvent(self,e):
  if not self.diagnostics:self.paint_canvas(e);return
  start=time.perf_counter()
  try:self.paint_canvas(e)
  finally:self.paint_times.append((time.perf_counter()-start)*1000)
 def paint_canvas(self,e):
  painter=G.QPainter(self);painter.fillRect(self.rect(),G.QColor('#000000'));s=self.snap
  if not s or not s['raw']:painter.setPen(G.QColor('white'));painter.drawText(self.rect(),C.Qt.AlignCenter,'Press Start to view the camera');return
  p=s.get('display_packet') or s['raw'];a=np.ascontiguousarray(p.image);h,w=a.shape[:2]
  fmt=G.QImage.Format_Grayscale8 if a.ndim==2 else G.QImage.Format_BGR888;im=G.QImage(a.data,w,h,a.strides[0],fmt)
  scale=min(self.width()/w,self.height()/h);ox=(self.width()-w*scale)/2;oy=(self.height()-h*scale)/2;self.image_rect=(ox,oy,w*scale,h*scale,w,h)
  painter.translate(ox,oy);painter.scale(scale,scale);painter.drawImage(0,0,im)
  state=s['state'];valid=s['valid'];arrows=[]
  g=self.active_geometry if p.source_mode=='LIVE_CAMERA' else state.get('geometry') if state else None
  if g and tuple(g['image_wh'])==p.original_resolution and self.layers.get('Mask',True):
   painter.setPen(G.QPen(G.QColor('#55dfd0'),2/scale));painter.drawPolygon(G.QPolygonF([C.QPointF(*v) for v in g['channel_polygon']]))
   painter.setBrush(G.QColor(30,30,40,80))
   for poly in g['excluded_polygons']:painter.drawPolygon(G.QPolygonF([C.QPointF(*v) for v in poly]))
   painter.setBrush(C.Qt.NoBrush)
  if state and s.get('fresh') and state.get('geometry') :
   marks,arrows=self.marks_for(state,w,h,scale,valid)
   painter.setOpacity(overlay_opacity(s.get('age_ms'),s.get('ttl_ms',250.)))
   painter.drawImage(C.QRectF(0,0,w,h),marks)
  painter.setOpacity(1.)
  if self.points:
   painter.setPen(G.QPen(G.QColor('#ffcc66'),3/scale));painter.drawPolyline(G.QPolygonF([C.QPointF(*v) for v in self.points]))
   for point in self.points:painter.drawEllipse(C.QPointF(*point),4/scale,4/scale)
  painter.resetTransform()
  if arrows:
   y=self.height()-54;reference=arrow_shape([0,0],[REFERENCE_SPEED,0]);shift=np.array([22,y])
   painter.fillRect(10,y-18,int(220+REFERENCE_SPEED*GAIN_SECONDS*scale),36,G.QColor(20,25,30,220))
   self.paint_arrows(painter,[tuple(v*scale+shift for v in reference)],1.)
   smoothed=display_is_smoothed(state,self.layers)
   painter.setPen(G.QColor('#f0f0f0'));painter.drawText(int(34+REFERENCE_SPEED*GAIN_SECONDS*scale),y+5,f'{REFERENCE_SPEED:g} px/s | '+('smoothed' if smoothed else 'raw'))
   velocities=np.asarray(display_vectors(state,smoothed)).reshape(-1,4)
   if np.any(np.linalg.norm(velocities[:,2:],axis=1)*GAIN_SECONDS>MAX_LENGTH_PX):painter.drawText(280,y+5,'Long arrows capped')
  # Keep this band in place instead of flashing it on/off at each exact match.
  if state and state.get('geometry'):
   delay=max(0.,(p.arrival_monotonic-state['arrival_monotonic'])*1000)
   painter.fillRect(10,self.height()-30,self.width()-20,24,G.QColor(20,25,30,210));painter.setPen(G.QColor('#e0c99a'));painter.drawText(20,self.height()-13,overlay_caption(delay) if s.get('fresh') else 'Waiting for fresh analysis')
  if p.source_mode=='LIVE_CAMERA' and (self.editing or self.active_geometry is None):
   painter.fillRect(10,10,self.width()-20,30,G.QColor(20,25,30,210));painter.setPen(G.QColor('#e0c99a'));painter.drawText(20,31,'Drawing mask — double-click / Enter to apply' if self.editing else 'No active mask — Draw LIVE mask to start analysis')
  elif not valid:
   painter.fillRect(10,10,self.width()-20,36,G.QColor(50,20,20,210));painter.setPen(G.QColor('#ffc27b'));painter.drawText(20,34,s['status']+' — motion / trend invalid')
  elif state and state.get('motion_reference')=='IMAGE_SPACE_UNREGISTERED':
   painter.fillRect(10,10,self.width()-20,30,G.QColor(20,25,30,210));painter.setPen(G.QColor('#e0c99a'));painter.drawText(20,31,'Image-space motion | camera drift not corrected')
 def marks_for(self,state,w,h,scale,valid):
  key=(w,h,scale,valid,tuple(self.layers.get(k,v) for k,v in DEFAULT_LAYERS.items()))
  if state is getattr(self,'_marks_state',None) and key==self._marks_key:return self._marks
  image=G.QImage(max(1,round(w*scale)),max(1,round(h*scale)),G.QImage.Format_ARGB32_Premultiplied);image.fill(C.Qt.transparent)
  painter=G.QPainter(image);painter.scale(scale,scale)
  candidates,arrows,arrow_paths=self.paths_for(state);x,y=state['geometry']['crop_xywh'][:2]
  try:
   if self.layers.get('Candidates',True):
    painter.setPen(G.QPen(G.QColor('#73e9fa'),1.2/scale));painter.drawPath(candidates)
   if self.layers.get('Validity',False):
    painter.setPen(G.QPen(G.QColor(245,195,90,90),1/scale))
    for (px,py),ok in zip(state['window_grid_xy'],state['window_valid']):
     if not ok:painter.drawPoint(C.QPointF(px+x,py+y))
   if (self.layers.get('Regions',False) or self.layers.get('Shape',False)) and state.get('aggregation'):
    painter.setPen(G.QPen(G.QColor('#e3a6d8'),1.3/scale))
    for cluster in state['aggregation']['clusters']:
     if self.layers.get('Regions',False):
      a,b,d,f=cluster['extent'];painter.drawRect(C.QRectF(a+x,b+y,d-a,f-b))
     if self.layers.get('Shape',False):
      px,py=cluster['center'];r=cluster['rg_px'];painter.drawEllipse(C.QPointF(px+x,py+y),r,r)
   if valid and self.layers.get('PIV',True):self.paint_arrows(painter,arrows,scale,arrow_paths)
   else:arrows=[]
  finally:painter.end()
  self._marks_state=state;self._marks_key=key;self._marks=(image,arrows)
  return self._marks
 def paths_for(self,state):
  # The published observation is immutable. Cache only its drawing paths, not
  # an old result to keep displaying after invalidation, a mask edit or expiry.
  smoothed=display_is_smoothed(state,self.layers)
  if state is not self._draw_state or smoothed!=getattr(self,'_draw_mode',None):
   x,y=state['geometry']['crop_xywh'][:2];candidates=G.QPainterPath()
   for px,py in state['candidates']:candidates.addEllipse(C.QPointF(px+x,py+y),4,4)
   arrows=piv_arrows(display_vectors(state,smoothed),(x,y),sample=not smoothed)
   self._draw_paths=(candidates,arrows,self.arrow_paths(arrows));self._draw_state=state;self._draw_mode=smoothed
  return self._draw_paths
 def arrow_paths(self,arrows):
  shafts=G.QPainterPath();heads=G.QPainterPath()
  for start,end,head in arrows:
   shafts.moveTo(C.QPointF(*start));shafts.lineTo(C.QPointF(*end))
   heads.moveTo(C.QPointF(*head[1]));heads.lineTo(C.QPointF(*head[0]));heads.lineTo(C.QPointF(*head[2]))
  return shafts,heads
 def paint_arrows(self,painter,arrows,scale,paths=None):
  painter.save();painter.setRenderHint(G.QPainter.Antialiasing,True)
  shafts,heads=paths if paths is not None else self.arrow_paths(arrows)
  painter.setBrush(C.Qt.NoBrush)
  painter.setPen(G.QPen(G.QColor(*OUTLINE_RGB),HALO_PX/scale,C.Qt.SolidLine,C.Qt.RoundCap,C.Qt.RoundJoin));painter.drawPath(shafts);painter.drawPath(heads)
  painter.setPen(G.QPen(G.QColor(*ARROW_RGB),STROKE_PX/scale,C.Qt.SolidLine,C.Qt.RoundCap,C.Qt.RoundJoin));painter.drawPath(shafts);painter.drawPath(heads)
  painter.restore()
class Curves(W.QWidget):
 def __init__(self):super().__init__();self.setMinimumHeight(140);self.rows=[];self.events=[]
 def paintEvent(self,e):
  p=G.QPainter(self);p.fillRect(self.rect(),G.QColor('#000000'));labels=['Candidate count','Group radius (px)','Radial motion (px/s)'];cols=['#73e9fa','#ffd180','#ef9a9a'];indices=[3,6,7]
  for j,name in enumerate(labels):
   x=j*self.width()/3;y=0;w=self.width()/3;h=self.height();p.setPen(G.QColor(cols[j]));p.drawText(int(x+8),int(y+18),name)
   vals=[r[indices[j]] for r in self.rows[-150:]]
   if len(vals)>2:
    finite=[v for v in vals if v is not None and np.isfinite(v)];lo=min(finite) if finite else 0;hi=max(finite) if finite else 1
    path=G.QPainterPath();started=False
    for k,v in enumerate(vals):
     if v is None:started=False;continue
     pt=C.QPointF(x+8+k*(w-16)/(len(vals)-1),y+h-10-(v-lo)/max(hi-lo,1e-6)*(h-38))
     if started:path.lineTo(pt)
     else:path.moveTo(pt);started=True
    p.drawPath(path);p.drawText(int(x+w-80),int(y+18),f'{vals[-1]:.2f}' if vals[-1] is not None else 'null')
    visible=self.rows[-150:];timeline=[r[2] for r in visible]
    if timeline[0] is not None and timeline[-1] is not None and timeline[-1]>timeline[0]:
     for event in self.events:
      et=event.get('source_pts')
      if event.get('session_id')==visible[-1][0] and et is not None and timeline[0]<=et<=timeline[-1]:
       ex=x+8+(et-timeline[0])/(timeline[-1]-timeline[0])*(w-16);p.setPen(G.QColor('#ffc27b'));p.drawLine(C.QPointF(ex,y+25),C.QPointF(ex,y+h-8));p.drawText(int(ex+2),int(y+35),event['state'])
  if self.events:p.setPen(G.QColor('#ffc27b'));p.drawText(10,self.height()-2,f'Event markers: {len(self.events)} (record only; host time)')
class ScreenRecorder:
 def __init__(self,path):self.path=path;self.q=queue.Queue(8);self.done=threading.Event();self.dropped=0;self.thread=threading.Thread(target=self.run,daemon=True);self.thread.start()
 def submit(self,image,t):
  a=image.convertToFormat(G.QImage.Format_RGB888);arr=np.frombuffer(a.bits(),np.uint8,count=a.sizeInBytes()).reshape(a.height(),a.bytesPerLine())[:,:a.width()*3].reshape(a.height(),a.width(),3).copy()
  try:self.q.put_nowait((arr,t))
  except queue.Full:self.dropped+=1
 def run(self):
  from fractions import Fraction
  out=av.open(str(self.path),'w');stream=None
  try:
   while not self.done.is_set() or not self.q.empty():
    try:a,t=self.q.get(timeout=.05)
    except queue.Empty:continue
    if stream is None:
     stream=out.add_stream('libx264',rate=20);stream.width=a.shape[1]//2*2;stream.height=a.shape[0]//2*2;stream.pix_fmt='yuv420p';stream.time_base=Fraction(1,1000000);stream.codec_context.time_base=Fraction(1,1000000);stream.options={'crf':'22','preset':'ultrafast','threads':'2'}
    f=av.VideoFrame.from_ndarray(a[:stream.height,:stream.width],format='rgb24');f.pts=round(t*1e6);f.time_base=Fraction(1,1000000)
    for packet in stream.encode(f):out.mux(packet)
   if stream:
    for packet in stream.encode():out.mux(packet)
  finally:out.close()
 def close(self):self.done.set();self.thread.join(15)
class Window(W.QMainWindow):
 def __init__(self,cfg,output,auto=False,seconds=None,capture=None):
  super().__init__();self.cfg=cfg;self.output=Path(output);self.pipe=None;self.record=None;self.started=time.monotonic();self.elapsed_limit=seconds;self.capture=capture;self.metrics_last=None;self.saved_preview=False
  self.setWindowTitle('F10Y · Chronos collective observation');self.resize(1480,1000)
  self.setStyleSheet('QWidget{background:#000000;color:#e9eff4;font-size:13px;} QPushButton,QComboBox,QLineEdit{background:#121212;border:1px solid #4c5962;border-radius:4px;padding:7px;} QPushButton:hover{background:#262626;} QGroupBox{border:1px solid #3d454c;margin-top:12px;padding-top:12px;} QGroupBox::title{subcontrol-origin:margin;left:10px;} QTabBar::tab{background:#101010;padding:6px;border:1px solid #3d454c;} QTabBar::tab:selected{background:#242424;} QToolTip{background:#161616;color:#e9eff4;}')
  central=W.QWidget();self.setCentralWidget(central);grid=W.QGridLayout(central);self.canvas=Canvas();self.canvas.diagnostics=bool(cfg.get('diagnostics',False));grid.addWidget(self.canvas,0,0);side=W.QWidget();side.setFixedWidth(450);form=W.QVBoxLayout(side);grid.addWidget(side,0,1);self.curves=Curves();grid.addWidget(self.curves,1,0,1,2);grid.setRowStretch(0,1)
  title=W.QLabel('CHRONOS | OBSERVATION');title.setWordWrap(True);title.setStyleSheet('font-size:15px;font-weight:600');form.addWidget(title)
  group=W.QGroupBox('1  Source and session');g=W.QGridLayout(group);g.setColumnStretch(0,1);g.setColumnStretch(1,1);self.mode=W.QComboBox();self.mode.addItems(['REALTIME_REPLAY','LIVE_CAMERA']);self.mode.setCurrentText(cfg['source_mode']);g.addWidget(self.mode,0,0,1,2)
  self.speed=W.QComboBox();self.speed.addItems(['0.5','1.0','2.0']);self.speed.setCurrentText(str(cfg['speed']));self.loop=W.QCheckBox('Loop (reset session)');self.loop.setChecked(cfg['loop']);g.addWidget(self.speed,1,0);g.addWidget(self.loop,1,1)
  self.buttons={}
  for i,(name,fn) in enumerate([('Start',self.start),('Stop',self.stop),('Pause',self.pause),('Resume',self.resume),('Reconnect',self.reconnect),('Open records',self.open_records)]):
   b=W.QPushButton(name);b.clicked.connect(fn);g.addWidget(b,2+i//2,i%2);self.buttons[name]=b
  form.addWidget(group)
  tabs=W.QTabWidget();form.addWidget(tabs,1);observation=W.QWidget();observation_form=W.QVBoxLayout(observation);observation_form.setContentsMargins(4,4,4,4)
  observation_scroll=W.QScrollArea();observation_scroll.setWidgetResizable(True);observation_scroll.setWidget(observation);tabs.addTab(observation_scroll,'Observation')
  self.pzt_panel=PZTPanel(cfg);pzt_scroll=W.QScrollArea();pzt_scroll.setWidgetResizable(True);pzt_scroll.setWidget(self.pzt_panel);tabs.addTab(pzt_scroll,'PZT status / layout');tabs.setCurrentIndex(1);self.side_tabs=tabs
  form=observation_form
  group=W.QGroupBox('2  Image layers and channel mask');g=W.QGridLayout(group);g.setColumnStretch(0,1);g.setColumnStretch(1,1);self.checks={}
  for i,name in enumerate(['Candidates','PIV','Mask','Regions','Shape','Curves','Measured frame','Smooth PIV']):
   b=W.QCheckBox(name);b.setChecked(DEFAULT_LAYERS.get(name,name=='Curves'));self.checks[name]=b;g.addWidget(b,i//2,i%2)
  self.checks['Smooth PIV'].setToolTip('Stabilize the displayed field on a fixed grid. Scientific velocities and recorded raw vectors are unchanged. Uncheck to inspect raw arrows.')
  self.geometry_buttons={}
  for col,(name,fn) in enumerate([('Draw LIVE mask',self.draw_mask),('Clear current mask',self.clear_mask)]):
   b=W.QPushButton(name);b.clicked.connect(fn);g.addWidget(b,4,col);self.geometry_buttons[name]=b
  b=W.QPushButton('Load saved mask (optional)');b.clicked.connect(self.load_mask);g.addWidget(b,5,0,1,2);self.geometry_buttons['Load saved mask']=b
  self.mask_help=W.QLabel('Click boundary points; double-click or Enter to apply.\nBackspace: undo point. Esc: cancel.');self.mask_help.setWordWrap(True);self.mask_help.setStyleSheet('color:#b8cedc;font-size:12px;');g.addWidget(self.mask_help,6,0,1,2)
  self.canvas.mask_finished.connect(self.apply_drawn_mask);self.canvas.mask_cancelled.connect(lambda:self.mask_help.setText('Drawing cancelled. The previous mask is unchanged.'))
  if cfg.get('live_geometry'):
   try:self.canvas.active_geometry=json.loads(Path(cfg['live_geometry']).read_text())
   except (OSError,ValueError):self.mask_help.setText('Saved mask unavailable. Draw a new mask on the live image.')
  form.addWidget(group)
  group=W.QGroupBox('3  PZT event log (record only)');g=W.QGridLayout(group);g.setColumnStretch(0,1);g.setColumnStretch(1,1);self.event_source=W.QComboBox();self.event_source.addItems(['SIMULATED_EVENT','MANUAL_ANNOTATION']);g.addWidget(self.event_source,0,0,1,2);self.channel=W.QComboBox();self.channel.addItems(['ALL','PZT1','PZT2','PZT3','PZT4']);g.addWidget(self.channel,1,0,1,2);self.note=W.QLineEdit();self.note.setPlaceholderText('Conditions / event note');g.addWidget(self.note,2,0,1,2)
  for i,state in enumerate(['ON','OFF']):
   b=W.QPushButton('Mark '+state);b.clicked.connect(lambda _,s=state:self.mark_event(s));g.addWidget(b,3,i)
  form.addWidget(group);self.info=W.QLabel('Stopped');self.info.setWordWrap(True);self.info.setMinimumHeight(110);form.addWidget(self.info);form.addStretch()
  self.video_status=W.QLabel('Start records the experiment. Stop saves both videos.');self.statusBar().addPermanentWidget(self.video_status,1)
  for check in self.checks.values():check.toggled.connect(self.update_video_layers)
  self.mode.currentTextChanged.connect(self.update_source_controls);self.update_source_controls()
  self.timer=C.QTimer(self);self.timer.timeout.connect(self.tick);self.timer.start(round(1000/cfg['ui_hz']))
  if auto:C.QTimer.singleShot(100,self.start)
 def update_source_controls(self,*args):
  replay=self.mode.currentText()=='REALTIME_REPLAY'
  for widget in [self.speed,self.loop,self.buttons['Pause'],self.buttons['Resume'],self.event_source]:widget.setVisible(replay)
  self.event_source.setCurrentText('SIMULATED_EVENT' if replay else 'MANUAL_ANNOTATION')
 def start(self):
  if self.pipe and self.pipe.worker and self.pipe.worker.is_alive():return
  if self.pipe:self.pipe.stop()
  self.cfg.update(source_mode=self.mode.currentText(),speed=float(self.speed.currentText()),loop=self.loop.isChecked())
  self.cfg['video_layers']={k:v.isChecked() for k,v in self.checks.items()}
  if self.cfg['source_mode']=='REALTIME_REPLAY':self.cfg['source_sha256']=sha(self.cfg['video'])
  self.output.mkdir(parents=True,exist_ok=True);self.pipe=Pipeline(self.cfg,self.output/f'session_{time.time_ns()}');self.pipe.start();self.started=time.monotonic()
  if self.capture and self.record is None:self.record=ScreenRecorder(self.capture)
 def stop(self):
  if self.pipe:
   if self.pipe.video_record and not self.pipe.video_record.closed:
    self.video_status.setText('Finishing video files. Please wait...');W.QApplication.processEvents(C.QEventLoop.ExcludeUserInputEvents)
   self.pipe.stop();self.update_video_status()
 def pause(self):
  if self.pipe and self.cfg['source_mode']=='REALTIME_REPLAY':self.pipe.source.pause()
 def resume(self):
  if self.pipe:self.pipe.source.resume()
 def reconnect(self):self.stop();self.start()
 def open_records(self):
  directory=self.pipe.output if self.pipe else self.output
  try:ok,message=open_directory(directory)
  except OSError as e:ok,message=False,str(e)
  if ok:self.statusBar().showMessage(message,6000)
  else:
   box=W.QMessageBox(self);box.setWindowTitle('Records folder');box.setText('Unable to open the folder. Copy the path and open it in your file manager.');box.setInformativeText(str(directory.resolve()));box.setDetailedText(message)
   copy=box.addButton('Copy path',W.QMessageBox.ActionRole);box.addButton(W.QMessageBox.Close);box.exec()
   if box.clickedButton()==copy:W.QApplication.clipboard().setText(str(directory.resolve()))
 def update_video_layers(self,*args):
  if self.pipe:self.pipe.set_video_layers({k:v.isChecked() for k,v in self.checks.items()})
 def update_video_status(self):
  if not self.pipe or not self.pipe.video_record:return
  r=self.pipe.video_record.status();state=r['state'];dropped=r['dropped_frames']
  if r['error']:text='Recording failed. Check disk space and session records.'
  elif state=='SAVED':text='Videos saved. Select Open records to view them.'
  elif state=='NO_FRAMES':text='No images received. No video was created.'
  elif state=='WAITING_FOR_FRAMES':text='Recording ready. Waiting for images.'
  else:text='REC | Saving camera images and detection overlays'
  if dropped:text+=f' | Warning: {dropped} frames not recorded'
  self.video_status.setText(text);self.video_status.setStyleSheet('color:#ff967f;' if r['error'] or dropped else 'color:#b8dfce;')
 def mark_event(self,state):
  if self.pipe:self.pipe.event(state,self.event_source.currentText(),self.channel.currentText(),self.note.text())
 def live_packet(self):
  p=self.pipe.raw if self.pipe else None
  if p is None or p.source_mode!='LIVE_CAMERA' or self.mode.currentText()!='LIVE_CAMERA':
   self.mask_help.setText('Select LIVE_CAMERA, press Start, and wait for an image.');return None
  return p
 def draw_mask(self):
  p=self.live_packet()
  if p is None:return
  self.drawing_resolution=p.original_resolution;self.canvas.editing=True;self.canvas.points=[];self.canvas.setCursor(C.Qt.CrossCursor);self.canvas.setFocus();self.canvas.update()
  self.mask_help.setText('Click along the channel. Double-click, Enter, or click the first point to apply.\nBackspace: undo point. Esc: cancel. Changes save automatically.')
 def apply_drawn_mask(self):
  p=self.live_packet()
  if p is None:return
  try:
   if p.original_resolution!=self.drawing_resolution:raise ValueError('Camera resolution changed. Click Draw LIVE mask and redraw.')
   self.apply_geometry(new_geometry(self.canvas.points,p.original_resolution))
  except (ValueError,OSError) as e:self.mask_help.setText(str(e))
 def apply_geometry(self,g):
  # Immutable per-change snapshots make existing observations reproducible.
  path=self.output/'geometry'/f'live_mask_{time.time_ns()}.json';dump(path,g)
  self.cfg['live_geometry']=str(path);self.pipe.set_live_geometry(path)
  self.canvas.active_geometry=g;self.canvas.editing=False;self.canvas.points=[];self.canvas.unsetCursor();self.checks['Mask'].setChecked(True);self.canvas.update();self.curves.rows=[]
  self.mask_help.setText('Mask applied and saved. Draw again to replace it; Clear to remove it.')
 def clear_mask(self):
  if self.mode.currentText()!='LIVE_CAMERA':self.mask_help.setText('The replay mask is fixed. Clear is available in LIVE mode.');return
  self.cfg['live_geometry']=None
  if self.pipe and self.pipe.cfg['source_mode']=='LIVE_CAMERA':self.pipe.set_live_geometry(None)
  self.canvas.active_geometry=None;self.canvas.editing=False;self.canvas.points=[];self.canvas.unsetCursor();self.canvas.update();self.curves.rows=[]
  self.mask_help.setText('Mask cleared. Select Draw LIVE mask to draw a new region.')
 def load_mask(self):
  p=self.live_packet()
  if p is None:return
  path,_=W.QFileDialog.getOpenFileName(self,'Load saved LIVE mask',str(self.output),'JSON (*.json)')
  if path:
   try:self.apply_geometry(validate_geometry(json.loads(Path(path).read_text()),p.original_resolution))
   except (ValueError,OSError,KeyError,TypeError) as e:self.mask_help.setText('Mask not loaded: '+str(e))
 def tick(self):
  t=time.perf_counter()
  if self.pipe:
   self.update_video_status()
   s=self.pipe.snapshot();s['display_packet']=s['analyzed_packet'] if s.get('fresh') and self.checks['Measured frame'].isChecked() else s['raw'];self.canvas.snap=s;self.canvas.layers={k:v.isChecked() for k,v in self.checks.items()};self.curves.setVisible(self.checks['Curves'].isChecked());self.canvas.update()
   self.info.setText(observation_text(s))
   self.curves.rows=list(self.pipe.metric_history);self.curves.events=self.pipe.events;self.curves.update()
   if self.cfg.get('diagnostics',False):self.pipe.ui_times.append((time.perf_counter()-t)*1000);self.pipe.counter['ui_ticks']+=1
   if self.record:self.record.submit(self.grab().toImage(),time.monotonic()-self.started)
   if self.cfg.get('diagnostics',False) and not self.saved_preview and time.monotonic()-self.started>4:
    self.grab().save(str(self.output/'ui_actual.png'));self.saved_preview=True
   if self.elapsed_limit and time.monotonic()-self.started>=self.elapsed_limit:self.close()
 def closeEvent(self,e):
  self.timer.stop();self.pzt_panel.shutdown();self.stop()
  if self.record:self.record.close()
  summary=dict(closed=True,source_mode=self.cfg['source_mode'],last_session=str(self.pipe.output) if self.pipe else None)
  if self.cfg.get('diagnostics',False):summary.update(canvas_paint_ms=stats(self.canvas.paint_times),ui_state_update_ms=stats(self.pipe.ui_times if self.pipe else []))
  dump(self.output/'ui_summary.json',summary)
  e.accept()
