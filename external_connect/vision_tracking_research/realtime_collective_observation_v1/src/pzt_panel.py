"""Fixed experimental apparatus with explicit AFG-channel-to-PZT assignment."""
import json
import math
import time
from pathlib import Path

from PySide6 import QtCore as C, QtGui as G, QtWidgets as W
from pzt_status import (StatusReader, EXPERIMENT_MODEL, ENDPOINTS, PHYSICAL_PZTS,
                        experiment_layout, endpoint_view, physical_views, save_layout, validate_layout)
from pzt_video import PZTVideoSource


class TopView(W.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(220)
        self.setSizePolicy(W.QSizePolicy.Expanding, W.QSizePolicy.Expanding)
        self.layout_data = experiment_layout()
        self.views = physical_views({}, self.layout_data)
        self.show_output = True
        self.setToolTip('XY top view of the current 32.5 × 9 × 9 mm PDMS model. Four 7 mm PZTs face inward.\n'
                        'Arrows show emitting normals. Top-port centres retain the model reference assumption.\n'
                        'Channel assignments are entered by the user, not sensed by the generator.')

    def paintEvent(self, event):
        p = G.QPainter(self)
        p.setRenderHint(G.QPainter.Antialiasing)
        p.fillRect(self.rect(), G.QColor('#000000'))
        p.setPen(G.QColor('#d9e6ee'))
        p.drawText(12, 20, 'Top view')
        p.drawText(self.width()-111, 20, '+X right · +Y up')
        model = EXPERIMENT_MODEL
        length, width = model['block_mm'][:2]
        scale = min(max(1, self.width()-130)/(length+2), max(1, self.height()-140)/(width+1))
        cx, cy = self.width()/2, (self.height()+14)/2
        def point(x, y):
            return C.QPointF(cx+x*scale, cy-y*scale)
        # Model footprint and the sealed water channel, using one XY scale.
        block = C.QRectF(cx-length*scale/2, cy-width*scale/2, length*scale, width*scale)
        p.setPen(G.QPen(G.QColor('#8199ac'), 1.2))
        p.setBrush(G.QColor(94, 121, 145, 45))
        p.drawRect(block)
        diameter = model['channel_diameter_mm']
        p.setBrush(G.QColor(91, 166, 190, 45))
        p.setPen(G.QPen(G.QColor('#71a0b3'), 1))
        p.drawRect(C.QRectF(block.left(), cy-diameter*scale/2, block.width(), diameter*scale))
        for port in model['ports']:
            p.setBrush(G.QColor('#243d4b'))
            p.setPen(G.QPen(G.QColor('#a3c0cf'), 1))
            r=port['diameter_mm']*scale/2
            p.drawEllipse(point(*port['center_xy_mm']), r, r)
        for name in PHYSICAL_PZTS:
            row=model['pzts'][name]
            x,y,_=row['center_mm']
            at=point(x,y)
            view=self.views.get(name,{})
            color=G.QColor(view.get('color','#a3afbd'))
            p.save()
            p.translate(at)
            p.rotate(-math.degrees(math.atan2(row['normal'][1],row['normal'][0])))
            disc_length=row['diameter_mm']*scale
            gradient=G.QLinearGradient(-6,0,0,0)
            gradient.setColorAt(0,G.QColor('#647c8e'))
            gradient.setColorAt(.65,G.QColor('#d4e0e7'))
            gradient.setColorAt(1,G.QColor('#8fa5b6'))
            p.setBrush(G.QBrush(gradient))
            p.setPen(G.QPen(color,1.2))
            p.drawRoundedRect(C.QRectF(-6,-disc_length/2,6,disc_length),1.5,1.5)
            p.setPen(G.QPen(color,1.3))
            p.drawLine(C.QPointF(3,0),C.QPointF(20,0))
            p.drawLine(C.QPointF(20,0),C.QPointF(15,-3.5))
            p.drawLine(C.QPointF(20,0),C.QPointF(15,3.5))
            p.restore()
            text=name+'\n'+(view.get('endpoint') or '—')
            if view.get('endpoint') and self.show_output:
                text+='\n'+view.get('output','UNKNOWN')
            if name=='X-': rect=C.QRectF(1,cy-28,max(58,block.left()-9),56)
            elif name=='X+': rect=C.QRectF(block.right()+8,cy-28,self.width()-block.right()-9,56)
            elif name=='Y+': rect=C.QRectF(cx-65,block.top()-63,130,54)
            else: rect=C.QRectF(cx-65,block.bottom()+9,130,54)
            p.setPen(color)
            font=p.font();font.setPixelSize(12);p.setFont(font)
            p.drawText(rect,C.Qt.AlignCenter,text)
        p.setPen(G.QPen(G.QColor('#b9cbd7'),1.2))
        p.drawLine(C.QPointF(14,self.height()-16),C.QPointF(14+5*scale,self.height()-16))
        p.drawText(int(20+5*scale),self.height()-12,'5 mm')


class ConnectionEditor(W.QDialog):
    def __init__(self,data,parent=None):
        super().__init__(parent)
        self.setWindowTitle('PZT connections')
        self.resize(650,620)
        outer=W.QVBoxLayout(self)
        self.preview=TopView()
        self.preview.show_output=False
        self.preview.setMinimumHeight(240)
        outer.addWidget(self.preview,1)
        grid=W.QGridLayout()
        grid.setColumnStretch(3,1)
        outer.addLayout(grid)
        for col,title in enumerate(('Generator','USB serial','Channel','Connected PZT')):
            grid.addWidget(W.QLabel(title),0,col)
        self.fields={}
        for i,(endpoint,row) in enumerate(data['assignments'].items(),1):
            for col,text in enumerate((row['device'],row['serial'],f"CH{row['channel']}")):
                grid.addWidget(W.QLabel(text),i,col)
            combo=W.QComboBox()
            combo.addItem('Unassigned',None)
            for pzt in PHYSICAL_PZTS:combo.addItem(pzt,pzt)
            combo.setCurrentIndex(combo.findData(row['pzt']))
            grid.addWidget(combo,i,3)
            self.fields[endpoint]=combo
        self.error=W.QLabel()
        self.error.setStyleSheet('color:#edaf91;')
        self.error.setMinimumHeight(20)
        outer.addWidget(self.error)
        buttons=W.QDialogButtonBox(W.QDialogButtonBox.Save|W.QDialogButtonBox.Cancel)
        self.save_button=buttons.button(W.QDialogButtonBox.Save)
        self.save_button.setText('Save connections')
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)
        for combo in self.fields.values():combo.currentIndexChanged.connect(self.validate_selection)
        self.validate_selection()

    def value(self):
        result=experiment_layout()
        for endpoint,combo in self.fields.items():result['assignments'][endpoint]['pzt']=combo.currentData()
        return validate_layout(result)

    def validate_selection(self,*args):
        try:
            data=self.value()
        except ValueError as exc:
            self.error.setText(str(exc));self.save_button.setEnabled(False)
            return
        self.error.setText('');self.save_button.setEnabled(True)
        self.preview.layout_data=data
        self.preview.views=physical_views({},data)
        self.preview.update()


class PZTPanel(W.QWidget):
    def __init__(self,cfg,parent=None,reader_factory=StatusReader):
        super().__init__(parent)
        self.layout_path=Path(cfg.get('pzt_layout_path',Path(__file__).resolve().parents[1]/'configs/pzt_layout.json')).expanduser().resolve()
        self.status_path=Path(cfg.get('pzt_status_path','/mnt/f/Tektronix_control/observation_status.json'))
        self.layout_data=experiment_layout()
        load_error=None
        if self.layout_path.exists():
            try:self.layout_data=validate_layout(json.loads(self.layout_path.read_text()))
            except (OSError,ValueError,KeyError,TypeError) as exc:load_error=f'Confirm connections: {exc}'
        self.reader=reader_factory(self.status_path)
        self.video_source=PZTVideoSource(self.reader,self.layout_data)
        self.data={}
        self.pending_request=None
        outer=W.QVBoxLayout(self);outer.setContentsMargins(0,0,0,0)
        self.hint=W.QLabel(load_error or '');self.hint.setWordWrap(True);self.hint.setVisible(bool(load_error));outer.addWidget(self.hint)
        cards=W.QGridLayout();cards.setColumnStretch(0,1);cards.setColumnStretch(1,1)
        self.cards={}
        for i,endpoint in enumerate(ENDPOINTS):
            box=W.QGroupBox(endpoint);box.setStyleSheet('QGroupBox{margin-top:10px;padding-top:5px;}');box.setMinimumWidth(185)
            layout=W.QVBoxLayout(box);layout.setContentsMargins(7,5,7,5)
            text=W.QLabel();text.setTextFormat(C.Qt.PlainText);text.setWordWrap(True);text.setMinimumHeight(58)
            layout.addWidget(text);cards.addWidget(box,i//2,i%2);self.cards[endpoint]=(box,text)
        outer.addLayout(cards)
        buttons=W.QGridLayout();buttons.setColumnStretch(0,1);buttons.setColumnStretch(1,1)
        self.edit_button=W.QPushButton('PZT connections');self.edit_button.clicked.connect(self.edit_connections)
        self.read_button=W.QPushButton('Read PZT status');self.read_button.clicked.connect(self.read_status)
        self.read_button.setToolTip('Query the four generator channels through the Windows controller. Does not change outputs.')
        buttons.addWidget(self.read_button,0,0);buttons.addWidget(self.edit_button,0,1);outer.addLayout(buttons)
        self.top=TopView();outer.addWidget(self.top,1)
        self.timer=C.QTimer(self);self.timer.timeout.connect(self.refresh);self.timer.start(50);self.refresh()

    def read_status(self):
        try:
            request_id=self.reader.request_refresh()
            self.pending_request=(request_id,time.monotonic())
            self.read_button.setEnabled(False);self.read_button.setText('Reading…')
            self.hint.clear();self.hint.hide()
        except (OSError,RuntimeError,ValueError) as exc:
            self.hint.setText(str(exc));self.hint.show()

    def edit_connections(self):
        dialog=ConnectionEditor(self.layout_data,self)
        if dialog.exec()==W.QDialog.Accepted:
            try:
                self.layout_data=save_layout(self.layout_path,dialog.value())
                self.video_source.set_layout(self.layout_data)
                self.hint.clear();self.hint.hide();self.refresh()
            except (OSError,ValueError,KeyError,TypeError) as exc:W.QMessageBox.warning(self,'Connections not saved',str(exc))

    def refresh(self):
        self.data=self.reader.snapshot()
        if self.pending_request:
            request_id,started=self.pending_request
            result=self.data.get('refresh') or {}
            failure=self.data.get('_request_error') or {}
            message=None
            if failure.get('request_id')==request_id:
                message='Read request failed. Check the Windows controller.'
            elif result.get('request_id')==request_id and result.get('state') in ('COMPLETE','PARTIAL','FAILED'):
                message='Status updated.' if result['state']=='COMPLETE' else 'Read incomplete. See channel status.'
            elif time.monotonic()-started>20:
                message='Read timed out. Check the AFG controller.'
            if message:
                self.pending_request=None;self.read_button.setEnabled(True);self.read_button.setText('Read PZT status')
                self.hint.setText(message);self.hint.show()
                if message=='Status updated.':C.QTimer.singleShot(3000,self.hint.hide)
        for endpoint,(box,label) in self.cards.items():
            view=endpoint_view(self.data,endpoint)
            target=self.layout_data['assignments'][endpoint]['pzt'] or 'Unassigned'
            title=f'{endpoint} → {target}'
            if box.title()!=title:box.setTitle(title)
            params=view['parameters']
            lines=[f"{view['connection']} | {view['output']}"]
            if params:
                mode='Burst' if params.get('burst_enabled') else str(params.get('frequency_mode','?'))
                lines += [f"{params['frequency_hz']/1e6:g} MHz · {params['amplitude_vpp']:g} Vpp",
                          f"Phase {params['phase_deg']:.1f}° · {params.get('waveform','?')} / {mode}"]
            else:lines+=['Parameters unavailable']
            text='\n'.join(lines)
            if label.text()!=text:label.setText(text)
            if label.property('status_color')!=view['color']:
                label.setStyleSheet(f"color:{view['color']};font-size:12px;")
                label.setProperty('status_color',view['color'])
            label.setToolTip('Generator readback; Vpp is not acoustic pressure. ON is not physical cable sensing.\n'+
                json.dumps(dict(identity=ENDPOINTS[endpoint],last_read_s=view['age_s'],parameters=params,error=view['error']),indent=2))
        views=physical_views(self.data,self.layout_data)
        drawing_key=tuple((n,v.get('endpoint'),v['output'],v['color']) for n,v in views.items())
        self.top.layout_data=self.layout_data;self.top.views=views
        if drawing_key!=getattr(self,'_drawing_key',None):
            self._drawing_key=drawing_key;self.top.update()

    def shutdown(self):
        self.timer.stop();self.reader.close()
