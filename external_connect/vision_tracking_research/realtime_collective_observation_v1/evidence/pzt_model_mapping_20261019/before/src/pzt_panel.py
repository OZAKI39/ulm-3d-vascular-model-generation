"""PZT telemetry cards and an explicitly user-defined XY top view."""
import json
import math
from pathlib import Path
import time

from PySide6 import QtCore as C, QtGui as G, QtWidgets as W
from pzt_status import NAMES, StatusReader, blank_layout, channel_view, save_layout, validate_layout
from desktop import open_directory


class LayoutEditor(W.QDialog):
    def __init__(self, data, parent=None):
        super().__init__(parent)
        self.setWindowTitle('PZT layout — manual lab coordinates')
        self.resize(880, 380)
        outer = W.QVBoxLayout(self)
        help_text = W.QLabel('View from +Z toward XY. +X is right; +Y is up. Coordinates are in mm.\n'
            'Normal angle: 0° points +X; 90° points +Y. Set the emitting direction.\n'
            'Enable only known positions. Wiring is your manual cable check, not a device measurement.')
        help_text.setWordWrap(True)
        outer.addWidget(help_text)
        grid = W.QGridLayout()
        outer.addLayout(grid)
        for col, text in enumerate(('PZT', 'Position set', 'X (mm)', 'Y (mm)', 'Normal (°)', 'Diameter (mm)', 'Wiring check')):
            grid.addWidget(W.QLabel(text), 0, col)
            if col > 1:
                grid.setColumnStretch(col, 1)
        self.fields = {}
        for i, name in enumerate(NAMES, 1):
            row = data['pzts'][name]
            grid.addWidget(W.QLabel(name), i, 0)
            placed = W.QCheckBox()
            placed.setChecked(row['placed'])
            grid.addWidget(placed, i, 1)
            fields = dict(placed=placed)
            for col, key in enumerate(('x_mm', 'y_mm', 'normal_deg', 'diameter_mm'), 2):
                spin = W.QDoubleSpinBox()
                spin.setDecimals(3)
                spin.setRange(.01 if key == 'diameter_mm' else -10000, 1000 if key == 'diameter_mm' else 10000)
                spin.setSingleStep(.1 if key != 'normal_deg' else 5.)
                spin.setValue(row[key])
                spin.setMinimumWidth(100)
                grid.addWidget(spin, i, col)
                fields[key] = spin
            wiring = W.QComboBox()
            wiring.addItems(['Not checked', 'Connected', 'Disconnected'])
            wiring.setCurrentText(row['wiring'])
            grid.addWidget(wiring, i, 6)
            fields['wiring'] = wiring
            self.fields[name] = fields
        outer.addStretch()
        buttons = W.QDialogButtonBox(W.QDialogButtonBox.Save | W.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)

    def value(self):
        result = blank_layout()
        for name, fields in self.fields.items():
            result['pzts'][name] = {k: (w.isChecked() if k == 'placed' else w.currentText() if k == 'wiring' else w.value())
                                    for k, w in fields.items()}
        return validate_layout(result)


class TopView(W.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(180)
        self.setSizePolicy(W.QSizePolicy.Expanding, W.QSizePolicy.Expanding)
        self.layout_data = blank_layout()
        self.views = {}
        self.setToolTip('Schematic XY top view. Bars represent vertical PZT discs; arrows are user-entered emitting normals.\n'
                        'No automatic camera registration, acoustic force or pressure field is inferred.')

    def paintEvent(self, event):
        p = G.QPainter(self)
        p.setRenderHint(G.QPainter.Antialiasing)
        p.fillRect(self.rect(), G.QColor('#101b26'))
        p.setPen(G.QColor('#d9e6ee'))
        p.drawText(12, 20, 'TOP VIEW · XY (mm) · user layout')
        rows = {n: r for n, r in self.layout_data['pzts'].items() if r['placed']}
        if not rows:
            p.setPen(G.QColor('#c6d5df'))
            p.drawText(self.rect().adjusted(20, 35, -20, -20), C.Qt.AlignCenter | C.Qt.TextWordWrap,
                       'Positions are not set.\nSelect Edit layout and enter the actual PZT positions.')
            return
        xs, ys = [0.], [0.]
        for row in rows.values():
            half = row['diameter_mm'] * .65
            xs += [row['x_mm'] - half, row['x_mm'] + half]
            ys += [row['y_mm'] - half, row['y_mm'] + half]
        lo_x, hi_x, lo_y, hi_y = min(xs), max(xs), min(ys), max(ys)
        # One isotropic scale: never distort relative positions or diameters.
        scale = min(max(1, self.width()-115) / max(hi_x-lo_x, 1.),
                    max(1, self.height()-100) / max(hi_y-lo_y, 1.))
        center_x, center_y = (hi_x+lo_x)/2, (hi_y+lo_y)/2
        def point(x, y):
            return C.QPointF(self.width()/2+(x-center_x)*scale, (self.height()+20)/2-(y-center_y)*scale)
        p.setPen(G.QPen(G.QColor('#4b6376'), 1, C.Qt.DashLine))
        p.drawLine(point(lo_x, 0), point(hi_x, 0))
        p.drawLine(point(0, lo_y), point(0, hi_y))
        p.setPen(G.QColor('#b6c9d6'))
        p.drawText(self.width()-115, 20, '+X right · +Y up')
        p.drawText(point(0, 0)+C.QPointF(4, 14), '0')
        for name, row in rows.items():
            at = point(row['x_mm'], row['y_mm'])
            view = self.views.get(name, {})
            color = G.QColor(view.get('color', '#a3afbd'))
            p.save()
            p.translate(at)
            p.rotate(-row['normal_deg'])
            length = row['diameter_mm'] * scale
            gradient = G.QLinearGradient(-4, 0, 4, 0)
            gradient.setColorAt(0, G.QColor('#596c7e'))
            gradient.setColorAt(.65, G.QColor('#ccd8e2'))
            gradient.setColorAt(1, G.QColor('#8298ab'))
            p.setBrush(G.QBrush(gradient))
            p.setPen(G.QPen(color, 1.5))
            p.drawRoundedRect(C.QRectF(-4, -length/2, 8, length), 2, 2)
            p.setPen(G.QPen(color, 1.6))
            p.drawLine(C.QPointF(6, 0), C.QPointF(25, 0))
            p.drawLine(C.QPointF(25, 0), C.QPointF(20, -4))
            p.drawLine(C.QPointF(25, 0), C.QPointF(20, 4))
            p.restore()
            p.setPen(color)
            # Put each label radially outside the group to keep the centre open.
            dx, dy = row['x_mm']-center_x, row['y_mm']-center_y
            label = f"{name} {view.get('output', 'UNKNOWN')}"
            rect = C.QRectF(at.x()-48, at.y()+max(12, length/2+5) if dy <= 0 else at.y()-max(32, length/2+25), 96, 20)
            if abs(dx) > abs(dy):
                rect = C.QRectF(at.x()+12 if dx > 0 else at.x()-108, at.y()-10, 96, 20)
            # Keep all labels visible even for layouts with strongly unequal spans.
            rect.moveLeft(max(3, min(rect.left(), self.width()-99)))
            rect.moveTop(max(26, min(rect.top(), self.height()-42)))
            p.fillRect(rect, G.QColor(16, 27, 38, 225))
            p.drawText(rect, C.Qt.AlignCenter, label)
        bar_mm = 10 ** math.floor(math.log10(max((hi_x-lo_x)/5, .01)))
        p.setPen(G.QPen(G.QColor('#c1d1db'), 1.4))
        p.drawLine(C.QPointF(14, self.height()-17), C.QPointF(14+bar_mm*scale, self.height()-17))
        p.drawText(int(20+bar_mm*scale), self.height()-12, f'{bar_mm:g} mm')
        p.drawText(self.width()-160, self.height()-12, 'Arrows: emitting normal')


class PZTPanel(W.QWidget):
    def __init__(self, cfg, parent=None, reader_factory=StatusReader):
        super().__init__(parent)
        self.layout_path = Path(cfg.get('pzt_layout_path', Path(__file__).resolve().parents[1]/'configs/pzt_layout.json'))
        self.status_path = Path(cfg.get('pzt_status_path', '/mnt/f/Tektronix_control/observation_status.json'))
        self.layout_data = blank_layout()
        load_error = None
        if self.layout_path.exists():
            try:
                self.layout_data = validate_layout(json.loads(self.layout_path.read_text()))
            except (OSError, ValueError, KeyError, TypeError) as exc:
                load_error = f'Layout not loaded: {exc}'
        self.reader = reader_factory(self.status_path)
        self.data = {}
        outer = W.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        self.hint = W.QLabel(load_error or 'Read-only. AFG controller → Connect Devices (Read Only).')
        self.hint.setWordWrap(True)
        outer.addWidget(self.hint)
        cards = W.QGridLayout()
        cards.setColumnStretch(0, 1)
        cards.setColumnStretch(1, 1)
        self.cards = {}
        for i, name in enumerate(NAMES):
            box = W.QGroupBox(name)
            box.setStyleSheet('QGroupBox{margin-top:10px;padding-top:5px;}')
            box.setMinimumWidth(185)
            layout = W.QVBoxLayout(box)
            layout.setContentsMargins(7, 5, 7, 5)
            text = W.QLabel()
            text.setTextFormat(C.Qt.PlainText)
            text.setWordWrap(True)
            text.setMinimumHeight(72)
            layout.addWidget(text)
            cards.addWidget(box, i//2, i%2)
            self.cards[name] = (box, text)
        outer.addLayout(cards)
        self.top = TopView()
        buttons = W.QGridLayout()
        buttons.setColumnStretch(0, 1)
        buttons.setColumnStretch(1, 1)
        self.edit_button = W.QPushButton('Edit layout')
        self.edit_button.clicked.connect(self.edit_layout)
        folder = W.QPushButton('AFG controller folder')
        folder.clicked.connect(self.open_controller)
        buttons.addWidget(self.edit_button, 0, 0)
        buttons.addWidget(folder, 0, 1)
        outer.addLayout(buttons)
        outer.addWidget(self.top, 1)
        note = W.QLabel('Vpp: generator readback, not acoustic Pa.\nWiring: manual check; not sensed. ON: output enabled, not proof of emitted sound.')
        note.setWordWrap(True)
        note.setStyleSheet('color:#b8cedc;font-size:12px;')
        outer.addWidget(note)
        self.timer = C.QTimer(self)
        self.timer.timeout.connect(self.refresh)
        self.timer.start(500)
        self.refresh()

    def open_controller(self):
        try:
            ok, message = open_directory(self.status_path.parent)
            self.hint.setText('Open start_gui.cmd, then Connect Devices (Read Only).' if ok else message)
        except OSError as exc:
            self.hint.setText(str(exc))

    def edit_layout(self):
        dialog = LayoutEditor(self.layout_data, self)
        if dialog.exec() == W.QDialog.Accepted:
            try:
                self.layout_data = save_layout(self.layout_path, dialog.value())
                self.hint.setText('Layout saved. Positions and wiring checks are user-entered.')
                self.refresh()
            except (OSError, ValueError, KeyError, TypeError) as exc:
                W.QMessageBox.warning(self, 'Layout not saved', str(exc))

    def refresh(self):
        self.data = self.reader.snapshot()
        views = {}
        for name, (box, label) in self.cards.items():
            view = channel_view(self.data, name)
            views[name] = view
            box.setTitle(f"{name} · {view['mapping']}" if view['mapping'] else name)
            params = view['parameters']
            lines = [f"{view['connection']} | {view['output']}"]
            if params:
                mode = 'Burst' if params.get('burst_enabled') else str(params.get('frequency_mode', '?'))
                lines += [f"{params['frequency_hz']/1e6:g} MHz  ·  {params['amplitude_vpp']:g} Vpp",
                          f"Phase {params['phase_deg']:.1f}° · {params.get('waveform', '?')} / {mode}"]
            else:
                lines += ['Parameters: unavailable', 'No fresh parameter readback']
            age = view['age_s']
            if age is not None:
                lines.append(f'Last read: {max(age, 0):.1f} s ago')
            else:
                lines.append('Last read: none')
            lines.append('Wire: ' + self.layout_data['pzts'][name]['wiring'])
            label.setText('\n'.join(lines))
            label.setStyleSheet(f"color:{view['color']};font-size:12px;")
            detail = dict(parameters=params, error=view['error'], serial=self.data.get('channels', {}).get(name, {}).get('serial'))
            label.setToolTip(json.dumps(detail, indent=2, ensure_ascii=False))
        self.top.layout_data = self.layout_data
        self.top.views = views
        self.top.update()

    def shutdown(self):
        self.timer.stop()
        self.reader.close()
