"""Frame-associated, read-only PZT context and a camera-video sidebar (no Qt capture)."""
from collections import deque
import copy
import hashlib
import json
import threading
import time

import cv2
import numpy as np

from pzt_status import ENDPOINTS, EXPERIMENT_MODEL, PHYSICAL_PZTS, endpoint_view, validate_layout
from video_text import draw_text


class PZTVideoSource:
    """Share the monitor's existing reader; no additional file or VISA polling."""
    def __init__(self, reader, layout):
        self.reader = reader
        self.lock = threading.Lock()
        self.layouts = deque(maxlen=128)
        self.set_layout(layout)

    def set_layout(self, layout):
        layout = validate_layout(layout)
        with self.lock:
            revision = self.layouts[-1][1] + 1 if self.layouts else 0
            self.layouts.append((time.monotonic(), revision, layout))

    def capture(self, packet):
        # The local frame-arrival time is the association boundary, not encoder
        # time or camera exposure time. Older preview clips have no AFG history.
        frame_wall = time.time() + packet.arrival_monotonic - time.monotonic()
        with self.lock:
            layout = next((r for r in reversed(self.layouts) if r[0] <= packet.arrival_monotonic), None)
            layout = copy.deepcopy(layout)
        sample = dict(scope='NOT_RECORDED', layout=layout[2] if layout else None,
                      layout_revision=layout[1] if layout else None, observation=None)
        if packet.source_mode == 'LIVE_CAMERA':
            observation = self.reader.snapshot_at(packet.arrival_monotonic)
            sample.update(scope='LIVE_HOST_READBACK', observation=observation)
            if observation:
                frame_wall = observation['observed_unix_s'] + packet.arrival_monotonic - observation['observed_monotonic']
        else:
            sample['scope'] = 'REPLAY_NO_SYNCHRONIZED_PZT'
        return dict(frame_unix_s=frame_wall, sample=sample)


def sample_id(sample):
    return hashlib.sha256(json.dumps(sample, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()[:24]


def frame_views(context):
    sample = (context or {}).get('sample') or {}
    observation = sample.get('observation') or {}
    data = observation.get('data') or {}
    if sample.get('scope') != 'LIVE_HOST_READBACK':
        data = {}
    return {key: endpoint_view(data, key, (context or {}).get('frame_unix_s', 0.)) for key in ENDPOINTS}


def _bgr(hex_color):
    return tuple(int(hex_color[i:i+2], 16) for i in (5, 3, 1))


class PZTVideoRenderer:
    """Cache the sidebar between readback changes; preserve native camera pixels."""
    def __init__(self):
        self.key = None
        self.panel = None

    def compose(self, camera, context):
        sample = (context or {}).get('sample') or dict(scope='NOT_RECORDED')
        views = frame_views(context)
        # Expiry must be reevaluated at each frame, even if the file stopped updating.
        display = {key: {k: v for k, v in row.items() if k not in ('age_s', 'error')} for key, row in views.items()}
        key = (json.dumps(sample.get('layout'), sort_keys=True), sample.get('scope'), json.dumps(display, sort_keys=True))
        if key != self.key:
            self.panel = self.draw(sample, views)
            self.key = key
        scale = max(1., min(2., camera.shape[0]/720.))
        panel = self.panel if scale == 1. else cv2.resize(self.panel, None, fx=scale, fy=scale, interpolation=cv2.INTER_LINEAR)
        h, w = camera.shape[:2]
        canvas = np.zeros((max(h, panel.shape[0]), w+panel.shape[1], 3), dtype=np.uint8)
        canvas[(canvas.shape[0]-h)//2:(canvas.shape[0]-h)//2+h, :w] = camera
        canvas[:panel.shape[0], w:] = panel
        return canvas

    @staticmethod
    def draw(sample, views):
        image = np.zeros((720, 480, 3), dtype=np.uint8)
        white = (232, 226, 216)
        muted = (185, 172, 155)
        def text(s, x, y, scale=.50, color=white):
            draw_text(image, s, x, y, max(12,round(scale*36)),color)
        text('PZT', 20, 34, .65)
        live = sample.get('scope') == 'LIVE_HOST_READBACK'
        layout = sample.get('layout') or {}
        assignments = layout.get('assignments', {})
        for i, endpoint in enumerate(ENDPOINTS):
            x, y = 16+(i%2)*232, 72+(i//2)*130
            view = views[endpoint]
            color = _bgr(view['color'])
            target = assignments.get(endpoint, {}).get('pzt') or 'Unassigned'
            cv2.rectangle(image, (x,y), (x+216,y+116), (83,74,63), 1)
            text(f'{endpoint} → {target}', x+8,y+20,.49)
            state = view['output'] if live else 'NOT RECORDED'
            text(state, x+8,y+42,.50,color)
            params = view['parameters']
            if params:
                text(f"{params['frequency_hz']/1e6:g} MHz | {params['amplitude_vpp']:g} Vpp", x+8,y+63,.46)
                text(f"Phase {params['phase_deg']:.1f}°", x+8,y+83,.46)
                mode = 'Burst' if params.get('burst_enabled') else params.get('frequency_mode','?')
                text(f"{params.get('waveform','?')} / {mode}", x+8,y+103,.43,muted)
            else:
                if live:text(view['connection'], x+8,y+70,.43,muted)
        scale, cx, cy = 10., 240, 505
        length, width = EXPERIMENT_MODEL['block_mm'][:2]
        def point(x,y): return (round(cx+x*scale),round(cy-y*scale))
        a, b = point(-length/2,width/2), point(length/2,-width/2)
        cv2.rectangle(image,a,b,(25,21,18),-1)
        cv2.rectangle(image,a,b,(172,153,129),1)
        d = EXPERIMENT_MODEL['channel_diameter_mm']
        cv2.rectangle(image,point(-length/2,d/2),point(length/2,-d/2),(50,33,24),-1)
        cv2.rectangle(image,point(-length/2,d/2),point(length/2,-d/2),(179,160,113),1)
        for port in EXPERIMENT_MODEL['ports']:
            cv2.circle(image,point(*port['center_xy_mm']),round(port['diameter_mm']*scale/2),(207,192,163),1,cv2.LINE_AA)
        for name in PHYSICAL_PZTS:
            pzt = EXPERIMENT_MODEL['pzts'][name]
            x,y = point(*pzt['center_mm'][:2])
            nx,ny = pzt['normal'][:2];ny = -ny
            endpoint = next((key for key,row in assignments.items() if row['pzt']==name), None)
            view = views.get(endpoint, {})
            color = _bgr(view.get('color','#a3afbd'))
            r = round(pzt['diameter_mm']*scale/2)
            for k,shade in enumerate((105,155,220,190,140)):
                off = k-2
                if nx:cv2.line(image,(x+off,y-r),(x+off,y+r),(shade+min(20,255-shade),shade,shade-15),1)
                else:cv2.line(image,(x-r,y+off),(x+r,y+off),(shade+min(20,255-shade),shade,shade-15),1)
            cv2.arrowedLine(image,(x+nx*6,y+ny*6),(x+nx*23,y+ny*23),color,1,cv2.LINE_AA,tipLength=.25)
            tx,ty = (10,492) if name=='X-' else (410,492) if name=='X+' else (204,393) if name=='Y+' else (204,585)
            text(name,tx,ty,.48,color)
            text(endpoint or '--',tx,ty+19,.40,color)
            text(view.get('output','--') if live else 'N/A',tx,ty+38,.39,color)
        cv2.line(image,(20,651),(70,651),muted,1)
        text('5 mm',78,656,.42,muted)
        return image
