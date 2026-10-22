"""Read-only Chronos focus aid: raw crop, relative edge score and optional display contrast.

Uses HTTP snapshots so the existing 60 fps RTSP preview can continue.
Never changes camera settings. Saved images are unmodified decoded frames.
Run only one HTTP screenshot/focus tool at a time on this legacy camera.
"""
from collections import deque
from datetime import datetime, timezone
from pathlib import Path
import argparse
import json
import threading
import time
import urllib.request

import cv2
import numpy as np

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seconds',type=float,default=360)
    parser.add_argument('--output',default=str(Path(__file__).resolve().parent/'data'/datetime.now().strftime('focus_%Y%m%d_%H%M%S')))
    args=parser.parse_args()
    folder=Path(args.output);folder.mkdir(parents=True,exist_ok=True)
    stop=threading.Event();lock=threading.Lock();latest=[None];errors=[]
    def reader():
        opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
        seq=0
        while not stop.is_set():
            try:
                with opener.open('http://192.168.0.10/cgi-bin/screenCap',timeout=4) as response:
                    data=response.read(4_000_000)
                raw=cv2.imdecode(np.frombuffer(data,np.uint8),cv2.IMREAD_COLOR)
                if raw is None:raise ValueError('Invalid camera image')
                seq+=1
                with lock:latest[0]=(seq,time.monotonic(),raw)
            except Exception as exc:
                errors.append(str(exc));del errors[:-10]
            stop.wait(.10)
    worker=threading.Thread(target=reader,daemon=True);worker.start()
    title='Chronos FOCUS | left: overview | right: 2x crop | C: display contrast | Q: finish'
    cv2.namedWindow(title,cv2.WINDOW_NORMAL);cv2.resizeWindow(title,1280,610)
    center=[670,515];reset=[False]
    def click(event,x,y,flags,param):
        if event==cv2.EVENT_LBUTTONDOWN and 0<=x<640 and 0<=y<512:
            center[:]=[int(x*2),int(y*2)];reset[0]=True
    cv2.setMouseCallback(title,click)
    history=deque(maxlen=5);best=0;baseline=None;last_seq=-1;contrast=False;count=0
    started=time.monotonic();last_write=0;last_log=0;last_status={};panel=None;records=[]
    try:
        while time.monotonic()-started<args.seconds:
            with lock:packet=latest[0]
            if packet and packet[0]!=last_seq:
                last_seq,received,raw=packet
                h,w=raw.shape[:2];x=max(0,min(w-320,center[0]-160));y=max(0,min(h-256,center[1]-128))
                if reset[0]:history.clear();best=0;baseline=None;reset[0]=False
                crop=raw[y:y+256,x:x+320]
                gray=cv2.cvtColor(crop,cv2.COLOR_BGR2GRAY).astype(np.float32)
                smooth=cv2.GaussianBlur(gray,(3,3),0)
                gx=cv2.Sobel(smooth,cv2.CV_32F,1,0,ksize=3)
                gy=cv2.Sobel(smooth,cv2.CV_32F,0,1,ksize=3)
                score=float(np.mean(gx*gx+gy*gy)/max(float(np.mean(gray))**2,1)*1000)
                history.append(score);stable=float(np.median(history))
                if baseline is None and len(history)==5:
                    baseline=stable;cv2.imwrite(str(folder/'focus_start.png'),raw)
                if len(history)==5 and stable>best:
                    best=stable;cv2.imwrite(str(folder/'focus_best.png'),raw)
                cv2.imwrite(str(folder/'focus_latest.png'),raw)
                overview=cv2.resize(raw,(640,512),interpolation=cv2.INTER_AREA)
                detail=crop.copy()
                if contrast:
                    # Fixed display mapping; never used by scoring or saved images.
                    detail=np.clip((detail.astype(np.float32)-95)*(255/100),0,255).astype(np.uint8)
                detail=cv2.resize(detail,(640,512),interpolation=cv2.INTER_NEAREST)
                cv2.rectangle(overview,(x//2,y//2),((x+320)//2,(y+256)//2),(0,255,255),2)
                panel=np.zeros((610,1280,3),np.uint8);panel[:512]=np.hstack([overview,detail])
                lines=[f'Focus score: {stable:.2f}   Best: {best:.2f}   Baseline: {(baseline or 0):.2f}',
                       'Keep lighting/position fixed. Turn focus slowly; pause 1 second. Compare channel edges.',
                       f'C: contrast aid {"ON (display only)" if contrast else "OFF"} | click left to move ROI | B: reset best | S: save | Q: finish']
                for i,line in enumerate(lines):cv2.putText(panel,line,(14,538+i*26),cv2.FONT_HERSHEY_SIMPLEX,.60,(240,240,240),1,cv2.LINE_AA)
                last_status={'elapsed_s':round(time.monotonic()-started,2),'score':stable,'best':best,'baseline':baseline,
                             'roi_xywh':[x,y,320,256],'contrast_display_only':contrast,'errors':errors[-5:]}
                records.append(last_status.copy())
                temporary=folder/'focus_status.tmp'
                temporary.write_text(json.dumps(last_status,indent=2),encoding='utf-8');temporary.replace(folder/'focus_status.json')
                if time.monotonic()-last_log>10:
                    print(json.dumps(last_status),flush=True);last_log=time.monotonic()
            if panel is not None:cv2.imshow(title,panel)
            key=cv2.waitKey(20)&0xff
            if key in (27,ord('q'),ord('Q')):break
            if key in (ord('c'),ord('C')):contrast=not contrast
            if key in (ord('b'),ord('B')):history.clear();best=0;baseline=None
            if key in (ord('s'),ord('S')) and packet:
                count+=1;cv2.imwrite(str(folder/f'manual_{count:02d}.png'),packet[2])
            if cv2.getWindowProperty(title,cv2.WND_PROP_VISIBLE)<1:break
            if packet is None and time.monotonic()-started>15:raise RuntimeError('No camera frames: '+str(errors))
    finally:
        stop.set();worker.join(5);cv2.destroyAllWindows()
        (folder/'focus_history.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
        print('Focus aid finished. Camera settings were not changed.',flush=True)

if __name__=='__main__':main()
