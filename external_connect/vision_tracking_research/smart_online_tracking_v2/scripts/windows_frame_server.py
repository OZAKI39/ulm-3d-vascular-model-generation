"""Request/response IPC adapter around the ORIGINAL Windows Chronos receiver.

No new RTSP implementation, no dependency installation in the camera environment.
At most one request/response is in flight; the camera receiver retains one frame.
"""
import argparse,json,os,sys,struct,time,threading,queue,uuid,socket
from pathlib import Path

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',default='F:/CHRONOS_control')
    p.add_argument('--seconds',type=float,default=90);p.add_argument('--summary',required=True)
    p.add_argument('--roi',type=int,nargs=4);a=p.parse_args()
    if not 1<=a.seconds<=180:raise ValueError('Finite 1..180 s only')
    sys.path.insert(0,a.source)
    from chronos_common import load_config,video_options,LatestFrameStream
    cfg=load_config(Path(a.source)/'chronos_config.json')
    if os.environ.get('OPENCV_FFMPEG_CAPTURE_OPTIONS')!=video_options(cfg):
        raise RuntimeError('Parent must pass the original capture options at process creation')
    if os.name=='nt':
        import msvcrt
        msvcrt.setmode(sys.stdout.fileno(),os.O_BINARY)
        msvcrt.setmode(sys.stdin.fileno(),os.O_BINARY)
    for port in (cfg['udp_port'],cfg['udp_port']+1):
        with socket.socket(socket.AF_INET,socket.SOCK_DGRAM) as s:
            if hasattr(socket,'SO_EXCLUSIVEADDRUSE'):s.setsockopt(socket.SOL_SOCKET,socket.SO_EXCLUSIVEADDRUSE,1)
            s.bind(('0.0.0.0',port))
    messages=queue.Queue(maxsize=1);stop=threading.Event()
    def read_commands():
        for line in sys.stdin.buffer:
            messages.put(line)
            if line.strip()==b'QUIT':break
        stop.set()
    threading.Thread(target=read_commands,daemon=True).start()
    # Only this adapter process is subject to the finite emergency deadline.
    watchdog=threading.Timer(a.seconds+15,lambda:os._exit(124));watchdog.daemon=True;watchdog.start()
    stream=LatestFrameStream(cfg).start();session=uuid.uuid4().hex
    start=time.monotonic();cpu_start=time.process_time();records=[];previous=0
    try:
        while time.monotonic()-start<a.seconds and not stop.is_set():
            try:command=messages.get(timeout=.1).strip()
            except queue.Empty:continue
            if command==b'QUIT':break
            if command!=b'GET':raise ValueError('Unexpected IPC command')
            deadline=time.monotonic()+.05
            packet=stream.snapshot()
            while packet is not None and packet[0]==previous and time.monotonic()<deadline:
                time.sleep(.001);packet=stream.snapshot()
            stats=stream.stats();snapshot_at=time.monotonic()
            header={'session':f'{session}:{stats["connections"]}','source':'LIVE_CAMERA',
                    'camera_ip':cfg['camera_ip'],'received_frames':stats['received_frames'],
                    'connections':stats['connections'],'failures':stats['failures'],'receiver_state':stats['state'],
                    'input_fps':stats['receive_fps_last_5s'],'exposure_timestamp':None,
                    'timestamp_status':'EXPOSURE_TIMESTAMP_UNKNOWN','camera_config':cfg,
                    'status':'NO_NEW_FRAME','shape':None,'payload_bytes':0,'seq':previous}
            payload=b''
            if packet is not None and packet[0]!=previous:
                seq,decoded,utc,frame=packet;previous=seq
                h,w=frame.shape[:2];roi=a.roi or [0,0,w,h];x,y,rw,rh=roi
                if x<0 or y<0 or rw<1 or rh<1 or x+rw>w or y+rh>h:raise ValueError('ROI out of bounds')
                data=frame[y:y+rh,x:x+rw]
                payload=data.tobytes(order='C')
                header.update(status='FRAME',seq=seq,shape=list(data.shape),roi_xywh=roi,
                              original_resolution_wh=[w,h],payload_bytes=len(payload),
                              sender_decode_mono=decoded,sender_snapshot_mono=snapshot_at,
                              sender_frame_age_ms=(snapshot_at-decoded)*1000,decode_utc=utc,
                              prepare_payload_ms=(time.monotonic()-snapshot_at)*1000)
            raw=json.dumps(header,separators=(',',':')).encode()
            sending=time.monotonic()
            sys.stdout.buffer.write(b'CV2FRAME'+struct.pack('<I',len(raw))+raw+payload);sys.stdout.buffer.flush()
            records.append({'seq':header['seq'],'status':header['status'],'bytes':len(payload),
                            'send_write_ms':(time.monotonic()-sending)*1000,
                            'sender_frame_age_ms':header.get('sender_frame_age_ms')})
    finally:
        stream.stop();watchdog.cancel()
        summary={'receiver':stream.stats(),'wall_s':time.monotonic()-start,
                 'cpu_seconds':time.process_time()-cpu_start,'messages':records,
                 'source_reused':str(Path(a.source)/'chronos_common.py'),'camera_config':cfg,
                 'hardware_control':False,'max_inflight_responses':1}
        Path(a.summary).write_text(json.dumps(summary,indent=2))

if __name__=='__main__':
    sys.dont_write_bytecode=True
    main()
