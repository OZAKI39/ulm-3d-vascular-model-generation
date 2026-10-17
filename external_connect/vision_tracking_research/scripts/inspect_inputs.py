from pathlib import Path
import av,json,cv2,time,socket,platform,subprocess
from PIL import Image,ImageDraw
ROOT=Path(__file__).resolve().parents[1]
video=Path(json.loads((ROOT/'evidence/input_identity.json').read_text())['path'])
c=av.open(str(video)); s=c.streams.video[0]
meta={'codec':s.codec_context.name,'width':s.width,'height':s.height,'average_rate':str(s.average_rate),'base_rate':str(s.base_rate),'time_base':str(s.time_base),'frames':s.frames,'duration_s':float(c.duration/av.time_base) if c.duration else None,'stream_metadata':dict(s.metadata),'container_metadata':dict(c.metadata),'capture_fps':'UNKNOWN','pixel_scale_um_per_px':'UNKNOWN','interpretation':'encoded playback clock, not verified camera exposure clock'}
(ROOT/'evidence/video_metadata.json').write_text(json.dumps(meta,indent=2)); print(meta)
ids=[0,30,120,240,480,720,960,1200,1500,1800,2100,2400]; thumbs=[]; selected=[]; pts=[]
for i,frame in enumerate(c.decode(s)):
 pts.append(float(frame.time))
 if i in ids:
  bgr=frame.to_ndarray(format='bgr24'); cv2.imwrite(str(ROOT/f'evidence/frame_{i:05d}.png'),bgr)
  im=Image.fromarray(cv2.cvtColor(bgr,cv2.COLOR_BGR2RGB)); im.thumbnail((560,330)); panel=Image.new('RGB',(560,360),'white'); panel.paste(im,(0,24)); ImageDraw.Draw(panel).text((8,5),f'Frame {i} | PTS {frame.time:.3f} s',fill='black'); thumbs.append(panel); selected.append(i)
c.close(); sheet=Image.new('RGB',(560*3,360*((len(thumbs)+2)//3)),'#eeeeee')
for j,im in enumerate(thumbs):sheet.paste(im,((j%3)*560,(j//3)*360))
sheet.save(ROOT/'evidence/contact_sheet.png')
import numpy as np
(ROOT/'evidence/video_pts.json').write_text(json.dumps({'decoded_frames':len(pts),'first_pts_s':pts[0],'last_pts_s':pts[-1],'delta_unique_s':np.unique(np.round(np.diff(pts),8)).tolist(),'nonincreasing':int(np.sum(np.diff(pts)<=0)),'inspection_frames':selected},indent=2))
probe=[]
for port in [80,554]:
 t=time.monotonic()
 try:
  with socket.create_connection(('192.168.0.10',port),timeout=2):status='TCP_REACHABLE'
 except OSError as e:status=f'{type(e).__name__}: {e}'
 probe.append({'ip':'192.168.0.10','port':port,'result':status,'elapsed_s':time.monotonic()-t})
(ROOT/'evidence/live_probe.json').write_text(json.dumps({'probe':'TCP connect only; no instrument or camera setting commands','results':probe,'live_status':'LIVE_CAMERA_NOT_AVAILABLE' if all('TCP_REACHABLE'!=p['result'] for p in probe) else 'READ_ONLY_STREAM_TEST_NEEDED'},indent=2)); print(probe)
