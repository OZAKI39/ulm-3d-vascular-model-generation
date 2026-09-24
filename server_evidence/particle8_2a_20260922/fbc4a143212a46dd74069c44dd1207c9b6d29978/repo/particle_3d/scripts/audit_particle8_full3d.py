#!/usr/bin/env python3
"""Complete MP4/GIF decode plus contact sheets from actual encoded frames."""
from pathlib import Path
import sys,argparse,hashlib
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import imageio_ffmpeg
from PIL import Image,ImageDraw
from particle_3d.particle8_replay import read,write,digest
from particle_3d.particle8_full3d_data import OUTPUT,CASES,LABELS
from particle_3d.particle8_full3d_visuals import text
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=OUTPUT);a=p.parse_args();root=a.output
rows=[];overview=Image.new('RGB',(1800,2100),'white');overview_draw=ImageDraw.Draw(overview)
text(overview_draw,(25,12),'Particle-8 | Full 3D replay — decoded first / middle / last frames',26,bold=True)
for row,case in enumerate(CASES):
 name='particle8_full3d_'+case;manifest=read(root/'frames'/(name+'_manifest.json'));count=manifest['frame_count'];indices=[0,count//2,count-1]
 reader=imageio_ffmpeg.read_frames(str(root/'mp4'/(name+'.mp4')),pix_fmt='rgb24');meta=next(reader)
 sheet=Image.new('RGB',(1800,450),'white');draw=ImageDraw.Draw(sheet);text(draw,(15,8),case+' | decoded MP4',22,bold=True)
 hashes=set();found=[]
 for i,pixels in enumerate(reader):
  hashes.add(hashlib.sha256(pixels).hexdigest())
  if i in indices:
   col=indices.index(i);im=Image.frombytes('RGB',meta['size'],pixels);im.save(root/'keyframes'/f'{name}_decoded_{i:04d}.png')
   sheet.paste(im.resize((600,375)),(col*600,43));frame=manifest['frames'][i]
   label=f'Frame {i} | t={frame["physical_time_s"]:.9f} s | az={frame["camera_main"]["azimuth_deg"]:.1f}°'
   text(draw,(col*600+8,422),label,16);found.append(i)
   overview.paste(im.resize((580,363)),(col*600+10,row*400+65))
 gif=Image.open(root/'gif'/(name+'.gif'))
 for k in range(gif.n_frames):gif.seek(k);gif.convert('RGB').load()
 sheet.save(root/'storyboard'/(name+'_storyboard.png'))
 rows.append(dict(case=case,mp4=str(Path('mp4')/(name+'.mp4')),gif=str(Path('gif')/(name+'.gif')),
  mp4_sha256=digest(root/'mp4'/(name+'.mp4')),gif_sha256=digest(root/'gif'/(name+'.gif')),
  decoded_mp4_frames=i+1,decoded_gif_frames=gif.n_frames,distinct_video_frames=len(hashes),
  frame_size=meta['size'],fps=meta['fps'],codec=meta['codec'],duration_s=(i+1)/meta['fps'],
  sampled_frame_indices=found,pass_checks=(i+1)==count and found==indices and len(hashes)>40 and gif.n_frames>=30))
 print(case,i+1,'MP4 frames;',gif.n_frames,'GIF frames',flush=True)
overview.save(root/'storyboard/full3d_storyboard.png')
summary=Image.new('RGB',(1600,1000),'white');d=ImageDraw.Draw(summary)
text(d,(60,45),'Particle-8 | Scientific scope and display choices',32,bold=True)
items=[('REAL FROZEN GEOMETRY','Full original wall + inlet + 3 outlets. No synthetic particles placed in it.'),
 ('REAL SINGLE MB','Saved P6.5 via P7: 1 ms / 0.440609 µm. Not full-vessel transit.'),
 ('REAL MIXED SMOKE','One admitted RBC surrogate held STATIC; pending is counted, never drawn.'),
 ('SYNTHETIC CONTROL','10 nm open section; finite shapes straddle caps. Bookkeeping only.'),
 ('CAMERA','Orthographic azimuth 35° → 115° (or fixed 35°); elevation 28°.'),
 ('TAILS','Recent saved polyline; only clipped display endpoints are interpolated.'),
 ('DISPLAY vs PHYSICS','Synthetic z ×20,000 is DISPLAY ONLY. Physical coordinates stay in SI.'),
 ('CORE LIMITATION','Real RBC continuous passage remains unresolved. No full suspension / PK.')]
for i,(title,body) in enumerate(items):
 y=125+i*82;text(d,(60,y),title,19,'#168eb1' if i<3 else '#c95444',True);text(d,(60,y+28),body,19)
for i,label in enumerate(LABELS[2:]):text(d,(60,835+i*35),label,14,'#8a4036')
summary.save(root/'storyboard/scientific_scope_summary.png')
write(root/'data/media_audit.json',dict(animations=rows,all_pass=all(r['pass_checks'] for r in rows),
 decoder=imageio_ffmpeg.get_ffmpeg_version(),method='EVERY_FRAME_DECODED; FIRST_MIDDLE_LAST_EXTRACTED_FROM_MP4'))
