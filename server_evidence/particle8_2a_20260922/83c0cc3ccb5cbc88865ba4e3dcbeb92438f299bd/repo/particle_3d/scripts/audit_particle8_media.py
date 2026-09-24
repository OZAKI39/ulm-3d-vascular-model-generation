#!/usr/bin/env python3
"""Extract inspection sheets from actual encoded MP4 streams, not render buffers."""
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import imageio_ffmpeg
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from particle_3d.particle8_replay import DEFAULT_OUTPUT,read,write,digest
from particle_3d.particle8_visuals import ANIMATIONS
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT);a=p.parse_args();root=a.output
out=root/'inspection';out.mkdir(exist_ok=True)
font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',18)
records=[]
for name in ANIMATIONS:
 path=root/'animations'/(name+'.mp4');m=read(root/'data'/(name+'_frames.json'));n=m['frame_count']
 picks=sorted(set(np.linspace(0,n-1,9,dtype=int).tolist()))
 sheet=Image.new('RGB',(1440,1010),'white');draw=ImageDraw.Draw(sheet)
 draw.text((12,7),name+' | decoded MP4 frames',fill='#183449',font=font)
 reader=imageio_ffmpeg.read_frames(str(path),pix_fmt='rgb24');meta=next(reader);k=0
 for i,pixels in enumerate(reader):
  if i in picks:
   picture=Image.frombytes('RGB',meta['size'],pixels)
   row,col=divmod(k,3);x,y=col*480,40+row*322
   sheet.paste(picture.resize((480,300)),(x,y))
   t=m['frames'][i]['physical_or_display_time_s'];draw.text((x+8,y+298),f'frame {i:03d} | clock {t:.9f} s',fill='#183449',font=font);k+=1
  if i in m['keyframe_indices']:
   Image.frombytes('RGB',meta['size'],pixels).save(out/f'{name}_decoded_{i:04d}.png')
 target=out/(name+'_contact_sheet.png');sheet.save(target)
 records.append(dict(animation=str(path.relative_to(root)),animation_sha256=digest(path),
  contact_sheet=str(target.relative_to(root)),contact_sheet_sha256=digest(target),decoded_frame_indices=picks,
  total_decoded_frames=i+1,expected_frames=n,complete=(i+1)==n))
 print(name,i+1,'decoded',flush=True)
write(root/'data/decoded_contact_sheets.json',dict(role='ACTUAL_MP4_DECODE_NOT_RENDER_PREVIEW',records=records))
