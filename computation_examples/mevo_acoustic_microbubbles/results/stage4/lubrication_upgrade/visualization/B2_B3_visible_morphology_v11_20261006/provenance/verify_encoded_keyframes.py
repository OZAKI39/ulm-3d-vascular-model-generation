"""Check actual H.264 keyframes against rendered blue spheres, not only JSON geometry."""
import gzip
import json
from pathlib import Path
import subprocess
import sys
import imageio_ffmpeg
import numpy as np
from PIL import Image

out=Path(__file__).resolve().parents[1]
case=sys.argv[1]
path=out/f'stage4L_{case}_full_vessel_morphology_x24_4k60.mp4'
command=[imageio_ffmpeg.get_ffmpeg_exe(),'-v','error','-threads','2','-i',str(path),
         '-vf',"select='eq(n,0)+eq(n,600)+eq(n,1200)'",'-vsync','0','-f','rawvideo','-pix_fmt','rgb24','pipe:1']
video=subprocess.run(command,capture_output=True,check=True)
frames=np.frombuffer(video.stdout,dtype=np.uint8).reshape(3,2160,3840,3)
projections=json.load(gzip.open(out/f'{case}_SPHERE_PROJECTIONS.json.gz','rt'))
rows=[]
for frame,k in zip(frames,[0,150,300],strict=True):
    expected=np.asarray(Image.open(out/'keyframes'/f'{case}_{k:03d}.png').convert('RGB'))
    difference=abs(frame.astype(float)-expected.astype(float))
    assert difference.mean()<1.5
    mask=np.zeros((2160,3840),bool)
    p=projections[str(k)]
    for center,radius in zip(p['centers_px'],p['radii_px'],strict=True):
        x,y=center
        x0,x1=int(np.floor(x-radius)),int(np.ceil(x+radius))
        y0,y1=int(np.floor(y-radius)),int(np.ceil(y+radius))
        X,Y=np.meshgrid(np.arange(x0,x1)+.5,np.arange(y0,y1)+.5)
        mask[y0:y1,x0:x1]|=(X-x)**2+(Y-y)**2<=(radius*.85)**2
    def blue(im):
        c=im.astype(float)
        return (c[:,:,2]>1.18*c[:,:,1])&(c[:,:,2]>1.4*c[:,:,0])&(c[:,:,2]>60)
    visible=mask&blue(expected)
    assert visible.sum()>100
    recall=float(np.count_nonzero(visible&blue(frame))/visible.sum())
    assert recall>.95,recall
    rows.append(dict(video_frame=k*4,source_state=k,mean_absolute_rgb_error=float(difference.mean()),
                     expected_blue_sphere_pixels=int(visible.sum()),encoded_blue_sphere_pixel_recall=recall))
(out/f'{case}_ENCODED_FRAME_CHECK.json').write_text(json.dumps(dict(status='PASS',frames=rows),indent=2)+'\n')
print(case,rows)
