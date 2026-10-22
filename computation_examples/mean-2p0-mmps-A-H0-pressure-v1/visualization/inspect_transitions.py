"""Extract the largest annotation crossfade from each final MP4 for inspection."""
from pathlib import Path
import json
import numpy as np
import imageio_ffmpeg
from PIL import Image, ImageDraw, ImageFont, ImageOps

out=Path(__file__).resolve().parent/'results'
reports={}
for key,movie,kind in [('full_vessel','streamlines_full_vessel_enlarged.mp4','port'),('branch_vectors','velocity_vectors_branch_detail_enlarged.mp4','axis_title')]:
    trace=json.loads((out/f'{key}_annotations.json').read_text())['frames']
    choices=[]
    for frame in trace:
        rows=[r for r in frame['labels'] if r['kind']==kind and r['opacity']>.08]
        for j,a in enumerate(rows):
            for b in rows[j+1:]:
                if a['id'].split(':')[0]==b['id'].split(':')[0]:
                    score=a['opacity']*b['opacity']*np.linalg.norm(np.array(a['center_px'])-b['center_px'])
                    choices.append((score,frame['frame'],a['text']))
    _,center,name=max(choices)
    indexes=[(center+d)%432 for d in [-12,-7,-3,0,5,12]]
    boxes=[]
    for i in indexes:
        for r in trace[i]['labels']:
            if r['text']==name and r['opacity']>.001:
                boxes.append(r['bbox_px'])
                if 'arrow_tip_px' in r:
                    x,y=r['arrow_tip_px'];boxes.append([x,y,x,y])
    boxes=np.asarray(boxes)
    crop=(max(0,int(boxes[:,0].min()-60)),max(0,int(boxes[:,1].min()-60)),
          min(1728,int(boxes[:,2].max()+60)),min(1000,int(boxes[:,3].max()+60)))
    selected={}
    reader=imageio_ffmpeg.read_frames(str(out/'animations'/movie));next(reader)
    for i,raw in enumerate(reader):
        if i in indexes:
            selected[i]=Image.fromarray(np.frombuffer(raw,np.uint8).reshape(1080,1920,3).copy())
    sheet=Image.new('RGB',(1920,840),(8,10,14))
    draw=ImageDraw.Draw(sheet);font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',20)
    for k,i in enumerate(indexes):
        tile=ImageOps.contain(selected[i].crop(crop),(620,365))
        x,y=(k%3)*640,(k//3)*420
        sheet.paste(tile,(x+(640-tile.width)//2,y+48+(365-tile.height)//2))
        draw.text((x+16,y+12),f'{name} | frame {i} | {i/24:.3f} s',fill=(230,237,245),font=font)
    target=out/'inspection'/f'{key}_text_transition.png';sheet.save(target)
    reports[key]=dict(file=str(target.relative_to(out)),text=name,frames=indexes,source_crop_px=crop,
                      reason='Largest visible crossfade between two locations of one label')
(out/'TEXT_TRANSITION_INSPECTION.json').write_text(json.dumps(reports,indent=2)+'\n')
print(json.dumps(reports,indent=2))
