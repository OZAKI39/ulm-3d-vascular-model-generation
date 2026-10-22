"""Same saved real frames and frozen estimates; selected clip, no interpolation."""
import sys,json,time,functools,importlib.util
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
import cv2
from PIL import Image,ImageDraw
spec=importlib.util.spec_from_file_location('video_helpers',ROOT.parent/'collective_motion_v2/scripts/videos.py');vv=importlib.util.module_from_spec(spec);spec.loader.exec_module(vv);vv.OUT=OUT;vv.font=functools.lru_cache(maxsize=32)(vv.font)
def main():
    record=json.loads((OUT/'results/real_comparison_clip.json').read_text());begin,end=record['protocol']['selected_real_clip'];frames=old.read_frames(range(begin,end),VIDEO);names=['PIV','LRAFT_ST','CoTracker3','POT_unbalanced'];labels=['Frozen PIV','LRAFT-ST','CoTracker3 online','Unbalanced OT'];rows={m:{r['frame_id']:r for r in record['methods'][m]} for m in names};x,y,w,h=P['roi_xywh'];writer=vv.Writer(OUT/'videos/Video_F10W_3_pretrained_motion_comparison.mp4');start=time.perf_counter()
    for fid,(bgr,pts) in sorted(frames.items()):
        gray=cv2.cvtColor(bgr[y:y+h,x:x+w],cv2.COLOR_BGR2GRAY);im,d=vv.canvas('Chronos microscopy | pretrained motion comparison','SELECTED CLIP  |  identical source frames  |  real six-parameter accuracy NOT_VERIFIED')
        for k,m in enumerate(names):
            row=rows[m][fid];e=row['published_estimate'];ok=all(e['valid'].values());xx=30+k*395;vv.text(d,(xx,115),labels[k],24,bold=True);scale,_,_=vv.microscopy(im,gray,row['objects'],x=xx,y=158,w=350,h=530)
            if ok:
                a=np.array(im)
                for px,py,dx,dy,good in row['vectors']:
                    if not good or any(v is None for v in [px,py,dx,dy]):continue
                    p=(round(xx+px*scale),round(158+py*scale));q=(round(p[0]+dx*scale*8),round(p[1]+dy*scale*8));cv2.arrowedLine(a,p,q,(99,225,165),1,tipLength=.2)
                im=Image.fromarray(a)
            d=ImageDraw.Draw(im);vv.text(d,(xx,710),'QUALITY ACCEPTED' if ok else 'NOT OBSERVABLE',20,color='#6bd7a2' if ok else '#ffac80')
            for j,n in enumerate(NAMES):
                val=e['values'][n];name=['vx','vy','div','vort','s1','s2'][j];vv.text(d,(xx+(j%2)*175,750+(j//2)*45),name+': '+('null' if val is None else f'{val:+.3f}'),20)
            vv.text(d,(xx,900),f'Adapter {row["timing"]["adapter_total_ms"]:.1f} ms',18,color='#aab8ca')
        vv.text(d,(30,950),f'Source frame {fid} | PTS {pts:.4f} s | vx,vy: px/s; others: 1/s | arrows ×8',21);writer.add(im,pts-frames[begin][1])
    output=vv.validate(writer.close());dump(OUT/'videos/comparison_video_manifest.json',dict(identity=identity(),video=output,wall_s=time.perf_counter()-start,source_mapping=[dict(output_index=i,source_frame=fid,source_pts_s=frames[fid][1]) for i,fid in enumerate(sorted(frames))],data_sha256=sha(OUT/'results/real_comparison_clip.json')))
if __name__=='__main__':main()
