"""Two complete recordings; every source PTS retained, saved estimates only."""
import sys,time,json,importlib.util,functools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from wcommon import *
from full_records import full_records
import av,cv2,pyarrow.parquet as pq
from PIL import Image,ImageDraw
spec=importlib.util.spec_from_file_location('f10v_video_helpers',ROOT.parent/'collective_motion_v2/scripts/videos.py');vv=importlib.util.module_from_spec(spec);spec.loader.exec_module(vv)
vv.OUT=OUT;vv.font=functools.lru_cache(maxsize=32)(vv.font)
text=vv.text
LABELS=['vx (px/s)','vy (px/s)','divergence (1/s)','vorticity (1/s)','shear 1 (1/s)','shear 2 (1/s)']

def arrows(image,row,origin,scale,zoom=8):
    a=np.array(image);ox,oy=origin
    for x,y,dx,dy,valid in row['methods']['candidate_3'].get('vectors',[]):
        if any(v is None for v in [x,y,dx,dy]):continue
        p=(round(ox+x*scale),round(oy+y*scale))
        if valid:cv2.arrowedLine(a,p,(round(p[0]+dx*scale*zoom),round(p[1]+dy*scale*zoom)),(107,229,157),1,tipLength=.25)
        else:cv2.drawMarker(a,p,(244,124,119),cv2.MARKER_TILTED_CROSS,4,1)
    return Image.fromarray(a)

def full_panel(im,bgr,row,x=30,y=120,width=1080):
    scale=width/bgr.shape[1];height=round(bgr.shape[0]*scale);rgb=cv2.cvtColor(bgr,cv2.COLOR_BGR2RGB);im.paste(Image.fromarray(rgb).resize((width,height)),(x,y));d=ImageDraw.Draw(im)
    rx,ry,rw,rh=row['roi']
    # Image-only visual guides fixed from frame 0, not calibrated channel segmentation.
    for bx in [450,1330]:d.line((x+bx*scale,y,x+bx*scale,y+height),fill='#c2b58c',width=1)
    for o in row['full_frame_candidates']:
        xx,yy=o['xy'];cx=x+xx*scale;cy=y+yy*scale
        if not (rx<=xx<rx+rw and ry<=yy<ry+rh):d.ellipse((cx-1,cy-1,cx+1,cy+1),fill='#70b9cb')
    d.rectangle((x+rx*scale,y+ry*scale,x+(rx+rw)*scale,y+(ry+rh)*scale),outline='#39c9ff',width=3)
    for o in row['objects']:
        xx,yy=o['xy'];cx=x+(rx+xx)*scale;cy=y+(ry+yy)*scale;d.ellipse((cx-2,cy-2,cx+2,cy+2),outline='#39c9ff')
    center=row['distribution']['group_center_xy']
    if center:
        cx=x+(rx+center[0])*scale;cy=y+(ry+center[1])*scale;d.line((cx-9,cy,cx+9,cy),fill='#ffcf75',width=2);d.line((cx,cy-9,cx,cy+9),fill='#ffcf75',width=2)
    return arrows(im,row,(x+rx*scale,y+ry*scale),scale),height

def state_panel(d,row,x=1140,y=124):
    e=row['methods']['candidate_3']['estimate'];valid=all(e['valid'].values());text(d,(x,y),'Fixed historical ROI',24,bold=True)
    text(d,(x,y+39),f'Frame {row["frame_id"]} / 936',23);text(d,(x,y+72),f'PTS {row["pts_s"]:.4f} s',23)
    text(d,(x,y+119),'VALID GLOBAL ESTIMATE' if valid else row['processing_status'],20,color='#6bd7a2' if valid else '#ffac80')
    reason=next((r for r in e['reasons'].values() if r),'Quality gate only; not truth verified')
    reason=reason.replace('_',' ')
    import textwrap
    for k,line in enumerate(textwrap.wrap(reason,width=32)):text(d,(x,y+150+26*k),line,18,color='#aab8ca')
    for j,n in enumerate(NAMES):
        yy=y+240+j*66;value=e['values'][n];text(d,(x,yy),LABELS[j],19,color='#9bb0c9');text(d,(x+240,yy),'null' if value is None else f'{value:+.3f}',23)
    text(d,(x,y+664),f'Candidates: {len(row["objects"])}',21)
    text(d,(x,y+698),f'Valid vectors: {row["methods"]["candidate_3"].get("valid_vector_count",0)}',21)
    text(d,(x,y+733),'Arrows: displacement x8',18,color='#aab8ca')

def make_curves():
    rows=pq.read_table(OUT/'results/full_video_piv_per_frame.parquet').to_pylist();rows=[r for r in rows if r['method']=='candidate_3'];t=np.array([r['pts_s'] for r in rows]);vals=np.array([[r['values'][n] if r['values'][n] is not None else np.nan for n in NAMES] for r in rows]);lims=[]
    for j in range(6):
        finite=vals[:,j][np.isfinite(vals[:,j])];lo=float(finite.min()) if len(finite) else -1;hi=float(finite.max()) if len(finite) else 1;pad=max((hi-lo)*.1,.01);lims.append((lo-pad,hi+pad))
    return t,vals,lims
def curves(d,fid,t,vals,lims):
    for j in range(6):
        x=690;y=140+j*113;w=835;h=78;lo,hi=lims[j]
        text(d,(x,y-27),LABELS[j],18,color='#b7c8dc');d.line((x,y+h,x+w,y+h),fill='#536579');d.line((x,y,x,y+h),fill='#536579')
        text(d,(x-52,y-4),f'{hi:.2g}',14,color='#9bb0c9');text(d,(x-52,y+h-15),f'{lo:.2g}',14,color='#9bb0c9')
        points=[]
        for k in range(fid+1):
            if np.isfinite(vals[k,j]):points.append((x+w*t[k]/t[-1],y+h*(hi-vals[k,j])/(hi-lo)))
            else:
                if len(points)>1:d.line(points,fill='#5dc9ef',width=2)
                points=[]
        if len(points)>1:d.line(points,fill='#5dc9ef',width=2)
        xx=x+w*t[fid]/t[-1];d.line((xx,y,xx,y+h),fill='#ffcf75',width=1)
        if j==5:text(d,(x,y+h+5),'0',15);text(d,(x+w-65,y+h+5),f'{t[-1]:.2f} s',15)

def main():
    start=time.perf_counter();writers=[vv.Writer(OUT/'videos'/f'Video_F10W_{n}_{title}.mp4') for n,title in [(1,'full_duration_PIV'),(2,'full_duration_PIV_diagnostics')]];t,vals,lims=make_curves();count=0
    with av.open(str(VIDEO)) as c:
        for fid,(frame,row) in enumerate(zip(c.decode(video=0),full_records(),strict=True)):
            assert row['frame_id']==fid and abs(float(frame.time)-row['pts_s'])<1e-8;bgr=frame.to_ndarray(format='bgr24')
            im,d=vv.canvas('Full Chronos recording | frozen single-pass PIV','OFFLINE FULL-FRAME COVERAGE  |  causal 3-frame endpoints  |  REAL ACCURACY: NOT_VERIFIED')
            im,height=full_panel(im,bgr,row);d=ImageDraw.Draw(im);state_panel(d,row)
            text(d,(30,955),'Cyan: candidates / ROI   Gold: visible center / visual channel guide   Green: PIV   Red: rejected vectors',17,color='#b7c8dc');writers[0].add(im,row['pts_s'])
            im,d=vv.canvas('Full Chronos recording | motion diagnostics','ALL SOURCE FRAMES  |  displayed curves use saved past estimates only  |  no interpolated frames')
            im,height=full_panel(im,bgr,row,x=30,y=118,width=590);d=ImageDraw.Draw(im);text(d,(30,574),f'Frame {fid}  |  t = {row["pts_s"]:.4f} s  |  {row["processing_status"]}',19)
            x,y,w,h=row['roi'];gray=cv2.cvtColor(bgr[y:y+h,x:x+w],cv2.COLOR_BGR2GRAY);enh=preprocess(gray,'candidate',row['objects']);enh=np.uint8(np.clip(enh*5,0,255))
            scale,_,_=vv.microscopy(im,gray,row['objects'],x=30,y=635,w=230,h=310);im=arrows(im,row,(30,635),scale)
            vv.microscopy(im,enh,[],x=315,y=635,w=230,h=310);d=ImageDraw.Draw(im);text(d,(30,607),'ROI + local vectors',18);text(d,(315,607),'Candidate enhancement',18)
            curves(d,fid,t,vals,lims);e=row['methods']['candidate_3']['estimate'];reason=next((r for r in e['reasons'].values() if r),'QUALITY GATE PASSED')
            text(d,(690,874),f'Candidates {len(row["objects"])} | valid vectors {row["methods"]["candidate_3"].get("valid_vector_count",0)}',22)
            text(d,(690,912),reason.replace('_',' '),18,color='#ffcf75');text(d,(690,951),'Invalid parameters remain null; vector display scale x8.',18,color='#aab8ca');writers[1].add(im,row['pts_s']);count+=1
            if fid%200==0:print('rendered',fid,flush=True)
    audit=json.loads((OUT/'results/full_video_source_audit.json').read_text());assert count==audit['decoded_frames']
    records=[vv.validate(writer.close()) for writer in writers];dump(OUT/'videos/full_video_manifest.json',dict(identity=identity(),source_sha256=sha(VIDEO),videos=records,wall_s=time.perf_counter()-start,source_decoded_frames=count,display_channel_guide='frame-0 visual guide x=450 and 1330, not a calibrated channel mask',visual_review='PENDING'))
if __name__=='__main__':main()
