"""Actual contiguous selected clip, all backends independently computed."""
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from render_videos import vv,text
from xcommon import *
import av,cv2
from PIL import Image,ImageDraw
def main():
 start=time.perf_counter();d=json.loads((ROOT/'results/diagnostics_v1/backend_clip.json').read_text());rows={(r['backend'],r['frame_id']):r for r in d['rows']};writer=vv.Writer(ROOT/'videos/Video_F10X_3_backend_and_accuracy_comparison.mp4');mapping=[]
 for fid,pts,g in frames(range(360,420)):
  im,dr=vv.canvas('F10X | CPU / GPU measured displacement comparison','SELECTED CLIP: source frames 360-419 | same historical ROI and local enhancement | no physical truth')
  for i,n in enumerate(['B0','B1','B2','B3','B4']):
   r=rows[n,fid];px=25+(i%3)*520;py=118+(i//3)*430;roi=r['roi'];x,y,w,h=roi;rgb=cv2.cvtColor(g[y:y+h,x:x+w],cv2.COLOR_GRAY2RGB);good=np.array(r['valid']);u=np.array(r['u'],float);v=np.array(r['v'],float);xx=np.array(r['x']);yy=np.array(r['y'])
   for xx0,yy0,uu,vv0 in zip(xx[good],yy[good],u[good],v[good]):cv2.arrowedLine(rgb,(round(xx0),round(yy0)),(round(xx0+6*uu),round(yy0+6*vv0)),(30,185,100),1,tipLength=.25)
   im.paste(Image.fromarray(rgb).resize((235,352)),(px,py+30));dr=ImageDraw.Draw(im);text(dr,(px,py),['B0 Old OpenPIV','B1 SciPy','B2 Rust','B3 TorchPIV','B4 CuPy'][i],22,bold=True);text(dr,(px+248,py+55),f'{r["total_ms"]:.2f} ms',21);text(dr,(px+248,py+92),f'Valid: {int(good.sum())}',19)
   ref=rows['B1',fid];common=good&np.array(ref['valid']);delta=np.hypot(u-np.array(ref['u'],float),v-np.array(ref['v'],float));err=float(np.mean(delta[common])) if common.any() else None;text(dr,(px+248,py+139),'vs SciPy (shared)',17);text(dr,(px+248,py+170),'null' if err is None else f'{err:.4g} px',20)
   e=r['estimate'];text(dr,(px+248,py+214),'Six: accepted' if all(e['valid'].values()) else 'Six: rejected',18,color='#8fceb1');text(dr,(px+248,py+254),'Arrows x6',17);text(dr,(px+248,py+288),'Measured time',16,color='#b3bed1');text(dr,(px+248,py+314),'not video FPS',16,color='#b3bed1')
  text(dr,(1065,568),f'Frame {fid} | source t={pts:.4f} s',21);text(dr,(1065,615),'48-px window / 24-px overlap',18);text(dr,(1065,653),'Interval = 66.6 ms (gap 2)',18);text(dr,(1065,707),'Algorithms have different support.',17);text(dr,(1065,741),'Agreement is not physical accuracy.',17);text(dr,(1065,796),'Rust / SciPy: same correlation',17);text(dr,(1065,828),'CuPy: padded FFT, different result',17)
  outpts=pts-360*.0333;mapping.append(dict(frame=fid,source_pts_s=pts,output_pts_s=outpts));writer.add(im,outpts)
 record=vv.validate(writer.close());dump(ROOT/'videos/backend_video_manifest.json',dict(identity=identity(),video=record,source_sha256=sha(VIDEO),input_sha256=sha(ROOT/'results/diagnostics_v1/backend_clip.json'),mapping=mapping,wall_s=time.perf_counter()-start,scope='60 contiguous actual source frames, selected clip; no invented GPU result',visual_review='PENDING'));print('backend video verified',flush=True)
if __name__=='__main__':main()
