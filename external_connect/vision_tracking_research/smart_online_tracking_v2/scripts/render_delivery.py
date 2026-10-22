"""Actual microscope images, logged detections and timing. No synthetic trajectories."""
from pathlib import Path
from fractions import Fraction
import sys,json,hashlib
import cv2,numpy as np,av
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from runtime import identity
FIG=ROOT/'figures';VID=ROOT/'videos'
plt.rcParams.update({'font.size':11,'axes.titlesize':13,'axes.labelsize':11,'figure.dpi':150,'savefig.dpi':220,
                     'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})

def loadlines(p):return [json.loads(s) for s in p.read_text().splitlines() if s]
def rgb(a):return cv2.cvtColor(a,cv2.COLOR_BGR2RGB)
def text_at(canvas,txt,xy,size=.6,color=(220,225,235)):
    cv2.putText(canvas,txt,tuple(map(int,xy)),cv2.FONT_HERSHEY_SIMPLEX,size,color,1,cv2.LINE_AA)

def overlay(frame,msg,scale=2):
    result=cv2.resize(frame,None,fx=scale,fy=scale,interpolation=cv2.INTER_NEAREST)
    origin=np.array(msg['roi_xywh'][:2]);
    def p(xy):return np.round((np.asarray(xy)-origin)*scale).astype(np.int32)
    for sh in msg.get('unreviewed_region_shapes',[]):
        color=(75,175,250) if sh['candidate_class']=='AGGREGATE_CANDIDATE' else (240,205,70)
        pts=p(sh['contour_xy_px'])
        if len(pts)>1:cv2.polylines(result,[pts.reshape(-1,1,2)],True,color,1,cv2.LINE_AA)
    for obj in msg.get('candidates',[]):
        if 'xy_px' in obj:cv2.circle(result,tuple(p(obj['xy_px'])),3,(240,205,70),1,cv2.LINE_AA)
    for t in msg.get('targets',[]):
        if t['xy_px'] is None:continue
        hist=t['short_history'];pts=p([q[1:] for q in hist])
        if len(pts)>1:cv2.polylines(result,[pts.reshape(-1,1,2)],False,(170,240,120),1,cv2.LINE_AA)
        q=p(t['xy_px']);text_at(result,t['id'].split(':')[-1],q+[4,-4],.35)
    return result

def static_msg(out,frame,route):
    meta=json.loads((out/f'static_{frame:06d}_ROI_{route}.json').read_text());labs=np.load(out/f'static_{frame:06d}_ROI_{route}_labels.npy')
    from vision_engine import region_shape
    x,y,_,_=meta['roi_xywh'];shapes=[];objs=[]
    for o in meta['objects']:
        sh=region_shape(labs==o['label'],(x,y));cl='SINGLE_BUBBLE_CANDIDATE' if o['area']<=450 else 'AGGREGATE_CANDIDATE'
        shapes.append({**sh,'candidate_class':cl});objs.append({**o,'xy_px':(np.array(o['xy_local'])+[x,y]).tolist()})
    return {'roi_xywh':meta['roi_xywh'],'candidates':objs,'targets':[],'unreviewed_region_shapes':shapes}

def savefig(fig,name):
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight',facecolor='white');fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(fig)

def encode(path,frames,rate=60):
    first_time=None;prev=-1;count=0;times=[]
    with av.open(str(path),'w') as container:
        stream=container.add_stream('libx264',rate=rate,options={'crf':'21','preset':'veryfast','bf':'0'})
        stream.width=1280;stream.height=960;stream.pix_fmt='yuv420p';stream.time_base=Fraction(1,90000)
        stream.codec_context.time_base=Fraction(1,90000)
        for timestamp,bgr in frames:
            if first_time is None:first_time=timestamp
            pts=round((timestamp-first_time)*90000)
            if pts<=prev:raise ValueError('Non-monotonic source time; do not fabricate timestamps')
            prev=pts;f=av.VideoFrame.from_ndarray(bgr,format='bgr24');f.pts=pts;f.time_base=Fraction(1,90000)
            for packet in stream.encode(f):container.mux(packet)
            count+=1;times.append(timestamp)
        for packet in stream.encode():container.mux(packet)
    return {'frames_encoded':count,'source_first_time_s':times[0],'source_last_time_s':times[-1],
            'no_interpolation':True,'time_base':'1/90000','sha256':hashlib.sha256(path.read_bytes()).hexdigest()}

def main():
    FIG.mkdir(exist_ok=True);VID.mkdir(exist_ok=True)
    pair=ROOT/'results/paired_comparison';paired=loadlines(pair/'records.jsonl')
    replay=[r for r in paired if r['kind']=='CAUSAL_REPLAY_OFFLINE'];bysrc={}
    for r in replay:bysrc.setdefault(r['frame'],{})[r['route']]=r
    live=ROOT/'results/live_opencv_20s';frames=loadlines(live/'recorded/frames.jsonl')
    A=json.loads((live/'summary.json').read_text());B=json.loads((ROOT/'results/live_cellpose_20s_resumed/summary.json').read_text())
    audit=json.loads((ROOT/'results/annotation_audit_current.json').read_text())
    figure_frame=780;full=cv2.imread(str(pair/f'images/static_{figure_frame:06d}.png'));x,y,w,h=[650,350,180,350];crop=full[y:y+h,x:x+w]
    fig,axes=plt.subplots(2,3,figsize=(14,9),layout='constrained')
    axes[0,0].imshow(rgb(full));axes[0,0].add_patch(Rectangle((x,y),w,h,fill=False,edgecolor='#21b4cf',linewidth=2));axes[0,0].set_title('(a) Recorded Chronos image')
    axes[0,1].imshow(rgb(crop));axes[0,1].set_title('(b) Original ROI pixels')
    axes[0,2].axis('off');axes[0,2].text(.04,.9,'Expert reference',fontsize=15,weight='bold');axes[0,2].text(.04,.74,f'{audit["expert_confirmed_frames"]} confirmed frames\n\nEXPERT_GROUND_TRUTH_PENDING\n\nNo ground-truth masks shown.\nNo FP / FN claims assigned.',va='top',linespacing=1.6)
    confirmed=[r for r in audit['frames'] if r['expert_confirmed']]
    if confirmed:
        reviewed=confirmed[0]
        records=json.loads((ROOT/'data/candidate_manifest.json').read_text())['rows']
        record=next(r for r in records if r['id']==reviewed['id'])
        annotated=cv2.imread(str(ROOT/record['image']))
        with np.load(ROOT/'data/annotations'/reviewed['id']/reviewed['mask_file'],allow_pickle=False) as masks:
            for cls,color in [('SINGLE_BUBBLE',(230,180,30)),('AGGREGATE',(55,170,250)),('BACKGROUND_ARTIFACT',(210,120,180)),('UNCERTAIN',(80,235,240)),('OUT_OF_FOCUS',(180,180,180))]:
                for lab in np.unique(masks[cls]):
                    if not lab:continue
                    contours,_=cv2.findContours(np.uint8(masks[cls]==lab),cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
                    cv2.drawContours(annotated,contours,-1,color,2)
        axes[0,2].clear();axes[0,2].imshow(rgb(annotated));axes[0,2].axis('off')
        axes[0,2].set_title(f'(c) Human-confirmed {reviewed["split"]} frame\n{reviewed["id"]}: {reviewed["counts"]["SINGLE_BUBBLE"]} singles, {reviewed["counts"]["AGGREGATE"]} groups')
    for ax,route in zip(axes[1,:2],['opencv','cellpose']):
        msg=static_msg(pair,figure_frame,route);ax.imshow(rgb(overlay(crop,msg)))
        ax.set_title(f'{route.title()} candidates: {len(msg["candidates"])}\nUNREVIEWED')
    live_sample=cv2.imread(str(ROOT/'data/candidates/live_20261016_000.png'))
    axes[1,2].imshow(rgb(live_sample));axes[1,2].set_title('(f) Separate live capture\nScene / focus not expert classified')
    for ax in axes.flat:
        ax.set_xticks([]);ax.set_yticks([])
    fig.suptitle('Vision-2 | Real images and two unreviewed detection routes',fontsize=17)
    savefig(fig,'Fig_VisionV2_1_detection_and_labels')

    selected=max(frames,key=lambda r:len(r['message']['targets']))
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    image=cv2.imread(str(live/'recorded'/selected['image']));axes[0,0].imshow(rgb(overlay(image,selected['message'])));axes[0,0].axis('off')
    axes[0,0].set_title('(a) LIVE_CAMERA | short candidate tracks\nIDs do not establish physical identity')
    log=loadlines(live/'frames.jsonl');t=np.array([r['source_time_s'] for r in log]);t-=t[0]
    axes[0,1].plot(t,[r['observed_targets'] for r in log],label='Observed software tracks',lw=1,color='#316aa0')
    axes[0,1].plot(t,[r['lost_targets'] for r in log],label='Lost software tracks',lw=.8,color='#dd8c32')
    axes[0,1].set(xlabel='Host decode time since first frame (s)',ylabel='Count',title='(b) Short tracker state');axes[0,1].legend(fontsize=9)
    ex=bysrc[800]['opencv'];im=cv2.imread(str(pair/'images/replay_000800.png'));axes[1,0].imshow(rgb(overlay(im,ex['message'])));axes[1,0].axis('off')
    axes[1,0].set_title('(c) Recorded-video causal replay\nOnly past and current frames used')
    axes[1,1].axis('off');axes[1,1].text(.02,.95,'Not yet admitted',fontsize=15,weight='bold',va='top')
    axes[1,1].text(.02,.82,'True single-bubble identities: NOT_EVALUATED\nAssociation errors / turns: NOT_EVALUATED\nVerified group masks: unavailable\nGroup neighborhood: null\nIncorporation events: NOT_EVALUATED\n\nLost targets carry no old coordinates.\nNo fixed-flow or F10S identity prior.\nNo claims of physical fusion.',va='top',linespacing=1.7)
    fig.suptitle('Vision-2 | Observed short tracks and explicit uncertainty',fontsize=17);savefig(fig,'Fig_VisionV2_2_online_tracking')

    fig,axes=plt.subplots(1,2,figsize=(14,7),gridspec_kw={'width_ratios':[1.65,1]},layout='constrained')
    metrics=['Single-bubble precision / recall / F1','False positives / negatives per image','Single-bubble center error','Group detection / Dice / IoU','Group area / Rg / tensor / aspect error','Short association accuracy / ID breaks','Sudden turns / focus loss / entry / exit']
    axes[0].axis('off');table=axes[0].table(cellText=[[m,'NOT_EVALUATED','NOT_EVALUATED'] for m in metrics],colLabels=['Metric','OpenCV','Cellpose'],colWidths=[.56,.22,.22],cellLoc='left',loc='center')
    table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,2.7)
    for (rr,cc),cell in table.get_celld().items():
        cell.set_edgecolor('#d4d7dd');cell.set_facecolor('#e9eef5' if rr==0 else '#f7f8fa')
    axes[0].set_title('(a) Accuracy cannot be inferred from model output')
    axes[1].bar(['Candidate\nimages','Draft\nframes','Confirmed\nframes'],[audit['candidates'],audit['draft_frames'],audit['expert_confirmed_frames']],color=['#98a7bf','#dfb16c','#62a3ad'])
    axes[1].set(ylabel='Frame count',ylim=(0,100),title='(b) Human annotation status')
    for pos,value in enumerate([audit['candidates'],audit['draft_frames'],audit['expert_confirmed_frames']]):
        axes[1].text(pos,value+2,str(value),ha='center',fontsize=12)
    axes[1].text(.98,.93,'No model predictions\nused as truth.\n\nGroup shape not admitted.\n\nCross-session generalization\nnot verified.',transform=axes[1].transAxes,va='top',ha='right',fontsize=10)
    fig.suptitle('Vision-2 | Accuracy admission is pending human labels',fontsize=17);savefig(fig,'Fig_VisionV2_3_accuracy')

    fig,axes=plt.subplots(2,2,figsize=(14,9),layout='constrained')
    live_logs=[log,loadlines(ROOT/'results/live_cellpose_20s_resumed/frames.jsonl')]
    axes[0,0].boxplot([[r['software_ms'] for r in q] for q in live_logs],tick_labels=['OpenCV','Cellpose'],showfliers=True)
    axes[0,0].axhline(1000/60,color='#c53d48',ls='--',label='16.7 ms reference');axes[0,0].set_yscale('log')
    axes[0,0].set(ylabel='Software processing (ms, log scale)',title='(a) All processed live frames');axes[0,0].legend(fontsize=9)
    keys=['detect_ms','track_ms','morphology_ms','flow_ms','serialize_send_ms'];xx=np.arange(len(keys))
    for i,(d,col) in enumerate([(A,'#357bac'),(B,'#c98e4b')]):axes[0,1].bar(xx+(i-.5)*.35,[d['stage_ms'][k]['p95'] for k in keys],width=.35,label=d['route'],color=col)
    axes[0,1].set_xticks(xx,['Detection','Tracking','Morphology','Sparse LK','Send'],rotation=18)
    axes[0,1].set_yscale('log');axes[0,1].set(ylabel='Component p95 (ms)',title='(b) Measured processing components');axes[0,1].legend()
    kinds=['Host input','Processed','Fresh candidate','Validated bubble'];xx=np.arange(4)
    for i,(d,col) in enumerate([(A,'#357bac'),(B,'#c98e4b')]):
        axes[1,0].bar(xx+(i-.5)*.35,[d[k] for k in ['input_host_decode_hz','processed_hz','fresh_candidate_hz','validated_bubble_hz']],width=.35,label=d['route'],color=col)
    axes[1,0].set_xticks(xx,kinds,rotation=12);axes[1,0].set(ylabel='Updates / measured wall second',title='(c) Throughput and semantic admission');axes[1,0].legend()
    axes[1,1].axis('off');table_rows=[]
    for name,d in [('OpenCV',A),('Cellpose',B)]:
        table_rows.append([name,f'{d["stage_ms"]["software_ms"]["p95"]:.2f}',str(d['unprocessed_source_frames_in_span']),str(d['stale_result_rejected']),f'{d["resources"]["rss_mb"]["maximum"]:.0f}',f'{d["resources"]["cpu_pct_one_core_basis"]["p95"]:.0f}',f'{d["resources"]["gpu_device_util_pct"]["p95"]:.0f}'])
    tab=axes[1,1].table(cellText=table_rows,colLabels=['Route','p95 ms','Skipped','Stale','RAM MB','CPU %','GPU %'],cellLoc='center',bbox=[0,.60,1,.30]);tab.auto_set_font_size(False);tab.set_fontsize(9)
    axes[1,1].text(.01,.52,'CPU: process p95 (one-core basis); GPU: whole-device p95.\nSequential LIVE_CAMERA trials; input scenes may differ.\nPaired saved-video tests are reported separately.\n\nEXPOSURE_TIMESTAMP_UNKNOWN\nCAMERA_TO_CONTROLLER_END_TO_END_LATENCY\n= NOT_MEASURED\n\nHardware control disabled; no LED test.',va='top',fontsize=9,linespacing=1.5)
    fig.suptitle('Vision-2 | Actual camera processing, not end-to-end latency',fontsize=17);savefig(fig,'Fig_VisionV2_4_live_latency')
    if '--figures-only' in sys.argv:
        (ROOT/'evidence/figure_render_identity.json').write_text(json.dumps(identity(ROOT),indent=2))
        return

    def live_video():
        for r in frames:
            canvas=np.full((960,1280,3),(20,16,13),np.uint8);img=cv2.imread(str(live/'recorded'/r['image']))
            over=overlay(img,r['message']);canvas[105:105+over.shape[0],55:55+over.shape[1]]=over
            text_at(canvas,'LIVE_CAMERA | OpenCV + Norfair',(45,40),.9)
            text_at(canvas,'UNREVIEWED CANDIDATES | HARDWARE CONTROL DISABLED',(45,76),.65,(110,210,250))
            x=630;yy=155
            lines=[f'Source frame: {r["frame_id"]}',f'Host decode time: {r["source_relative_s"]:.3f} s',
                   f'Software processing: {r["timing"]["software_ms"]:.2f} ms',f'Whole-run processed: {A["processed_hz"]:.2f} Hz',
                   f'Fresh candidates: {A["fresh_candidate_hz"]:.2f} Hz','Validated bubble output: 0 Hz','',
                   f'Current candidates: {len(r["message"]["candidates"])}','Cyan: small image regions','Orange: large image regions','Green: last 2-5 observed points',
                   'IDs are short software associations.','','group_shape = null','GROUP_SHAPE_NOT_ADMITTED','',
                   'Original ROI: x496 y568, 256 x 384','Display enlarged 2x; inference 1x',
                   'Source decode timing retained (VFR).','No interpolated frames.','','Exposure time unknown.','Camera end-to-end latency: NOT_MEASURED']
            for line in lines:text_at(canvas,line,(x,yy),.53);yy+=29
            text_at(canvas,'Source: Chronos Windows receiver -> WSL vision -> mock only',(45,925),.6)
            yield r['source_time_s'],canvas
    def comparison_video():
        for frame,methods in sorted(bysrc.items()):
            raw=cv2.imread(str(pair/f'images/replay_{frame:06d}.png'));canvas=np.full((960,1280,3),(20,16,13),np.uint8)
            text_at(canvas,'CAUSAL_REPLAY_OFFLINE | same recorded frame and ROI',(40,35),.82)
            text_at(canvas,'UNREVIEWED | No expert reference yet | Control disabled',(40,70),.65,(110,210,250))
            for j,route in enumerate(['opencv','cellpose']):
                r=methods[route];img=overlay(raw,r['message']);left=140+j*620
                canvas[140:840,left:left+360]=img
                text_at(canvas,route.upper()+' + NORFAIR',(left-40,115),.75)
                text_at(canvas,f'Candidates: {r["count"]} | engine: {r["engine_ms"]:.1f} ms',(left-40,875),.59)
            r=methods['opencv'];text_at(canvas,f'Video frame {frame} | source PTS {r["source_time_s"]:.4f} s | 1x encoded source time',(45,917),.59)
            text_at(canvas,'Sequential causal processing; playback speed does not imply real-time inference.',(45,948),.55)
            yield r['source_time_s'],canvas
    reports={}
    reports['Video_VisionV2_live_tracking.mp4']=encode(VID/'Video_VisionV2_live_tracking.mp4',live_video())
    reports['Video_VisionV2_method_comparison.mp4']=encode(VID/'Video_VisionV2_method_comparison.mp4',comparison_video(),30)
    for name,report in reports.items():
        with av.open(str(VID/name)) as c:
            st=c.streams.video[0];decoded=list(c.decode(st));pts=[float(f.time) for f in decoded]
            assert len(decoded)==report['frames_encoded'] and all(b>a for a,b in zip(pts,pts[1:]))
            report.update(decoded_frames=len(decoded),codec=st.codec_context.name,resolution_wh=[st.width,st.height],
                          average_rate=str(st.average_rate),duration_s=float(st.duration*st.time_base),monotonic_pts=True)
            for i in [0,len(decoded)//2,len(decoded)-1]:
                cv2.imwrite(str(VID/f'{Path(name).stem}_check_{i:04d}.png'),decoded[i].to_ndarray(format='bgr24'))
    (ROOT/'evidence/video_validation.json').write_text(json.dumps(reports,indent=2))
    (ROOT/'evidence/render_identity.json').write_text(json.dumps(identity(ROOT),indent=2))
    print(json.dumps(reports,indent=2))
if __name__=='__main__':main()
