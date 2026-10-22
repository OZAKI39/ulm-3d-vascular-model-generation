"""Render frozen human ROI evaluation and actual causal replay, without new inference."""
from pathlib import Path
import sys, json, hashlib
import numpy as np, cv2, av
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from runtime import identity
from annotation_store import sha
from render_delivery import overlay as basic_overlay, encode, text_at, rgb

TAG = 'supervised_20261016_v1'
OUT = ROOT/'results'/TAG
FIG = ROOT/'figures'/TAG
VID = ROOT/'videos'/TAG
FREEZE = ROOT/'results/expert_freezes/vision2_supervised_20261016_v1'
COLORS = ['#147eab','#7550a3','#b66431']
NAMES = ['OpenCV','Cellpose pretrained','Cellpose tuned (rejected)']
MODELS = ['opencv','cellpose_pretrained','cellpose_tuned']
plt.rcParams.update({'font.size':10,'axes.titlesize':12,'axes.labelsize':10,
                    'pdf.fonttype':42,'savefig.dpi':240,'axes.spines.top':False,'axes.spines.right':False})


def save(fig,name):
    for ext in ['png','pdf']:
        fig.savefig(FIG/f'{name}.{ext}',facecolor='white',bbox_inches='tight')
    plt.close(fig)


def overlay(frame,msg,scale=1):
    # Dark outline keeps software IDs readable against the pale real image.
    image=basic_overlay(frame,{**msg,'targets':[]},scale)
    origin=np.asarray(msg['roi_xywh'][:2])
    def point(xy):return np.round((np.asarray(xy)-origin)*scale).astype(np.int32)
    for target in msg['targets']:
        if target['xy_px'] is None:continue
        pts=point([h[1:] for h in target['short_history']])
        if len(pts)>1:cv2.polylines(image,[pts.reshape(-1,1,2)],False,(65,130,35),1,cv2.LINE_AA)
        pos=tuple(point(target['xy_px'])+[4,-4]);label=target['id'].split(':')[-1]
        cv2.putText(image,label,pos,cv2.FONT_HERSHEY_SIMPLEX,.30,(30,30,30),2,cv2.LINE_AA)
        cv2.putText(image,label,pos,cv2.FONT_HERSHEY_SIMPLEX,.30,(245,245,245),1,cv2.LINE_AA)
    return image


def contours(ax,labels,ids,color,lw=.8):
    for i in ids:
        curves,_=cv2.findContours(np.uint8(labels==i),cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
        for curve in curves:
            points=curve[:,0]
            if len(points)>1:
                points=np.vstack([points,points[:1]])
                ax.plot(points[:,0],points[:,1],color=color,lw=lw)


def comparison(ax,image,record,arrays):
    ax.imshow(rgb(image));m=record['metrics']
    contours(ax,arrays['labels'],[r['pred'] for r in m['matches']],'#05a2c2')
    contours(ax,arrays['labels'],m['false_positive_ids'],'#d64938')
    for i in m['false_negative_ids']:
        yy,xx=np.where(arrays['truth']==i)
        ax.plot(xx.mean(),yy.mean(),'+',color='#e6a800',ms=9,mew=1.2)
    ax.set_xticks([]);ax.set_yticks([])


def figures():
    summary=json.loads((OUT/'summary.json').read_text())
    replay=json.loads((OUT/'replay/summary.json').read_text())
    paced=json.loads((OUT/'paced_replay/summary.json').read_text())
    rows=[json.loads(x) for x in (OUT/'replay/records.jsonl').read_text().splitlines()]
    frozen=json.loads((FREEZE/'frozen.json').read_text())
    human=next(r for r in frozen['rows'] if r['record']['id']=='video_000360')
    full=cv2.imread(str(FREEZE/human['frozen_image_file']))
    specs=[('opencv','ROI'),('cellpose_pretrained','ROI'),('cellpose_tuned','ROI'),('cellpose_pretrained','FULL_FRAME')]
    records=[]
    for model,domain in specs:
        p=OUT/'predictions'/f'{model}_{domain}_video_000360_roi00'
        records.append((json.loads(p.with_suffix('.json').read_text()),np.load(p.with_suffix('.npz'))))
    y0,x0,y1,x1=records[0][0]['box_yxyx'];crop=full[y0:y1,x0:x1]
    fig,axes=plt.subplots(2,3,figsize=(13,7.7),layout='constrained')
    axes[0,0].imshow(rgb(crop));axes[0,0].set_title('(a) Raw validation ROI | frame 360')
    axes[0,1].imshow(rgb(crop));truth=records[0][1]['truth']
    contours(axes[0,1],truth,np.unique(truth)[1:],'#9543b3',1)
    axes[0,1].set_title('(b) Human masks | 21 eligible / 23 drawn')
    titles=['(c) OpenCV','(d) Cellpose pretrained','(e) Fine-tuned candidate: rejected','(f) Full-frame inference: failure']
    for ax,(record,arrays),title in zip([axes[0,2],*axes[1]],records,titles):
        comparison(ax,crop,record,arrays);m=record['metrics']['center']
        ax.set_title(f'{title}\nTP {m["tp"]} / FP {m["fp"]} / FN {m["fn"]}')
    for ax in axes.flat:
        ax.set_xticks([]);ax.set_yticks([])
    axes[0,0].plot([8,58],[145,145],c='white',lw=3);axes[0,0].text(33,138,'50 px',ha='center',color='white')
    fig.legend(handles=[Line2D([],[],c='#9543b3',label='Human mask'),Line2D([],[],c='#05a2c2',label='Center-matched prediction'),
                        Line2D([],[],c='#d64938',label='Unmatched prediction (FP)'),Line2D([],[],c='#e6a800',marker='+',ls='',label='Missed human object (FN)')],
               loc='outside lower center',ncol=4,frameon=False)
    fig.suptitle('Human-reviewed ROI benchmark | first reviewed validation ROI\n6 px center matching; border/clipped labels excluded; unreviewed pixels are not ground truth',fontsize=13)
    save(fig,'Fig_VisionV2_1_detection_and_labels')

    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for ax,route,title in zip(axes[0],['opencv','cellpose'],NAMES):
        row=next(r for r in rows if r['frame']==390 and r['route']==route)
        image=cv2.imread(str(OUT/'replay/images'/row['image']))
        ax.imshow(rgb(overlay(image,row['message'],scale=1)))
        ax.set_title(f'{title} | frame 390, PTS {row["source_pts_s"]:.3f} s')
        ax.set_xticks([]);ax.set_yticks([])
    for route,color,name in zip(['opencv','cellpose'],COLORS,NAMES):
        a=[r for r in rows if r['route']==route];t=np.array([r['source_pts_s'] for r in a]);t-=t[0]
        axes[1,0].plot(t,[len(r['message']['candidates']) for r in a],c=color,label=name)
        axes[1,1].plot(t,[sum(z['xy_px'] is None for z in r['message']['targets']) for r in a],c=color,label=name)
    axes[1,0].set(xlabel='Recorded time since frame 360 (s)',ylabel='Image candidates per frame',title='(c) Candidate counts are not verified bubble counts')
    axes[1,1].set(xlabel='Recorded time since frame 360 (s)',ylabel='Unobserved software tracks per frame',title='(d) Missing observations: no old coordinates sent')
    for ax in axes[1]:ax.legend(fontsize=9);ax.grid(alpha=.15)
    fig.suptitle('CAUSAL_REPLAY | unchanged Norfair short association\nNo human track IDs: identity, turn handling and ID-break accuracy NOT_EVALUATED',fontsize=13)
    save(fig,'Fig_VisionV2_2_online_tracking')

    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    m=[summary['metrics'][f'{model}/validation/ROI'] for model in MODELS]
    x=np.arange(3)
    for k,(key,label) in enumerate([('precision','Precision'),('recall','Recall'),('f1','Center F1')]):
        bars=axes[0,0].bar(x+(k-1)*.24,[a['center'][key] for a in m],.24,label=label,color=['#2a8ab0','#7fb6b0','#67539c'][k])
        axes[0,0].bar_label(bars,fmt='%.2f',fontsize=8,padding=2)
    axes[0,0].set(xticks=x,xticklabels=['OpenCV','Pretrained','Tuned\n(rejected)'],ylim=(0,1.19),ylabel='Center detection score',title='(a) 129 eligible GT objects in 5 reviewed ROIs')
    axes[0,0].legend(ncol=3,fontsize=8,loc='upper center')
    axes[0,1].bar(x-.17,[a['center']['fp'] for a in m],.34,color='#d07869',label='False positives')
    axes[0,1].bar(x+.17,[a['center']['fn'] for a in m],.34,color='#dbc368',label='False negatives')
    axes[0,1].set(xticks=x,xticklabels=['OpenCV','Pretrained','Tuned'],ylabel='Count',title='(b) Fewer predictions do not mean better detection');axes[0,1].legend(fontsize=9)
    bars=axes[1,0].bar(x,[a['mask']['f1'] for a in m],color=COLORS)
    axes[1,0].bar_label(bars,fmt='%.3f',padding=3)
    axes[1,0].set(xticks=x,xticklabels=['OpenCV','Pretrained','Tuned'],ylim=(0,.65),ylabel='Mask F1 at IoU >= 0.50',title='(c) Boundary agreement differs from center agreement')
    for model,name,color in zip(MODELS,NAMES,COLORS):
        vals=[summary['metrics'][f'{model}/validation/{d}']['center']['recall'] for d in ['ROI','FULL_FRAME']]
        axes[1,1].plot([0,1],vals,'o-',label=name,c=color)
    axes[1,1].set(xticks=[0,1],xticklabels=['ROI inference','Full-frame inference'],ylim=(-.04,1.05),ylabel='Recall within reviewed ROIs',title='(d) Full-frame Cellpose inference fails here')
    axes[1,1].legend(fontsize=8)
    fig.suptitle('Pilot validation | 2 source frames, 133 drawn / 129 eligible single bubbles\nSame video; validation selected the model; held-out test, groups and tracking remain unassessed',fontsize=13)
    save(fig,'Fig_VisionV2_3_accuracy')

    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    stages=['detect_ms','track_ms','morphology_ms','flow_ms','serialize_send_ms','software_ms']
    xx=np.arange(len(stages))
    for route,name,color,delta in zip(['opencv','cellpose'],NAMES,COLORS,[-.18,.18]):
        values=[paced[route]['stage_ms'][k]['p95'] for k in stages]
        axes[0,0].bar(xx+delta,values,.36,label=name,color=color)
    axes[0,0].axhline(1000/60,color='#b13d3d',ls='--',lw=1,label='60 Hz budget (16.7 ms)')
    axes[0,0].set(yscale='log',xticks=xx,xticklabels=['Detect','Track','Shape','Flow','Send','Total'],ylabel='p95 time (ms; log scale)',title='(a) Paced replay | 180 x 350 px, native pixels')
    axes[0,0].legend(fontsize=8)
    for route,name,color in zip(['opencv','cellpose'],NAMES,COLORS):
        vals=[paced[route]['processed_hz'],paced[route]['fresh_candidate_hz'],0]
        offset=-.18 if route=='opencv' else .18
        bars=axes[0,1].bar(np.arange(3)+offset,vals,.36,label=name,color=color)
        axes[0,1].bar_label(bars,fmt='%.1f',fontsize=9,padding=3)
    axes[0,1].set(xticks=range(3),xticklabels=['Processed','Fresh candidate','Validated bubble'],ylabel='States per wall-clock second',ylim=(0,36),title='(b) Recorded input ~30 fps; this is not live 60 Hz')
    axes[0,1].legend(fontsize=8)
    values=[replay['routes'][r]['timing_ms']['engine_ms']['p95'] for r in ['opencv','cellpose']]
    bars=axes[1,0].bar(['OpenCV','Cellpose pretrained'],values,color=COLORS[:2]);axes[1,0].bar_label(bars,fmt='%.1f ms',padding=3)
    axes[1,0].set(ylabel='Engine p95 (ms)',ylim=(0,1250),title='(c) Larger 550 x 560 px ROI | 60 ordered frames')
    axes[1,1].axis('off')
    axes[1,1].text(0,.98,'Current camera run: NO LIVE FRAMES',va='top',weight='bold',fontsize=12)
    axes[1,1].text(0,.84,'Camera disconnected (user confirmed).\n\nLIVE_CAMERA_PASS = false\nLive update rate / latency: NOT_MEASURED\nEnd-to-end exposure latency: NOT_MEASURED\n\nGroup shape: null | control permission: false\nNo LED test; no AFG or PZT command.\n\nGPU timings include the existing annotation workload.',va='top',linespacing=1.5)
    fig.suptitle('Measured offline software performance and admission limits\nLatest-slot replay includes detection, tracking, morphology and mock-only messaging',fontsize=13)
    save(fig,'Fig_VisionV2_4_live_latency')


def human_overlay(image,index,freeze):
    """Only draw frozen, reviewed masks at their actual recorded frame."""
    result=image.copy();x,y,_,_=json.loads((OUT/'replay/summary.json').read_text())['roi_xywh']
    row=next((r for r in freeze['rows'] if r['record']['source_frame_index']==index),None)
    if row is None:return result,False
    # Load evaluated ROI truth: excludes all unknown parts of the frame.
    for p in sorted((OUT/'predictions').glob(f'opencv_ROI_video_{index:06d}_roi*.json')):
        meta=json.loads(p.read_text());y0,x0,y1,x1=meta['box_yxyx']
        truth=np.load(p.with_suffix('.npz'))['truth']
        for lab in np.unique(truth)[1:]:
            cc,_=cv2.findContours(np.uint8(truth==lab),cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_SIMPLE)
            for c in cc:
                c=c+np.array([x0-x,y0-y]);cv2.drawContours(result,[c],-1,(205,80,190),1,cv2.LINE_AA)
        cv2.rectangle(result,(x0-x,y0-y),(x1-x-1,y1-y-1),(205,80,190),1)
    return result,True


def videos():
    rows=[json.loads(x) for x in (OUT/'replay/records.jsonl').read_text().splitlines()]
    byframe={}
    for row in rows:byframe.setdefault(row['frame'],{})[row['route']]=row
    freeze=json.loads((FREEZE/'frozen.json').read_text())
    def frames(comparison=False):
        for index,pair in sorted(byframe.items()):
            a=pair['opencv'];raw=cv2.imread(str(OUT/'replay/images'/a['image']))
            canvas=np.full((960,1280,3),(23,18,15),np.uint8)
            text_at(canvas,'CAUSAL_REPLAY | recorded microscope video',(36,48),.9)
            text_at(canvas,f'Source frame {index} | recorded PTS {a["source_pts_s"]:.3f} s | no frame interpolation',(36,89),.65)
            if comparison:
                for route,left,title in [('opencv',50,'OpenCV + Norfair'),('cellpose',680,'Cellpose pretrained + Norfair')]:
                    r=pair[route];im=overlay(raw,r['message'],scale=1);im,hasgt=human_overlay(im,index,freeze)
                    canvas[185:745,left:left+550]=im
                    text_at(canvas,title,(left,142),.67)
                    text_at(canvas,f'Processing {r["timing"]["engine_ms"]:.1f} ms | candidates {len(r["message"]["candidates"])}',(left,171),.5)
                text_at(canvas,'Magenta: human masks in reviewed regions on this frame' if hasgt else 'This frame has no frozen human reference; candidate overlays are UNREVIEWED',(36,805),.65)
                text_at(canvas,'Matched source frames and crop; these are offline results, not live processing speed.',(36,850),.62)
            else:
                im=overlay(raw,a['message'],scale=1.2)
                canvas[150:822,40:700]=im
                text_at(canvas,'OpenCV + Norfair',(750,190),.8)
                text_at(canvas,f'Engine: {a["timing"]["engine_ms"]:.1f} ms',(750,245),.67)
                text_at(canvas,f'Candidates: {len(a["message"]["candidates"])}',(750,290),.67)
                text_at(canvas,'Cyan: image candidates',(750,370),.61)
                text_at(canvas,'Green: short observed paths',(750,410),.57)
                text_at(canvas,'IDs are unverified.',(750,480),.65)
                text_at(canvas,'No old locations for lost tracks.',(750,520),.5)
                text_at(canvas,'Group shape: NOT_ADMITTED',(750,590),.57)
                text_at(canvas,'Control: DISABLED',(750,630),.65)
                text_at(canvas,'Recorded-time playback does not demonstrate real-time throughput.',(36,867),.7)
            text_at(canvas,'Pixel coordinates only | no exposure-time or physical-size calibration',(36,920),.66)
            yield a['source_pts_s'],canvas
    outputs={}
    for stem,compare in [('Video_VisionV2_CAUSAL_REPLAY_tracking',False),('Video_VisionV2_CAUSAL_REPLAY_method_comparison',True)]:
        path=VID/f'{stem}.mp4'
        outputs[path.name]=encode(path,frames(compare),rate=30)
    (OUT/'video_manifest.json').write_text(json.dumps(outputs,indent=2))


def main():
    FIG.mkdir(exist_ok=False);VID.mkdir(exist_ok=False)
    figures();videos()
    inputs=[OUT/'summary.json',OUT/'all_rows.json',OUT/'replay/summary.json',OUT/'replay/records.jsonl',
            OUT/'paced_replay/summary.json',FREEZE/'frozen.json']
    (OUT/'render_identity.json').write_text(json.dumps(dict(execution=identity(ROOT),
        inputs={str(p.relative_to(ROOT)):sha(p) for p in inputs},
        selection='Example: first reviewed validation ROI; replay panel: midpoint frame 390; all metrics retained',
        output_files={str(p.relative_to(ROOT)):sha(p) for directory in [FIG,VID] for p in sorted(directory.iterdir())}),indent=2))
    print('Rendered 4 PNG/PDF figures and 2 causal-replay MP4 videos.')


if __name__=='__main__':main()
