"""Extract real full-resolution frames and register temporal groups before inference."""
from pathlib import Path
import hashlib,json,sys,datetime
import av,cv2,numpy as np
ROOT=Path(__file__).resolve().parents[1]
VIDEO=ROOT.parent.parent/'sample_videos/Recording 2026-10-15 185901 截取视频.mp4'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    out=ROOT/'data/candidates'
    manifest=ROOT/'data/candidate_manifest.json'
    if manifest.exists():raise RuntimeError('Candidate manifest already exists; refusing overwrite')
    source_sha=sha(VIDEO)
    selections={}
    blocks={'train':(0,270,120),'validation':(360,570,450),'test':(660,936,780)}
    for split,(start,end,seq_start) in blocks.items():
        seq=set(range(seq_start,seq_start+8))
        spaced=np.linspace(start,end,16,dtype=int)
        for i in spaced:
            while int(i) in seq:i+=1
            selections[int(i)]=(split,'independent_time_sample',None)
        for i in seq:selections[i]=(split,'short_sequence',f'video_{split}_sequence')
    rows=[]
    with av.open(str(VIDEO)) as c:
        stream=c.streams.video[0]
        stream.codec_context.thread_count=1
        for i,frame in enumerate(c.decode(stream)):
            if i not in selections:continue
            split,kind,sequence=selections[i]
            image=frame.to_ndarray(format='bgr24')
            name=f'video_{i:06d}'
            path=out/f'{name}.png'
            assert not path.exists()
            cv2.imwrite(str(path),image,[cv2.IMWRITE_PNG_COMPRESSION,3])
            rows.append({'id':name,'image':path.relative_to(ROOT).as_posix(),
                         'image_sha256':sha(path),'source_kind':'SAVED_VIDEO',
                         'source_file':str(VIDEO),'source_sha256':source_sha,'source_frame_index':i,
                         'source_pts':frame.pts,'source_time_base':str(frame.time_base),
                         'source_pts_s':float(frame.time),'resolution_wh':[image.shape[1],image.shape[0]],
                         'saved_roi_xywh':[0,0,image.shape[1],image.shape[0]],'resized':False,
                         'split':split,'time_block':f'video_{split}','selection_kind':kind,'sequence_id':sequence,
                         'review_status':'UNREVIEWED','camera_configuration':None,
                         'capture_fps':None,'objective_magnification':None,'adapter':None,'focal_plane':None,
                         'scene_class':'EXPERT_SCENE_CONFIRMATION_PENDING',
                         'proposed_roi_xywh':[640,80,620,960]})
    assert len(rows)==72 and len({r['id'] for r in rows})==72
    manifest.write_text(json.dumps({'version':1,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                   'rows':rows,'planned_live_candidates':12,'expert_confirmed_frames':0},indent=2))
    protocol={'status':'PREREGISTERED_PENDING_EXPERT_LABELS','source_sha256':source_sha,
              'split_ranges_inclusive':{s:[a,b] for s,(a,b,_) in blocks.items()},
              'sequence_length':8,'no_random_adjacent_split':True,'cross_session_generalization':'NOT_VERIFIED',
              'live_candidates':'QUARANTINE_UNASSIGNED; do not train or evaluate until provenance/review frozen',
              'expert_label_classes':['SINGLE_BUBBLE','AGGREGATE','BACKGROUND_ARTIFACT','UNCERTAIN','OUT_OF_FOCUS'],
              'default_scope':'FULL_FRAME; ROI-only reviews must be explicitly recorded and excluded from full-frame metrics',
              'evaluation':{'instance_match_iou':0.5,'single_center_match_radius_px':6,
                            'assignment':'scipy.optimize.linear_sum_assignment, one-to-one with explicit unmatched objects',
                            'ignore_classes':['UNCERTAIN','OUT_OF_FOCUS'],'ignore_prediction_overlap_fraction':0.5,
                            'group_metrics':['IoU','Dice','centroid_error_px','area_error_px2','Rg_error_px','Q_error_px2','aspect_ratio_error'],
                            'confidence_selection':'training/validation only; no test tuning',
                            'short_track_metrics':'Only on human-matched 8-frame sequences; physical identity not inferred from direction'},
              'inference_budget':{'cellpose_models':1,'model':'Cellpose 3 cyto3 lightweight alternative',
                                  'finetune':'NOT_RUN until expert training labels and fixed validation split exist',
                                  'first_model_test_max_frames':12,'causal_comparison_window_s':3,
                                  'live_run_s_each':20,'live_roi_xywh':'register after observing current camera view, before comparing models'},
              'no_sam2':True,'no_led_test':True,'hardware_control':False}
    (ROOT/'protocol.json').write_text(json.dumps(protocol,indent=2))
    print(json.dumps({'frames':len(rows),'splits':{s:sum(r['split']==s for r in rows) for s in blocks},'manifest':str(manifest)}))

if __name__=='__main__':main()
