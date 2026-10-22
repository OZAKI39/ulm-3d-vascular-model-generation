"""Prepare an isolated existing-workbench snapshot from the supplied microscope clip.

No labels, reviewed channel mask or detection performance are manufactured.
Overlay geometry was visually checked before choosing the source-pixel crop.
"""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess

import av
import cv2
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
VIDEO=ROOT.parent.parent/'sample_videos/experiment_overlay 截取视频.mp4'
OUT=ROOT/'microscope_sessions/microscope_20261019_v1'
SCHEDULE={'train':[5,15,25,35,45,55,65,75,85,105,115,125],
          'validation':[140,145,150,155,160,165], 'test':[180,183,186,189,192,195]}
CROP=[0,40,1280,984]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def dump(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2))


def main():
    if OUT.exists():raise RuntimeError('Session already exists; refusing to replace data or labels')
    source_sha=sha(VIDEO)
    for d in ['data/candidates','data/annotations','data/auto_training','src','scripts','results','evidence','logs']:
        (OUT/d).mkdir(parents=True,exist_ok=False)
    # Reuse the current tested annotator and its installed environments/weights.
    copied=[]
    for pattern in ['annotation_*.py','sam_*.py']:
        for src in (ROOT/'src').glob(pattern):
            shutil.copy2(src,OUT/'src'/src.name);copied.append(src.relative_to(ROOT).as_posix())
    for name in ['annotate.py','sam_worker.py','preannotation_worker.py','sam_auto_train.py','sam_auto_service.py']:
        shutil.copy2(ROOT/'scripts'/name,OUT/'scripts'/name);copied.append('scripts/'+name)
    shutil.copy2(ROOT/'start_annotation.sh',OUT/'start_annotation.sh')
    for name in ['.venv','.venv_annotation','.venv_preannotation','models']:
        (OUT/name).symlink_to((ROOT/name).resolve(),target_is_directory=True)
    protocol=dict(scope='NEW_MICROSCOPE_ANNOTATION_ONLY',source=str(VIDEO),source_sha256=source_sha,
        selection_times_s=SCHEDULE,split_policy='Fixed separated time blocks; single recording only; no cross-session validation claim',
        source_crop_xywh=CROP,camera_area_in_source_xywh=[0,0,1280,1024],
        excluded='Right PZT panel and top timestamp strip; no inpainting, resizing or optical enhancement',
        channel_mask_status='NOT_DRAWN_NOT_HUMAN_VERIFIED',annotation_status='NO_HUMAN_LABELS',
        model_training_status='NOT_RUN_NO_CONFIRMED_LABELS',realtime_model_connected=False,
        coordinate_transform='source_x=annotation_x; source_y=annotation_y+40',
        provenance='Frames decoded from a compressed overlay export, not raw sensor data')
    dump(OUT/'protocol.json',protocol) # Freeze selection before extraction or model inference.
    rows=[];previews=[]
    with av.open(str(VIDEO)) as container:
        stream=container.streams.video[0]
        assert (stream.width,stream.height)==(1964,1024),'Unexpected overlay layout; inspect before cropping'
        for split,times in SCHEDULE.items():
            for target in times:
                container.seek(int(target/stream.time_base),stream=stream,backward=True)
                frame=next(f for f in container.decode(stream) if float(f.time)>=target)
                full=frame.to_ndarray(format='bgr24');x,y,w,h=CROP;image=full[y:y+h,x:x+w].copy()
                assert image.shape==(984,1280,3)
                name=f'microscope_20261019_{split}_{target:03d}s';path=OUT/'data/candidates'/f'{name}.png'
                assert cv2.imwrite(str(path),image,[cv2.IMWRITE_PNG_COMPRESSION,3])
                assert np.array_equal(cv2.imread(str(path)),image),'Extracted PNG pixels changed'
                rows.append(dict(id=name,image=path.relative_to(OUT).as_posix(),image_sha256=sha(path),
                    source_kind='CROPPED_OVERLAY_VIDEO',source_file=str(VIDEO),source_sha256=source_sha,
                    source_pts=frame.pts,source_time_base=str(frame.time_base),source_pts_s=float(frame.time),
                    requested_time_s=target,source_resolution_wh=[1964,1024],resolution_wh=[1280,984],
                    saved_roi_xywh=CROP,resized=False,split=split,time_block=f'new_microscope_{split}',
                    selection_kind='fixed_time_sample',sequence_id=None,review_status='UNREVIEWED',
                    camera_configuration=None,capture_fps=None,objective_magnification=None,adapter=None,
                    focal_plane=None,scene_class='HUMAN_REVIEW_REQUIRED',channel_mask_status='NOT_DRAWN',
                    proposed_roi_xywh=[600,300,256,256]))
                thumb=cv2.resize(image,(320,246));cv2.rectangle(thumb,(0,0),(320,24),(0,0,0),-1)
                cv2.putText(thumb,f'{len(rows):02d}  {split}  {float(frame.time):.1f}s',(6,17),cv2.FONT_HERSHEY_SIMPLEX,.44,(240,240,240),1,cv2.LINE_AA)
                previews.append(thumb)
    dump(OUT/'data/candidate_manifest.json',dict(version=1,session_id=OUT.name,rows=rows,expert_confirmed_frames=0))
    # Open a later training frame first, with fewer visible objects, without labels.
    dump(OUT/'data/annotations/session.json',dict(index=10,label_opacity_pct=18,bulk_diameter_filter='REFERENCE_ONLY'))
    dump(OUT/'data/auto_training/config.json',dict(enabled=True))
    cv2.imwrite(str(OUT/'evidence/frame_contact_sheet.jpg'),np.vstack([np.hstack(previews[i:i+4]) for i in range(0,24,4)]))
    dump(OUT/'EXECUTION.json',dict(base=str(ROOT),base_git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        reused_sources={p:sha(ROOT/p) for p in copied},prepare_script_sha256=sha(__file__),source_sha256=source_sha,
        shared_dependencies=['.venv','.venv_annotation','.venv_preannotation','models'],
        isolated=['data','results','evidence','logs'],frame_count=len(rows)))
    print(json.dumps(dict(workspace=str(OUT),frames=len(rows),split_counts={k:len(v) for k,v in SCHEDULE.items()},source_sha256=source_sha)))


if __name__=='__main__':main()
