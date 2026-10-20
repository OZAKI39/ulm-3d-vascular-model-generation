"""One real stage image to review video rendering while the H protocol runs."""
from common import *
import render_videos as rv
from PIL import Image, ImageDraw, ImageFont

ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'BIPHASIC_FIELD.npz')
evidence=json.loads((O/'TIME_2P550_MS_NEAREST_CASE_COMPARISON_READBACK.json').read_text())
states={m:dict(path=O/evidence[m]['source_path'],state=load_npz(O/evidence[m]['source_path'])) for m in CASES}
registries={}
for mode in CASES:
    s=states[mode]['state'];lab,sz,att=rv.graph.topology(ref,s['d']);registries[mode]={}
    for label in np.flatnonzero(~att):
        ids=tuple((np.flatnonzero(lab==label)+1).tolist());seed=min(ids)
        registries[mode][ids]=(seed,seed)
points=np.concatenate([ref['X'],states['uniform']['state']['x'],states['biphasic']['state']['x'],load_npz(O/'cases/uniform/keys/final.npz')['x']])*1e6
lo=points.min(0);hi=points.max(0);pad=np.maximum((hi-lo)*.06,25);lo-=pad;hi+=pad
bounds=np.array([lo[0],hi[0],lo[1],hi[1],lo[2],hi[2]]);center=(lo+hi)/2;span=max(hi-lo)
camera=[center+np.array([1.4,-2.9,1.0])*span,center,[0,0,1]]
rv.F=F/'previews/stage_video_2p55';rv.STATUS={'uniform':'COMPLETED','biphasic':'RUNNING'}
original_overlay=rv.overlay
def stage_overlay(*args,**kwargs):
    im=Image.fromarray(original_overlay(*args,**kwargs));d=ImageDraw.Draw(im)
    font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
    d.rectangle((0,0,1920,55),fill='#fcfcfc')
    d.text((40,14),'STAGE PREVIEW · Biphasic damage and cohesive graph components',font=ImageFont.truetype(font,32),fill='#243d4e')
    d.rectangle((0,1025,1920,1057),fill='#fcfcfc')
    d.text((40,1029),'Gray: fixed anchors. Tan: 1-2 node fines. Other colors: current true components. Full lineage animation is pending.',font=ImageFont.truetype(font,21),fill='#536675')
    return np.array(im)
rv.overlay=stage_overlay
row=dict(requested_time_s=.00255,**states)
rv.render(2,[row],ref,field,bounds,camera,registries,only_frame=True,smoke=True)
dump(O/'STAGE_VIDEO_PREVIEW_EVIDENCE.json',dict(scope='Single actual stage image for rendering review, not final video or final smoke acceptance.',sources={m:{k:evidence[m][k] for k in ['source_path','source_sha256','time_s']} for m in CASES},displacement_scale=1,bounds_um=bounds,camera=camera,color_rule='Components >=3 nodes use their minimum original atom ID as color seed; 1–2 node fines use muted tan. No lineage claim in this preview.',parallel_scale_um=rv.fitted_parallel_scale(bounds,camera),final_video_visual_review_satisfied=False))
print('STAGE_VIDEO_PREVIEW_DONE',flush=True)
