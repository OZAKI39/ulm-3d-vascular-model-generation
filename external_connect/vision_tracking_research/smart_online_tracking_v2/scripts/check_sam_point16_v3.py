"""Bounded follow-up to the retained poor default32 result, TRAIN prompts only."""
from pathlib import Path
import sys,json,contextlib,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from sam_worker import SamWorker
BASE=ROOT/'evidence/annotation_repair_v3_20261016/sam_comparison';OUT=BASE/'point16_followup';OUT.mkdir(exist_ok=False)
plan=json.loads((BASE/'plan.json').read_text());(OUT/'plan.json').write_text(json.dumps({'reference_plan':str(BASE/'plan.json'),'target_width_px':16,'reason':'32 px point-only boxes often include broad background; eight fixed training instances only; explicit user box remains preferred'},indent=2))
p=ROOT/'evidence/annotation_repair_v3_20261016/labels_backup/video_000000/latest.json';meta=json.loads(p.read_text())
with np.load(p.parent/meta['mask_file']) as z:a=z['SINGLE_BUBBLE']
worker=SamWorker();rows=[]
for instance in plan['instance_ids']:
 req=json.loads((BASE/f'{instance:02d}_local_point_default32_request.json').read_text());req.update(target_width_px=16,token=f'point16_{instance}',output=str(OUT/f'{instance:02d}.npz'))
 (OUT/f'{instance:02d}_request.json').write_text(json.dumps(req,indent=2))
 with (OUT/f'{instance:02d}.log').open('w') as f,contextlib.redirect_stdout(f):worker.run(req)
 with np.load(req['output']) as z:m=z['mask'];d=json.loads(str(z['details']))
 t=a==instance;row={'instance':instance,'allowed':d['proposals'][0]['allowed'],'area_px':int(m.sum()),'iou_to_training_reference':float((m&t).sum()/max((m|t).sum(),1)),'wall_s':d['wall_s'],'cold_model':d['cold_model']};rows.append(row);print(row,flush=True)
summary={'scope':'TRAIN_RETROSPECTIVE_NOT_GENERALIZATION','rows':rows,'allowed':sum(v['allowed'] for v in rows),'median_iou':float(np.median([r['iou_to_training_reference'] for r in rows])),'warm_median_s':float(np.median([r['wall_s'] for r in rows if not r['cold_model']]))}
(OUT/'summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
