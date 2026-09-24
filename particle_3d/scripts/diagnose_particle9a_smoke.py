"""Inspect recorded failures without new trajectory integrations or model changes."""
from pathlib import Path
import sys,json,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.particle81_simulation import environment,dump
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_gap import wall_gap
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block
from particle_3d.validation_boundary import ValidationBoundaryClassifier
OUT=ROOT/'particle_3d/outputs/particle9a_2mmps';REPORT=ROOT/'particle_3d/reports/particle9a_2mmps'
env=environment();cap=ValidationBoundaryClassifier({'INLET':env.boundaries['INLET']});rows=[]
for path in sorted((OUT/'trajectories').glob('mb_*.json')):
    if '.receipt.' in path.name:continue
    meta=json.loads(path.read_text());samples=np.load(OUT/meta['samples_path'])['samples']
    p=samples[0,1:4];a=meta['radius_m'];s=env.field.sample(p);gap=wall_gap(Sphere(p,a),env.wall)
    _,_,d=planar_wall_affine_block(a,env.mu,gap.gap_m,gap.normal_inward,s.velocity_gradient_s_inv,np.r_[s.velocity_m_s,.5*s.vorticity_s_inv])
    inlet=env.boundaries['INLET'];tri=np.asarray(inlet.points)[inlet.faces.reshape(-1,4)[:,1:]]
    ns=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);n=ns.sum(0);n/=np.linalg.norm(n)
    if n@s.velocity_m_s<0:n=-n
    speed=samples[1,4:7] if len(samples)>1 else np.full(3,np.nan)
    end=samples[-1,1:4];end_sample=env.field.sample(end)
    crossing=cap.first_event(samples[0,1:4],samples[1,1:4]) if len(samples)>1 else None
    rows.append(dict(particle_id=meta['particle_id'],end_reason=meta['end_reason'],failure_detail=meta['failure_detail'],
        first_gap_ratio=d['gap_ratio'],first_wall_weight=d['wall_weight'],
        first_bulk_speed_m_s=float(np.linalg.norm(s.velocity_m_s)),
        first_accepted_speed_m_s=float(np.linalg.norm(speed)),
        first_inward_inlet_velocity_m_s=float(speed@n),
        first_shear_magnitude_s_inv=float(np.linalg.norm(d['shear_vector_xyz'])),
        first_inlet_crossing_detected=bool(crossing is not None),last_center_inside=end_sample.inside_lumen,
        minimum_gap_m=meta['minimum_original_wall_gap_m'],last_elapsed_s=float(samples[-1,0]),
        last_64_accepted_progress_s=float(samples[-1,0]-samples[max(0,len(samples)-64),0]),
        minimum_handoff_margin_m=float(samples[:,15].min()),
        max_depth=meta['maximum_refinement_depth'],accepted_steps=meta['accepted_steps'],
        rejected_trials=meta['rejected_trials']))
dump(REPORT/'data/smoke_diagnosis.json',dict(scope='POSTHOC_RECORDED_SMOKE_ONLY_NO_NEW_INTEGRATION',rows=rows,
    conclusion='FAIL: two initial upstream cap escapes, four deterministic progress guards near continuum handoff, one residence horizon; no sampled WALL penetration.',
    scientific_boundary='Do not force velocity inward, move births, lower the guard, add curvature/lift, or select successful IDs to claim production pass.'))
with (REPORT/'data/smoke_diagnosis.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
print(json.dumps(rows,indent=2))
