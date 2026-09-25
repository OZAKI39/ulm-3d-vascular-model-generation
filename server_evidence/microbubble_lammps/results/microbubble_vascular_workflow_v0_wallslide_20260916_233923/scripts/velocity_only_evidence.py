"""Analytic concave curved-wall counterexample, before enabling position corrections."""
from pathlib import Path
import json,numpy as np
S=Path(__file__).resolve().parents[1];rows=[]
R=5e-6;a=1e-6;m=1e-10;r=R-a-m;V=1e-4
for shape in ['cylinder','sphere']:
 for dt in [1e-4,5e-5,2.5e-5]:
  x=np.array([r,0.,0.]);n=np.array([-1.,0.,0.]);raw=np.array([V,V,0.]);used=raw-min(raw@n,0)*n
  mid=x+dt/2*used
  violation=np.linalg.norm(mid[:2] if shape=='cylinder' else mid)-r
  rows.append(dict(shape=shape,dt=dt,midpoint_violation_m=float(violation),raw=raw.tolist(),used=used.tolist(),safe=bool(violation<=0)))
assert all(r['midpoint_violation_m']>0 for r in rows)
result={'status':'VELOCITY_ONLY_CURVED_GEOMETRY_TANGENTIAL_STALL_DEMONSTRATED','reason':'At exact concave contact every nonzero straight tangent midpoint leaves the offset cylinder/sphere. Violation is O(dt^2), so halving alone does not admit a strictly safe step.','rows':rows}
(S/'validation/VELOCITY_ONLY_CURVED_EVIDENCE.json').write_text(json.dumps(result,indent=2)+'\n')
p=S/'contracts/WALL_SLIDING_CONTRACT.json';c=json.loads(p.read_text());c.update(position_correction_policy='BOUNDED_GEOMETRIC_OFFSET_PROJECTION enabled only after at least one actual wall retry, only if this stage projected inward velocity, only endpoint deficit <= cap. Recheck corrected entire segment. Reject deeper proposals and halve; never accept deep penetration.',position_correction_cap_m=m/4,roundoff_cushion='64 * machine epsilon * max(norm(endpoint),radius+margin); included in cap',correction_target='gap = inherited margin + roundoff cushion',maximum_corrections_per_stage=1,evidence='validation/VELOCITY_ONLY_CURVED_EVIDENCE.json');p.write_text(json.dumps(c,indent=2)+'\n')
print(json.dumps(result,indent=2))
