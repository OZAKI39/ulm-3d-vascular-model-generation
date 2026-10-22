"""Full-window numerical admission and action distinguishability; no target relabelling."""
import numpy as np,h5py
from .common import *
from .analyze import numerical

def numerical_admission():
 refs=read(BASE/'discovery/INDEX.json');others=read(BASE/'refinement/INDEX.json');checks=[]
 for b in others:
  a=next((r for r in refs if r['state_id']==b.get('state_id') and r.get('action')==b.get('action')),None)
  if a is None:continue
  if a['status']!='PASS' or b['status']!='PASS':r=dict(reference=a['run_id'],comparison=b['run_id'],original=dict(admitted=False,reason='INCOMPLETE_MATCHED_WINDOW'),shape_uncertainty={})
  else:r=numerical(a,b)
  r.update(state_id=a['state_id'],action=a['action'],variant=b.get('variant'));checks.append(r)
  write(BASE/'analysis/NUMERICAL_COMPARISONS.json',checks)
 admissions=[]
 for a in refs:
  cc=[r for r in checks if r['reference']==a['run_id']]
  if not cc:continue
  exact=len(cc)==2 and all(r['original']['admitted'] for r in cc)
  noise={k:max([r['shape_uncertainty'].get(k) or 0 for r in cc],default=0) for k in ['similarity_rmsd_3d_m','deviatoric_stretch_3d','shear_norm_3d','directional_xy_stretch_3d']}
  m=a['final_shape'];resolved=exact and m['similarity_rmsd_3d_m']>max(1e-9,5*noise['similarity_rmsd_3d_m'])
  aff=bool(m['affine_valid_3d']);anis=resolved and aff and m['deviatoric_stretch_3d']>max(1e-4,5*noise['deviatoric_stretch_3d'])
  admissions.append(dict(run_id=a['run_id'],state_id=a['state_id'],action=a['action'],numerical_pass=exact,
    NUMERICALLY_RESOLVED=resolved,affine_anisotropy_resolved=anis,noise=noise,affine_valid=aff,
    topology_agrees=all(r['original'].get('event_topology_agrees',False) for r in cc),
    classification='NUMERICALLY_RESOLVED' if resolved else 'SHAPE_RESPONSE_BELOW_NUMERICAL_RESOLUTION',
    scope='Fixed-ID full window. Non-affine cases retain covariance/distance residuals; no invalid polar axes.'))
 write(BASE/'analysis/ADMISSION.json',admissions)
 print([(r['state_id'],r['action'],r['numerical_pass'],r['NUMERICALLY_RESOLVED'],r['affine_anisotropy_resolved']) for r in admissions])

def confirmation():
 refs=read(BASE/'discovery/INDEX.json');repeated=read(BASE/'confirmation/INDEX.json');checks=[]
 for b in repeated:
  if b.get('variant') not in ['DT05','FINE']:continue
  a=next(r for r in repeated if r['state_id']==b['state_id'] and r['action']==b['action'] and r['variant']=='DT1')
  check=numerical(a,b);check.update(state_id=a['state_id'],action=a['action'],variant=b['variant']);checks.append(check)
 write(BASE/'analysis/CONFIRMATION_NUMERICAL.json',checks)
 rows=[r for r in refs+repeated if r['status']=='PASS' and r['variant']=='DT1'];groups=[]
 for shape in ['S2','S4','S6']:
  for action in ['A5_X180','A7_X210']:
   rr=sorted([r for r in rows if r['state_id'].startswith(shape+'_') and r['action']==action],key=lambda r:r['state_id']);values=[]
   for r in rr:
    a=r['initial_shape'];m=r['final_shape'];values.append(dict(run_id=r['run_id'],state_id=r['state_id'],aspect_fraction=m['aspect_3d']/a['aspect_3d']-1,
     signed_xy=m['directional_xy_stretch_3d'],similarity_nm=m['similarity_rmsd_3d_m']*1e9,
     com_nm=(np.array(m['center_3d_m'])-a['center_3d_m'])*1e9,components=m['proximity_components'],largest=m['largest_component']))
   signs=[np.sign(v['aspect_fraction']) for v in values];signed=[v['signed_xy'] for v in values]
   groups.append(dict(shape=shape,action=action,values=values,replicate_count=len(rr),
    aspect_direction_repeated=len(rr)==3 and len(set(signs))==1,
    signed_xy_direction_repeated=len(rr)==3 and all(v is not None for v in signed) and len(set(np.sign(signed)))==1,
    all_connected=all(v['components']==1 for v in values),
    full_seed_numerical_confirmation=shape=='S2' and action=='A5_X180' and len(checks)==4 and all(c['original']['admitted'] for c in checks),
    scope='r0 discovery; r1/r2 held-out from action selection. Frames are not independent experiments.'))
 write(BASE/'analysis/REPLICATION.json',dict(groups=groups,identity=identity()))
 print([(g['shape'],g['action'],g['aspect_direction_repeated'],g['all_connected']) for g in groups])

if __name__=='__main__':numerical_admission()
