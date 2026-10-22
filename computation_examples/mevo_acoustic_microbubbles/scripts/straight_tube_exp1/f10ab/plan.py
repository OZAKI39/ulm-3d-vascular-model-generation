"""Deterministic finite selection; hold-out seeds are never used for discovery."""
import argparse,numpy as np
from .common import *

def task(stage,state,action,duration=.002,variant='DT1',schedule=None):
 out=dict(stage=stage,run_id=f'{stage}__{state}__{action}__{round(duration*1e6)}us__{variant}',state_id=state,action=action,duration_s=duration,variant=variant)
 if schedule is not None:out['schedule']=schedule
 return out

def freeze(stage,tasks,budget,reason):
 p=BASE/'plans'/f'{stage}.json';assert not p.exists()
 write(p,dict(stage=stage,stage_budget_s=budget,workers=2,tasks=tasks,selection_reason=reason,identity=identity(),pre_execution=True))
 print(stage,len(tasks),reason)

def pilot():
 data=read(BASE/'static/SUMMARY.json')['records'];names=read(BASE/'ACTIONS.json');scores={}
 for action in names:
  rows=[r for r in data if r['action']==action]
  legal=all(r['status']=='PASS' and all(r.get('single_step_gates',{}).values()) for r in rows)
  if legal:scores[action]=dict(deviatoric=float(np.mean([r['free_deviatoric_rate_s1'] for r in rows])),signed=float(np.mean([r['free_signed_xy_rate_s1'] for r in rows])))
 xnames=[k for k in scores if k.startswith('A') and k not in ['A0_OFF','A8_YM','A9_YP','A10_Y90','A11_FOUR']]
 high=max(xnames,key=lambda k:scores[k]['deviatoric']);low=min(xnames,key=lambda k:scores[k]['signed'])
 chosen=list(dict.fromkeys(['A7_X210',high,low]))
 # Keep exact squared-amplitude counterpart for the strongest eligible pair.
 match={'A3_X0':'B3_X0_L2','A4_X90':'B4_X90_L2','A7_X210':'B7_X210_L2'}.get(high,'B7_X210_L2')
 chosen=list(dict.fromkeys(chosen+[match]))
 chosen=[a for a in chosen if a in scores]
 write(BASE/'ACTION_SELECTION.json',dict(scores=scores,selected=chosen,diagnostic_only=[a for a in names if a not in chosen],scope='Static eligibility; no dynamic control claim'))
 states=[f'S{k}_r0' for k in range(7)]
 freeze('pilot',[task('pilot',s,a,.0005) for s in states for a in ['A0_OFF']+chosen],450,'Frozen static ranking; all seven discovery states; no held-out seed inspected')

def discovery():
 pilot=read(BASE/'pilot/INDEX.json');ok={(r['state_id'],r.get('action')):r for r in pilot if r['status']=='PASS'}
 tasks=[]
 for k in range(7):
  for action in ['A0_OFF','A7_X210','A5_X180']:
   state=f'S{k}_r0'
   if (state,action) in ok:tasks.append(task('discovery',state,action))
 for state in ['S2_r0','S6_r0']:
  if (state,'B3_X0_L2') in ok:tasks.append(task('discovery',state,'B3_X0_L2'))
  # Single-source comparator is legal in static diagnostics, cost forecast
  # borrowed from the corresponding L2-normalized two-source pilot.
  tasks.append(task('discovery',state,'A1_XM'))
 freeze('discovery',tasks,700,'A7/A5 complete all pilots; same-state OFF; two L2-normalized pair and single-source comparisons. A3 failed dense states and is not extended. Full-window numerical confirmation remains required.')

def refinement():
 records=read(BASE/'discovery/INDEX.json');ok={(r['state_id'],r.get('action')) for r in records if r['status']=='PASS'}
 pairs=[('S0_r0','A7_X210')]+[(s,a) for s in ['S2_r0','S4_r0','S6_r0'] for a in ['A7_X210','A5_X180']]
 tasks=[task('refinement',s,a,.002,v) for s,a in pairs if (s,a) in ok for v in ['DT05','FINE']]
 freeze('refinement',tasks,600,'Full 2ms dt/2 and fine-cache controls, selected before inspecting numerical differences; original thresholds unchanged.')

def confirmation():
 admission=read(BASE/'analysis/ADMISSION.json');valid={(r['state_id'].split('_')[0],r['action']) for r in admission if r['numerical_pass'] and r['NUMERICALLY_RESOLVED']}
 tasks=[]
 for shape in ['S2','S4','S6']:
  actions=[a for a in ['A7_X210','A5_X180'] if (shape,a) in valid]
  if not actions:continue
  for rep in [1,2]:
   for action in ['A0_OFF']+actions:tasks.append(task('confirmation',f'{shape}_r{rep}',action))
 # Repeat the promising anisotropic branch with dt/cache on each new S2 realization.
 if ('S2','A5_X180') in valid:
  for rep in [1,2]:
   for variant in ['DT05','FINE']:tasks.append(task('confirmation',f'S2_r{rep}','A5_X180',.002,variant))
 freeze('confirmation',tasks,500,'Unseen seed realizations only; eligible full-window action/shape pairs from original numerical gate; S2 A5 held-out seeds receive independent dt/cache controls.')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage');a=ap.parse_args()
 if a.stage=='pilot':pilot()
 elif a.stage=='discovery':discovery()
 elif a.stage=='refinement':refinement()
 elif a.stage=='confirmation':confirmation()
if __name__=='__main__':main()
