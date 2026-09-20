"""Compact metrics and adaptive timestep history from actual artifacts only."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3q';rows=[];adaptive={}
for p in sorted(R.glob('*_acceptance.json')):
 d=json.loads(p.read_text());r=d.get('reuse') or {};events=(d.get('profile') or {}).get('events',{});reason=d.get('observed_KSP_reasons',[])
 counts=[x['iterations'] for x in reason if x.get('iterations') is not None]
 rows.append(dict(name=d['name'],candidate=d['candidate'],mode=d['mode'],status=d['status'],wall_time_s=d['wall_time_s'],timesteps=d['steps'],KSP_attempts=len(reason),logical_solves=d['statistics'].get('KSP_solves'),total_iterations=sum(counts),mean_iterations=sum(counts)/len(counts) if counts else None,max_iterations=max(counts) if counts else None,ILU_rebuild_count=r.get('ILU_rebuild_count'),ILU_reuse_count=r.get('ILU_reuse_count'),recovery_count=r.get('recovery_count'),MatLUFactorNum=events.get('MatLUFactorNum'),PCSetUp=events.get('PCSetUp'),PCSetUpOnBlocks=events.get('PCSetUpOnBlocks'),PCApply=events.get('PCApply'),KSPSolve=events.get('KSPSolve'),MatMult=events.get('MatMult'),observed_GPU_memory_MiB=d['sampled_device_memory_max_MiB'],CPU_fallback_warning_lines=d['CPU_fallback_warning_lines'],events_nested_not_additive=True))
 if d['candidate']=='RA':
  bystep={}
  for a in r.get('trace',[]):
   e=a['outcome'];bystep.setdefault(str(a['step']),[]).append(dict(logical_solve=a['logical'],attempt=a['attempt'],age_since_rebuild=a['age'],I_ref_before=a['ref'],I_ref_after=e['ref'],I_current=e['iterations'],trigger_ratio=e['ratio'],rebuild=not bool(a['reuse']),rebuild_reason=a['rebuild_reason'],KSP_reason=e['reason'],recovery=e['recovery'],rebuild_required_next=bool(e['pending'])))
  adaptive[d['name']]=bystep
(R/'runtime_metrics.json').write_text(json.dumps(dict(status='PASS',cases=rows),indent=2)+'\n')
(R/'adaptive_timestep_history.json').write_text(json.dumps(dict(status='PASS',rule='Strict >1.5 * first healthy fresh-solve iterations, or age>=5',cases=adaptive),indent=2)+'\n')
print('Summarized',len(rows),'native cases; RA trace retained by timestep.')
