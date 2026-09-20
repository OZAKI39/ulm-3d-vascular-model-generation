"""Interpret the two small diagnostic runs without treating nested times as additive."""
import json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3o'
d=json.loads((R/'profile_runs.json').read_text());summaries={}
for name in set((d['baseline'],d['winner'])):
 p=ROOT/'outputs/sv1_3o'/name/'petsc_profile.txt';text=p.read_text();events={};by_stage={};stages={};stage=None
 assert 'with 1 process ' in text, 'This summation is defined for the frozen single-rank path only'
 for line in text.splitlines():
  s=re.match(r'^\s*\d+:\s+(.+?):\s+([0-9.eE+-]+)\s+([0-9.]+)%',line)
  if s:stages[s[1]]=dict(time_s=float(s[2]),percent_total=float(s[3]),raw=line.strip())
  s=re.match(r'^--- Event Stage \d+: (.+)$',line)
  if s:stage=s[1];by_stage[stage]={};continue
  m=re.match(r'^([A-Za-z][A-Za-z0-9_: ]*?)\s{2,}(\d+)\s+([0-9.]+)\s+([0-9.eE+-]+)\s+([0-9.]+)\s+',line)
  if m and stage is not None:
   event=m[1].rstrip();row=dict(count=int(m[2]),time_s=float(m[4]),raw=line.strip())
   by_stage[stage][event]=row
   total=events.setdefault(event,dict(count=0,time_s=0.,stage_rows=[]))
   total['count']+=row['count'];total['time_s']+=row['time_s'];total['stage_rows'].append(dict(stage=stage,**row))
 assert 'KSPSolve' in events and 'PCApply' in events and 'PCSetUp' in events
 a=json.loads((R/(name+'_acceptance.json')).read_text())
 summaries[name]=dict(wall_time_s=a['wall_time_s'],events=events,events_by_stage=by_stage,stages=stages,event_aggregation='Same event summed over recorded stages for one MPI rank; distinct nested events must not be summed',profile_path=str(p.relative_to(ROOT)),factor_packages=a['factor_packages'],fallback_warning_lines=a['CPU_fallback_warning_lines'],FEM_assembly='NOT_SEPARATELY_MEASURED',VTU_restart_IO='NOT_SEPARATELY_MEASURED',GPU_residency='NOT_MEASURED',nested_events_are_not_additive=True)
out=dict(status='PASS',runs=d,profiles=summaries,interpretation='Inclusive PETSc event durations from diagnostic windows, never candidate ranking times. MatAssembly is PETSc matrix finalization, not total FEM assembly. Unattributed time includes application assembly, startup, I/O and instrumentation.',official_sources=['https://petsc.org/release/manualpages/PC/PCSetUpOnBlocks/','https://petsc.org/release/manual/profiling/','https://petsc.org/release/manualpages/PC/PCILU/','https://petsc.org/release/manualpages/KSP/KSPGMRESSetRestart/'])
(R/'profile_summary.json').write_text(json.dumps(out,indent=2)+'\n');print('Baseline/winner lightweight profiles parsed; nested timings not added together.')
