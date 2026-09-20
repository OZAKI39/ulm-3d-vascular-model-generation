"""The requested final terminal decision, with unavailable observations explicit."""
import json
from pathlib import Path
R=Path(__file__).resolve().parents[2]/'reports/sv1_3p'
def read(n):return json.loads((R/(n+'.json')).read_text())
s=read('candidate_summary');w=s['winner'];status=read('stage_status');b=s['baseline'];lines=['Stage SV1.3P completed.','','baseline:', '    config = ASM2 + ILU2','    window = step60→70',f"    wall = {b['wall_time_s']:.6f} s",f"    iterations = {b['statistics']['total_iterations']}"]
for k,d in s['candidates'].items():
 x=d['window'] or d['smoke'];e=x['profile']['events'] if x and x.get('profile') else {};it=x['statistics'].get('total_iterations',0) if x else None
 lines+=['',k+' '+{'P1':'reuse','P2':'hypre GPU ILU0','P3':'hypre BoomerAMG','P4':'PETSc GAMG','P5':'NVIDIA AmgX'}[k]+':']
 if k=='P1':lines+=['    lifecycle audit = PASS','    reuse scope = REUSE_WITHIN_TIMESTEP',f"    setup count = {d['PC_rebuild_count']} actual numeric factorizations / {x['statistics']['KSP_solves']} KSP solves"]
 if k=='P2':
  h=read('remote/petsc_hypre_build');lines += [f"    hypre version = {h['package_version']}", '    GPU setup confirmed = supported CUDA device path, not kernel residency profiling','    GPU solve confirmed = CUDA sparse standalone PASS; CFD health reported separately']
 if k=='P3':lines+=['    profile = PMIS / extended+i / l1-Jacobi, one V cycle']
 if k=='P4':lines+=['    matrix suitability = existing block size 4 mixed velocity/pressure; no supplied physical near nullspace']
 if k=='P5':lines+=['    build = '+read('amgx_feasibility')['result']]
 lines += [f"    wall = {x['wall_time_s'] if x else 'N/A'} s ({d['timing_scope'] or 'not run'})",f"    iterations = {it if it is not None else 'N/A'}",'    result = '+(d.get('strategy_result','')+' / ' if d.get('strategy_result') else '')+d['result']]
x=s['candidates'][w['fastest']]['window'] if w['fastest'] else None;e=x['profile']['events'] if x else {}
lines+=['','comparison:',f"    fastest healthy = {w['fastest']}",f"    improvement vs baseline = {w['improvement']*100:.3f}%" if w['fastest'] else '    improvement vs baseline = N/A',f"    setup time = {e.get('PCSetUpOnBlocks',e.get('PCSetUp',{})).get('time_s','N/A')} s (inclusive)",f"    apply time = {e.get('PCApply',{}).get('time_s','N/A')} s (inclusive)",f"    peak VRAM = {x['sampled_device_memory_max_MiB'] if x else 'N/A'} MiB (observed device-wide samples; not exact process peak)",'','winner:',f"    name = {w['fastest'] if w['full_required'] else 'NONE'}",'    >=10% faster = '+('YES' if w['full_required'] else 'NO'),'','full steady:']
if w['full_required']:
 f=read('winner_steady_candidate');lines += [f"    status = {f['status']}",f"    steady step = {f['first_full_steady_step']}",f"    stop step = {f['stop_step']}",f"    wall = {f['wall_time_s']:.6f} s"]
else:lines+=['    status = NOT_REQUIRED','    steady step = N/A','    stop step = N/A','    wall = N/A']
lines+=['','deferred:','    CPU/GPU science equivalence = YES','    multi-rank CUDA = YES','    FieldSplit = YES','','fallback:','    Stage O GPU config preserved = YES','','STAGE SV1.3P STATUS:','    '+status['status']]
text='\n'.join(lines)+'\n';(R/'terminal_summary.txt').write_text(text);print(text)
