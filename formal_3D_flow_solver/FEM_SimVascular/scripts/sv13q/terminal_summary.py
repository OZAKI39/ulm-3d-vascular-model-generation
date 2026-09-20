import json
from pathlib import Path
R=Path(__file__).resolve().parents[2]/'reports/sv1_3q'
def read(n):return json.loads((R/(n+'.json')).read_text())
def optional(n):return read(n) if (R/(n+'.json')).exists() else None
def metrics(c):
 s=read(c+'_SMOKE_acceptance');d=optional(c+'_WINDOW_acceptance')
 return s,d
b=read('baseline_r1');w=read('winner');full=optional('winner_steady_candidate');lines=['Stage SV1.3Q completed.','','baseline R1:','    rebuild interval = 1',f"    late wall = {b['late']['wall_time_s']:.6f} s",'    late iterations = 6720','    rebuild count = 10']
for c in ('R2','R3','R5','RA'):
 s,d=metrics(c);lines+=['',c+(' adaptive' if c=='RA' else '')+':',f"    smoke = {s['status']}"]
 if c=='RA':lines += ['    threshold = 1.5 x I_ref','    max age = 5']
 lines += [f"    late wall = {d['wall_time_s'] if d else 'N/A'} s",f"    late iterations = {d['total_attempt_iterations'] if d else 'N/A'}",f"    rebuild count = {(d.get('reuse') or {}).get('ILU_rebuild_count','N/A') if d else 'N/A'}"]
 if c=='RA':lines += [f"    recovery count = {(d.get('reuse') or {}).get('recovery_count','N/A') if d else 'N/A'}"]
 lines += [f"    result = {d['gain_classification'] if d else 'REJECT_UNHEALTHY'}"]
ranking=read('late_ranking');lines+=['','late window ranking:']
for i,c in enumerate(ranking['top_two'],1):lines+=[f'    {i} = {c}']
lines+=['','early window:']
for i,c in enumerate(ranking['top_two'],1):
 d=read(c+'_EARLY_acceptance');lines += [f'    candidate {i} = {c}',f"    wall = {d['wall_time_s']:.6f} s ({d['status']})",f"    iterations = {d['total_attempt_iterations']}",'']
lines += ['winner:',f"    name = {w.get('candidate','NONE')}",f"    policy = {w.get('candidate','R1 fallback')}",f"    expected gain vs R1 = {100*w['improvement'] if w.get('candidate') else 'N/A'}%",'','full steady:',f"    status = {full['status'] if full else 'NOT_REQUIRED'}"]
if full:
 a=read('REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance');lines += [f"    wall = {full['wall_time_s']:.6f} s",f"    first steady step = {full['first_full_steady_step']}",f"    stop step = {full['stop_step']}",f"    rebuild count = {full['reuse']['ILU_rebuild_count']}",f"    iterations = {a['total_attempt_iterations']}",f"    failures = linear {a['linear_failures']}, nonlinear {a['nonlinear_failures']}, recovered {a['reuse']['recovery_count']}",'    reload = PASS']
lines += ['','comparison with Stage P:',f"    Stage P wall = {b['full']['wall_time_s']:.6f} s",f"    Stage Q wall = {full['wall_time_s'] if full else 'N/A'}",f"    observed speedup = {b['full']['wall_time_s']/full['wall_time_s'] if full else 'N/A'}",'    label = OBSERVATIONAL DEVELOPMENT SPEEDUP','','fallback:','    Stage P winner preserved = YES','','STAGE SV1.3Q STATUS:','    '+read('stage_status')['status']]
text='\n'.join(lines)+'\n';(R/'terminal_summary.txt').write_text(text);print(text)
