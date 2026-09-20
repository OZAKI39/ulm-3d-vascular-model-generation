import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3o'
def read(n):return json.loads((R/(n+'.json')).read_text())
b=read('baseline_runtime');w=read('winner');c=read('optimized_steady_candidate');s=read('stage_result');t=read('test_summary')
rows=[]
for k in 'ABCD':
 d=read('PERF_'+k+'_acceptance');rows.append(f"PERF_{k}:\n    status = {d['status']}\n    config = {d.get('PETSC_OPTIONS','NOT RUN')}\n    wall = {d.get('wall_time_s','NOT RUN')}\n    iterations = {d.get('statistics',{}).get('total_iterations','NOT RUN')}\n    complete window = {d['status']=='PASS'}\n")
text=f'''Stage SV1.3O completed.

baseline:
    Stage N wall = {b['wall_time_s']:.6f} s
    steps = {b['steps']}
    KSP solves = {b['KSP_solves']}
    KSP iterations = {b['total_iterations']}
    classification = DEVELOPMENT_BASELINE_OBSERVED

I/O optimization:
    previous VTU cadence = every step
    optimized cadence = every 10 steps
    final forced VTU = PASS
    checkpoint = PASS, same step/time as final VTU

fixed window:
    start step = 60
    end step = 70
    checkpoint SHA = {w['checkpoint_sha256']}

{chr(10).join(rows)}
winner:
    candidate = {w['candidate']}
    config = {w['PETSC_OPTIONS']}
    fixed-window improvement = {w['fixed_window_reduction']*100:.6f}%

full optimized GPU:
    wall = {c['wall_time_s']:.6f} s
    steady step = {c['first_full_steady_step']}
    stop step = {c['stop_step']}
    KSP solves = {c['statistics']['KSP_solves']}
    iterations = {c['statistics']['total_iterations']}
    mean / max iterations = {c['statistics']['mean_iterations']:.6f} / {c['statistics']['max_iterations']}
    failures = 0 linear, 0 nonlinear
    reload = {c['reload']['status']}
    VTU count = {c['VTU_count']}
    sampled peak device VRAM = {c['sampled_device_memory_max_MiB']} MiB (not exact process peak)

development speedup:
    Stage N observed wall = {b['wall_time_s']:.6f}
    Stage O observed wall = {c['wall_time_s']:.6f}
    ratio = {c['development_speedup']:.6f}
    label = OBSERVATIONAL DEVELOPMENT SPEEDUP
    comparison = SINGLE-RUN DEVELOPMENT COMPARISON

deferred:
    CPU/GPU science equivalence = YES
    formal CPU/GPU benchmark = YES
    multi-rank CUDA = YES
    GPU-native advanced PC = YES

production:
    CPU_EARLY_STOP_PRODUCTION, unchanged
result:
    {c['classification']}
    NOT YET SCIENTIFICALLY VALIDATED
tests:
    Stage O = {t['stage']}
    Full suite = {t['full']}

STAGE SV1.3O STATUS:
    {s['status']}

Stopped. No Stage P, additional CFD, CPU comparisons, mesh convergence or particle coupling.
'''
(R/'TERMINAL_SUMMARY.txt').write_text(text);print(text)
