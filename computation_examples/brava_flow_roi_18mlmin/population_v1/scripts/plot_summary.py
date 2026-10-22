"""English scientific figures from completed records, no prescribed outcomes."""
from pathlib import Path
import argparse,csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ap=argparse.ArgumentParser();ap.add_argument('--run-name',required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];out=root/'figures';out.mkdir(exist_ok=True);analysis=root/'analysis'/args.run_name
rows=list(csv.DictReader((analysis/'trajectory_catalog.csv').open()));timeline=list(csv.DictReader((analysis/'population_timeline.csv').open()))
audit=json.loads((analysis/'INDEPENDENT_AUDIT.json').read_text());s=audit['summary']
colors={'OUTLET_01':'#9c71ee','OUTLET_02':'#1fa9bd','OUTLET_03':'#e58b2b','SOLVER_FAILURE':'#e34949','OBSERVATION_CENSORED':'#de60a0','SUPPORTED_STATIONARY':'#778799'}
fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
t=np.array([float(r['physical_time_s']) for r in timeline]);active=[int(r['active_count']) for r in timeline]
axes[0,0].plot(t,active,color='#167ca0');axes[0,0].axvline(s['infusion_stop_time_s'],ls='--',color='k',label='Source stopped: all outlets >= 20')
axes[0,0].set(xlabel='Absolute physical time (s)',ylabel='Simultaneously active particles',title='Shared population clock');axes[0,0].legend(fontsize=9)
axes[0,0].annotate(f"{s['censored_count']} observation-censored\n(not outlet crossings)",
 xy=(s['physical_time_s'],s['censored_count']),xytext=(.52,.28),textcoords='axes fraction',
 fontsize=9,color=colors['OBSERVATION_CENSORED'],arrowprops=dict(arrowstyle='->',color=colors['OBSERVATION_CENSORED']))
for o in s['outlets']:axes[0,1].plot(t,[int(r[o]) for r in timeline],color=colors[o],label=o.replace('OUTLET_0','O'))
axes[0,1].axhline(20,color='k',ls=':',label='Coverage target');axes[0,1].set(xlabel='Absolute physical time (s)',ylabel='Actual cumulative center crossings',title='Observed outlet coverage');axes[0,1].legend()
keys=sorted(set(r['outlet'] or r['status'] for r in rows));counts=[sum((r['outlet'] or r['status'])==k for r in rows) for k in keys]
labels=[k.replace('OUTLET_0','O').replace('SOLVER_FAILURE','Solver stop').replace('OBSERVATION_CENSORED','Censored').replace('SUPPORTED_STATIONARY','Supported stop') for k in keys]
bars=axes[1,0].bar(labels,counts,color=[colors[k] for k in keys]);axes[1,0].bar_label(bars);axes[1,0].tick_params(axis='x',labelrotation=18)
axes[1,0].set(ylabel='Accepted microbubbles',title=f'All {len(rows)} accepted births: separate terminal classes')
for k,label in zip(keys,labels):
 r=[float(q['residence_time_s']) for q in rows if (q['outlet'] or q['status'])==k]
 axes[1,1].hist(r,bins=np.linspace(0,max(float(q['residence_time_s']) for q in rows),35),histtype='step',lw=1.6,label=label,color=colors[k])
axes[1,1].set(xlabel='Recorded residence / observation time (s)',ylabel='Count',title='Censored/failed times are not completed transit times');axes[1,1].legend(fontsize=8)
fig.savefig(out/'03_population_accounting.png',dpi=300);fig.savefig(out/'03_population_accounting.pdf');plt.close(fig)
gpu=json.loads((root/'validation/phase9_benchmark_gpu.json').read_text())['rows'];cpu=json.loads((root/'validation/phase9_benchmark_cpu.json').read_text())['rows']
fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
for records,label,color in [(cpu,'CPU component reference (1 core)','#555555'),(gpu,'GPU hybrid (RTX 4090)','#167ca0')]:
 axes[0].loglog([r['N_active'] for r in records],[r['particle_updates_s'] for r in records],'o-',label=label,color=color)
axes[0].set(xlabel='Active particles in synthetic benchmark',ylabel='Certified particle updates / s',title='Full 0.5 ms step; IO excluded');axes[0].legend(fontsize=9)
g=gpu[-1];names=['FEM sample','Wall query','Neighbors','Resistance','Contact solve','Wall path','Transfer + IO'];values=[g[k] for k in ['field_sampling_s','wall_query_s','neighbor_build_and_exact_gap_s','resistance_assembly_s','contact_solve_s','continuous_wall_certificate_s','batched_state_transfer_and_compressed_IO_s']]
axes[1].barh(names,np.array(values)*1000,color='#167ca0');axes[1].set(xlabel='Wall time / CUDA event duration (ms)',title='N = 10,000; isolated stage measurements')
fig.suptitle('Benchmark fixtures are not the physical infusion population; single full-step timings')
fig.savefig(out/'04_performance.png',dpi=300);fig.savefig(out/'04_performance.pdf');plt.close(fig)
with (analysis/'performance_summary.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=['N_active','CPU_wall_s','GPU_wall_s','CPU_updates_s','GPU_updates_s','CPU_over_GPU_wall_ratio']);w.writeheader()
 for c,g in zip(cpu,gpu):w.writerow(dict(N_active=c['N_active'],CPU_wall_s=c['wall_s'],GPU_wall_s=g['wall_s'],CPU_updates_s=c['particle_updates_s'],GPU_updates_s=g['particle_updates_s'],CPU_over_GPU_wall_ratio=c['wall_s']/g['wall_s']))
print('FIGURES_COMPLETE',flush=True)
