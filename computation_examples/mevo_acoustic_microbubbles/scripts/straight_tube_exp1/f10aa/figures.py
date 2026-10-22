"""Eight publication figures drawn only from saved F10AA measurements."""
import numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import *
COL={'B0':'#60666d','B1':'#2466a1','PACKED':'#27877e','B3_CuPy':'#7956a6','B2_t2':'#b27635','B2_t4':'#a94557'}
LABEL={'B0':'Reference CPU','B1':'Matrix-free CPU (1 thread)','PACKED':'CPU + packed fields','B3_CuPy':'CuPy FP64 hybrid','B2_t2':'CPU (2 threads)','B2_t4':'CPU (4 threads)'}
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,'pdf.fonttype':42,'ps.fonttype':42,'figure.facecolor':'white','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})

def main():
    out=BASE/'figures';out.mkdir(exist_ok=True);perf=read(BASE/'analysis/PERFORMANCE.json');b=BASE/'benchmark_v3'
    fixed=read(b/'FIXED_BENCHMARK.json');geom=read(b/'GEOMETRY.json');linear=read(b/'LINEAR_ALGEBRA.json');step=read(b/'FULL_STEP.json')
    profiles=read(BASE/'baseline/PROFILING.json');audits=read(BASE/'analysis/AUDITS.json');cmp=read(BASE/'full/COMPARISONS.json');extra=read(BASE/'analysis/EXTRA_COMPARISONS.json')
    runs=read(BASE/'analysis/RUNS.json');parallel=read(BASE/'parallel/BATCHES.json');transfer=read(BASE/'profile_extra/TRANSFER_ALLOCATION.json');manifest=[]
    def save(fig,name,caption):
        fig.suptitle(caption,fontsize=14,x=.06,ha='left')
        fig.text(.06,.015,'EXPLORATORY | Uncalibrated physics | Homogeneous-water proxy | training_ready=false',fontsize=8,color='#4b4b4b')
        fig.tight_layout(rect=(0,.04,1,.94))
        for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext))
        manifest.append(dict(name=name,caption=caption,png_sha256=sha(out/(name+'.png')),pdf_sha256=sha(out/(name+'.pdf'))));plt.close(fig)
    def bars(ax,labels,values,ylabel,colors=None):
        ax.bar(np.arange(len(values)),values,color=colors or '#2466a1');ax.set_xticks(np.arange(len(values)),labels);ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.15)
    def groups(ax,data,names,ylabel):
        x=np.arange(len(names));width=.8/len(data)
        for k,(label,val) in enumerate(data.items()):ax.bar(x+(k-(len(data)-1)/2)*width,val,width,label=LABEL.get(label,label),color=COL.get(label))
        ax.set_xticks(x,names);ax.set_ylabel(ylabel);ax.grid(axis='y',alpha=.15)
    # Figure 1: nested costs explicitly not stacked.
    fig,axs=plt.subplots(1,3,figsize=(14,5))
    keys=['dense_contact_geometry','resistance','Cholesky_SPD','linear_solve','active_set_QP','SQP_total','fields_total','wall_total','D2H_wait_and_copy','H2D_upload']
    for ax,n,state in [(axs[0],64,'historical_end'),(axs[1],128,'dispersed')]:
        p=next(p for p in profiles if p['case']==f'N{n}_{state}')
        vals=[p['instrumented_inclusive'].get(k,[0])[0]*1e3 for k in keys];ax.barh(keys[::-1],vals[::-1],color='#60666d');ax.set_xlabel('Instrumented inclusive time (ms)');ax.set_title(f'{n} bubbles: {state.replace("_"," ")}\nNested times; do not add')
    ks=next(r for r in step if r['N']==128 and r['backend']=='B0')['kernels'];names=['Acoustic + APECSS','F2 pair + reduction','Two wall checks']
    vals=[sum(k['elapsed_ms'] for k in ks if 'cubic' in k['name']),sum(k['elapsed_ms'] for k in ks if 'pair_' in k['name']),sum(k['elapsed_ms'] for k in ks if 'cavity' in k['name'])]
    bars(axs[2],names,vals,'CUDA event kernel time (ms)');axs[2].tick_params(axis='x',rotation=25);axs[2].set_title('128 bubbles: actual GPU events\nField interpolation and APECSS are fused')
    save(fig,'Fig_F10AA_1_solver_profiling','Original solver: measured computational bottlenecks')
    fig,axs=plt.subplots(2,2,figsize=(11,8));ns=[32,64,128]
    groups(axs[0,0],{k:[perf[str(n)][k]['step']['p50']*1e3 for n in ns] for k in COL if k in perf['32']},[str(n) for n in ns],'Full-trajectory step p50 (ms)');axs[0,0].legend(fontsize=8)
    groups(axs[0,1],{'All pair constraints':[n*(n-1)//2 for n in ns]},list(map(str,ns)),'Number of constraints')
    gd=[next(r for r in geom if r['case']==f'N{n}_dispersed_ON') for n in ns]
    groups(axs[1,0],{'Dense B':[r['dense_bytes']/2**20 for r in gd],'Matrix-free pairs':[r['matrix_free_bytes']/2**20 for r in gd]},list(map(str,ns)),'B storage (MiB)');axs[1,0].set_yscale('log');axs[1,0].legend()
    groups(axs[1,1],{k:[perf[str(n)][k]['speedup'] for n in ns] for k in ['B1','PACKED','B3_CuPy','B2_t2','B2_t4']},list(map(str,ns)),'Complete-output trajectory speedup');axs[1,1].axhline(1,color='k',lw=.7)
    for ax in axs.flat:ax.set_xlabel('Number of bubbles')
    save(fig,'Fig_F10AA_2_bubble_number_scaling','Scaling at unchanged dt, physics and per-step outputs')
    fig,axs=plt.subplots(2,2,figsize=(11,8))
    groups(axs[0,0],{'Dense B':[r['dense_bytes']/2**20 for r in gd],'Matrix-free':[r['matrix_free_bytes']/2**20 for r in gd]},list(map(str,ns)),'Constraint storage (MiB)');axs[0,0].set_yscale('log');axs[0,0].legend()
    groups(axs[0,1],{'Dense B':[r['times']['dense_B']['p50']*1e3 for r in gd],'Matrix-free':[r['times']['matrix_free_B']['p50']*1e3 for r in gd]},list(map(str,ns)),'Geometry p50 (ms)');axs[0,1].set_yscale('log');axs[0,1].legend()
    selected=[r for r in fixed if r['backend']=='B1' and r['case'].endswith('_ON')];x=np.arange(len(selected));axs[1,0].bar(x,[r['active'] for r in selected],color='#2466a1');axs[1,0].set_xticks(x,[r['case'].replace('_ON','').replace('historical_end','history') for r in selected],rotation=75,fontsize=8);axs[1,0].set_ylabel('Active constraints (all pairs checked)')
    errs=[r['max_fixed_ID_m']*1e9 for r in selected];bars(axs[1,1],['Maximum fixed-ID difference'],[max(errs)],'Reference difference (nm)');axs[1,1].set_title('All frozen ON states\nIncludes 532 simultaneous FCC contacts')
    save(fig,'Fig_F10AA_3_dense_vs_matrix_free','All-pair matrix-free constraints: storage, cost and equivalence')
    fig,axs=plt.subplots(2,2,figsize=(11,8));cases=['N32_dispersed_ON','N64_historical_end_ON','N128_historical_end_ON']
    for ax,key in [(axs[0,0],'factor_s'),(axs[0,1],'solve_s')]:
        vals={label:[next(r for r in linear if r['case']==case and r['backend']==back and r['threads']==t)[key]['p50']*1e3 for case in cases] for label,back,t in [('CPU 1 thread','cpu',1),('CPU 2 threads','cpu',2),('CuPy FP64','cupy',1)]}
        groups(ax,vals,list(map(str,ns)),('Cholesky + upload' if key=='factor_s' else '8-RHS solve + transfers')+' p50 (ms)');ax.legend(fontsize=8)
    groups(axs[1,0],{k:[perf[str(n)][k]['step']['p50']*1e3 for n in ns] for k in ['B1','B3_CuPy']},list(map(str,ns)),'Complete step p50 (ms)');axs[1,0].legend()
    vals={k:[sum(next(r for r in transfer if r['N']==n and r['backend']==k)['cupy_transfer'][d] for d in ['H2D_bytes','D2H_bytes'])/1024 for n in ns] for k in ['B1','B3_CuPy']}
    groups(axs[1,1],vals,list(map(str,ns)),'Extra CuPy transfers per step (KiB)');axs[1,1].set_title('Additional to unchanged Warp transfers')
    for ax in axs.flat:ax.set_xlabel('Number of bubbles')
    save(fig,'Fig_F10AA_4_cpu_vs_gpu','CPU and CuPy FP64: factorization, solve and complete-step costs')
    fig,axs=plt.subplots(2,2,figsize=(11,8))
    for ax,key,label in [(axs[0,0],'fixed_id_m','Fixed-ID difference (nm)'),(axs[0,1],'COM_m','COM difference (nm)'),(axs[1,0],'Rg_m','Rg difference (nm)')]:
        groups(ax,{b:[cmp[f'{n}_{b}']['errors'][key]*1e9 for n in ns] for b in ['B1','PACKED','B3_CuPy']},list(map(str,ns)),label);ax.set_xlabel('Number of bubbles')
    vals=[];names=[]
    for n in ns:
        for variant in ['DT05','FINE']:
            vals.append(cmp[f'{n}_{variant}']['errors']['fixed_id_m']*1e9);names.append(f'{n}\n{variant}')
    bars(axs[1,1],names,vals,'Refinement fixed-ID difference (nm)');axs[1,1].set_title('50 nm original fixed-ID limit\nAll accepted frames, not endpoints only')
    axs[0,0].legend(fontsize=8)
    save(fig,'Fig_F10AA_5_numerical_equivalence','Complete trajectories: fixed IDs, COM, Rg and refinement error')
    fig,axs=plt.subplots(2,2,figsize=(11,8));difficult=[r for r in fixed if r['backend']=='B1' and any(x in r['case'] for x in ['cubic_ON','fcc_ON'])];names=[r['case'].replace('_ON','') for r in difficult]
    bars(axs[0,0],names,[r['active'] for r in difficult],'Active constraints');axs[0,0].tick_params(axis='x',rotation=30)
    bars(axs[0,1],names,[r['iterations'] for r in difficult],'QP working-set iterations');axs[0,1].tick_params(axis='x',rotation=30)
    bars(axs[1,0],names,[r['complementarity_m']/policy()['actual_geometry_complementarity_m'] for r in difficult],'True-gap complementarity / original limit');axs[1,0].tick_params(axis='x',rotation=30)
    keys=list(audits);vals=[audits[k]['force_gate_fraction_max'] for k in keys];axs[1,1].plot(np.arange(len(vals)),vals,'o',ms=4,color='#2466a1');axs[1,1].set_xlabel('Complete trajectory index');axs[1,1].set_ylabel('Independent force residual / original limit');axs[1,1].set_title(f'{sum(a["checked_steps"] for a in audits.values()):,} accepted steps independently checked')
    save(fig,'Fig_F10AA_6_contact_validation','Contact validation: manufactured hard cases and acoustic trajectories')
    fig,axs=plt.subplots(2,2,figsize=(11,8));names=['32\n2 ms','64\n2 ms','128\n1 ms']
    groups(axs[0,0],{k:[perf[str(n)][k]['wall_s'] for n in ns] for k in ['B0','B1','PACKED','B3_CuPy','B2_t2','B2_t4']},names,'Complete-output wall time (s)');axs[0,0].set_ylim(0,65);axs[0,0].legend(fontsize=8,loc='upper left')
    groups(axs[0,1],{k:[perf[str(n)][k]['speedup'] for n in ns] for k in ['B1','PACKED','B3_CuPy','B2_t2','B2_t4']},names,'End-to-end speedup');axs[0,1].axhline(1,color='k',lw=.7)
    groups(axs[1,0],{k:[perf[str(n)][k]['save_s'] for n in ns] for k in ['B0','B1']},names,'Full output writing time (s)');axs[1,0].set_title('Positions, forces and reactions every step')
    er=read(BASE/'extension/RUNS.json') if (BASE/'extension/RUNS.json').exists() else []
    if er:bars(axs[1,1],[r['backend']+'\n'+(r['variant'] if r['status']=='PASS' else 'BUDGET STOP') for r in er],[r['wall_s'] for r in er],'128-bubble extension wall time (s)');axs[1,1].set_title('2 ms optimized; reference stopped at 1.659 ms')
    else:axs[1,1].text(.5,.5,'Longer duration: NOT_RUN',ha='center')
    save(fig,'Fig_F10AA_7_full_trajectory_speedup','Fair end-to-end comparison: identical outputs and physical duration')
    fig,axs=plt.subplots(2,2,figsize=(11,8));workers=[1,2,4]
    for phase,marker in [('cold','s'),('warm','o')]:
        rr=[next(r for r in parallel if r['workers']==w and r['round']==phase) for w in workers]
        axs[0,0].plot(workers,[r['throughput_physical_s_per_wall_s']*1e3 for r in rr],marker+'-',label=phase)
        axs[0,1].plot(workers,[np.mean([j['result']['wall_s'] for j in r['jobs']]) for r in rr],marker+'-',label=phase)
    axs[0,0].set_ylabel('Batch throughput (physical ms / wall s)');axs[0,1].set_ylabel('Mean per-trajectory wall time (s)')
    for ax in axs[0]:ax.set_xlabel('Independent processes (1 BLAS thread each)');ax.set_xticks(workers);ax.legend()
    resource=read(BASE/'full/resources.json');rt=[r['elapsed_s'] for r in resource]
    axs[1,0].plot(rt,[r['cpu_percent']/100 for r in resource],label='CPU core equivalent');axs[1,0].plot(rt,[r.get('gpu_percent',0)/100 for r in resource],label='Whole-device GPU fraction');axs[1,0].set_xlabel('Full campaign elapsed (s)');axs[1,0].set_ylabel('Utilization');axs[1,0].legend(fontsize=8)
    pr=read(BASE/'parallel/resources.json');axs[1,1].plot([r['elapsed_s'] for r in pr],[r.get('gpu_MiB',0) for r in pr],color='#7956a6');axs[1,1].set_xlabel('Parallel campaign elapsed (s)');axs[1,1].set_ylabel('Whole-device VRAM (MiB)');axs[1,1].set_title('GPU capacity used by independent processes')
    save(fig,'Fig_F10AA_8_resource_efficiency','Task throughput is distinct from single-trajectory latency')
    write(out/'MANIFEST.json',dict(figures=manifest,source_sha256=sha(Path(__file__)),input_sha256={str(p.relative_to(BASE)):sha(p) for p in [BASE/'analysis/PERFORMANCE.json',BASE/'analysis/AUDITS.json',BASE/'full/COMPARISONS.json',BASE/'parallel/BATCHES.json',b/'FIXED_BENCHMARK.json',b/'GEOMETRY.json',b/'FULL_STEP.json',b/'LINEAR_ALGEBRA.json',BASE/'profile_extra/TRANSFER_ALLOCATION.json',BASE/'baseline/PROFILING.json',BASE/'full/resources.json',BASE/'parallel/resources.json',BASE/'extension/RUNS.json']},identity=identity()))
if __name__=='__main__':main()
