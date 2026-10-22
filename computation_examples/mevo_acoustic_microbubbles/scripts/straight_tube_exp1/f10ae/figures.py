"""Twelve publication figure groups from measured trajectories and radial solves."""
import time,numpy as np,h5py,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import *
COLORS=['#375e9b','#8768a5','#a56745','#39867b']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'#20252a','axes.labelcolor':'#20252a','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':220})
def traces(r):
    p=BASE/r['trajectory_path']
    if not p.exists():p=BASE/'dataset'/(r['id']+'.h5')
    with h5py.File(p) as h:return h['time_s'][:],h['positions_m'][:]
def metrics(r):return dict(np.load(BASE/'runs'/r['id']/'METRICS.npz'))
def main():
    start=time.perf_counter();out=BASE/'figures';out.mkdir(exist_ok=True)
    rows=read(BASE/'ALL_RUNS.json');by={r['id']:r for r in rows};sel=read(BASE/'SELECTED_FROZEN.json')['selected'];chosen={k:by[v['source']] for k,v in sel.items()};validation=read(BASE/'VALIDATION.json');interpretation=read(BASE/'INTERPRETATION.json');manifest=[]
    def save(fig,n,title,sources):
        fig.suptitle(title,fontsize=14,fontweight='bold');fig.text(.012,.009,'EXPLORATORY | Homogeneous-water control proxy | Shell, contact and electrical-to-acoustic mapping uncalibrated',fontsize=8,color='#3f4852');fig.tight_layout(rect=[0,.035,1,.94])
        stem=f'Fig_F10AE_{n:02d}_'+title.lower().replace(' ','_').replace('/','_')
        for ext in ['png','pdf']:fig.savefig(out/(stem+'.'+ext),bbox_inches='tight')
        manifest.append(dict(number=n,title=title,png=stem+'.png',pdf=stem+'.pdf',sources=sources,sha256={e:sha(out/(stem+'.'+e)) for e in ['png','pdf']}));plt.close(fig)
    # 1: independent regression evidence, including all three historical frequencies.
    reg=read(BASE/'BASELINE_VALIDATION.json');fig,ax=plt.subplots(1,3,figsize=(12,4.3))
    rr=[by[f'baseline_elite{i}'] for i in range(3)]
    for i,r in enumerate(rr):
        m=metrics(r);ax[0].plot(m['time_s']*1e3,(m['aspect_3d']/m['aspect_3d'][0]-1)*100,label=f'Elite {i}, {r["task"].get("frequency",2.5e6)/1e6:g} MHz',color=COLORS[i])
    ax[0].set(xlabel='Physical time (ms)',ylabel='3D aspect change (%)');ax[0].legend(fontsize=8)
    ax[1].bar(['E0','E1','E2'],[z['errors']['fixed_id_m']*1e9 for z in reg['regressions']],color=COLORS[:3]);ax[1].set(ylabel='Maximum fixed-ID difference (nm)',ylim=(0,1));ax[1].text(.05,.85,'All tested historical positions\nbitwise identical',transform=ax[1].transAxes)
    ax[2].axis('off');ax[2].text(.02,.92,'Original F10AA mechanics retained\n\nN = 64; R = 1.78 um; zero flow\n\nCPU resistance + contact\nWarp GPU acoustics\n\nContinuous vs identical segments:\n'+('PASS' if reg['schedule']['admitted'] else 'FAIL')+'\nOFF-stage motion: '+f'{reg["schedule"]["OFF_stage_displacement_m"]*1e9:g} nm',va='top')
    save(fig,1,'Historical baseline and regression',['BASELINE_VALIDATION.json','BASELINE_RUNS.json'])
    # 2: measured radial and force scales; failed lowest refinement remains visible.
    low=read(BASE/'LOW_PRESSURE_ANALYSIS.json');admit=read(BASE/'LOW_PRESSURE_ADMISSION.json');raw=read(BASE/'LOW_PRESSURE_RAW.json');valid=[r for r in low['checks'] if r.get('radial_checks')];p=np.array([r['pressure_Pa'] for r in valid]);fig,ax=plt.subplots(2,3,figsize=(12,7.2))
    ax[0,0].loglog(p,[r['volume_amplitude_m3'] for r in valid],'o-');ax[0,0].axvline(.1,color=COLORS[2],ls='--',label='Historical lower limit');ax[0,0].set(xlabel='Pressure (Pa)',ylabel='Fundamental volume amplitude (m3)');ax[0,0].legend(fontsize=8)
    ax[0,1].semilogx(p,np.degrees([r['phase_rad'] for r in valid]),'o-');ax[0,1].set(xlabel='Pressure (Pa)',ylabel='Complex-response phase (deg)');ax[0,1].ticklabel_format(axis='y',style='plain',useOffset=False)
    ax[0,2].loglog(p,[r['fine_finest_complex'] for r in valid],'o-',label='Complex response');ax[0,2].loglog(p,[r['F2_energy_relative'] for r in valid],'s-',label='F2 derivative energy');ax[0,2].axhline(.01,c='#a04d46',ls='--');ax[0,2].set(xlabel='Pressure (Pa)',ylabel='Relative refinement difference');ax[0,2].legend(fontsize=8)
    checks=[r for r in low['interpolation'] if r.get('valid')];ax[1,0].loglog([r['pressure_Pa'] for r in checks],[r['complex_relative'] for r in checks],'o-');ax[1,0].axhline(.01,c='#a04d46',ls='--');ax[1,0].set(xlabel='Independent midpoint pressure (Pa)',ylabel='Complex interpolation error')
    forces=read(BASE/'LOW_PRESSURE_FINAL.json')['actual_force_checks']
    for key,label in [('F1_rms_N','Primary F1'),('F2_rms_N','Secondary F2')]:ax[1,1].loglog([r['pressure_min_Pa'] for r in forces],[r[key] for r in forces],'o-',label=label)
    ax[1,1].set(xlabel='Minimum incident pressure (Pa)',ylabel='Force RMS per component (N)');ax[1,1].legend(fontsize=8)
    ax[1,2].axis('off');ax[1,2].text(.02,.95,f'Official APECSS: {len(raw)} attempts\nLowest finest solve: TIMEOUT\n\n0.0001 Pa: NOT ADMITTED\n0.001 Pa extension: {"PASS" if admit["production_extension_admitted"] else "FAIL"}\n\nHistorical join difference:\n{low["junction_relative_error"]:.3g}\n\nOriginal files and >=0.1 Pa queries\nremain unchanged.',va='top',fontsize=10)
    save(fig,2,'Low pressure APECSS validation',['LOW_PRESSURE_RAW.json','LOW_PRESSURE_ADMISSION.json','LOW_PRESSURE_FINAL.json'])
    # 3: target is a metric interval, not an invented particle trajectory.
    fig=plt.figure(figsize=(12,4.6));a=fig.add_subplot(121,projection='3d');_,x=traces(chosen['T1']);a.scatter(*x[0].T*1e6,s=24,c=COLORS[0],depthshade=True);a.set(xlabel='X (um)',ylabel='Y (um)',zlabel='Z (um)',title='Actual S2 discovery initial state');a.set(xlim=(-28,28),ylim=(-12,12),zlim=(-12,12));a.set_box_aspect((56,24,24));a.view_init(26,-60)
    a=fig.add_subplot(122);a.axis('off');table=a.table(cellText=[['T1','Aspect +3% +/-0.5 pp; COM <=20 nm'],['T2','Aspect +0.4% +/-0.2 pp; connected endpoint'],['T3','Absolute RMS axis > max(1 nm, 5 x noise)'],['T4-R1/R2','Rg / all axes within 0.5%; >=50% recovery'],['T4-R3/R4','Aligned Q error <=1%; fixed-ID RMS <=5 nm'],['T4-R5','Neighbor Jaccard >=0.95; >=50% recovery']],colLabels=['Target','Frozen criterion'],loc='center',cellLoc='left',colWidths=[.18,.82]);table.auto_set_font_size(False);table.set_fontsize(9);table.scale(1,2.1)
    save(fig,3,'Target shape definitions',['INITIAL_STATES.h5','PREREGISTRATION.json'])
    # 4: failed trials have no endpoint loss.
    comp=read(BASE/'OPTIMIZER_COMPARISON.json');fig,ax=plt.subplots(1,3,figsize=(12,4.2))
    for i,(file,label) in enumerate([('SEARCH_SOBOL.json','Sobol'),('SEARCH_AX.json','Ax / BoTorch')]):
        d=read(BASE/file);best=np.inf;history=[]
        for r in d['rows']:
            if r['status']=='PASS' and r['targets']['COM_TARGET_PASS']:best=min(best,r['targets']['T1_loss'])
            history.append(best if np.isfinite(best) else np.nan)
        ax[0].plot(np.arange(1,len(history)+1),history,label=label,color=COLORS[i]);ax[1].bar(i,d['wall_s'],color=COLORS[i]);ax[2].bar(i,d['optimizer_s'],color=COLORS[i])
    ax[0].axhline(1,c='#a04d46',ls='--',label='T1 error tolerance');ax[0].set(xlabel='Charged physical observations',ylabel='Best COM-feasible squared target loss',yscale='log');ax[0].legend(fontsize=8)
    for a,label in zip(ax[1:],['Charged method wall (s)','Optimizer ask/tell wall (s)']):a.set(xticks=[0,1],xticklabels=['Sobol','Ax / BoTorch'],ylabel=label)
    save(fig,4,'Constrained optimizer comparison',['OPTIMIZER_COMPARISON.json','SEARCH_SOBOL.json','SEARCH_AX.json'])
    # 5: true particle centers and trajectories for selected targets.
    fig=plt.figure(figsize=(12,7))
    for i,k in enumerate(['T1','T2','T3']):
        r=chosen[k];t,x=traces(r);a=fig.add_subplot(2,3,i+1,projection='3d');a.scatter(*x[0].T*1e6,s=12,c='#aab2bb',alpha=.5,label='Initial');a.scatter(*x[-1].T*1e6,s=20,c=COLORS[i],label='Actual endpoint');a.set(xlabel='X (um)',ylabel='Y (um)',zlabel='Z (um)',title=k+' | '+r['id']+'\nIndependent target: '+('PASS' if interpretation['targets'][k]['confirmed'] else 'NOT CONFIRMED'));a.set(xlim=(-28,28),ylim=(-12,12),zlim=(-12,12));a.set_box_aspect((56,24,24));a.view_init(25,-60)
        if i==0:a.legend(fontsize=7)
        m=metrics(r);a=fig.add_subplot(2,3,i+4);a.plot(m['time_s']*1e3,(m['aspect_3d']/m['aspect_3d'][0]-1)*100,color=COLORS[i]);a.set(xlabel='Physical time (ms)',ylabel='3D aspect change (%)')
        if k in ['T1','T2']:
            cfg=read(CONFIG)[k];a.axhspan((cfg['aspect_gain']-cfg['tolerance'])*100,(cfg['aspect_gain']+cfg['tolerance'])*100,color='#b5c8d8',alpha=.35)
    save(fig,5,'Targeted shape trajectories',[v['source'] for v in sel.values()])
    # 6: retain the incidental seed-specific extension, without claiming cross-seed control.
    fig,ax=plt.subplots(1,3,figsize=(12,4.6))
    examples=[('extension_028','Discovery state: contraction'),('validate_T1_S2_validation_1_DT1','New validation state: split distribution'),('final_T1_0','Final test state: refinement NOT_RUN')]
    for a,(rid,title) in zip(ax,examples):
        r=by[rid];m=metrics(r);length=m['principal_rms_lengths_3d_m']
        for j in range(3):a.plot(m['time_s']*1e3,(length[:,j]-length[0,j])*1e9,label=f'L{j+1}',color=COLORS[j])
        a.axhline(0,c='#444444',lw=.7);a.set(xlabel='Physical time (ms)',ylabel='Absolute RMS axis change (nm)',title=title);a.legend(fontsize=8)
        a.text(.03,.05,f'Final L1 change: {r["targets"]["axis_delta_m"][0]*1e9:+.3f} nm\nComponents: {r["targets"]["components"]}',transform=a.transAxes,fontsize=9)
    save(fig,6,'Absolute principal axis extension',['SEED_SPECIFIC_EXTENSION.json','VALIDATION.json','FINAL_TEST_RUNS.json'])
    # 7: real per-channel actions and clock.
    r=chosen['T4'];fig,ax=plt.subplots(2,2,figsize=(12,6.8));segments=r['task']['segments'];tt=np.array([0,.002,.002,.004])*1e3
    for j,name in enumerate(['Y-','Y+','X-','X+']):
        ax[0,0].plot(tt,[segments[0]['amplitudes'][j]]*2+[segments[1]['amplitudes'][j]]*2,label=name,color=COLORS[j]);ax[0,1].plot(tt,np.degrees([segments[0]['phases'][j]]*2+[segments[1]['phases'][j]]*2),label=name,color=COLORS[j])
    ax[0,0].set(ylabel='Amplitude (Pa-equivalent)');ax[0,0].legend(ncol=4);ax[0,1].set(ylabel='Phase offset (deg)')
    m=metrics(r);ax[1,0].plot(m['time_s']*1e3,m['Rg_3d_m']*1e6);ax[1,0].set(ylabel='Rg (um)');ax[1,1].plot(m['time_s']*1e3,np.linalg.norm(m['center_3d_m']-m['center_3d_m'][0],axis=1)*1e9);ax[1,1].set(ylabel='COM displacement (nm)')
    for a in ax.flat:a.axvline(2,c='#888888',ls='--',lw=.8);a.set(xlabel='Physical time (ms)')
    save(fig,7,'Two stage acoustic action sequence',['RECOVERY_PLAN.json',r['trajectory_path']])
    # 8: actual forward simulations; no reverse playback.
    fig,ax=plt.subplots(2,2,figsize=(12,6.8))
    for r,label,c in [(chosen['T4'],'Selected active second action',COLORS[0]),(by['schedule_C'],'First action then OFF',COLORS[2])]:
        m=metrics(r);t=m['time_s']*1e3;ax[0,0].plot(t,(m['Rg_3d_m']/m['Rg_3d_m'][0]-1)*100,label=label,c=c);ax[0,1].plot(t,m['kabsch_rmsd_3d_m']*1e9,c=c);ax[1,0].plot(t,m['pair_distance_rms_change_3d_m']*1e9,c=c);rot=m['rigid_rotation_3d'];q=m['covariance_3d_m2'];aligned=np.swapaxes(rot,-1,-2)@q@rot;err=np.linalg.norm(aligned-q[0],axis=(1,2))/np.linalg.norm(q[0]);ax[1,1].plot(t,err*100,c=c)
    for a,label in zip(ax.flat,['Rg change (%)','Fixed-ID Kabsch RMS (nm)','Fixed-ID pair-distance RMS change (nm)','Aligned covariance error (%)']):a.axvline(2,c='#888888',ls='--');a.set(xlabel='Physical time (ms)',ylabel=label)
    ax[0,0].legend(fontsize=8);save(fig,8,'Deformation and attempted recovery',['RECOVERY_RUNS.json','BASELINE_RUNS.json'])
    # 9: global distribution vs actual geometric components.
    fig,ax=plt.subplots(1,3,figsize=(12,4.4))
    for i,(k,r) in enumerate(chosen.items()):
        m=metrics(r);ax[0].step(m['time_s']*1e3,m['proximity_components'],where='post',label=k,c=COLORS[i]);ax[1].plot(m['time_s']*1e3,m['largest_component'],label=k,c=COLORS[i]);d=r['component_diagnostic'];ax[2].bar(i,d['within_component_Rg_squared_m2']*1e12,color=COLORS[0],label='Within groups' if i==0 else None);ax[2].bar(i,d['between_component_Rg_squared_m2']*1e12,bottom=d['within_component_Rg_squared_m2']*1e12,color=COLORS[2],label='Between groups' if i==0 else None)
    ax[0].set(xlabel='Physical time (ms)',ylabel='Components at 2.1 R');ax[0].legend();ax[1].set(xlabel='Physical time (ms)',ylabel='Largest component size');ax[2].set(xticks=range(len(chosen)),xticklabels=list(chosen),ylabel='Rg squared decomposition (um2)');ax[2].legend(fontsize=8)
    save(fig,9,'Connectivity and splitting',['SELECTED_FROZEN.json','VALIDATION.json'])
    # 10: numerical admission is explicit per frequency.
    stress=read(BASE/'ROBUSTNESS_RUNS.json');cc=read(BASE/'ROBUSTNESS_CONVERGENCE.json');fig,ax=plt.subplots(1,3,figsize=(12,4.2))
    for i,k in enumerate(['T1','T2']):
        rrs=sorted([r for r in stress if r['id'].startswith('freq_'+k) and r['status']=='PASS'],key=lambda r:r['task']['frequency'])
        for a,key,mult in zip(ax,['aspect_gain','Rg_gain','COM_m'],[100,100,1e9]):
            for r in rrs:
                numeric=all(c['admitted'] for c in cc if c.get('reference')==r['id']);a.scatter(r['task']['frequency']/1e6,r['targets'][key]*mult,marker='o' if numeric else 'x',c=COLORS[i],label=k if r is rrs[0] else None)
    for a,y in zip(ax,['3D aspect change (%)','Rg change (%)','COM displacement (nm)']):a.set(xlabel='Frequency (MHz)',ylabel=y,xticks=[2.4,2.5,2.6]);a.legend(fontsize=8)
    save(fig,10,'Frequency robustness',['ROBUSTNESS_RUNS.json','ROBUSTNESS_CONVERGENCE.json'])
    # 11: shell and spatial offsets including explicit Y-pair node diagnostic.
    fig,ax=plt.subplots(1,3,figsize=(12,4.2))
    for i,k in enumerate(['T1','T2']):
        rr=sorted([r for r in stress if (r['id'].startswith('shell_'+k) or r['id']==f'freq_{k}_2500000') and r['status']=='PASS'],key=lambda r:r['task'].get('shell',1.));ax[0].plot([r['task'].get('shell',1.) for r in rr],[r['targets']['aspect_gain']*100 for r in rr],'o-',label=k,c=COLORS[i])
        rr=sorted([r for r in stress if r['id'].startswith('position_'+k) and r['status']=='PASS'],key=lambda r:r['task']['offset_m'][1]);ax[1].plot([r['task']['offset_m'][1]*1e6 for r in rr],[r['targets']['aspect_gain']*100 for r in rr],'o-',label=k,c=COLORS[i])
    rr=sorted([r for r in stress if r['id'].startswith('Ynode') and r['status']=='PASS'],key=lambda r:r['task']['offset_m'][1]);ax[2].plot([r['task']['offset_m'][1]*1e6 for r in rr],[r['targets']['aspect_gain']*100 for r in rr],'o-',c=COLORS[2]);ax[2].set_title('Y pair, 700 Pa-eq, 90 deg')
    ax[0].set(xlabel='Shell elasticity multiplier',ylabel='3D aspect change (%)');ax[1].set(xlabel='Initial Y offset (um)',ylabel='3D aspect change (%)');ax[2].set(xlabel='Initial Y offset (um)',ylabel='3D aspect change (%)');ax[0].legend();ax[1].legend()
    save(fig,11,'Shell and spatial sensitivity',['ROBUSTNESS_RUNS.json'])
    # 12: all independent comparisons and real cost.
    fig,ax=plt.subplots(2,2,figsize=(12,7));comparisons=validation['comparisons']+cc
    for key,color in [('DT05',COLORS[0]),('FINE',COLORS[1])]:
        rr=[r for r in comparisons if key in r.get('comparison','') and 'errors' in r];ax[0,0].plot([r['errors']['fixed_id_m']*1e9 for r in rr],'.',label=key,c=color)
    ax[0,0].axhline(50,c='#a04d46',ls='--');ax[0,0].set(xlabel='Comparison index',ylabel='Maximum fixed-ID difference (nm)');ax[0,0].legend()
    audit=read(BASE/'AUDIT.json');a=[r['audit'] for r in audit['rows'] if r['scientific_status']=='VALID_COMPLETE_TRAJECTORY'];ax[0,1].plot([r['force_gate_fraction_max'] for r in a],'.',c=COLORS[0]);ax[0,1].axhline(1,c='#a04d46',ls='--');ax[0,1].set(xlabel='Complete trajectory index',ylabel='Independent force residual / limit');ax[0,1].text(.03,.9,f'Max overlap: {max(z["overlap_m"] for z in a)*1e9:.3g} nm\nMax complementarity distance: {max(z["actual_complementarity_m"] for z in a)*1e9:.3g} nm',transform=ax[0,1].transAxes,fontsize=9)
    groups=validation['groups'];ax[1,0].bar([g['target'] for g in groups],[sum(z['confirmed'] for z in g['replicates'][1:]) for g in groups],color=COLORS);ax[1,0].set(ylabel='Independent validation confirmations / 3',ylim=(0,3.4))
    ax[1,1].axis('off');valid=[r for r in rows if r['status']=='PASS'];ax[1,1].text(.02,.95,f'Complete trajectories: {len(valid)} / {len(rows)}\nIndependently audited steps: {audit["checked_steps"]:,}\n\nMedian complete solver wall: {np.median([r["wall_s"] for r in valid]):.2f} s\nScientific workflow wall: {read(BASE/"WORKFLOW_COMPLETE.json")["wall_s"]:.1f} s\n\nNumerical admission is separate from target success.\nNo experimental calibration; training_ready=false.',va='top')
    save(fig,12,'Numerical and scientific admission',['VALIDATION.json','AUDIT.json','WORKFLOW_COMPLETE.json'])
    write(out/'MANIFEST.json',dict(figures=manifest,wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
