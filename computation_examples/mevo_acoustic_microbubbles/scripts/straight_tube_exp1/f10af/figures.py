"""Twelve data-backed publication figure groups; every source is indexed."""
import time
import numpy as np,h5py,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import *
C=['#426b9c','#83629f','#bd7958','#478f86']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'#20252a','axes.labelcolor':'#20252a','axes.edgecolor':'#40454b','xtick.color':'#30353b','ytick.color':'#30353b','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':220})
def series(r):return dict(np.load(BASE/r['series_path']))
def main():
    begin=time.perf_counter();out=BASE/'figures';out.mkdir(exist_ok=True);sci=read(BASE/'SCIENTIFIC_RESULTS.json');rows=sci['rows'];by={r['id']:r for r in rows};manifest=[]
    def save(fig,i,title,source):
        fig.suptitle(title,fontsize=15,fontweight='bold');fig.text(.012,.008,'EXPLORATORY | Homogeneous water | No flow or streaming | Uncalibrated physics and electrical mapping',fontsize=8,color='#46505c');fig.tight_layout(rect=(0,.04,1,.94));name=f'Fig_F10AF_{i:02d}_'+title.lower().replace(' ','_').replace('/','_')
        for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext))
        manifest.append(dict(number=i,title=title,sources=source,files={ext:dict(path=name+'.'+ext,sha256=sha(out/(name+'.'+ext))) for ext in ['png','pdf']}));plt.close(fig)
    def panel_text(ax,txt):ax.axis('off');ax.text(.02,.98,txt,va='top',linespacing=1.6,fontsize=9,transform=ax.transAxes)
    fig,ax=plt.subplots(1,3,figsize=(12.8,4.7));reg=read(BASE/'HISTORICAL_REPLAY.json');panel_text(ax[0],'Protected F10AE baseline\n\nHistorical trajectory replay:\n'+('BITWISE IDENTICAL' if reg['bitwise_positions_equal'] else 'FAIL')+'\n\nOriginal F10AA mechanics retained.\nCPU contact + Warp GPU forces.')
    with h5py.File(BASE/'INITIAL_STATES.h5') as h:
        for name,col in [('N64_wide_r0',C[0]),('N64_compact_reference',C[2])]:x=h[name][:];ax[1].scatter(x[:,0]*1e6,x[:,1]*1e6,s=11,c=col,label='Dispersed cloud' if 'wide' in name else 'Compact reference')
    ax[1].set(xlabel='X (µm)',ylabel='Y (µm)',title='Actual initial particle centers');ax[1].set_aspect('equal');ax[1].legend(fontsize=8)
    panel_text(ax[2],'New question: spatial organization\n\nP1: stripes\nP2: lattice / composite pattern\nP3: collective rotation\nP4: switching (conditional)\n\nR = 1.78 µm; monodisperse\nCamera observability UNVERIFIED')
    save(fig,1,'Baseline and observable-scale scope',['HISTORICAL_REPLAY.json','INITIAL_STATES.h5'])
    fig,g2=plt.subplots(2,4,figsize=(14,8.2));ax=g2[0];fld=BASE/'fields/2500000/DIAGNOSTICS.h5'
    with h5py.File(fld) as h:
        a=h['xy_axis_m'][:]*1e6
        for j,name in enumerate(['X_PAIR','XY_INPHASE','ROT_PLUS']):
            p=abs(h[name+'/pressure_Pa'][:]).reshape(len(a),len(a));im=ax[j].imshow(p.T,origin='lower',extent=[a[0],a[-1],a[0],a[-1]],cmap='magma',aspect='equal');ax[j].set(title={'X_PAIR':'X pair','XY_INPHASE':'XY in-phase','XY_QUADRATURE':'XY quadrature','ROT_PLUS':'Positive phase order'}.get(name,name),xlabel='X (µm)',ylabel='Y (µm)');fig.colorbar(im,ax=ax[j],shrink=.7,label='Pressure amplitude (Pa)')
    with h5py.File(fld) as h:
        for j,name in enumerate(['X_PAIR','XY_INPHASE','ROT_PLUS']):
            g=np.linalg.norm(h[name+'/gradient_Pa_m'][:],axis=1).reshape(len(a),len(a))/1e6;im=g2[1,j].imshow(g.T,origin='lower',extent=[a[0],a[-1],a[0],a[-1]],cmap='magma',aspect='equal');g2[1,j].set(title='Complex-gradient norm',xlabel='X (µm)',ylabel='Y (µm)');fig.colorbar(im,ax=g2[1,j],shrink=.7,label='Norm of gradient P (MPa/m)')
    panel_text(g2[1,3],'Cached central water domain\n\nX: -620 to 620 µm\nY: -520 to 520 µm\nZ: -100 to 100 µm\n\nCoarse / fine: 10 / 5 µm\nShown plane: Z = 0\nFour independent source bases\nNo extrapolation outside cache')
    for k,col in [('coarse_gradient',C[0]),('fine_gradient',C[1])]:ax[3].plot([2.4,2.5,2.6],[100*read(BASE/'fields'/str(f)/'FIELD_GATE.json')['checks'][k] for f in [2400000,2500000,2600000]],'o-',label=k.replace('_gradient',''),c=col)
    ax[3].set(xlabel='Frequency (MHz)',ylabel='Relative gradient difference (%)');ax[3].legend();save(fig,2,'Extended acoustic fields and convergence',['fields/*/FIELD_GATE.json','fields/2500000/DIAGNOSTICS.h5'])
    fig,ax=plt.subplots(1,3,figsize=(12.8,4.5))
    with h5py.File(fld) as h:
        for j,deg in enumerate([0,90,180]):
            x=h[f'line_{deg}/x_m'][:]*1e6;ax[0].plot(x,abs(h[f'line_{deg}/pressure_Pa'][:]),label=f'{deg} deg',c=C[j]);ax[1].plot(x,h[f'line_{deg}/force_N'][:,0]*1e12,c=C[j])
    ax[0].set(xlabel='X (µm)',ylabel='Pressure amplitude (Pa)');ax[0].legend();ax[1].axhline(0,c='#888',lw=.7);ax[1].set(xlabel='X (µm)',ylabel='Primary force Fx (pN)')
    for j,r in enumerate([z for z in read(BASE/'fields/2500000/STATIC_FEASIBILITY.json')['stripes'] if z['phase_deg'] in [0,90,180]]):
        ax[1].scatter([q['x_m']*1e6 for q in r['stable_X_roots']],np.zeros(len(r['stable_X_roots'])),marker='o',s=24,facecolors='white',edgecolors=C[j])
    ax[1].text(.02,.98,'Open circles: X-restoring zeros',va='top',transform=ax[1].transAxes,fontsize=8)
    for f,col in zip([2400000,2500000,2600000],C):
        d=read(BASE/'fields'/str(f)/'STATIC_FEASIBILITY.json')['stripes'];ax[2].plot([z['phase_deg'] for z in d],[np.mean(z['spacing_m'])*1e6 for z in d],'o-',label=f'{f/1e6:g} MHz',c=col)
    ax[2].axhline(300,c='#666',ls='--',label='Ideal half wavelength (2.5 MHz)');ax[2].set(xlabel='X+ phase (deg)',ylabel='Adjacent stable-X spacing (µm)');ax[2].legend(fontsize=8,loc='center',bbox_to_anchor=(.5,.7));save(fig,3,'Force-based stripe stability',['fields/*/STATIC_FEASIBILITY.json','fields/2500000/DIAGNOSTICS.h5'])
    fig,ax=plt.subplots(1,3,figsize=(12.8,4.5));extra=read(BASE/'STATIC_EXTRA.json')
    with h5py.File(fld) as h:
        a=h['xy_axis_m'][:]*1e6
        for j,name in enumerate(['XY_INPHASE','XY_QUADRATURE','ROT_PLUS']):
            P=abs(h[name+'/pressure_Pa'][:]).reshape(101,101);F=h[name+'/primary_force_N'][:].reshape(101,101,3);im=ax[j].imshow(P.T,origin='lower',extent=[-500,500,-500,500],cmap='magma',alpha=.75);fig.colorbar(im,ax=ax[j],shrink=.65,label='Pressure amplitude (Pa)');ax[j].quiver(a[::10],a[::10],F[::10,::10,0].T,F[::10,::10,1].T,color='white',width=.004)
            roots=next(r['roots'] for r in extra['equilibria'] if r['frequency']==2500000 and r['action']==name)
            for q in roots:ax[j].plot(q['position_m'][0]*1e6,q['position_m'][1]*1e6,'o' if q['stable_3D'] else 'x',c='#77f1c1' if q['stable_3D'] else '#eee',ms=8)
            ax[j].set(title={'X_PAIR':'X pair','XY_INPHASE':'XY in-phase','XY_QUADRATURE':'XY quadrature','ROT_PLUS':'Positive phase order'}.get(name,name),xlabel='X (µm)',ylabel='Y (µm)',aspect='equal')
    fig.text(.4,.045,'Circle: locally stable in 3D | Cross: not stable in 3D | Arrows: primary force (panel-wise scale)',ha='center',fontsize=8)
    save(fig,4,'XY equilibria are not particle lattices',['STATIC_EXTRA.json','fields/2500000/DIAGNOSTICS.h5'])
    fig,grid=plt.subplots(2,3,figsize=(12.8,7.8));ax=grid[0];source=read(BASE/'SOURCE_ORDER.json')
    q=np.array([v['center_local_m'][:2] for v in source['elements']])*1e3;ax[0].scatter(q[:,0],q[:,1],s=180,c=C)
    for i in range(4):ax[0].annotate(source['elements'][i]['id']+f'\n{np.degrees(source["positive_phase"][i]):g} deg',q[i],xytext=(4,5),textcoords='offset points',fontsize=9)
    norm=np.linalg.norm(q,axis=1);ax[0].quiver(q[:,0],q[:,1],-3*q[:,0]/norm,-3*q[:,1]/norm,angles='xy',scale_units='xy',scale=1,color='#555',width=.006);ax[0].set(xlim=(-22,22),ylim=(-10,10),xlabel='X (mm)',ylabel='Y (mm)',title='7 mm PZT discs: centers and normals');ax[0].set_aspect('equal')
    for name,col in [('ROT_PLUS',C[0]),('ROT_MINUS',C[1])]:
        rr=[r for r in extra['winding'] if r['frequency']==2500000 and r['action']==name];ax[1].plot([r['radius_m']*1e6 for r in rr],[r['phase_winding'] for r in rr],'o-',label=name,c=col);ax[2].plot([r['radius_m']*1e6 for r in rr],[r['mean_tangential_F1_N']*1e12 for r in rr],'o-',label=name,c=col)
    ax[1].set(xlabel='Diagnostic ring radius (µm)',ylabel='Pressure phase winding');ax[2].set(xlabel='Diagnostic ring radius (µm)',ylabel='Mean tangential F1 (pN)');ax[1].legend();ax[2].axhline(0,c='#888',lw=.7)
    for name,col in [('ROT_PLUS',C[0]),('ROT_MINUS',C[1])]:
        candidates=[r for r in rows if r['status']=='PASS' and r['task'].get('action')==name and r['task']['state']=='N64_wide_r0' and not r['task'].get('variant') and not r['task'].get('offset_m')];r=max(candidates,key=lambda z:z['completed_s']);m=series(r)
        for a,key in zip(grid[1],['Omega_z_s','cumulative_angle_deg','rigid_fraction']):a.plot(m['time_s']*1e3,m[key],label=name,c=col)
    for a,label in zip(grid[1],['Actual Omega_z (rad/s)','Cumulative center motion (deg)','Rigid rotation energy fraction']):a.set(xlabel='Physical time (ms)',ylabel=label);a.axhline(0,c='#aaa',lw=.6)
    grid[1,0].legend(fontsize=8);grid[1,1].text(.04,.93,'Registered rotation threshold: 1 deg',va='top',transform=grid[1,1].transAxes,fontsize=8);grid[1,2].axhline(.8,c=C[2],ls='--');grid[1,2].set_ylim(0,1)
    save(fig,5,'Phase order and bounded rotation diagnostics',['SOURCE_ORDER.json','STATIC_EXTRA.json','series/'])
    fig,g6=plt.subplots(2,3,figsize=(12.8,7.8));ax=g6[0];null=read(BASE/'METRIC_NULL_CALIBRATION.json')['rows']
    ax[0].boxplot([r['null_peak_S'] for r in null],tick_labels=[str(r['N']) for r in null]);ax[0].axhline(12,c=C[2],ls='--',label='Registered threshold');ax[0].set(xlabel='Random particle count',ylabel='Maximum structure factor');ax[0].legend(fontsize=8)
    panel_text(ax[1],'Metric unit tests: PASS\n\nUniform random clouds\nSynthetic stripes and grid\nSynthetic ring\nRigid rotation + translation\nPure shear negative control\n\nSynthetic cases test analysis only.')
    panel_text(ax[2],'F10Y distribution() reused unchanged\n\nCOM, Rg, covariance, aspect\nNearest neighbor and local density\n\nOmega_z: center-motion rotation\nomega_affine = A21 - A12\nPlanar rigid motion: omega = 2 Omega\n\nNo single-bubble spin state.')
    for name,col in [('X_PAIR',C[0]),('XY_QUADRATURE',C[1]),('ROT_MINUS',C[2])]:
        r=max([r for r in rows if r['status']=='PASS' and r['task'].get('action')==name and r['task']['state']=='N64_wide_r0' and not r['task'].get('variant')],key=lambda r:r['completed_s']);m=series(r)
        for a,key,factor in zip(g6[1],['Rg_m','aspect_xy','local_components'],[1e6,1,1]):a.plot(m['time_s']*1e3,m[key]*factor,label=name,c=col)
    for a,label in zip(g6[1],['Actual Rg (µm)','Actual XY aspect ratio','Local components (20 µm cutoff)']):a.set(xlabel='Physical time (ms)',ylabel=label)
    g6[1,0].legend(fontsize=8)
    save(fig,6,'Spatial descriptors and false-positive controls',['METRIC_NULL_CALIBRATION.json','tests/f10af/test_metrics.py'])
    fig=plt.figure(figsize=(12.8,7.8));index=read(BASE/'INITIAL_STATE_INDEX.json')['states'];names=['N32_wide_r0','N64_wide_r0','N200_wide_r0']
    with h5py.File(BASE/'INITIAL_STATES.h5') as h:
        for j,name in enumerate(names):
            a=fig.add_subplot(2,3,j+1,projection='3d');x=h[name][:]*1e6;a.scatter(*x.T,s=9,c=C[j]);a.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',title=f'{len(x)} dispersed bubble centers');a.set_box_aspect((1,.8,.2));a.set_zlim(-100,100)
            a=fig.add_subplot(2,3,j+4);a.scatter(x[:,0],x[:,1],s=10,c=C[j]);a.set(xlabel='X (µm)',ylabel='Y (µm)',title='XY projection',xlim=(-550,550),ylim=(-450,450),aspect='equal')
    save(fig,7,'Frozen legal dispersed initial states',['INITIAL_STATES.h5','INITIAL_STATE_INDEX.json'])
    fig,g8=plt.subplots(2,3,figsize=(12.8,7.8));ax=g8[0];pp=[p for p in sci['paired'] if p.get('max_fixed_ID_m') is not None]
    ax[0].plot([p['duration_s']*1e3 for p in pp],[p['max_fixed_ID_m']*1e9 for p in pp],'o',ms=4);ax[0].set_yscale('symlog',linthresh=.001);fig.text(.3,.043,'Symmetric-log error axis: linear below 0.001 nm; exact zeros retained',fontsize=8);ax[0].axhline(50,c=C[2],ls='--');ax[0].set(xlabel='Same physical duration (ms)',ylabel='Maximum fixed-ID difference (nm)')
    ax[1].scatter([p['HF_wall_s'] for p in pp],[p['LF_wall_s'] for p in pp],c=C[0]);lim=max(p['HF_wall_s'] for p in pp);ax[1].plot([0,lim],[0,lim],'--',c='#999');ax[1].set(xlabel='HF wall time (s)',ylabel='LF wall time (s)')
    panel_text(ax[2],'All original F1 and F2 retained\nFar resistance: exactly A = I\nAdaptive SciPy RK45\nConservative geometric interval guard\nOriginal solver near contact\n\nCompact failure reproduced at 92 us\nNo relaxed overlap / wall limits\n\nLF scope is action- and time-specific.')
    for a,key,factor,label in zip(g8[1],['Rg_error_m','structure_peak_error','HF_fallback_steps'],[1e9,1,1],['Endpoint Rg error (nm)','Endpoint structure-factor error','Original HF fallback steps']):
        a.scatter([p['duration_s']*1e3 for p in pp],[p.get(key,0)*factor for p in pp],s=14,c=C[0]);a.set(xlabel='Same physical duration (ms)',ylabel=label)
    save(fig,8,'Paired slow-time validation',['PAIRED_VALIDATION.json','EXTENSION_RESULTS.json','GENERALIZATION_RESULTS.json','LONG_TIME_RESULTS.json'])
    fig,ax=plt.subplots(1,3,figsize=(12.8,4.5));rr=[r for r in rows if r['status']=='PASS' and r['id'].startswith(('extend_','long')) and not r['task'].get('variant')]
    for name,col in [('X_PAIR',C[0]),('ROT_PLUS',C[1]),('ROT_MINUS',C[2])]:
        a=sorted([r for r in rr if r['task']['action']==name],key=lambda r:r['completed_s']);ax[0].plot([r['completed_s']*1e3 for r in a],[r['max_fixed_ID_displacement_m']*1e6 for r in a],'o-',label=name,c=col);ax[1].plot([r['completed_s']*1e3 for r in a],[r['wall_s'] for r in a],'o-',c=col)
    ax[0].set(xlabel='Completed HF physical time (ms)',ylabel='Maximum displacement (µm)');ax[0].legend(fontsize=8);ax[1].set(xlabel='Completed HF physical time (ms)',ylabel='HF task wall time (s)');r=[by[f'scale_N{n}'] if f'scale_N{n}' in by else by[f'cost_N{n}'] for n in [32,128,200]];ax[2].bar([str(z['N']) for z in r],[z['wall_s'] for z in r],color=C[:3]);ax[2].set(xlabel='Particle count (2 ms attempt)',ylabel='Task wall time (s)');save(fig,9,'Physical time and measured compute cost',['SCIENTIFIC_RESULTS.json'])
    fig,ax=plt.subplots(1,3,figsize=(12.8,4.5));api=read(BASE/'OPTIMIZER_API_REPLAY_V2.json');methods=api['methods'];ax[0].bar(['SMAC3','BoTorch\nqMFKG','Sobol'],[methods[k]['actual_overhead_s'] for k in ['SMAC3','BoTorch_qMFKG','Sobol_replay']],color=C[:3]);ax[0].set(ylabel='Offline API replay overhead (s)',title='Not physical compute savings')
    so=[r for r in rows if r['id'].startswith('sobol_')];ax[1].bar([r['id'].split('_')[-1] for r in so],[r['stripe_order_gain'] for r in so],color=C[0]);ax[1].set(xlabel='Bounded full-physics Sobol trial',ylabel='Actual stripe-order gain')
    panel_text(ax[2],'Formal multi-fidelity search: NOT_RUN\nNo HF-positive pattern yet\nFalse-negative sensitivity unidentifiable\n\nOfficial APIs executed on frozen data\nNo new physics from replay\nNo aggressive early stopping\n\nCMA-MAE superiority: NOT ESTABLISHED')
    save(fig,10,'Optimizer integration and scientific admission',['OPTIMIZER_API_REPLAY_V2.json','SOBOL_RESULTS.json'])
    ref=max([r for r in rows if r['status']=='PASS' and r['task'].get('action')=='X_PAIR' and not r['task'].get('variant')],key=lambda r:r['completed_s']);m=series(ref);fig,ax=plt.subplots(1,3,figsize=(12.8,4.5))
    for j,k in enumerate([0,-1]):ax[j].scatter(m['positions_m'][k,:,0]*1e6,m['positions_m'][k,:,1]*1e6,s=15,c=C[0]);ax[j].set(xlim=(-550,550),ylim=(-450,450),xlabel='X (µm)',ylabel='Y (µm)',title='Initial' if k==0 else f'Actual endpoint: {ref["completed_s"]*1e3:g} ms');ax[j].set_aspect('equal')
    ax[2].plot(m['time_s']*1e3,m['structure_peak_S'],c=C[0],label='Measured peak S');ax[2].axhline(12,c=C[2],ls='--',label='Threshold (also needs gain/periods)');ax[2].set(xlabel='Physical time (ms)',ylabel='Structure factor peak');ax[2].legend(fontsize=8);save(fig,11,'Actual particle pattern evolution',[ref['series_path'],ref['id']])
    fig,ax=plt.subplots(1,3,figsize=(12.8,4.5));conv=sci['convergence'];ax[0].bar(range(len(conv)),[r.get('max_fixed_ID_m',np.nan)*1e9 for r in conv],color=C[0]);ax[0].axhline(50,c=C[2],ls='--');ax[0].set(xlabel='dt / cache refinement comparison',ylabel='Maximum fixed-ID difference (nm)')
    ax[1].bar(['Stripe','Lattice','Rotation'],[len(sci[k]) for k in ['stripe_candidates','lattice_candidates','rotation_candidates']],color=C[:3]);ax[1].set(ylabel='HF candidates meeting criteria',ylim=(0,max(1,max(len(sci[k]) for k in ['stripe_candidates','lattice_candidates','rotation_candidates']))+1))
    panel_text(ax[2],f'HF attempts: {sci["HF_attempts"]}\nComplete: {sci["HF_complete"]}\nFailures retained: {len(sci["failures"])}\n\nStreaming: NOT MODELED\nSingle-bubble spin: NOT MODELED\nReal PDMS field: UNCALIBRATED\nChronos observability: UNVERIFIED\nWorld model / MPC: NOT TRAINED\ntraining_ready = false')
    save(fig,12,'Numerical quality and scientific gates',['SCIENTIFIC_RESULTS.json','AUDIT.json'])
    write(out/'MANIFEST.json',dict(figures=manifest,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
