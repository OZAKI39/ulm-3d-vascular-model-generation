"""Six publication figures from saved audit arrays; no simulations or fitted data."""
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'verification/baqp_near_boundary_kinematic_energy_audit_v1'
FIG = ROOT/'results/baqp_near_boundary_kinematic_energy_audit_v1/figures'
BLUE, ORANGE, GREEN, RED, GRAY = '#176b9b', '#c46522', '#388368', '#a73c50', '#596777'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10, 'axes.titlesize':11,
    'axes.spines.top':False, 'axes.spines.right':False, 'axes.grid':True, 'grid.alpha':.18,
    'pdf.fonttype':42, 'ps.fonttype':42, 'legend.frameon':False, 'savefig.dpi':600})
SAMPLING = {}
LAYOUT_QA = {}


def read(name): return json.loads((OUT/name).read_text())
def rows(name):
    with (OUT/name).open() as s: return list(csv.DictReader(s))
def select(r, **kw): return [v for v in r if all(v[k] == value for k,value in kw.items())]
def col(r,k): return np.array([float(v[k]) for v in r])
def canvas(title):
    fig, ax = plt.subplots(2,2,figsize=(11,8.5))
    fig.subplots_adjust(left=.10, right=.97, bottom=.16, top=.89, hspace=.48, wspace=.33)
    fig.suptitle(title, y=.967, fontsize=14)
    return fig,ax
def save(fig,name,caption):
    fig.text(.10,.045,caption,fontsize=8.7,color=GRAY,va='bottom',linespacing=1.5)
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer();width,height=fig.canvas.get_width_height()
    outside=[]
    texts=list(fig.texts)
    for legend in fig.legends:texts += legend.get_texts()
    for a in fig.axes:
        texts += [a.title,a.xaxis.label,a.yaxis.label]+list(a.texts)
        if a.get_legend():texts += a.get_legend().get_texts()
        for axis,lim in [(a.xaxis,a.get_xlim()),(a.yaxis,a.get_ylim())]:
            texts.append(axis.get_offset_text())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if min(lim)<=tick.get_loc()<=max(lim):texts += [tick.label1,tick.label2]
    for t in texts:
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        if b.x0 < -1 or b.y0 < -1 or b.x1 > width+1 or b.y1 > height+1:outside.append(t.get_text())
    LAYOUT_QA[name]={'canvas_pixels_at_check':[width,height],'text_outside_canvas':outside}
    if outside:raise RuntimeError(f'{name}: text outside canvas: {outside}')
    fig.savefig(FIG/(name+'.png'),dpi=600,facecolor='white')
    fig.savefig(FIG/(name+'.pdf'),facecolor='white')
    plt.close(fig)
    print('FIGURE',name,flush=True)
def sample(n,key,limit=1500):
    stride=max(1,int(np.ceil(n/limit))); ids=np.arange(0,n,stride)
    SAMPLING[key]={'population':n,'indices':'arange(0, population, stride)','stride':stride,'shown':len(ids)}
    return ids
def depth_label(a): a.set_xlabel(r'Nearest-face depth $d/\delta$')
def norm(a): return np.sqrt(np.sum(a*a,axis=tuple(range(1,a.ndim))))


def fbar():
    r=rows('FBAR_DEPTH_PROFILES.csv'); fig,ax=canvas('Third moments explain the near-boundary nodal gradient bias')
    for field,c,marker in [('shear',BLUE,'o'),('axial',ORANGE,'s')]:
        q=select(r,field=field); x=col(q,'depth_over_delta'); label=field.capitalize()
        for a,key in [(ax[0,0],'actual_max_scaled'),(ax[0,1],'predicted_max_scaled')]:
            a.plot(x,col(q,key),color=c,lw=1.4,marker=marker,ms=5,mfc='none',label=label)
        ax[1,0].semilogy(x,col(q,'residual_max_scaled'),color=c,marker=marker,ms=4,mfc='none',label=label+' max')
        ax[1,0].semilogy(x,col(q,'residual_RMS_scaled'),color=c,ls=':',label=label+' RMS')
    ax[0,0].set(title=r'(a) Actual $\|\bar F-F_{exact}\|$',ylabel='Maximum error / 0.004')
    ax[0,1].set(title=r'(b) Independent $M_3 K^{-1}$ prediction',ylabel='Maximum error / 0.004')
    ax[1,0].axhline(1e-10,c=RED,ls='--',lw=1,label='Max gate');ax[1,0].axhline(1e-11,c=GRAY,ls='--',lw=1,label='RMS gate')
    ax[1,0].set(title='(c) Moment-identity residual',ylabel='Residual / 0.004',ylim=(1e-15,2e-10))
    q=select(r,field='shear')
    ax[1,1].plot(col(q,'depth_over_delta'),col(q,'M3_zz_norm_max_m6')/1e-24,'o-',c=GREEN,ms=4)
    ax[1,1].set(title=r'(d) Shared geometry: $\|M_{3,zz:}\|$',ylabel=r'Maximum moment [$10^{-24}$ m$^6$]')
    ax[0,0].legend();ax[0,1].legend();ax[1,0].legend(fontsize=8,ncol=2,loc='lower right')
    for a in ax.flat:depth_label(a)
    save(fig,'third_moment_fbar_causality',
         'All 6,859 physical nodes; depth-bin maxima. Axial z-row and shear x-row map to the same geometric predictor.\n'
         'RMS uses physical nodal volumes. Overlapping curves are expected; no fitted correction is applied.')


def variants():
    m=np.load(OUT/'MOMENT_TENSORS.npz');s=read('SOURCE_DECOMPOSITION.json');l=read('FULL_LATTICE_REFERENCE.json')
    fig,ax=canvas('Separate support truncation from physical volume weighting')
    names=['Continuous\nfull ball','Full\nlattice','Truncated\nuniform V','Physical\ntrapezoid V']
    values=[0,l['M3_norm'],float(norm(m['M3_uniform']).max()),float(norm(m['M3']).max())]
    errs=[0]+[s['variants'][k]['max_scaled'] for k in ['R1_full','R2_support_uniform','R3_physical']]
    for a,v,title,ylabel,threshold in [(ax[0,0],values,r'(a) Third moments',r'Maximum $\|M_3\|$ [m$^6$]',1e-39),
                                        (ax[0,1],errs,'(b) Moment-predicted gradient bias','Maximum error / 0.004',1e-19)]:
        a.bar(np.arange(4),v,color=[GRAY,BLUE,GREEN,ORANGE],width=.55)
        a.set_xticks(np.arange(4),names,fontsize=8.5);a.set_yscale('symlog',linthresh=threshold)
        a.set(title=title,ylabel=ylabel);a.margins(y=.25)
        for k,value in enumerate(v):a.annotate('exact 0' if value==0 else f'{value:.2e}',(k,value),xytext=(0,5),textcoords='offset points',ha='center',fontsize=8)
    ax[0,0].set_ylim(-1e-40,1e-21);ax[0,0].set_yticks([0,1e-38,1e-34,1e-30,1e-26,1e-22])
    ax[0,1].set_ylim(-1e-20,.3);ax[0,1].set_yticks([0,1e-18,1e-14,1e-10,1e-6,1e-2])
    d=np.rint(m['distance_face_m']/6.25e-5).astype(int);levels=np.unique(d);x=levels/4.02
    for key,gkey,c,label in [('M3_uniform','G_uniform',GREEN,'R2: support only'),('M3','G_physical',ORANGE,'R3: physical volumes')]:
        yy=norm(m[key]);gg=norm(m[gkey])*.001/.001125/.004
        ax[1,0].plot(x,[yy[d==z].max()/1e-24 for z in levels],'o-',c=c,ms=4,label=label)
        ax[1,1].plot(x,[gg[d==z].max() for z in levels],'o-',c=c,ms=4,label=label)
    ax[1,0].set(title='(c) Third moment vs physical depth',ylabel=r'Maximum $\|M_3\|$ [$10^{-24}$ m$^6$]')
    ax[1,1].set(title='(d) Predicted bias vs physical depth',ylabel='Maximum error / 0.004')
    for a in ax[1]:depth_label(a);a.legend(fontsize=9)
    save(fig,'full_vs_truncated_moment_comparison',
         'Upper panels use a symmetric-log scale so analytical zero and lattice roundoff remain distinguishable.\n'
         'R2 and R3 use identical retained neighbors. The vector volume contribution is R3 minus R2; maxima are not additive.')


def ftilde():
    r=rows('FTILDE_PROJECTION_CLOSURE.csv');o=select(rows('BOND_ORIENTATION_PROFILES.csv'),field='shear')
    b=np.load(OUT/'SHEAR_BOND_DIAGNOSTICS.npz');fig,ax=canvas('Bond-associated bias is a transverse projection of nodal bias')
    ids=sample(len(b['actual_F_error']),'ftilde_projection_bias:panel_a')
    ax[0,0].scatter(b['predicted_F_error'][ids]/.004,b['actual_F_error'][ids]/.004,s=9,c=BLUE,alpha=.4,edgecolors='none')
    lim=b['actual_F_error'].max()/.004*1.03;ax[0,0].plot([0,lim],[0,lim],c=GRAY,ls='--',lw=1)
    ax[0,0].set(title='(a) Shear: actual vs projected error',xlabel='Predicted norm / 0.004',ylabel='Actual norm / 0.004')
    for field,c in [('shear',BLUE),('axial',ORANGE)]:
        q=select(r,field=field);x=col(q,'depth_over_delta')
        ax[0,1].semilogy(x,col(q,'residual_max_scaled'),'o-',c=c,ms=3,label=field.capitalize()+' max')
        ax[0,1].semilogy(x,col(q,'residual_RMS_scaled'),':',c=c,label=field.capitalize()+' RMS')
    ax[0,1].axhline(1e-10,c=RED,ls='--',lw=1);ax[0,1].axhline(1e-11,c=GRAY,ls='--',lw=1)
    ax[0,1].set(title='(b) Exact projection-identity residual',xlabel=r'Midpoint depth $d/\delta$',ylabel='Residual / 0.004',ylim=(1e-16,2e-10))
    ax[0,1].legend(fontsize=8,ncol=2)
    bins=['[0,.25)','[.25,.5)','[.5,.75)','[.75,.95)','[.95,1]'];x=np.arange(5)
    ax[1,0].plot(x,col(o,'actual_F_max_scaled'),'o-',c=BLUE,label='Actual max')
    ax[1,0].plot(x,col(o,'predicted_F_max_scaled'),'s',mfc='none',c=ORANGE,label='Projected max')
    ax[1,0].plot(x,col(o,'actual_F_RMS_scaled'),'--',c=GREEN,label='Actual weighted RMS')
    ax[1,0].set_xticks(x,bins,fontsize=8)
    ax[1,0].set(title=r'(c) Shear: orientation bins',xlabel=r'Bond orientation $|n_z|$',ylabel='Gradient error / 0.004')
    ax[1,0].legend(fontsize=8.5)
    q=select(r,field='shear');x=col(q,'depth_over_delta')
    ax[1,1].plot(x,col(q,'actual_max_scaled'),'o-',ms=4,c=BLUE,label='Actual max')
    ax[1,1].plot(x,col(q,'predicted_max_scaled'),'s',ms=4,mfc='none',c=ORANGE,label='Projected max')
    ax[1,1].plot(x,col(q,'actual_RMS')/.004,'--',c=GREEN,label='Actual weighted RMS')
    ax[1,1].set(title='(d) Shear: midpoint-depth profile',xlabel=r'Midpoint depth $d/\delta$',ylabel='Gradient error / 0.004');ax[1,1].legend(fontsize=8.5)
    save(fig,'ftilde_projection_bias',
         r'Prediction: $\Delta\widetilde F=D(I-n\otimes n)$, with $D=(\Delta\bar F_i+\Delta\bar F_j)/2$.'+'\n'
         'All 690,020 bonds enter gates and profiles; panel (a) uses deterministic stride sampling. RMS uses BAQP energy weights.')


def positive():
    r=rows('PSI_BIAS_RECONSTRUCTION.csv');summary=read('CAUSAL_IDENTITY_SUMMARY.json');ref=summary[0]['Psi_ref_Pa']
    fig,ax=canvas('The unchanged nonlinear spectral split reproduces the positive-energy bias')
    for j,(kind,file,c) in enumerate([('Fbar','SHEAR_NODE_DIAGNOSTICS.npz',BLUE),('Ftilde','SHEAR_BOND_DIAGNOSTICS.npz',ORANGE)]):
        q=select(r,field='shear',sample=kind);x=col(q,'depth_over_delta')
        ax[0,j].plot(x,col(q,'actual_max_scaled'),'o-',c=c,ms=4,label='Actual max absolute bias')
        ax[0,j].plot(x,col(q,'predicted_max_scaled'),'s',c=GRAY,mfc='none',ms=5,label='Moment-predicted max')
        ax[0,j].set(title=f'({"ab"[j]}) Shear: {kind} energy-density bias',xlabel=r'Depth $d/\delta$',ylabel=r'Maximum $|\Delta\Psi^+|/\Psi_{ref}$')
        ax[0,j].legend(fontsize=8.5)
        a=np.load(OUT/file);ids=sample(len(a['Psi_actual_bias']),'psi_plus_bias_causal_reconstruction:'+kind)
        xp=a['Psi_predicted_bias'][ids]/ref;yp=a['Psi_actual_bias'][ids]/ref
        ax[1,j].scatter(xp,yp,s=10,c=c,alpha=.4,edgecolors='none')
        lo=min(xp.min(),yp.min());hi=max(xp.max(),yp.max());ax[1,j].plot([lo,hi],[lo,hi],c=GRAY,ls='--',lw=1)
        ax[1,j].set(title=f'({"cd"[j]}) {kind}: signed bias, sampled',xlabel=r'Predicted $\Delta\Psi^+/\Psi_{ref}$',ylabel=r'Observed $\Delta\Psi^+/\Psi_{ref}$')
    save(fig,'psi_plus_bias_causal_reconstruction',
         r'Shear $\Psi_{ref}=0.004016032032$ Pa, common to nodes and bonds. Prediction starts from independent $M_3 K^{-1}$ errors.'+'\n'
         'Full nonlinear spectral evaluation; no linearized sensitivity. All samples enter gates; scatter uses deterministic strides.')


def energy():
    e=read('ENERGY_SOURCE_DECOMPOSITION.json');fig,ax=canvas('Global energy: kinematic bias and weight / quadrature discrepancy')
    fig.subplots_adjust(top=.85)
    for j,field in enumerate(['shear','axial']):
        for shift,channel,c,label in [(-.18,'intact',BLUE,'Intact W'),(.18,'positive',ORANGE,r'Positive $\Psi^+$')]:
            v=e[field][channel];x=np.arange(3)+shift
            val=np.array([v[k+'_J'] for k in ['U0','U1','U2']])/1e-12
            ax[j,0].bar(x,val,width=.34,color=c,label=label)
            for xx,yy in zip(x,val):ax[j,0].text(xx,yy+.04*val.max(),f'{yy:.5f}',ha='center',va='bottom',fontsize=7.5,rotation=90)
            keys=['kinematic_bias_J','weight_quadrature_bias_J','total_bias_J'];delta=np.array([v[k] for k in keys])/1e-14
            ax[j,1].bar(x,delta,width=.34,color=c,label=label)
        ax[j,0].set_xticks(np.arange(3),[r'$U_0$',r'$U_1$',r'$U_2$']);ax[j,0].set_ylim(0,ax[j,0].get_ylim()[1]*1.33)
        ax[j,0].set(title=f'({"ac"[j]}) {field.capitalize()}: total energies',ylabel=r'Energy [$10^{-12}$ J]')
        ax[j,1].set_xticks(np.arange(3),['Kinematics\n'+r'$U_0-U_1$','Weights / quadrature\n'+r'$U_1-U_2$','Total\n'+r'$U_0-U_2$'],fontsize=8.5)
        ax[j,1].set(title=f'({"bd"[j]}) {field.capitalize()}: signed discrepancies',ylabel=r'Energy difference [$10^{-14}$ J]')
        ax[j,1].axhline(0,c=GRAY,lw=.8)
    handles,labels=ax[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.52,.935),ncol=2,fontsize=10)
    save(fig,'energy_error_source_decomposition',
         r'$U_0$: current BAQP. $U_1$: exact midpoint kinematics with unchanged bond weights. $U_2$: continuum Gauss integral.'+'\n'
         r'$U_1$ is diagnostic only. Absolute joule scales are shown; axial intact and positive energies nearly coincide.')


def invariance():
    a=np.load(OUT/'AXIAL_NODE_DIAGNOSTICS.npz');s=np.load(OUT/'SHEAR_NODE_DIAGNOSTICS.npz')
    ar=a['actual_delta_F'][:,2,:];sr=s['actual_delta_F'][:,0,:];pr=s['predicted_delta_F'][:,0,:]
    d=np.rint(s['depth_over_delta']*4.02).astype(int);levels=np.unique(d);x=levels/4.02
    fig,ax=canvas('Axial and shear biases share one geometric predictor')
    for v,c,marker,label in [(sr,BLUE,'o','Shear x-row'),(ar,ORANGE,'s','Axial z-row'),(pr,GREEN,'^','Moment prediction')]:
        y=np.linalg.norm(v,axis=1)/.004
        ax[0,0].plot(x,[y[d==z].max() for z in levels],marker=marker,ms=5,c=c,mfc='none',lw=1,label=label)
    ax[0,0].set(title='(a) Mapped row-norm depth profile',xlabel=r'Nearest-face depth $d/\delta$',ylabel='Maximum row error / 0.004');ax[0,0].legend(fontsize=9)
    X=s['X'];ids=np.flatnonzero((abs(X[:,0])<1e-12)&(abs(X[:,1])<1e-12));ids=ids[np.argsort(X[ids,2])]
    for v,c,marker,label in [(sr,BLUE,'o','Shear xz'),(ar,ORANGE,'s','Axial zz'),(pr,GREEN,'^','Moment prediction')]:
        ax[0,1].plot(X[ids,2]*1e6,v[ids,2]/.004,marker=marker,c=c,ms=5,mfc='none',lw=1,label=label)
    ax[0,1].set(title='(b) Signed component along the central Z line',xlabel=r'$Z$ [$\mu$m]',ylabel='Signed component error / 0.004');ax[0,1].legend(fontsize=9)
    ids=sample(len(sr),'axial_shear_geometric_bias_invariance:panel_c',1200)
    for k,c,label in [(0,BLUE,'x'),(1,ORANGE,'y'),(2,GREEN,'z')]:
        ax[1,0].scatter(sr[ids,k]/.004,ar[ids,k]/.004,s=8,c=c,alpha=.35,edgecolors='none',label='Input '+label)
    lim=max(abs(sr).max(),abs(ar).max())/.004*1.03;ax[1,0].plot([-lim,lim],[-lim,lim],c=GRAY,ls='--',lw=1)
    ax[1,0].set(title='(c) Corresponding signed row components',xlabel='Shear x-row component / 0.004',ylabel='Axial z-row component / 0.004');ax[1,0].legend(fontsize=8)
    diff=np.linalg.norm(sr-ar,axis=1)/.004
    ax[1,1].semilogy(x,[diff[d==z].max() for z in levels],'o-',c=BLUE,ms=4,label='Actual mapped difference')
    ax[1,1].axhline(1e-10,c=RED,ls='--',label=r'Gate: $10^{-10}$')
    ax[1,1].set(title='(d) Difference after output-row mapping',xlabel=r'Nearest-face depth $d/\delta$',ylabel='Maximum difference / 0.004',ylim=(1e-13,2e-10));ax[1,1].legend(fontsize=9)
    save(fig,'axial_shear_geometric_bias_invariance',
         r'Both fields have $H_{azz}=2\beta/L$: $a=x$ for shear and $a=z$ for axial. The same $M_{3,zz:}K^{-1}$ drives both.'+'\n'
         'Panel (b) retains signs instead of averaging opposite faces. Full-population checks pass; panel (c) is stride sampled.')


if __name__=='__main__':
    FIG.mkdir(parents=True,exist_ok=True)
    for plot in [fbar,variants,ftilde,positive,energy,invariance]: plot()
    (OUT/'FIGURE_SAMPLING.json').write_text(json.dumps(SAMPLING,indent=2)+'\n')
    (OUT/'FIGURE_LAYOUT_QA.json').write_text(json.dumps(LAYOUT_QA,indent=2)+'\n')
