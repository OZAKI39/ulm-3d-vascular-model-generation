"""Six requested figures from saved shear audit data; no synthetic samples."""
import argparse,csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/nonlocal_boundary_nonaffine_shear_localization_v1'
FIG=ROOT/'results/nonlocal_boundary_nonaffine_shear_localization_v1/figures'
BLUE='#176b9b';ORANGE='#c46522';GREEN='#388368';RED='#a73c50';GRAY='#596777'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,
    'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,
    'pdf.fonttype':42,'ps.fonttype':42,'legend.frameon':False,'savefig.dpi':600})


def read(n):return json.loads((OUT/n).read_text())
def csvread(n):
    with (OUT/n).open() as s:return list(csv.DictReader(s))
def col(rows,k):return np.array([float(r[k]) for r in rows])
def save(fig,name,caption):
    fig.text(.07,.018,caption,fontsize=8.5,color=GRAY,va='bottom')
    fig.savefig(FIG/(name+'.png'),dpi=600,facecolor='white')
    fig.savefig(FIG/(name+'.pdf'),facecolor='white')
    plt.close(fig);print('FIGURE',name,flush=True)


def references():
    n27=read('N27_SHEAR_DEPTH_CHECK.json')['points']
    with np.load(OUT/'NODEWISE_SHEAR_DECOMPOSITION.npz') as a:
        z,iz=np.unique(a['X'][:,2],return_inverse=True)
        fig,ax=plt.subplots(2,2,figsize=(11,8));fig.subplots_adjust(left=.09,right=.98,bottom=.14,top=.84,hspace=.49,wspace=.31)
        channels=[('f_FD','Finite-domain internal',GRAY),('f_bulk_target','Bulk BAQP target',BLUE),
                  ('f_body_local','Local body load',GREEN),('f_BC_localized','Localized conversion',ORANGE),('f_model','Bulk-model residual',RED)]
        for j,component in enumerate([0,2]):
            for key,label,color in channels:
                mean=np.array([a[key][iz==k,component].mean() for k in range(len(z))])
                ax[0,j].plot(z*1e6,mean*1e9,label=label,color=color,lw=1.7)
            ax[0,j].set(title=f'({"ab"[j]}) {"xz"[j]} component: plane-mean nodal forces',xlabel=r'$Z$ [$\mu$m]',ylabel='Nodal force [nN]')
            model=np.array([a['f_model'][iz==k,component].mean() for k in range(len(z))])
            corrected=a['f_FD']+a['f_BC_localized']+a['f_body_local']
            actual=np.array([corrected[iz==k,component].mean() for k in range(len(z))])
            ax[1,j].plot(z*1e6,model,color=RED,lw=1.8,label='Bulk-model residual')
            ax[1,j].plot(z*1e6,actual,'o',mfc='none',color=BLUE,ms=4,label='Corrected residual')
            ax[1,j].set(title=f'({"cd"[j]}) {"xz"[j]} component: retained model residual',xlabel=r'$Z$ [$\mu$m]',ylabel='Plane-mean residual [N]')
            ax[1,j].ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
        handles,labels=ax[0,0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.53,.985),ncol=3,fontsize=9)
        save(fig,'nonaffine_shear_residual_decomposition',
             'N19 physical trapezoid volumes; every point averages one Z plane. Bottom panels resolve the small model residual.\n'
             'Red lines: independently measured bulk-model load. Open blue circles: corrected finite-domain residual.\n'
             'The z model difference is at roundoff scale; it is shown without a near-zero component normalization.')
    rows=csvread('SHEAR_BOUNDARY_LOCALIZATION.csv');groups={n:[r for r in rows if r['correction']==n] for n in ['localized','naive']}
    fig,ax=plt.subplots(1,2,figsize=(11,5));fig.subplots_adjust(left=.09,right=.98,bottom=.25,top=.89,wspace=.30)
    for name,color,label in [('localized',ORANGE,'Localized'),('naive',BLUE,'Naive cancellation')]:
        r=groups[name];x=col(r,'depth_over_delta')
        ax[0].semilogy(x,col(r,'load_mean_N')*1e9,'o-',ms=4,color=color,label=label+' mean')
        ax[0].semilogy(x,col(r,'load_max_N')*1e9,'--',lw=1,color=color,label=label+' max')
        ax[1].semilogy([r['depth_over_delta'] for r in n27],[r[name+'_scaled'] for r in n27],'o-',ms=5,color=color,label=label)
    ax[0].set(title='(a) N19 nearest-face depth profile',xlabel=r'$d_{\min}/\delta$',ylabel='Correction magnitude [nN]')
    ax[1].set(title='(b) N27 lateral depth at Z = 0',xlabel=r'Lateral face distance / $\delta$',ylabel=r'Correction magnitude / $F_{\rm ref}$',xlim=(.9,3.4))
    ax[0].legend(fontsize=8.5,loc='lower left');ax[1].legend(fontsize=9,loc='center left')
    ax[1].scatter(n27[-1]['depth_over_delta'],n27[-1]['localized_scaled'],s=105,facecolors='none',edgecolors=GREEN,zorder=5)
    save(fig,'shear_naive_vs_localized_conversion',
         'N19 does not contain a fully closed force stencil. The circled N27 center is fully closed.\n'
         'Naive cancellation still contains the legitimate bulk shear force there; localized conversion does not.')
    fig,ax=plt.subplots(1,2,figsize=(11,5.8));fig.subplots_adjust(left=.085,right=.98,bottom=.33,top=.78,wspace=.28)
    fig.suptitle('CONVERSION SUPPORT, NOT PHYSICAL BC THICKNESS',y=.97,fontsize=12,weight='bold')
    for name,color,label in [('localized',ORANGE,'Localized'),('naive',BLUE,'Naive')]:
        r=groups[name];x=col(r,'depth_over_delta');y=col(r,'cumulative_L1_fraction')
        for a in ax:
            a.step(x,y,where='post',color=color,label=label,lw=1.8);a.plot(x,y,'o',color=color,ms=3)
    ax[0].set(title='(a) Cumulative nodal magnitude',xlabel=r'$d_{\min}/\delta$',ylabel='Cumulative magnitude fraction',ylim=(0,1.04))
    ax[1].set(title='(b) Tail detail',xlabel=r'$d_{\min}/\delta$',ylabel='Cumulative magnitude fraction',ylim=(.87,1.004),xlim=(.4,2.28))
    for value in [.9,.95,.99]:ax[1].axhline(value,color=GRAY,lw=.8,ls=':')
    ax[0].legend(loc='lower right');qs=read('SHEAR_SUPPORT_QUANTILES.json')
    tableax=fig.add_axes([.13,.11,.75,.12]);tableax.axis('off')
    cells=[[name.capitalize()]+[f'{qs[name][q]["depth_um"]:.1f}' for q in ['d50','d90','d95','d99']] for name in ['localized','naive']]
    tab=tableax.table(cellText=cells,colLabels=['Correction','d50 [µm]','d90 [µm]','d95 [µm]','d99 [µm]'],cellLoc='center',bbox=[0,0,1,1]);tab.auto_set_font_size(False);tab.set_fontsize(9)
    for cell in tab.get_celld().values():cell.set_edgecolor('white')
    save(fig,'shear_boundary_localization_support','Quantiles describe the distributed conversion load for this fixed N19 benchmark; no boundary-layer thickness is inferred.')
    fig,ax=plt.subplots(figsize=(8.5,5.2));fig.subplots_adjust(left=.13,right=.97,bottom=.25,top=.89)
    x=np.array([p['depth_over_delta'] for p in n27]);y=np.array([p['localized_scaled'] for p in n27])
    ax.semilogy(x,y,'o-',color=ORANGE,lw=1.8,ms=6,label='Localized conversion')
    ax.axhline(1e-10,color=RED,ls='--',lw=1.3,label=r'Closed-center gate: $10^{-10}$')
    ax.scatter(x[-1],y[-1],s=125,facecolors='none',edgecolors=GREEN,zorder=5,label='Fully closed center')
    ax.set(xlabel=r'Nearest lateral face distance / $\delta$',ylabel=r'$|f_{\rm BC,localized}|/F_{\rm ref}$',
           title='N27 nonaffine shear: lateral boundary decay',ylim=(2e-16,1),xlim=(.85,3.45))
    ax.set_xticks(x);ax.set_xticklabels(['4h\n0.995δ','8h\n1.990δ','12h\n2.985δ','13h\n3.234δ'])
    ax.legend(loc='upper right',fontsize=9)
    ax.annotate('Small at 12h,\nbut not fully closed',xy=(x[2],y[2]),xytext=(1.65,2e-13),fontsize=9,
                arrowprops=dict(arrowstyle='-',color=GRAY),ha='center')
    save(fig,'n27_nonaffine_shear_boundary_decay','Same polynomial and horizon, fixed Z = 0; an independent shifted N29 stencil supplies the bulk reference.\nExplicit D1/D2/D3 closure, not small force alone, identifies the fully closed center.')


def bias():
    rows=csvread('SHEAR_FIELD_BIAS.csv');summary=read('SHEAR_FIELD_BIAS_SUMMARY.json')
    fig,ax=plt.subplots(1,3,figsize=(13,5.4));fig.subplots_adjust(left=.065,right=.98,bottom=.31,top=.84,wspace=.35)
    for a,name,title in zip(ax[:2],['Fbar','Ftilde'],['(a) Nodal mean gradient','(b) Bond-associated gradient']):
        r=[v for v in rows if v['sample']==name];x=col(r,'depth_over_delta')
        a.semilogy(x,np.maximum(col(r,'F_error_max_scaled'),1e-16),'o-',ms=3,color=BLUE,label='Maximum')
        a.semilogy(x,np.maximum(col(r,'F_error_RMS_scaled'),1e-16),'s--',ms=2.7,color=GREEN,label='Weighted RMS')
        a.set(title=title,xlabel=r'$d_{\min}/\delta$',ylabel=r'Gradient error / 0.004',ylim=(1e-14,.07))
    for name,color,label in [('Fbar',BLUE,r'$\bar F$'),('Ftilde',ORANGE,r'$\widetilde F$')]:
        r=[v for v in rows if v['sample']==name];x=col(r,'depth_over_delta')
        ax[2].semilogy(x,np.maximum(col(r,'Psi_plus_error_max_Pa'),1e-18),'-',color=color,lw=1.8,label=label+' max')
        ax[2].semilogy(x,np.maximum(col(r,'Psi_plus_error_RMS_Pa'),1e-18),'--',color=color,lw=1.3,label=label+' RMS')
    ax[2].set(title=r'(c) Positive energy $\Psi^+$ bias',xlabel=r'$d_{\min}/\delta$',ylabel='Absolute energy-density error [Pa]',ylim=(1e-16,7e-4))
    ax[0].legend(loc='lower left',fontsize=9);ax[2].legend(loc='lower left',fontsize=9)
    fig.suptitle('Exact target positions: reconstruction and positive-energy bias',y=.97,fontsize=13)
    s=summary['samples'];fig.text(.065,.20,
        f'Max gradient bias / 0.004: Fbar {s["Fbar"]["F_error_max_scaled"]:.4f}; Ftilde {s["Ftilde"]["F_error_max_scaled"]:.4f}.  Descriptor: MATERIAL.',fontsize=10)
    fig.text(.065,.15,
        f'Psi+ max / common reference: Fbar {100*s["Fbar"]["Psi_plus_error_max_scaled"]:.2f}%; Ftilde {100*s["Ftilde"]["Psi_plus_error_max_scaled"]:.2f}%.  No fracture acceptance rating.',fontsize=10)
    save(fig,'nonaffine_shear_kinematic_energy_bias','Fbar is compared at nodes; Ftilde at bond midpoints. Depth uses the same respective locations.\n'
         'RMS weights: physical nodal volume, or BAQP bond energy weights. The unchanged spectral split is diagnostic only.')


def recovery():
    fig,ax=plt.subplots(1,2,figsize=(11,5.1));fig.subplots_adjust(left=.09,right=.98,bottom=.29,top=.87,wspace=.30)
    for a,kind,title in zip(ax,['NL','classical'],['(a) BAQP-manufactured body force','(b) Classical manufactured body force']):
        rows=csvread('SHEAR_'+kind+'_RECOVERY_PROFILE.csv');x=col(rows,'depth_over_delta')
        a.semilogy(x,col(rows,'F_recovery_max'),'o-',color=BLUE,ms=4,label='Maximum discrete-target difference')
        a.semilogy(x,col(rows,'F_recovery_RMS'),'s--',color=GREEN,ms=3,label='RMS discrete-target difference')
        a.axhline(1e-5,ls=':',color=GRAY,label='NL recovery gate (reference)')
        a.set(title=title,xlabel=r'$d_{\min}/\delta$',ylabel=r'$|\bar F(x)-\bar F(x_L)|/0.004$',ylim=(1e-12,3e-5))
        result=read('SHEAR_NL_MANUFACTURED_RECOVERY.json' if kind=='NL' else 'SHEAR_CLASSICAL_MANUFACTURED_RESULT.json')
        a.text(.04,.84,f'Force residual: {result["force_residual"]:.2e}',transform=a.transAxes,color=GRAY,fontsize=9)
    handles,labels=ax[0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.53,.19),ncol=3,fontsize=8.5)
    save(fig,'nonaffine_shear_manufactured_recovery','Both solves start from the same smooth nonzero perturbation and hold the boundary conversion fixed.\n'
         'The classical load retains the independently measured bulk-model difference. Analytic-gradient bias is reported separately.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=['references','recovery','all'],default='all');a=p.parse_args()
    if a.stage in ['references','all']:references();bias()
    if a.stage in ['recovery','all']:recovery()
