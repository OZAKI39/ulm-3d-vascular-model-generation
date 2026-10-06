"""Publication figures from saved numerical arrays; no generated/synthetic data."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'verification/nonlocal_boundary_nonaffine_localization_v1'
FIG=ROOT/'results/nonlocal_boundary_nonaffine_localization_v1/figures'
BLUE='#176b9b';ORANGE='#c46522';GREEN='#388368';RED='#a73c50';GRAY='#596777'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,
                     'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.18,
                     'pdf.fonttype':42,'ps.fonttype':42,'legend.frameon':False,'savefig.dpi':600})


def read(n):return json.loads((OUT/n).read_text())
def csvread(n):
    with (OUT/n).open() as s:return list(csv.DictReader(s))
def col(rows,k):return np.array([float(r[k]) for r in rows])
def save(fig,name,caption):
    fig.text(.06,.018,caption,fontsize=8.5,color=GRAY,va='bottom')
    fig.savefig(FIG/(name+'.png'),dpi=600,facecolor='white')
    fig.savefig(FIG/(name+'.pdf'),facecolor='white')
    plt.close(fig);print('FIGURE',name,flush=True)


def references():
    ref=read('EXACT_NONAFFINE_REFERENCE.json');n27=read('N27_DEPTH_CHECK.json')['points']
    with np.load(OUT/'NODEWISE_RESIDUAL_DECOMPOSITION.npz') as a:
        z,iz=np.unique(a['X'][:,2],return_inverse=True)
        fig,ax=plt.subplots(2,2,figsize=(11,7.3));fig.subplots_adjust(left=.085,right=.98,bottom=.15,top=.93,hspace=.48,wspace=.30)
        channels=[('f_FD','Finite-domain internal',GRAY),('f_bulk_target','Bulk BAQP target',BLUE),
                  ('f_body_local','Local body load',GREEN),('f_BC_localized','Localized conversion',ORANGE),('f_model','Bulk-model residual',RED)]
        for key,label,c in channels:
            mean=np.array([a[key][iz==k,2].mean() for k in range(len(z))])
            ax[0,0].plot(z*1e6,mean*1e9,label=label,color=c,lw=1.8)
        ax[0,0].set(title='(a) N19 plane-mean axial nodal forces',xlabel=r'$Z$ [$\mu$m]',ylabel='Axial nodal force [nN]')
        ax[0,0].legend(fontsize=8.5,ncol=1,loc='upper center')
        q=np.array([a['q_bulk'][iz==k,2].mean() for k in range(len(z))])
        b=np.array([a['b_local'][iz==k,2].mean() for k in range(len(z))])
        ax[0,1].plot(z*1e6,q,color=BLUE,lw=2,label=r'$q_{\rm bulk,z}$')
        ax[0,1].plot(z*1e6,-b,color=GREEN,lw=1.5,ls='--',label=r'$-b_{\rm local,z}=\mathrm{Div}\,P_z$')
        ax[0,1].set(title='(b) Bulk force densities nearly coincide',xlabel=r'$Z$ [$\mu$m]',ylabel=r'Force density [N/m$^3$]');ax[0,1].legend()
        ax[1,0].plot(z*1e6,(q+b)*1e5,'o-',color=RED,ms=3)
        ax[1,0].set(title='(c) Small nonzero bulk-model difference',xlabel=r'$Z$ [$\mu$m]',ylabel=r'$e_{\rm model,z}$ [$10^{-5}$ N/m$^3$]',ylim=(0,8.6))
        ax[1,0].ticklabel_format(axis='y',style='plain',useOffset=False)
        modelmean=np.array([a['f_model'][iz==k,2].mean() for k in range(len(z))])
        corrected=a['f_FD']+a['f_BC_localized']+a['f_body_local']
        actualmean=np.array([corrected[iz==k,2].mean() for k in range(len(z))])
        ax[1,1].plot(z*1e6,modelmean*1e18,color=RED,lw=2,label=r'$e_{\rm model}V_i$')
        ax[1,1].plot(z*1e6,actualmean*1e18,'o',mfc='none',color=BLUE,ms=4,label='Corrected residual')
        ax[1,1].set(title='(d) Correction preserves bulk-model residual',xlabel=r'$Z$ [$\mu$m]',ylabel='Plane-mean axial residual [aN]');ax[1,1].legend(fontsize=9)
        save(fig,'nonaffine_residual_decomposition','N19 physical trapezoid volumes; plane averages include boundary-node weights.\nThe small nonzero model residual is retained, not absorbed into the boundary load.')
    rows=csvread('BOUNDARY_LOCALIZATION.csv');groups={name:[r for r in rows if r['correction']==name] for name in ['localized','naive']}
    fig,ax=plt.subplots(1,2,figsize=(11,4.8));fig.subplots_adjust(left=.085,right=.98,bottom=.23,top=.89,wspace=.29)
    for name,color,label in [('localized',ORANGE,'Localized'),('naive',BLUE,'Naive cancellation')]:
        r=groups[name];x=col(r,'depth_over_delta')
        ax[0].semilogy(x,col(r,'load_mean_N')*1e9,'o-',ms=4,color=color,label=label+' mean')
        ax[0].semilogy(x,col(r,'load_max_N')*1e9,'--',lw=1,color=color,label=label+' max')
        ax[1].semilogy([r['depth_over_delta'] for r in n27],[r[name+'_scaled'] for r in n27],'o-',ms=5,color=color,label=label)
    ax[0].set(title='(a) N19: no fully closed center stencil',xlabel=r'Nearest-face distance $d_{\min}/\delta$',ylabel='Correction magnitude [nN]');ax[0].legend(fontsize=8.5)
    ax[1].set(title='(b) N27: lateral depth at Z = 0',xlabel=r'Nearest-face distance $d_{\min}/\delta$',ylabel=r'Correction magnitude / $F_{\rm ref}$')
    ax[1].axvspan(3.16,3.30,color=GREEN,alpha=.10);ax[1].set_xlim(.9,3.3);ax[1].legend(loc='center left',fontsize=9)
    ax[1].text(.97,.10,'Closed center',transform=ax[1].transAxes,ha='right',color=GREEN,fontsize=9)
    save(fig,'naive_vs_localized_boundary_conversion','Only the localized correction vanishes at the fully closed N27 center.\nNaive cancellation retains the legitimate nonzero bulk force. A small N19 center residual does not establish closure.')
    fig,ax=plt.subplots(1,2,figsize=(11,5.2));fig.subplots_adjust(left=.08,right=.98,bottom=.32,top=.89,wspace=.26)
    for name,color,label in [('localized',ORANGE,'Localized'),('naive',BLUE,'Naive')]:
        r=groups[name];x=col(r,'depth_over_delta');y=col(r,'cumulative_L1_fraction')
        for a in ax:a.step(x,y,where='post',color=color,label=label,lw=1.8);a.plot(x,y,'o',color=color,ms=3)
    ax[0].set(title='(a) Cumulative nodal force magnitude',xlabel=r'$d_{\min}/\delta$',ylabel='Cumulative magnitude fraction',ylim=(0,1.04))
    ax[1].set(title='(b) Outer-tail detail',xlabel=r'$d_{\min}/\delta$',ylabel='Cumulative magnitude fraction',ylim=(.87,1.004),xlim=(.4,2.28))
    for value in [.9,.95,.99]:ax[1].axhline(value,color=GRAY,lw=.8,ls=':')
    ax[0].legend(loc='lower right');qs=read('SUPPORT_QUANTILES.json')
    tableax=fig.add_axes([.14,.11,.73,.12]);tableax.axis('off')
    cells=[[name.capitalize()]+[f'{qs[name][q]["depth_um"]:.1f}' for q in ['d50','d90','d95','d99']] for name in ['localized','naive']]
    tab=tableax.table(cellText=cells,colLabels=['Correction','d50 [µm]','d90 [µm]','d95 [µm]','d99 [µm]'],cellLoc='center',bbox=[0,0,1,1])
    tab.auto_set_font_size(False);tab.set_fontsize(9)
    for cell in tab.get_celld().values():cell.set_edgecolor('white')
    save(fig,'boundary_localization_cumulative_support','CONVERSION SUPPORT — NOT PHYSICAL BC THICKNESS. Quantiles summarize the finite N19 benchmark only.')
    fig,ax=plt.subplots(figsize=(8,4.8));fig.subplots_adjust(left=.13,right=.97,bottom=.23,top=.89)
    x=np.array([p['depth_over_delta'] for p in n27]);y=np.array([p['localized_scaled'] for p in n27])
    ax.semilogy(x,y,'o-',color=ORANGE,lw=1.8,ms=6,label='Localized correction')
    ax.axhline(1e-10,color=RED,ls='--',lw=1.3,label=r'Closed-center gate: $10^{-10}$')
    ax.scatter(x[-1],y[-1],s=100,facecolors='none',edgecolors=GREEN,zorder=5)
    ax.set(xlabel=r'Nearest-face distance $d_{\min}/\delta$',ylabel=r'$|f_{\rm BC,localized}|/F_{\rm ref}$',
           title='N27 nonaffine boundary decay: unchanged BAQP operator',ylim=(5e-17,1e-3),xlim=(.85,3.4))
    ax.set_xticks(x);ax.set_xticklabels(['4h\n0.995δ','8h\n1.990δ','12h\n2.985δ','13h\n3.234δ'])
    ax.legend(loc='upper right',fontsize=9)
    ax.annotate('Numerically small,\nbut not fully closed',xy=(x[2],y[2]),xytext=(1.6,2e-12),fontsize=9,
                arrowprops=dict(arrowstyle='-',color=GRAY),ha='center')
    ax.annotate('Fully closed\ncenter',xy=(x[-1],y[-1]),xytext=(3.17,2e-16),fontsize=9,color=GREEN,ha='center')
    save(fig,'n27_nonaffine_boundary_decay','Lateral x-depth samples at Z = 0; independent N29 reference and explicit D1/D2/D3 closure checks.\nNumerical cancellation at 12h does not establish a sufficient physical boundary-layer thickness.')


def recovery():
    fig,ax=plt.subplots(1,2,figsize=(11,5));fig.subplots_adjust(left=.085,right=.98,bottom=.30,top=.87,wspace=.30)
    for a,kind,title in zip(ax,['NL','classical'],['(a) BAQP-manufactured body force','(b) Classical manufactured body force']):
        rows=csvread(kind+'_RECOVERY_PROFILE.csv');x=col(rows,'depth_over_delta')
        a.semilogy(x,col(rows,'F_recovery_max'),'o-',color=BLUE,ms=4,label=r'Max $|\bar F(x)-\bar F(x_L)|/0.004$')
        a.semilogy(x,col(rows,'F_recovery_RMS'),'s--',color=GREEN,ms=3,label='RMS recovery difference')
        a.axhline(1e-5,ls=':',color=GRAY,label='NL gate (reference)')
        a.set(title=title,xlabel=r'$d_{\min}/\delta$',ylabel='Discrete target-gradient difference',ylim=(1e-12,3e-5))
        result=read('NL_MANUFACTURED_RECOVERY.json' if kind=='NL' else 'CLASSICAL_MANUFACTURED_RESULT.json')
        a.text(.04,.84,f'Force residual: {result["force_residual"]:.2e}',transform=a.transAxes,color=GRAY,fontsize=9)
    handles,labels=ax[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.53,.195),ncol=3,fontsize=9)
    response=read('CLASSICAL_VS_NL_ACCOUNTING.json')
    save(fig,'manufactured_solution_recovery','Both solves start from the same nonzero displacement perturbation; the boundary conversion is frozen.\nErrors use the discrete target gradient. Its analytic-gradient reconstruction bias is reported separately.\n'+f'Force-response closure: {100*response["force_response_closure_relative_to_model"]:.3f}% of the model load; NL displacement error: {100*response["NL_recovery_error_over_classical_deviation"]:.2f}% of the classical shift.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['references','recovery','all'],default='all');args=ap.parse_args()
    if args.stage in ['references','all']:references()
    if args.stage in ['recovery','all']:recovery()
