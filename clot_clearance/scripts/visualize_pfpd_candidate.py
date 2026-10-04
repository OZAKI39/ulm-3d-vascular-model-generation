"""Seven review figures, black background / Arial regular / no arrows."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'verification/pfpd_candidate_v1'
os.environ.setdefault('MPLCONFIGDIR',str(O/'provenance/runtime_cache/matplotlib'))
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pd_clot.pfpd_candidate.degradation import energetic,kinematic
font=Path('/mnt/c/Windows/Fonts/arial.ttf')
if font.exists():font_manager.fontManager.addfont(str(font));family='Arial'
else:family='DejaVu Sans'
plt.rcParams.update({'font.family':family,'font.weight':'normal','axes.titleweight':'normal','font.size':13,'axes.titlesize':17,'axes.labelsize':14,'figure.facecolor':'black','axes.facecolor':'black','savefig.facecolor':'black','text.color':'white','axes.labelcolor':'white','xtick.color':'white','ytick.color':'white','axes.edgecolor':'#717171','grid.color':'#333333','legend.frameon':False,'pdf.fonttype':42,'figure.constrained_layout.use':True})
colors=['#ecae48','#45bde5','#dddddd','#e777aa'];figdir=O/'figures'

def read(p):return json.loads((O/p).read_text())
def finish(fig,name):
    if name in ['historical_elastic_energy_change','griffith_normalization_convergence','formula_consistency']:
        fig.set_layout_engine('constrained',rect=(0,.065,1,.935))
    for ax in fig.axes:
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.45);ax.set_axisbelow(True)
    for ext in ['png','pdf']:fig.savefig(figdir/f'{name}.{ext}',dpi=200)
    plt.close(fig)

rows=read('historical_compatibility/ALL_35_JUMPS.json');summary=read('historical_compatibility/SUMMARY.json');smax=summary['max_D']
fig,axs=plt.subplots(1,2,figsize=(12.8,5),width_ratios=[1.5,1])
for ax,s in [(axs[0],np.linspace(0,1,1001)),(axs[1],np.linspace(0,.004,301))]:
    ax.plot(s,energetic(s),color=colors[0],lw=3,label='Energetic g(s)')
    ax.plot(s,kinematic(s,.95),color=colors[1],lw=2,label='Kinematic h, sc = 0.95')
    ax.plot(s,kinematic(s,1),color=colors[2],ls='--',lw=1.6,label='Kinematic h, sc = 1')
    ax.plot(s,kinematic(s,0),color=colors[3],ls=':',lw=2,label='Kinematic h, sc = 0 (= g)')
    ax.axvspan(0,smax,color='#487a55',alpha=.45);ax.set_xlabel('Prescribed bond phase proxy s');ax.set_ylabel('Degradation factor')
axs[0].set_title('Energetic and kinematic degradation');axs[1].set_title('Current historical range');axs[0].legend(fontsize=11,loc='lower left');axs[1].set_ylim(.991,1.001)
axs[1].text(.03,.05,f'Historical s_proxy <= {smax:.6f}\nKinematic factors remain exactly 1',transform=axs[1].transAxes,fontsize=11)
finish(fig,'energetic_vs_kinematic_degradation')

constants=read('normalization/NORMALIZATION_CONSTANTS.json');fig,ax=plt.subplots(figsize=(9.5,5.3));xx=np.arange(3)
ax.bar(xx-.15,[r['analytic'] for r in constants],.3,color=colors[1],label='Analytic')
ax.bar(xx+.15,[r['c0'] for r in constants],.3,color=colors[0],label='Numerical quadrature')
ax.set_xticks(xx,['Constant','Linear','Current quadratic']);ax.set_ylabel('Normalization constant c0');ax.set_ylim(0,.45);ax.set_title('Kernel normalization');ax.legend(loc='upper right')
for i,r in enumerate(constants):ax.text(i,r['c0']+.016,f"{r['c0']:.3f}",ha='center')
finish(fig,'normalization_constant')

aff=read('baqp_verification/MECHANICS_GATE.json')['affine'];fig,ax=plt.subplots(figsize=(11,5.4));xx=np.arange(4)
for n,(key,label) in enumerate([('bulk_Fbar_error','Point F, bulk'),('boundary_Fbar_error','Point F, boundary'),('bulk_Ftilde_error','Bond F, bulk'),('boundary_Ftilde_error','Bond F, boundary')]):
    ax.bar(xx+(n-1.5)*.18,[r[key] for r in aff],.18,color=colors[n],label=label)
ax.set_yscale('log');ax.set_ylim(1e-19,1e-10);ax.axhline(5e-12,color='white',ls=':',label='Declared tolerance')
ax.set_xticks(xx,[r['name'].capitalize() for r in aff]);ax.set_ylabel('Maximum absolute gradient error');ax.set_title('Intact affine reproduction');ax.legend(fontsize=11,ncol=2,loc='upper right',bbox_to_anchor=(1,.80))
finish(fig,'baqp_affine_reproduction')

labels=['Current production','Separated g, h = 1','Delayed h, sc = 0.95','Diagnostic h = g'];x=np.arange(1,36)
fig,axs=plt.subplots(1,2,figsize=(13,5.4))
for ax,field,title in zip(axs,['F','K'],['Deformation gradient change','Shape tensor change']):
    for k,label,color,style in zip('ABCD',labels,colors,['-','-','--','-']):
        y=[r[k+'_max_relative_delta_'+field] for r in rows]
        ax.plot(x,y,color=color,ls=style,lw=2,label=label,marker='o' if k=='C' else None,markersize=3,markerfacecolor='black')
    ax.set_yscale('symlog',linthresh=1e-16);ax.set_xlabel('Historical jump index');ax.set_ylabel('Maximum relative change');ax.set_title(title)
axs[0].legend(fontsize=10,loc='center right');axs[1].text(.05,.09,'B and C: exact zero for all 35 jumps',transform=axs[1].transAxes,fontsize=11)
fig.suptitle('Kinematic response to prescribed damage',fontsize=18);finish(fig,'historical_KF_change')

fig,ax=plt.subplots(figsize=(12.5,5.9))
for k,label,color in zip('ABCD',labels,colors):
    for lo,hi in [(0,20),(20,30),(30,35)]:
        ax.plot(x[lo:hi],np.array([r[k+'_delta_U_J'] for r in rows[lo:hi]])/1e-18,color=color,lw=2,ls='--' if k=='C' else '-',marker='o' if k=='C' else None,markersize=4,markerfacecolor='black',label=label if lo==0 else None)
ax.axhline(0,color='#888888',lw=1);ax.set_xlabel('Historical jump index (20 + 10 + 5 saved jumps)');ax.set_ylabel('Elastic energy change (new - old), 10^-18 J');ax.set_title('Elastic response at fixed positions');ax.legend(fontsize=11,ncol=2,loc='lower left')
fig.text(.5,.005,'KINEMATIC / ELASTIC COMPATIBILITY TEST ONLY',ha='center',fontsize=12)
for cut in [20.5,30.5]:ax.axvline(cut,color='#444444',ls=':')
finish(fig,'historical_elastic_energy_change')

norm=read('normalization/GRIFFITH_STUDY.json');fig,axs=plt.subplots(1,2,figsize=(13,5.6))
for m,col in zip([3.01,5.01,8.01],colors):
    rr=[r for r in norm['coupon'] if r['horizon_ratio']==m]
    axs[0].plot([r['spacing_m']*1e6 for r in rr],[r['ratio'] for r in rr],'o-',color=col,label=f'Horizon / spacing = {m}')
    axs[0].axhline(rr[0]['bulk_lattice_ratio'],color=col,ls=':',alpha=.6)
rr=norm['bulk_lattice_refinement'];axs[1].plot([r['horizon_ratio'] for r in rr],[r['ratio'] for r in rr],'o-',color=colors[1])
for ax in axs:ax.axhline(1,color='white',ls='--',label='Continuum target');ax.set_ylabel('Crack energy / (Gc x crack area)')
axs[0].set_xlabel('Point spacing (micrometers)');axs[0].set_title('Finite coupon, free edges');axs[0].invert_xaxis();axs[0].legend(fontsize=10,loc='lower left')
axs[1].set_xlabel('Horizon / spacing');axs[1].set_title('Bulk lattice resolution');axs[1].text(.32,.07,'No surface correction\nGc = 1 J/m², demonstration only',transform=axs[1].transAxes,fontsize=11)
fig.suptitle('Griffith energy normalization',fontsize=18);fig.text(.5,.005,'NUMERICAL NORMALIZATION TEST',ha='center',fontsize=12);finish(fig,'griffith_normalization_convergence')

Y=np.linspace(0,5,801);at1=np.maximum(0,1-1/np.maximum(Y,np.finfo(float).tiny));eq=Y/(Y+1)
np.savetxt(O/'formula_consistency/FORMULA_COMPARISON_ONLY.csv',np.column_stack([Y,at1,eq]),delimiter=',',header='Y_over_Yc,s_AT1_constrained,s_Eq325',comments='')
fig,ax=plt.subplots(figsize=(10,6));ax.plot(Y,at1,color=colors[1],lw=2.5,label='AT1 stationary constrained solution');ax.plot(Y,eq,color=colors[0],lw=2.5,label='Eq. 3.25 relation');ax.axvline(1,color='#999999',ls=':',label='AT1 initiation threshold')
ax.set_xlabel('Crack driver / threshold (Y / Yc)');ax.set_ylabel('Bond phase parameter s');ax.set_title('Paper v1 formula consistency');ax.legend(fontsize=11,loc='lower right');ax.text(.02,.94,'PAPER v1 FORMULA CONSISTENCY CHECK\nNOT A MODEL SELECTION',transform=ax.transAxes,va='top',fontsize=12)
fig.text(.5,.005,'FORMULA COMPARISON ONLY / NOT VERIFIED PFPD EVOLUTION',ha='center',fontsize=11);finish(fig,'formula_consistency')
(O/'figures/STYLE.json').write_text(json.dumps(dict(background='black',font=family,font_weight='normal',arrows=False,formats=['png','pdf']),indent=2))
print('7 figures written as PNG + PDF')
