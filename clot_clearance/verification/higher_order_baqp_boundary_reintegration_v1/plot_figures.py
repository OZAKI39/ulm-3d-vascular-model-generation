"""Publication figures from saved numerical records; no simulations."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch

V=Path(__file__).resolve().parent;R=V.parents[1]/'results/higher_order_baqp_boundary_reintegration_v1';F=R/'figures'
F.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':10,'legend.fontsize':8.5,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
BLUE='#176A9A';ORANGE='#C96B23';GREEN='#18866E';GRAY='#656B73'
KINDS=['normal_affine','shear_affine','nonaffine_axial','nonaffine_shear']
LABELS=['Normal affine','Shear affine','Quadratic axial','Quadratic shear']
def js(n):return json.loads((V/n).read_text())
def rows(n):
    with (V/n).open() as f:return list(csv.DictReader(f))
def vals(rs,key):return np.array([float(r[key]) for r in rs])
def positive(x):return np.maximum(np.asarray(x),1e-18)
def save(fig,name):
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(F/(name+'.png'),dpi=600,bbox_inches='tight')
    plt.close(fig);print('FIGURE',name,flush=True)
def clean(ax):ax.grid(axis='y',alpha=.18);ax.set_axisbelow(True)

def overview():
    fig,ax=plt.subplots(figsize=(11.5,6.5));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
    ax.text(.02,.97,'PD2 boundary reintegration',fontsize=17,weight='bold',va='top')
    ax.text(.02,.895,'Same surrogate field, separate finite-domain and fully closed mechanics',color=GRAY,fontsize=11)
    def box(x,y,w,h,title,body,color):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012',fc=color,ec='none'))
        ax.text(x+.018,y+h-.035,title,weight='bold',va='top',fontsize=10.5)
        ax.text(x+.018,y+h-.09,body,va='top',fontsize=9.5,linespacing=1.4)
    def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=14,lw=1.3,color=GRAY,connectionstyle='arc3'))
    box(.025,.50,.225,.27,'Classical surrogate',r'$x_L=X+u_L$'+'\n'+r'$b_{local}=-\mathrm{Div}\,P$'+'\n'+r'$t=PN$ on all six faces','#F0F2F5')
    box(.31,.61,.29,.20,'Finite-domain PD2',r'Internal force $f_{FD}(x_L)$'+'\nN19; physical nodal volumes','#E3EFF6')
    box(.31,.34,.29,.20,'Closed PD2 bulk oracle',r'$q_{bulk}(X)=f_{closed}/h^3$'+'\nN27; full dependency graph','#E3EFF6')
    box(.665,.45,.31,.28,'External equivalent load',r'$f_{BC}=q_{bulk}V-f_{FD}$'+'\nFrozen during equilibrium\nConversion support only','#FBEDDF')
    box(.33,.025,.61,.235,'Existing equilibrium solver + PD2',r'$f_{int}(x)+f_{BC}+bV=0$ on free nodes'+'\nNL: '+r'$b=-q_{bulk}$'+'     Classical shear: '+r'$b=b_{local}$'+'\nSame prescribed z faces; deterministic nonzero perturbation','#E5F3EC')
    arrow((.25,.66),(.31,.71));arrow((.25,.56),(.31,.46));arrow((.60,.71),(.665,.65));arrow((.60,.46),(.665,.52));arrow((.82,.45),(.82,.27))
    ax.text(.025,.20,'Internal mechanics',color=BLUE,weight='bold');ax.text(.025,.14,'External load',color=ORANGE,weight='bold')
    save(fig,'pd2_boundary_reintegration_overview')

def affine():
    fig,axs=plt.subplots(1,2,figsize=(11,4.4),layout='constrained');t=js('PD2_TARGET_SUMMARIES.json')
    labels=['Resultant','Moment','Virtual\nwork',r'$\bar F$',r'$\widetilde F$',r'$\Psi^+$'];gate=np.array([1e-6]*3+[1e-5]*2+[1e-4])
    for ax,k,label in zip(axs,KINDS[:2],LABELS[:2]):
        s=js('PD2_'+k.upper()+'.json');c=t[k]['load_checks']
        y=[c['resultant_error'],c['moment_error'],c['virtual_work_error'],s['fields']['Fbar']['F_error_max_scaled'],s['fields']['Ftilde']['F_error_max_scaled'],max(r['Psi_plus_error_max_scaled'] for r in s['fields'].values())]
        ax.bar(np.arange(6),positive(y),color=BLUE,width=.58,label='Measured error');ax.scatter(np.arange(6),gate,marker='_',s=400,lw=2,color=ORANGE,label='Predeclared gate',zorder=4)
        ax.set_yscale('log');ax.set_ylim(1e-18,1e-3);ax.set_xticks(range(6),labels);ax.set_ylabel('Normalized error');ax.set_title(label);clean(ax)
    axs[0].legend(loc='upper left');save(fig,'pd2_affine_boundary_consistency')

def localization():
    pro=rows('PD2_BOUNDARY_SUPPORT.csv');points=js('PD2_CLOSED_CENTER_LOCALIZATION.json');sc=js('PD2_TARGET_SUMMARIES.json')
    fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    for ax,k,label in zip(axs,KINDS[2:],LABELS[2:]):
        r=[x for x in pro if x['case']==k and x['backend']=='PD2'];scale=sc[k]['scales']['F_ref_N']
        ax.semilogy(vals(r,'depth_over_delta'),positive(vals(r,'load_max_N')/scale),'o-',color=BLUE,label='N19 layer maximum')
        ax.semilogy(vals(r,'depth_over_delta'),positive(vals(r,'load_mean_N')/scale),'--',color=BLUE,alpha=.65,label='N19 layer mean')
        p=[x for x in points if x['case']==k];x=np.array([a['depth_h']/4.02 for a in p]);y=np.array([a['boundary_force_scaled'] for a in p])
        ax.semilogy(x,y,'s:',color=ORANGE,label='N27 lateral-depth samples');ax.scatter(x[-1],y[-1],s=80,facecolors='none',edgecolors=GREEN,lw=2,zorder=5,label='Fully closed center')
        ax.axhline(1e-10,color=GRAY,ls='--',lw=1,label='Closed-center gate')
        ax.set_ylim(1e-17,3);ax.set_xlabel(r'Nearest boundary depth $d_{min}/\delta$');ax.set_ylabel(r'$|f_{BC}|/F_{ref}$');ax.set_title(label);clean(ax)
    axs[0].legend(loc='lower left',fontsize=8);save(fig,'pd2_nonaffine_boundary_localization')

def support():
    pro=rows('PD2_BOUNDARY_SUPPORT.csv');quant=rows('C1_PD2_BOUNDARY_COMPARISON.csv')
    fig,axs=plt.subplots(2,4,figsize=(13,6.2),layout='constrained')
    fig.suptitle('CONVERSION SUPPORT / NOT PHYSICAL BC THICKNESS',fontsize=14,weight='bold')
    for j,(k,label) in enumerate(zip(KINDS,LABELS)):
        for backend,color,shift in [('C1',ORANGE,-.16),('PD2',BLUE,.16)]:
            r=[a for a in pro if a['case']==k and a['backend']==backend];q=next(a for a in quant if a['case']==k and a['backend']==backend)
            axs[0,j].plot(vals(r,'depth_over_delta'),vals(r,'cumulative_L1_fraction'),'.-',color=color,label=backend)
            axs[1,j].bar(np.arange(4)+shift,[float(q[f'd{n}_delta']) for n in [50,90,95,99]],width=.30,color=color,label=backend)
            if backend=='PD2' and float(q['d50_delta'])==0:
                axs[1,j].text(shift,.025,'0',ha='center',va='bottom',color=BLUE,fontsize=9)
        axs[0,j].set_title(label);axs[0,j].set_xlabel(r'$d_{min}/\delta$');axs[0,j].set_ylim(0,1.05)
        axs[1,j].set_xticks(range(4),['d50','d90','d95','d99']);axs[1,j].set_ylim(0,1.1)
        clean(axs[0,j]);clean(axs[1,j])
    axs[0,0].set_ylabel('Cumulative load magnitude');axs[1,0].set_ylabel(r'Support quantile / $\delta$');axs[0,0].legend(loc='lower right')
    save(fig,'c1_vs_pd2_boundary_load_support')

def fieldquality():
    data=rows('C1_PD2_FIELD_QUALITY.csv');fig,axs=plt.subplots(2,3,figsize=(12,6.6),layout='constrained')
    for i,(k,label) in enumerate(zip(KINDS[2:],LABELS[2:])):
        for j,sample in enumerate(['Fbar','Ftilde']):
            ax=axs[i,j]
            for backend,color in [('C1_TARGET',ORANGE),('PD2_SOLVED',BLUE)]:
                r=[a for a in data if a['case']==k and a['backend']==backend and a['sample']==sample]
                ax.semilogy(vals(r,'depth_over_delta'),positive(vals(r,'F_error_max_scaled')),linestyle='-' if sample=='Fbar' else '--',color=color,lw=1.7,label='C1 saved target' if backend=='C1_TARGET' else 'PD2 solved')
            ax.set_title(label+(' — nodal gradient' if j==0 else ' — bond gradient'));ax.set_ylabel(r'Max $\|F-F_{analytic}\|/0.004$');ax.set_ylim(1e-14,1e-1);clean(ax)
        ax=axs[i,2]
        for backend,color in [('C1_TARGET',ORANGE),('PD2_SOLVED',BLUE)]:
            for sample,style in [('Fbar','-'),('Ftilde','--')]:
                r=[a for a in data if a['case']==k and a['backend']==backend and a['sample']==sample]
                ax.semilogy(vals(r,'depth_over_delta'),positive(vals(r,'Psi_plus_error_max_scaled')),style,color=color,lw=1.5,label=('C1 target' if backend=='C1_TARGET' else 'PD2 solved')+(' node' if sample=='Fbar' else ' bond'))
        ax.set_title(label+' — tensile-energy diagnostic');ax.set_ylabel(r'Max $|\Delta\Psi^+|/\Psi_{ref}$');ax.set_ylim(1e-14,1);clean(ax)
        for ax in axs[i]:ax.set_xlabel(r'Nearest sample depth $d_{min}/\delta$')
    handles,labels=axs[0,2].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=4,frameon=False)
    save(fig,'c1_vs_pd2_boundary_field_quality')

def recovery():
    data=rows('C1_PD2_FIELD_QUALITY.csv');fig,axs=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for ax,k,label,letter in zip(axs.flat,KINDS,LABELS,'ABCD'):
        for backend,color in [('PD2_TARGET',GRAY),('PD2_SOLVED',BLUE)]:
            for sample,style in [('Fbar','-'),('Ftilde','--')]:
                r=[a for a in data if a['case']==k and a['backend']==backend and a['sample']==sample]
                ax.semilogy(vals(r,'depth_over_delta'),positive(vals(r,'F_error_max_scaled')),style,color=color,label=('Target' if backend=='PD2_TARGET' else 'Solved')+(' nodal' if sample=='Fbar' else ' bond'))
        r=js('PD2_'+k.upper()+'.json');ax.axhline(1e-5,color=ORANGE,ls=':',label='Recovery gate')
        ax.set_title(f'({letter}) {label}');ax.text(.03,.965,f"Mechanical residual: {r['force_residual']:.2e}",transform=ax.transAxes,va='top',fontsize=9)
        ax.set_xlabel(r'Nearest sample depth $d_{min}/\delta$');ax.set_ylabel('Max normalized gradient error');ax.set_ylim(1e-14,1e-4);clean(ax)
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=5,frameon=False,fontsize=8.5)
    save(fig,'pd2_manufactured_recovery')

def classical():
    r=rows('PD2_CLASSICAL_SHEAR_SHIFT_PROFILE.csv');s=js('PD2_CLASSICAL_VS_NONLOCAL_SHEAR.json');fig,axs=plt.subplots(1,3,figsize=(12,4.3),layout='constrained')
    for ax,metric,label in zip(axs,['displacement_error_max_m','Fbar_error_max','Psi_error_max_Pa'],['Displacement deviation [m]',r'Nodal $\|F-F_{local}\|$',r'Nodal $|\Psi^+-\Psi^+_{local}|$ [Pa]']):
        for prefix,color,name in [('NL',BLUE,'Nonlocal manufactured'),('classical',ORANGE,'Classical loaded')]:
            ax.plot(vals(r,'Z_um'),vals(r,prefix+'_'+metric),color=color,label=name)
        ax.axhline(0,color=GRAY,ls=':',label='Analytic local target');ax.set_xlabel(r'Reference $Z$ [$\mu$m]');ax.set_ylabel(label);ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));clean(ax)
    handles,labels=axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,-.02),ncol=3,frameon=False)
    fig.suptitle(r'Classical shear retains $e_{model}=q_{bulk}+b_{local}$; no compensating load',fontsize=12)
    save(fig,'pd2_classical_vs_nonlocal_shear')

if __name__=='__main__':
    for fn in [overview,affine,localization,support,fieldquality,recovery,classical]:fn()
