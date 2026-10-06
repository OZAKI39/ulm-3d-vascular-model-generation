"""Publication figures from saved physical arrays only; no additional mechanics."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
FIG=ROOT/'results/higher_order_baqp_energy_force_candidate_v1/figures'
BLUE,ORANGE,GREEN,PURPLE,GRAY='#176b9b','#c46522','#388368','#8e557f','#586779'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
    'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.16,'axes.axisbelow':True,
    'pdf.fonttype':42,'ps.fonttype':42,'legend.frameon':False})
QA={}
def read(n):return json.loads((OUT/n).read_text())
def rows(n):
    with (OUT/n).open() as f:return list(csv.DictReader(f))
def save(fig,name,caption):
    fig.text(.09,.025,caption,fontsize=8.5,color=GRAY,va='bottom',linespacing=1.5)
    fig.canvas.draw();renderer=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height();texts=list(fig.texts)
    for legend in fig.legends:texts+=legend.get_texts()
    for a in fig.axes:
        texts += [a.title,a.xaxis.label,a.yaxis.label]+list(a.texts)
        if a.get_legend():texts+=a.get_legend().get_texts()
        for axis,lim in [(a.xaxis,a.get_xlim()),(a.yaxis,a.get_ylim())]:
            texts.append(axis.get_offset_text())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if min(lim)<=tick.get_loc()<=max(lim):texts.extend([tick.label1,tick.label2])
    outside=[]
    for t in texts:
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        if b.x0 < -1 or b.y0 < -1 or b.x1>w+1 or b.y1>h+1:outside.append(t.get_text())
    QA[name]=dict(text_outside_canvas=outside,all_data_from_saved_physical_results=True,random_sampling=False)
    if outside:raise RuntimeError((name,outside))
    fig.savefig(FIG/(name+'.png'),dpi=600,facecolor='white');fig.savefig(FIG/(name+'.pdf'),facecolor='white');plt.close(fig)
    print(name,flush=True)

def c1():
    summary=read('C1_FORCE_REGRESSION_SUMMARY.json')
    fig,ax=plt.subplots(1,3,figsize=(12.5,5.4));fig.subplots_adjust(left=.085,right=.98,bottom=.24,top=.80,wspace=.33)
    with np.load(OUT/'C1_FORCE_REGRESSION.npz') as data:
        for k,(name,color) in enumerate(zip(['affine','quadratic_shear','seeded_smooth'],[BLUE,ORANGE,GREEN])):
            old=data[name+'_historical_force'].ravel()*1e9;new=data[name+'_variational_force'].ravel()*1e9
            ax[k].scatter(old,new,s=14,facecolors='none',edgecolors=color,linewidths=.65)
            lo,hi=min(old.min(),new.min()),max(old.max(),new.max());pad=.1*(hi-lo)
            ax[k].plot([lo-pad,hi+pad],[lo-pad,hi+pad],'--',c=GRAY,lw=1)
            ax[k].set(title=['(a) Finite affine','(b) Quadratic shear','(c) Seeded smooth'][k],
                xlabel='Existing C1 force component [nN]',ylabel='Variational C1 component [nN]')
            r=summary[name]
            ax[k].text(.04,.96,f"Max / Fref: {r['max_normalized']:.2e}\nRMS / Fref: {r['RMS_normalized']:.2e}",
                       transform=ax[k].transAxes,va='top',fontsize=9)
    fig.suptitle('The same variational backend reproduces the existing C1 force',y=.95,fontsize=14)
    save(fig,'c1_variational_force_regression',
        'All 343 nodes and all three force components in each small-cube state; no subsampling. Dashed line: y = x.\n'
        'Physical-volume weighted RMS; unchanged intact SVK energy and unordered-bond prefactor 2 wbar Vi Vj.')

def fd():
    data=rows('ENERGY_FORCE_DIRECTIONAL_DERIVATIVES.csv');fig=plt.figure(figsize=(12.5,8.7))
    gs=fig.add_gridspec(2,3,left=.085,right=.98,bottom=.24,top=.85,hspace=.53,wspace=.34,height_ratios=[1,1])
    ax=[fig.add_subplot(gs[0,k]) for k in range(3)];error=fig.add_subplot(gs[1,:])
    fields=['finite_affine_shear','shear','cross'];titles=['Finite affine shear','Quadratic shear','XZ cross term']
    for k,(field,color) in enumerate(zip(fields,[BLUE,ORANGE,GREEN])):
        rr=[r for r in data if r['state']==field]
        xs=np.array([float(r['D_force_N'])*1e9 for r in rr]);ys=np.array([float(r['D_fd_N'])*1e9 for r in rr])
        lo,hi=min(xs.min(),ys.min()),max(xs.max(),ys.max());pad=max((hi-lo)*.12,1e-5)
        ax[k].plot([lo-pad,hi+pad],[lo-pad,hi+pad],'--',c=GRAY,lw=1)
        for e,marker in [(1e-5,'o'),(3e-6,'s'),(1e-6,'+')]:
            r=[v for v in rr if float(v['eps_rel'])==e]
            ax[k].plot([float(v['D_force_N'])*1e9 for v in r],[float(v['D_fd_N'])*1e9 for v in r],
                linestyle='none',marker=marker,ms=8 if marker=='o' else 5,mfc='none',c=color,label=f'{e:g}')
        ax[k].set(title=f'({"abc"[k]}) '+titles[k],xlabel=r'$-f\cdot v$ [nN]',ylabel='Central FD derivative [nN]')
        for direction,style in [(0,'-'),(1,'--')]:
            r=[v for v in rr if int(v['direction'])==direction]
            error.loglog([float(v['eps_rel']) for v in r],[float(v['normalized_error']) for v in r],
                style,marker='o' if direction==0 else 's',ms=4,c=color,label=titles[k]+f' / direction {direction+1}')
    error.axhline(5e-6,c=PURPLE,ls=':',lw=1);error.invert_xaxis()
    error.set(xlabel=r'Fixed perturbation $\epsilon/h$',ylabel='Normalized derivative error',title='(d) Every value from all three predeclared steps')
    error.text(.02,.91,'Middle-step gate: 5e-6',transform=error.transAxes,color=PURPLE,fontsize=9)
    fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.54,.945),ncol=3,
        title='Fixed steps (all coincide at the scatter scale)',fontsize=9,title_fontsize=9)
    fig.legend(*error.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.54,.104),ncol=3,fontsize=8.5)
    fig.suptitle('PD2 force is the first variation of the same discrete energy',y=.995,fontsize=14)
    save(fig,'pd2_energy_force_directional_derivative',
        'Exactly two deterministic directions per state; all 18 finite-difference evaluations shown. No step tuning.\n'
        'Stable normalization uses max(|Dfd|, |Dforce|, Eref ||v|| / h); unfloored errors are also saved in the CSV.')

def objectivity():
    obj=read('OBJECTIVITY_RESULTS.json');data=rows('GLOBAL_BALANCE_RESULTS.csv');cfg=read('CANDIDATE_CONFIGURATION.json')
    fig,ax=plt.subplots(2,2,figsize=(12,10.5),gridspec_kw={'height_ratios':[.8,1.25]})
    fig.subplots_adjust(left=.23,right=.98,bottom=.14,top=.88,hspace=.41,wspace=.33)
    energy=[obj[n]['energy_over_characteristic'] for n in ['translation','rotation']]
    force=[obj[n]['force_over_characteristic'] for n in ['translation','rotation']]
    ax[0,0].bar([0,1],energy,color=[BLUE,ORANGE],width=.55);ax[0,0].set_yscale('log')
    ax[0,0].set_xticks([0,1],['Translation','37 deg rotation'])
    ax[0,0].set(title='(a) Rigid-state energy',ylabel=r'$U/(\mu\,0.004^2 V)$')
    ax[0,0].axhline(cfg['rigid_energy_absolute_gate_J']/(1000*.004**2*cfg['primary_L_m']**3),ls=':',c=PURPLE)
    ax[0,1].bar([0,1],force,color=[BLUE,ORANGE],width=.55);ax[0,1].set_yscale('log')
    ax[0,1].set_xticks([0,1],['Translation','37 deg rotation'])
    ax[0,1].set(title='(b) Rigid-state maximum force',ylabel=r'$\max_i |f_i|/(4\,Pa\,h^2)$')
    ax[0,1].axhline(obj['translation']['force_gate_N']/(4*cfg['h_m']**2),ls=':',c=PURPLE)
    labels=[]
    for r in data:
        name={'quadratic_shear':'Quadratic shear','seeded_smooth':'Smooth','affine':'Affine',
            'finite_affine_shear':'Finite shear','axial':'Quadratic axial','shear':'Quadratic shear',
            'cross':'XZ','affine_axial':'Affine axial','affine_shear':'Affine shear',
            'superposed_finite_shear':'Rotated finite shear','cross_selected_DOF':'Selected-DOF XZ'}[r['state']]
        labels.append(f"{r['model']} N{r['N']} / {name}")
    colors=[GRAY if r['model']=='C1' else GREEN for r in data]
    for a,key,title,label in [(ax[1,0],'force_residual','(c) Global internal force',r'$|\sum_i f_i|/F_{ref}$'),
                              (ax[1,1],'moment_residual','(d) Current-position moment',r'$|\sum_i (x_i-x_c)\times f_i|/(F_{ref}L)$')]:
        a.scatter([float(r[key]) for r in data],range(len(data)),c=colors,s=25)
        a.set_xscale('log');a.set_yticks(range(len(data)),labels if a is ax[1,0] else [])
        a.invert_yaxis();a.axvline(1e-10,c=PURPLE,ls=':');a.set(title=title,xlabel=label)
        a.tick_params(axis='y',labelsize=8.5)
    fig.suptitle('Rigid motions and global balances of the intact variational candidate',y=.96,fontsize=13)
    save(fig,'pd2_objectivity_and_global_balance',
        'All validation states are included. Moments use current positions. Dotted lines: predeclared gates.\n'
        'The rigid-force gate retains the historical 1e-15 N noise floor; measured rigid forces are much smaller.\n'
        'No equilibrium solve, boundary-load conversion, phase degradation or fracture is used.')

def energy():
    data=rows('QUADRATIC_ENERGY_DECOMPOSITION.csv');fig,ax=plt.subplots(1,2,figsize=(12,6.2))
    fig.subplots_adjust(left=.085,right=.98,bottom=.28,top=.79,wspace=.32)
    x=np.arange(3);width=.18
    series=[('U_C1_J','C1',BLUE),('U_PD2_J','PD2',ORANGE),('U_exactkin_currentweights_J','Exact kinematics / current weights',GREEN),('U_continuum_J','Continuum Gauss',PURPLE)]
    for k,(key,name,color) in enumerate(series):
        ax[0].bar(x+(k-1.5)*width,[float(r[key])*1e12 for r in data],width,label=name,color=color)
    ax[0].set(title='(a) Same physical energy, different reconstruction',ylabel='Intact energy [pJ]');ax[0].set_xticks(x,['Axial','Shear','XZ'])
    for k,(key,label,color) in enumerate([('kinematic_correction_J','PD2 - C1',ORANGE),
                                          ('remaining_weight_quadrature_difference_J','PD2 - continuum',PURPLE)]):
        ax[1].bar(x+(k-.5)*.3,[float(r[key])*1e15 for r in data],.3,label=label,color=color)
    ax[1].axhline(0,c=GRAY,lw=.8);ax[1].set(title='(b) Kinematic correction and remaining difference',ylabel='Energy difference [fJ]')
    ax[1].set_xticks(x,['Axial','Shear','XZ']);ax[1].legend(loc='upper left',fontsize=9)
    fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.53,.905),ncol=2,fontsize=9)
    fig.suptitle('PD2 matches exact kinematics with the unchanged bond weights',y=.98,fontsize=14)
    save(fig,'c1_pd2_energy_decomposition',
        'All 690,020 physical bonds. PD2-to-exact-current-weight normalized error <= 4.63e-14.\n'
        'A remaining weight / quadrature difference from the continuum reference is retained, not fitted away.\n'
        'Intact SVK only; this comparison does not validate physical boundary conditions.')

if __name__=='__main__':
    for f in [c1,fd,objectivity,energy]:f()
    (OUT/'FIGURE_LAYOUT_QA.json').write_text(json.dumps(QA,indent=2)+'\n')
