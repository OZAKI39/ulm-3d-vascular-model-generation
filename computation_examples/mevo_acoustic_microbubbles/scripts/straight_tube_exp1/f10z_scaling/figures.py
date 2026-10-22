"""Academic plots from accepted particle positions, with physical-radius XY projections."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import numpy as np
from .common import BASE, R

def draw(groups,runs,metrics,arrays):
    out=BASE/'figures'; out.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.labelcolor':'#242424',
        'xtick.color':'#242424','ytick.color':'#242424','axes.spines.top':False,'axes.spines.right':False,
        'savefig.dpi':240,'pdf.fonttype':42})
    colors={32:'#6355a4',64:'#2185b0',128:'#b5653c'}
    selected={32:'historical32'}
    for n in [64,128]:
        if 'result' in groups[str(n)]: selected[n]=groups[str(n)]['result']['name']
    fig,axs=plt.subplots(2,3,figsize=(12,8),layout='constrained',sharex=True,sharey=True)
    for col,n in enumerate([32,64,128]):
        if n not in selected:
            for ax in axs[:,col]: ax.text(.5,.5,'NOT_RUN',transform=ax.transAxes,ha='center')
            continue
        t,x=arrays[selected[n]]; m=metrics[selected[n]]
        for row,i in enumerate([0,-1]):
            ax=axs[row,col]; q=x[i]*1e6
            order=np.argsort(q[:,2])
            circles=PatchCollection([Circle(q[k,:2],R*1e6) for k in order],facecolor='#3886bd',edgecolor='#214c70',linewidth=.35,alpha=.65)
            ax.add_collection(circles);ax.set_xlim(-34,34);ax.set_ylim(-34,34);ax.set_aspect('equal')
            ax.grid(alpha=.15);ax.set_axisbelow(True)
            ax.set_title(f'N = {n} | t = {t[i]*1e3:g} ms\n{m["components"][i]} components; largest = {m["largest_component"][i]}')
            if row==1:ax.set_xlabel('X (µm)')
            if col==0:ax.set_ylabel(('Initial\n' if row==0 else 'Accepted endpoint\n')+'Y (µm)')
    fig.suptitle('Dispersed bubbles under the same X-pair control\nXY projection; physical sphere radii; equal nominal initial density',fontsize=14)
    for ext in ['png','pdf']:fig.savefig(out/('Fig_particle_count_morphology.'+ext))
    plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for n,key in selected.items():
        m=metrics[key];t=m['time_s']*1e3;c=colors[n]
        label=f'{n} bubbles'+(' (prior 2 ms baseline)' if n==32 else '')
        axs[0,0].plot(t,100*(m['rg_m']/m['rg_m'][0]-1),color=c,label=label)
        axs[0,1].plot(t,np.linalg.norm(m['centroid_m']-m['centroid_m'][0],axis=1)*1e6,color=c)
        axs[1,0].plot(t,m['components']/n,color=c)
    axs[0,0].axhline(0,color='#555555',linestyle='--',linewidth=1,label='OFF: no motion')
    axs[0,0].set_ylabel('Change in Rg (%)');axs[0,0].legend(fontsize=9)
    axs[0,1].set_ylabel('Centroid displacement (µm)')
    axs[1,0].set_ylabel('Connected components / N')
    for ax in [axs[0,0],axs[0,1],axs[1,0]]:ax.set_xlabel('Physical time (ms)');ax.grid(alpha=.18)
    nvals=[32,64,128];costs=[]
    for n in nvals:
        r=next(v for v in runs.values() if v['N']==n and v['action']=='PAIR' and abs(v['requested_s']-1e-4)<1e-12)
        costs.append(r['mean_step_s']*1e3)
    ax=axs[1,1];bars=ax.bar([str(n) for n in nvals],costs,color=[colors[n] for n in nvals])
    ax.bar_label(bars,labels=[f'{x:.2f}' for x in costs],padding=3);ax.set_ylim(0,max(costs)*1.18)
    ax.set_xlabel('Bubble count N');ax.set_ylabel('Measured time / step (ms)')
    ax.set_title('Fresh 0.1 ms probes; same solver and dt = 1 µs',fontsize=10)
    fig.suptitle('Particle-count feasibility and aggregation response\nEXPLORATORY · Homogeneous-water proxy · Physics uncalibrated',fontsize=13)
    for ext in ['png','pdf']:fig.savefig(out/('Fig_particle_count_response.'+ext))
    plt.close(fig)
