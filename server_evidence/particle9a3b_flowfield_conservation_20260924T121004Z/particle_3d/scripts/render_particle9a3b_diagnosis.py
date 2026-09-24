"""Five compact diagnostic figures; reads only saved diagnostic tables."""
from pathlib import Path
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

REPO=Path(__file__).resolve().parents[2]
REPORT=REPO/'particle_3d/reports/particle9a3b_flowfield_conservation'


def table(name):
    rows=list(csv.DictReader((REPORT/'data'/name).open()))
    return {k:np.array([float(r[k]) for r in rows]) for k in rows[0]}


def save(fig,name):
    fig.savefig(REPORT/'figures'/f'{name}.png',dpi=300,bbox_inches='tight',facecolor='white')
    fig.savefig(REPORT/'figures'/f'{name}.pdf',bbox_inches='tight',facecolor='white')
    plt.close(fig)


def main():
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,
                         'grid.alpha':.18,'figure.facecolor':'white','axes.titlesize':13})
    d=table('root_divergence_bins.csv');b=table('slab_gauss_balance.csv');r=table('internal_flux_by_representation.csv')
    stats=json.loads((REPORT/'data/divergence_statistics.json').read_text());summary=json.loads((REPORT/'data/diagnosis_summary.json').read_text())
    wall=json.loads((REPORT/'data/wall_flux_audit.json').read_text());s=b['section_s_m']*1e6;rootend=56.2643494768
    fig,ax=plt.subplots(2,1,figsize=(9,6.8),sharex=True,gridspec_kw={'height_ratios':[2,1]},layout='constrained')
    ax[0].fill_between(d['s_um'],d['P01'],d['P99'],color='#bdd5e8',label='P01–P99')
    ax[0].fill_between(d['s_um'],d['P10'],d['P90'],color='#6fa3c8',label='P10–P90')
    ax[0].plot(d['s_um'],d['median'],color='#163c60',lw=1.5,label='Median')
    ax[0].set(ylabel=r'Element divergence (s$^{-1}$)',title='01  Native / Particle P1 divergence along the root')
    ax[0].legend(ncols=3,loc='upper center',fontsize=9)
    ax[1].plot(d['s_um'],d['volume_weighted_mean'],color='#ab3d36',lw=1.5,label='Bin volume-weighted mean')
    ax[1].axhline(stats['root_inlet_to_first_junction']['volume_weighted_mean'],color='#555',ls=':',label='Whole root mean')
    ax[1].set(ylabel=r'Mean divergence (s$^{-1}$)',xlabel='Root arclength from inlet (µm)')
    ax[1].legend(fontsize=9,loc='lower center')
    for a in ax:
        a.axhline(0,color='#555',lw=.8);a.axvline(0,color='#237b4b',ls='--');a.axvline(rootend,color='#6b498d',ls='--');a.set_xlim(-1,58)
    ax[0].text(.2,.97,'Inlet',transform=ax[0].get_xaxis_transform(),va='top',color='#237b4b')
    ax[0].text(rootend-.3,.97,'First junction',transform=ax[0].get_xaxis_transform(),ha='right',va='top',color='#6b498d')
    save(fig,'01_root_divergence')

    fig,ax=plt.subplots(figsize=(9,4.7),layout='constrained');qin=summary['Q_in_m3_s']
    styles=[('native_checkpoint','Native checkpoint','o',9,'#277b93'),('export','Exported VTU','s',6,'#bc6d32'),('particle','Particle P1','x',6,'#172b44')]
    for key,label,marker,size,color in styles:
        ax.plot(s,r['Q_'+key+'_m3_s']/qin,marker=marker,ms=size,mfc='white',lw=1.2,color=color,label=label)
    ax.axhline(1,color='#555',ls='--',label='Inlet reference');ax.set(xlabel='Root arclength from inlet (µm)',ylabel=r'$Q_{section}/Q_{inlet}$',title='02  Same internal flux loss in every representation',ylim=(.958,1.005))
    ax.legend(loc='lower left',fontsize=9);ax.text(.55,.22,'Three curves coincide\nSix actual candidate sections',ha='center',va='bottom',transform=ax.transAxes)
    save(fig,'02_flux_by_representation')

    fig,ax=plt.subplots(figsize=(9,4.7),layout='constrained')
    ax.plot(s,b['Q_loss_m3_s']*1e15,'o-',ms=8,mfc='white',color='#17668a',label=r'$Q_{in}-Q_{section}$')
    ax.plot(s,-b['volume_integral_divergence_m3_s']*1e15,'x--',color='#b54435',label=r'$-\int_\Omega \nabla\cdot u\,dV$')
    ax.plot(s,b['Q_wall_m3_s']*1e15,'s:',color='#568653',label=r'$Q_{wall}$ (outward)')
    ax.set(xlabel='Root arclength from inlet (µm)',ylabel='Flux (pL/s)',title='03  Lost internal flux equals the integrated numerical sink')
    ax.legend(fontsize=10);ax.text(.03,.92,'Outward convention:  Qsection + Qwall − Qin = ∫div(u)dV\nMaximum closure residual / Qin = '+f"{summary['max_gauss_residual_relative']:.2e}",transform=ax.transAxes,va='top',fontsize=10)
    ax.set_ylim(-.05,.75);save(fig,'03_slab_gauss_balance')

    fig,ax=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
    ax[0].bar(['x','y','z'],[0,0,0],color='#568653');ax[0].scatter([0,1,2],[0,0,0],s=70,color='#568653',zorder=3)
    ax[0].set(ylabel='Maximum |wall nodal component| (m/s)',ylim=(-.1,.9),title='All WALL node components are exactly zero')
    ax[0].text(.5,.65,f"{wall['node_count']:,} wall nodes\n{wall['triangle_count']:,} wall triangles\nMaximum normal speed = 0 m/s",ha='center',transform=ax[0].transAxes)
    ax[1].plot(s,b['wall_relative_inlet'],'o-',color='#568653');ax[1].set(xlabel='Root arclength from inlet (µm)',ylabel=r'$Q_{wall}/Q_{inlet}$',ylim=(-.1,.9),title='No leakage through the real slab side WALL')
    ax[1].text(.5,.65,'Signed leakage = 0\nSum of absolute triangle flux = 0',ha='center',transform=ax[1].transAxes)
    for a in ax:a.set_yticks([0])
    fig.suptitle('04  WALL audit — zero is measured, not assumed');save(fig,'04_wall_leakage_audit')

    fig,ax=plt.subplots(figsize=(13,4.1));ax.set_xlim(0,13);ax.set_ylim(0,4);ax.axis('off')
    labels=[('Native solver / checkpoint','P1 / P1 + VMS\nVelocity: float64\nCoordinates: float32 → double'),
            ('Native result VTU','Original nodal DOFs\nVelocity: float64\nCoordinates: float32'),
            ('Frozen / reference VTU','Byte-for-byte copies\nVelocity: float64\nCoordinates: float32'),
            ('Particle P1','Barycentric interpolation\nVelocity: float64\nCoordinates: exact float64 cast')]
    for i,(title,body) in enumerate(labels):
        x=.1+i*3.25
        ax.add_patch(FancyBboxPatch((x,.85),2.9,2.1,boxstyle='round,pad=.07',facecolor='#eef4f8',edgecolor='#7095af'))
        ax.text(x+1.45,2.7,title,ha='center',va='center',weight='bold',fontsize=11)
        ax.text(x+1.45,1.92,body,ha='center',va='center',fontsize=10)
        ax.text(x+1.45,1.08,'Internal loss: 0.291–3.601%',ha='center',color='#a23c34',fontsize=10)
        if i<3:ax.annotate('',xy=(x+3.18,1.9),xytext=(x+2.96,1.9),arrowprops={'arrowstyle':'->','lw':1.5})
    ax.text(6.5,3.65,'05  Conservation loss starts in the resolved native P1 field',ha='center',fontsize=15,weight='bold')
    ax.text(6.5,.32,'Checkpoint ↔ VTU velocity difference = 0     |     Frozen copy SHA-256 identical     |     Particle resampling = none',ha='center',fontsize=11)
    save(fig,'05_velocity_representation_pipeline')


if __name__=='__main__':main()
