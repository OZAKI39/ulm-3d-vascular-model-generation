"""Academic English plots from server-computed audit evidence."""
from pathlib import Path
import json,textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data/remote_results';OUT=ROOT/'figures'
COLORS=['#3174a1','#c36d34','#4e9983','#8c6aaf']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.8,'savefig.dpi':220,'axes.labelsize':12,'legend.frameon':False,'figure.facecolor':'white','pdf.fonttype':42})
def save(fig,name,note=None):
    if note:fig.text(.02,.014,'\n'.join(textwrap.wrap(note,width=int(fig.get_figwidth()*14))),fontsize=9,color='#475569')
    fig.tight_layout(rect=(0,.065,1,.96));fig.savefig(OUT/(name+'.png'));fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
def ecdf(ax,values,label,color):
    x=np.sort(values[np.isfinite(values)]);n=len(x)
    if not n:return
    idx=np.unique(np.linspace(0,n-1,min(4000,n)).astype(int));positive=x[idx]>0
    ax.plot(x[idx][positive],(idx[positive]+1)/n,label=label,color=color,lw=1.8)
    ax.set_xscale('log');ax.set_ylim(0,1);ax.set_ylabel('Cumulative fraction of states')
    ax.grid(alpha=.18,which='major')
def scatter(ax,x,y,color='#3174a1',maxn=30000):
    valid=np.flatnonzero(np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0))
    if len(valid)>maxn:valid=valid[np.linspace(0,len(valid)-1,maxn).astype(int)]
    ax.scatter(x[valid],y[valid],s=3,alpha=.35,c=color,rasterized=True)
    ax.set_xscale('log');ax.set_yscale('log');ax.grid(alpha=.18)

def cards(name,title,items,footer):
    fig,ax=plt.subplots(figsize=(12,6.8));ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)
    ax.text(.02,.95,title,fontsize=21,weight='bold',color='#173c57')
    for i,(heading,body) in enumerate(items):
        y=.78-i*.205
        box=FancyBboxPatch((.02,y-.115),.96,.16,boxstyle='round,pad=0.012',fc='#eef3f7',ec='#c9d6e0')
        ax.add_patch(box);ax.text(.04,y+.007,heading,fontsize=14,weight='bold',color=COLORS[i%4]);ax.text(.04,y-.055,body,fontsize=12,va='top')
    fig.text(.04,.05,footer,fontsize=11,color='#963c27');fig.savefig(OUT/(name+'.png'));fig.savefig(OUT/(name+'.pdf'));plt.close(fig)

def main():
    OUT.mkdir(exist_ok=True)
    cards('00_stage_scope','Shear-lift magnitude audit | Scope',[
      ('Existing dynamics','Frozen FEM → Stokes self resistance → Normal lubrication → Nonpenetration → Velocity'),
      ('Read-only audit','Saved trajectories + local velocity gradients → slip, Reynolds numbers and force scales'),
      ('Candidate comparison','Rigid-sphere Saffman magnitude → drag / lubrication comparison → applicability checks'),
      ('Stage boundary','No added force · No trajectory recomputation · No automatic model development')],
      'Candidate Saffman magnitude only. Not included in trajectory dynamics.')
    cards('11_scientific_limitations','What this audit can — and cannot — establish',[
      ('Interface physics','Rigid no-slip sphere ≠ clean bubble ≠ validated lipid-shell SonoVue model'),
      ('Confinement','A small particle Reynolds number does not establish an unbounded-flow regime.'),
      ('Three-dimensional flow','A local strain-rate scalar and a candidate direction do not validate simple-shear theory.'),
      ('Evidence boundaries','Saved-state timing, incomplete paths, zero-slip limits and missing outlet routes remain explicit.')],
      'The new 2.0 mm/s field has no matched MB trajectories: slip-dependent lift is not evaluable.')
    if not (DATA/'statistics.json').exists():return
    with np.load(DATA/'all_scalar_states.npz') as f:d={k:f[k] for k in f.files}
    stats=json.loads((DATA/'statistics.json').read_text());groups=[('Far field','region_far_field'),('Near wall','region_near_wall'),('Branch ROI','region_branch')]
    fig,axs=plt.subplots(1,2,figsize=(11.8,5.2))
    for label,key in groups[:2]:
        i=0 if 'Far' in label else 1;m=d[key]
        ecdf(axs[0],d['endpoint_slip_speed_m_s'][m],label,COLORS[i]);ecdf(axs[1],d['slip_speed_m_s'][m],label,COLORS[i])
    for ax,title in zip(axs,['(a) Endpoint mismatch diagnostic','(b) Force-evaluation-aligned slip']):
        ax.set_title(title);ax.set_xlabel('Slip speed (m/s)');ax.legend(loc='lower right')
    save(fig,'02_slip_velocity_distribution','Same saved velocities; different sampling positions. Exact zeros are included in the CDF denominator.')
    fig,ax=plt.subplots(figsize=(9,5.5));scatter(ax,d['Re_G'],d['Re_p'])
    lo=max(1e-10,np.nanmin(d['Re_G'])*.5);hi=max(.1,np.nanmax(d['Re_G'])*2);x=np.geomspace(lo,hi,100)
    ax.plot(x,np.sqrt(x),'k--',lw=1,label=r'$Re_p=\sqrt{Re_G}$ (ordering boundary)')
    ax.plot(x,.1*np.sqrt(x),color='#a1a1aa',ls=':',label='0.1 × boundary (descriptive only)')
    ax.set_xlabel(r'Shear Reynolds number, $Re_G=a^2\dot\gamma_E/\nu$');ax.set_ylabel(r'Slip Reynolds number, $Re_p=a|u-v|/\nu$')
    ax.set_title('Radius convention | Asymptotic requirement: '+r'$Re_p\ll\sqrt{Re_G}\ll1$');ax.legend(loc='lower right',fontsize=10)
    save(fig,'03_dimensionless_regime_map','No verified validity region: all audited centers lie inside the shear disturbance length. Zero-slip points omitted from log axes.')
    fig,ax=plt.subplots(figsize=(9,5.2))
    for i,(label,key) in enumerate(groups):ecdf(ax,d['CANDIDATE_SAFFMAN_MAGNITUDE'][d[key]],label,COLORS[i])
    ax.set_xlabel('Candidate Saffman magnitude (N)');ax.set_title('Absolute magnitude | Local-shear scalar proxy');ax.legend(loc='lower right')
    save(fig,'04_candidate_lift_magnitude_distribution','Groups overlap; numerical roundoff values are retained. No claim of a validated physical bubble lift.')
    meaningful=d['lift_drag_status']==0
    fig,ax=plt.subplots(figsize=(9,5.5));scatter(ax,d['drag_N'][meaningful],d['CANDIDATE_SAFFMAN_MAGNITUDE'][meaningful])
    positive=d['drag_N'][meaningful];xs=np.geomspace(positive[positive>0].min()*.5,positive.max()*2,100)
    for ratio,style in [(1,'-'),(.1,'--'),(.01,':')]:ax.plot(xs,ratio*xs,style,color='#6b7280',lw=1,label=f'Lift / drag = {ratio:g}')
    ax.set_xlabel('Stokes drag magnitude (N)');ax.set_ylabel('Candidate Saffman magnitude (N)');ax.set_title('Aligned states with numerically resolved drag');ax.legend()
    save(fig,'05_lift_vs_drag',f'{np.sum(~meaningful):,} near-zero or unresolved denominators excluded; reference ratios are descriptive only.')
    reps=json.loads((DATA/'representatives.json').read_text());fig,axs=plt.subplots(2,2,figsize=(12,7.3));axs=axs.ravel()
    for ax,r in zip(axs,reps):
        with np.load(DATA/'representatives'/f"{r['dataset']}_{r['id']:06d}.npz") as f:
            x=f['time_s'];ld=f['lift_drag_ratio'];pref=6.46/(6*np.pi)*np.sqrt(f['Re_G'])
            ax.plot(x,pref,'--',color='#9ca3af',lw=1,label='Nonzero-slip limit (not 0/0)')
            ax.plot(x,ld,color='#3174a1',lw=1.5,label='Defined candidate / drag');ax.set_title(f"{r['dataset']} #{r['id']} → {r['outlet'].replace('_',' ')}",fontsize=11)
            ax.set_xlabel('Elapsed physical time (s)');ax.set_ylabel('Candidate / drag');ax.grid(alpha=.18);ax.legend(fontsize=8)
    for ax in axs[len(reps):]:ax.set_axis_off()
    save(fig,'06_lift_to_drag_along_representative_trajectories','Complete saved trajectories only. Gaps in the solid curve identify numerically zero drag; no Outlet 01 completion exists.')
    near=d['region_near_wall'];good=near&(d['lift_lubrication_status']==0)
    fig,axs=plt.subplots(1,2,figsize=(11.8,5.2));scatter(axs[0],d['lubrication_N'][good],d['CANDIDATE_SAFFMAN_MAGNITUDE'][good])
    xs=np.geomspace(d['lubrication_N'][good].min(),d['lubrication_N'][good].max(),100)
    for ratio,style in [(1,'-'),(.1,'--'),(.01,':')]:axs[0].plot(xs,xs*ratio,style,color='#6b7280',lw=1,label=f'{ratio:g}')
    axs[0].set_xlabel('Equivalent normal lubrication contribution (N)');axs[0].set_ylabel('Candidate Saffman magnitude (N)');axs[0].legend(title='Lift / lubrication')
    ecdf(axs[1],d['lift_lubrication_ratio'][good],'Near-wall states',COLORS[1]);axs[1].set_xlabel('Candidate / lubrication');axs[1].set_title('Defined ratios only')
    save(fig,'07_lift_vs_wall_lubrication',f'Resistance contribution, not independent force integration. {int(near.sum()-good.sum()):,} near-wall denominators unresolved.')
    fig,axs=plt.subplots(1,2,figsize=(11.8,5.2));scatter(axs[0],d['h_over_a'],d['CANDIDATE_SAFFMAN_MAGNITUDE'],color='#c36d34')
    axs[0].set_xlabel('Geometric gap / radius, h/a');axs[0].set_ylabel('Candidate Saffman magnitude (N)')
    scatter(axs[1],d['h_over_a'],d['wall_margin'],color='#c36d34');axs[1].axhline(1,color='k',ls='--',label='Wall at shear disturbance length')
    axs[1].set_xlabel('Geometric gap / radius, h/a');axs[1].set_ylabel(r'Center-wall distance / $\sqrt{\nu/\dot\gamma_E}$');axs[1].legend(fontsize=8)
    save(fig,'08_validity_by_wall_distance','All plotted states fail the unbounded ordering. h/a alone is not a Saffman validity test; no arbitrary gap-validity threshold.')
    table=np.genfromtxt(DATA/'HYPOTHETICAL_ONLY_sensitivity.csv',delimiter=',',names=True)
    fig,axs=plt.subplots(1,3,figsize=(13,5.2))
    for g,c in zip([10.,100.,1000.],COLORS):
        m=(table['assumed_slip_m_s']==1e-5)&(table['shear_s_inv']==g)&(table['h_over_a']==1)
        axs[0].plot(table['radius_m'][m]*1e6,table['CANDIDATE_SAFFMAN_MAGNITUDE'][m],'-o',color=c,label=f'{g:g} '+r's$^{-1}$')
    axs[0].set_yscale('log');axs[0].set_xlabel('Radius at D10–D90 (µm)');axs[0].set_ylabel('Candidate magnitude (N)');axs[0].legend(fontsize=9);axs[0].set_title('Radius squared')
    for g,c in zip([10.,100.,1000.],COLORS):
        m=(table['diameter_quantile']==50)&(table['shear_s_inv']==g)&(table['h_over_a']==1)&(table['assumed_slip_m_s']>0)
        axs[1].loglog(table['assumed_slip_m_s'][m],table['CANDIDATE_SAFFMAN_MAGNITUDE'][m],'-o',color=c)
    axs[1].set_xlabel('Assumed slip speed (m/s)');axs[1].set_ylabel('Candidate magnitude (N)');axs[1].set_title('Linear slip; square-root shear')
    for g,c in zip([10.,100.,1000.],COLORS):
        m=(table['diameter_quantile']==50)&(table['shear_s_inv']==g)&(table['assumed_slip_m_s']==1e-5)
        axs[2].loglog(table['h_over_a'][m],table['wall_margin'][m],'-o',color=c)
    axs[2].axhline(1,color='k',ls='--',lw=1);axs[2].set_xlabel('Assumed gap / radius');axs[2].set_ylabel('Wall-distance / shear-length');axs[2].set_title('Confinement, absent from formula')
    for ax in axs:ax.grid(alpha=.18)
    save(fig,'10_radius_sensitivity','HYPOTHETICAL ONLY. Radius panel assumes slip = 10 µm/s; other panels use D50. Radius changes do not create new trajectories.')

if __name__=='__main__':main()
