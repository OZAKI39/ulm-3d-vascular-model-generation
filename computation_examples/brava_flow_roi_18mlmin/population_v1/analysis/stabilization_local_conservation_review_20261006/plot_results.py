from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':12,'axes.titlesize':13,'axes.labelsize':12,'legend.fontsize':10,'savefig.dpi':300})
COLORS=['#222222','#1673b1','#d87819']

def main():
    z=np.load(HERE/'STABILIZATION_TERMS.npz');values=z['values'];proto=read(PREVIOUS/'EVENT_FIELD_PROTOCOL.json')['events'];meta={}
    fig,ax=plt.subplots(4,2,figsize=(17,16),constrained_layout=True)
    for row,pid in enumerate(PIDS):
        kind='rise' if pid==822 else 'accepted';left,right=ax[row];twin=right.twinx()
        for ni,name in enumerate(NAMES):
            ix=np.flatnonzero((z['field']==name)&(z['pid']==pid)&(z['kind']==kind));x=z['arc_m'][ix]*1e3;v=values[ix]
            if name=='M0':
                for j,(label,col) in enumerate(zip(['kT / total','kU / total','kS / total'],['#ba274a','#118d85','#8b67a7'])):left.plot(x,v[:,5+j]*100,'.',ms=2.3,color=col,label='M0 '+label,rasterized=True)
            else:left.plot(x,v[:,5]*100,'.',ms=1.4,color=COLORS[ni],label=name.replace('M1_fixedplc_','')+' kT / total',rasterized=True)
            right.plot(x,v[:,3],'.',ms=1.8,color=COLORS[ni],label=['M0','medium','fine'][ni]+' tau_M',rasterized=True)
            twin.plot(x,v[:,4],ls='none',marker='x',ms=1.3,color=COLORS[ni],alpha=.55,label=['M0','medium','fine'][ni]+' tau_C',rasterized=True)
        left.set_yscale('log');left.set_ylim(1e-7,180);left.set_ylabel('Contribution (%) · log scale');left.set_title(f'ID{pid} · '+('fixed A → B, 168.6564 μm' if pid==822 else 'original accepted event path'))
        left.grid(alpha=.2);left.set_xlabel('Original path arc length (mm)');right.set_xlabel('Original path arc length (mm)');right.set_yscale('log');twin.set_yscale('log');right.set_ylabel('tau_M (m³ s / kg) · dots');twin.set_ylabel('tau_C (Pa s) · crosses');right.grid(alpha=.2)
        right.set_title('Same positions, each field’s own TET4 metric')
        if row==0:
            left.legend(loc='lower right',ncol=1);handles,labels=right.get_legend_handles_labels();h,l=twin.get_legend_handles_labels();right.legend(handles+h,labels+l,ncol=2,loc='lower right')
        if pid!=822:
            path=np.load(PREVIOUS/f'path_{pid}.npz');times=path['original_time_s'];arc=path['original_arc_m']*1e3
            for tag in ['A','B']:
                ii=np.flatnonzero(times==proto[str(pid)][tag]['time_s'])[0]
                left.axvline(arc[ii],color='.6',lw=.7);left.text(arc[ii],120,tag,ha='center')
    fig.suptitle('Existing BraVa fields: stabilization term balance | production CFD dt = 0.01 s\nDiscrete samples; no smoothing. STOP C: modified CFD NOT_RUN.',fontsize=17)
    fig.savefig(HERE/'01_STABILIZATION_TERM_BALANCE.png');plt.close(fig)
    r=np.load(HERE/'RT_REFERENCE.npz');res=read(HERE/'RT_REFERENCE.json');t=r['time_s'];fields=[r['original_path_u'],r['RT0_no_div_path_u'],r['RT_path_u']];labels=['Original M0 P1','RT0, no div constraint','RT0, div = 0'];colors=['#222222','#d87819','#187fb1']
    fig,ax=plt.subplots(2,3,figsize=(18,10),constrained_layout=True)
    for u,label,c in zip(fields,labels,colors):
        ax[0,0].plot(t,np.linalg.norm(u,axis=1)*1e3,color=c,lw=.9,label=label)
        q=np.einsum('ij,ij->i',u,r['inward_normals'])*1e3;ax[0,1].plot(t,q,color=c,lw=.9,label=label)
    ax[0,0].set(title='ID1189 · same accepted positions',xlabel='Original physical time (s)',ylabel='Speed (mm/s)');ax[0,0].legend()
    ax[0,1].set(title='Slow / release · original inward normals',xlabel='Original physical time (s)',ylabel='Signed normal velocity (mm/s)',xlim=(t.min(),3.969));ax[0,1].axhline(0,color='.6',lw=.6)
    a=res['representatives']['A'];x=np.arange(3)
    for j,(key,color) in enumerate(zip(a,colors)):
        v=a[key];ax[0,2].bar(x+(j-1)*.25,[v['speed_mm_s'],v['inward_normal_mm_s'],v['tangent_mm_s']],width=.25,color=color,label=labels[j])
    ax[0,2].set(title='Exact A · representation changes are large',ylabel='mm/s',xticks=x,xticklabels=['Speed','Inward normal','Tangent']);ax[0,2].axhline(0,color='.6',lw=.6)
    cv=res['control_volumes'];x=np.arange(3)
    ax[1,0].bar(x-.17,[c['original_Qnet_m3_s']*1e9 for c in cv],width=.34,color=colors[0],label='Original M0')
    ax[1,0].bar(x+.17,[c['RT_Qnet_m3_s']*1e9 for c in cv],width=.34,color=colors[2],label='RT0 div=0')
    ax[1,0].set(title='Same closed control volumes · ALL faces',ylabel='Net outward flux (10⁻⁹ m³/s)',xticks=x,xticklabels=['Layer 0','Layer 1','Layer 2']);ax[1,0].legend();ax[1,0].text(.03,.08,'RT |Qnet| ≤ 2.28 × 10⁻²⁴ m³/s\nEnforced by construction',transform=ax[1,0].transAxes,fontsize=11)
    cut=r['cut_faces'];old=r['original_flux_m3_s'][cut]*1e9;new=r['RT_flux_m3_s'][cut]*1e9;lo=min(old.min(),new.min());hi=max(old.max(),new.max())
    ax[1,1].scatter(old,new,s=13,alpha=.6,color=colors[2]);ax[1,1].plot([lo,hi],[lo,hi],color='.5',lw=.7);ax[1,1].set(title='390 artificial cut faces · allowed to change',xlabel='Original outward flux (10⁻⁹ m³/s)',ylabel='RT outward flux (10⁻⁹ m³/s)')
    ax[1,2].axis('off');ax[1,2].text(0,1,'DIAGNOSTIC ONLY\n\nOriginal 468 tetrahedra; no new mesh.\n132 WALL faces: normal flux = 0.\n390 artificial cut faces: FREE flux.\n\nOriginal ROI net = −2.043 × 10⁻⁹ m³/s.\nAll original cut fluxes cannot be retained.\nCut flux L1 change = 4.637 × 10⁻⁸ m³/s.\n\nRT0 div=0 ⇒ constant velocity per tet.\nTangential no-slip is not preserved.\nG / curl / shear are not particle-compatible.\nNo trajectories; no modified CFD.',va='top',fontsize=12)
    for a0 in ax.flat:
        if a0.axison:a0.grid(axis='y',alpha=.15)
    fig.suptitle('Local conservative reference · original BraVa ID1189\nLibrary RT0/P0 constrained L2 projection; NOT a validated production field',fontsize=17)
    fig.savefig(HERE/'05_CONSERVATIVE_REFERENCE.png');plt.close(fig)
    put('FIGURE_METADATA.json',dict(dpi=300,plot_source='plot_results.py',figures={'01_STABILIZATION_TERM_BALANCE.png':['STABILIZATION_TERMS.npz','STABILIZATION_TERMS.json'],'05_CONSERVATIVE_REFERENCE.png':['RT_REFERENCE.npz','RT_REFERENCE.json']},smoothing=False,omitted_failed_samples=False,new_CFD=0,conditional_figures_not_created=['02_PIPE_BEFORE_AFTER.png','03_BRAVA_LOCAL_BEFORE_AFTER.png','04_CONTROL_VOLUME_BEFORE_AFTER.png']))
    print('FIGURES_COMPLETE')
if __name__=='__main__':main()
