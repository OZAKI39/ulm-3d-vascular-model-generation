"""Scientific figures from saved task diagnostics; no evaluations or solves."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from scripts.finite_domain_boundary_audit import OUT,RESULTS,save
FOLDER=RESULTS/'figures';FOLDER.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.22,'savefig.dpi':180})
COLORS={'C5':'#2166ac','C7':'#b2182b','L050':'#1b7837','L075':'#762a83'}
def data(name,what):return list(csv.DictReader((RESULTS/'tables'/f'{name}_{what}.csv').open()))
def col(rows,key):return np.array([float(r[key]) for r in rows])
def summary(name):return json.loads((OUT/name/'SUMMARY.json').read_text())
written=[]
def finish(fig,name):
    fig.tight_layout()
    for ext in ['png','pdf']:
        path=FOLDER/f'{name}.{ext}';fig.savefig(path,bbox_inches='tight');written.append(str(path.relative_to(OUT.parents[1])))
    plt.close(fig)

def main():
    fig,axs=plt.subplots(1,2,figsize=(12,4.5))
    for name in ['C5','C7']:
        d=data(name,'layers');z=col(d,'z_m')*1e6
        axs[0].plot(z,col(d,'layer_fz_over_R'),'o-',label=name,color=COLORS[name])
        axs[1].plot(z,col(d,'layer_L2_over_R'),'o-',label=name,color=COLORS[name])
    for ax in axs:
        ax.set_xlabel('Reference z [um]');ax.legend();ax.axhline(0,color='k',lw=.7)
    axs[0].set_ylabel(r'$\sum_{layer} f_z^{int}/R_{cont}$');axs[1].set_ylabel(r'$\|f^{int}_{layer}\|_2/R_{cont}$')
    axs[0].set_title('Signed axial force: affine is not equilibrium');axs[1].set_title('Layer L2 force')
    for name in ['C5','C7']:
        n=summary(name)['n'];d=data(name,'layers')
        for k,text in [(n-1,'surface'),(n-2,'first'),(n-3,'second')]:
            offset={('C5','surface'):(-65,-13),('C7','surface'):(-65,14),
                    ('C5','first'):(-55,-4),('C7','first'):(-4,14),
                    ('C5','second'):(-35,12),('C7','second'):(-45,-17)}[(name,text)]
            axs[0].annotate(name+' '+text,(float(d[k]['z_m'])*1e6,float(d[k]['layer_fz_over_R'])),xytext=offset,textcoords='offset points',fontsize=8,color=COLORS[name])
    finish(fig,'affine_internal_force_layers')

    fig,axs=plt.subplots(1,3,figsize=(14,4.5),sharey=True)
    for ax,name in zip(axs,['C5','C7','L075']):
        rows=data(name,'cumulative');s=summary(name);q=s['top_required']
        for side,style in [('top','-'),('bottom','--')]:
            d=[r for r in rows if r['side']==side]
            ax.plot(col(d,'depth_over_delta'),col(d,'cumulative_over_R'),'o'+style,label=side,color=COLORS[name] if side=='top' else '0.3')
        ax.axhline(1,color='k',lw=.8,label='continuum')
        qs=[q[f'd{k}_m']*1e6 for k in [50,90,95,99]]
        ax.text(.04,.07,'d50/90/95/99 [um]\n'+'/'.join(f'{v:g}' for v in qs)+'\n'+('monotone' if q['monotone'] else 'signed overshoot: first crossings'),transform=ax.transAxes,fontsize=9)
        ax.set_xlabel(r'Depth / $\delta$');ax.set_title(f'{name}: final = {q["ratio_to_continuum"]:.5f}');ax.legend(fontsize=8)
        ax.set_ylim(0,1.65)
    axs[0].set_ylabel(r'Cumulative required outward force / $R_{cont}$')
    fig.suptitle('Quantiles refer to observed half-reaction; curves retain the continuum normalization',y=1.03)
    finish(fig,'cumulative_nonlocal_boundary_force')

    fig,axs=plt.subplots(1,3,figsize=(14,4.5))
    names=['C5','C7'];x=np.arange(2)
    axs[0].bar(x-.18,[summary(n)['first_layer'][0]['direct_fraction'] for n in names],.36,label='DIRECT')
    axs[0].bar(x+.18,[summary(n)['first_layer'][0]['Z_fraction'] for n in names],.36,label='Z CORRECTION')
    axs[0].set_xticks(x,names);axs[0].set_title('Bottom first-interior center');axs[0].legend()
    groups=['loaded_surface','same_first_interior_layer','deeper_specimen']
    for ax,name in zip(axs[1:],names):
        d=[r for r in data(name,'bond_groups') if r['side']=='bottom' and r['grouping']=='neighbor_group'];d={r['group']:r for r in d};den=summary(name)['first_layer'][0]['nodal_fz_N']
        for offset,key,label in [(-.18,'direct_fz_N','DIRECT'),(.18,'Z_fz_N','Z CORRECTION')]:
            ax.bar(np.arange(3)+offset,[float(d[g][key])/den for g in groups],.36,label=label)
        ax.set_xticks(range(3),['loaded\nsurface','same first\nlayer','deeper']);ax.set_title(name+' incident bonds');ax.legend(fontsize=8)
    for ax in axs:ax.axhline(0,color='k',lw=.8);ax.set_ylabel('Signed contribution / net particle force')
    finish(fig,'first_layer_force_decomposition')

    fig,axs=plt.subplots(1,3,figsize=(14,4.5));norms=json.loads((OUT/'NORMALIZATION_CHANNELS.json').read_text())['cases']
    for ax,keys,labels,title in [(axs[0],['index_i_half','index_j_half'],['i half','j half'],'Arbitrary unordered bond endpoint index'),
                              (axs[1],['loaded_surface_endpoint','nonloaded_endpoint'],['loaded surface endpoint','nonloaded endpoint'],'Geometric endpoint classes')]:
        for k,lab,off in zip(keys,labels,[-.18,.18]):ax.bar(x+off,[summary(n)['first_layer'][0]['endpoint_signed_fractions'][k] for n in names],.36,label=lab)
        ax.set_title(title,fontsize=10);ax.legend(fontsize=8)
    for key,lab,off in [('uniform_wbar_component_signed_fraction','lattice-weight component',-.18),('normalization_excess_component_signed_fraction','normalization excess',.18)]:
        axs[2].bar(x+off,[next(c for c in norms if c['case']==n)['rows'][0][key] for n in names],.36,label=lab)
    axs[2].set_title('Exact algebraic split; no model change',fontsize=10);axs[2].legend(fontsize=8)
    for ax in axs:ax.set_xticks(x,names);ax.set_ylabel('Signed fraction of first-layer force');ax.axhline(0,color='k',lw=.8)
    fig.suptitle('All C5/C7 endpoints have truncated families; loaded surface is a narrower class',y=1.03)
    finish(fig,'endpoint_normalization_decomposition')

    fig,axs=plt.subplots(2,3,figsize=(14,8),sharey=True)
    for k,name in enumerate(['C5','L050','L075']):
        d=data(name,'particles');full=np.array([r['full_support']=='True' for r in d]);y=np.maximum(abs(col(d,'fz_over_R')),1e-18)
        for a,key in [(0,'d_over_delta'),(1,'loading_d_over_delta')]:
            ax=axs[a,k]
            for mask,color,label in [(~full,'#777777','truncated family'),(full,'#d73027','full family')]:
                if mask.any():ax.scatter(col(d,key)[mask],y[mask],s=14,alpha=.55,c=color,label=label)
            ax.axvline(1,color='k',ls='--',lw=.7);ax.axhline(1e-10,color='navy',ls=':',lw=.7,label='noise gate')
            ax.set_yscale('log');ax.set_ylim(1e-18,.2);ax.set_xlabel(r'$d_{min}/\delta$' if a==0 else r'$d_{loading}/\delta$');ax.set_ylabel(r'$|f_z^{int}|/R_{cont}$')
            ax.set_title(f'{name}: L={summary(name)["L_m"]*1e3:g} mm; full={full.sum()}');ax.legend(fontsize=7)
    finish(fig,'large_domain_affine_patch_test')

    d=data('L075','particles');xyz=np.column_stack([col(d,k) for k in ['X_m','Y_m','Z_m']])*1e6;fz=col(d,'required_external_z_N');lim=max(abs(fz));cn=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim)
    fig=plt.figure(figsize=(15,5));ax=fig.add_subplot(131,projection='3d');sc=ax.scatter(*xyz.T,c=fz,s=9,cmap='coolwarm',norm=cn,alpha=.55)
    ax.set(xlabel='X [um]',ylabel='Y [um]',zlabel='Z [um]',title='Required force: all points')
    inner=(.00075/2-DELTA)*1e6
    for a in [-inner,inner]:
        for b in [-inner,inner]:
            ax.plot([-inner,inner],[a,a],[b,b],color='k',lw=1)
            ax.plot([a,a],[-inner,inner],[b,b],color='k',lw=1)
            ax.plot([a,a],[b,b],[-inner,inner],color='k',lw=1)
    ax2=fig.add_subplot(132);mask=abs(xyz[:,1])<1e-6
    ax2.scatter(xyz[mask,0],xyz[mask,2],c=fz[mask],s=50,cmap='coolwarm',norm=cn)
    for line in [-inner,inner]:ax2.axhline(line,color='k',ls='--',lw=.8);ax2.axvline(line,color='k',ls='--',lw=.8)
    ax2.set(xlabel='X [um]',ylabel='Z [um]',title=r'$y=0$; dashed: $\delta$-collar / core',aspect='equal')
    ax3=fig.add_subplot(133);layers=data('L075','layers');ax3.bar(col(layers,'z_m')*1e6,-col(layers,'layer_fz_over_R'),width=42,color=COLORS['L075'])
    ax3.axvspan(-inner,inner,color='orange',alpha=.15,label='loading-plane depth >= delta');ax3.set(xlabel='Z [um]',ylabel=r'$\sum f_z^{ext,req}/R_{cont}$',title='Signed required force per layer');ax3.legend(fontsize=7)
    fig.suptitle('MANUFACTURED NONLOCAL BOUNDARY FORCE\nNOT PHYSICAL FRACTURE LOAD',y=1.03,fontsize=13)
    # Manual spacing handles the combined 3D/colorbar axes without tight_layout overlap.
    fig.subplots_adjust(bottom=.28,wspace=.35,top=.85)
    color_axis=fig.add_axes([.14,.075,.45,.03])
    fig.colorbar(sc,cax=color_axis,label='Required external fz [N]',orientation='horizontal')
    for ext in ['png','pdf']:
        path=FOLDER/f'manufactured_external_force_map.{ext}';fig.savefig(path,bbox_inches='tight');written.append(str(path.relative_to(OUT.parents[1])))
    plt.close(fig)
    save(OUT/'FIGURES.json',dict(files=written,PNG_count=6,PDF_count=6,boundary_treatment_comparison='NOT_RUN: bulk patch gate failed'))

if __name__=='__main__':
    from scripts.finite_domain_boundary_audit import DELTA
    main()
