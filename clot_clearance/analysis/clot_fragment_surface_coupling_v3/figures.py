"""Eight figure groups from actual coordinates and accepted scalar histories."""
from common import *
from analyze import read_history
from current_surface import CurrentSurface
import csv,argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection,Poly3DCollection

plt.rcParams.update({'font.size':12,'axes.titlesize':14,'axes.labelsize':12,'legend.fontsize':11,'pdf.fonttype':42,'axes.grid':True,'grid.color':'#9ca6ad','grid.alpha':.5})
OLD='#bd704d';NEW='#217e9a';GREEN='#32977d';GRAY='#647582'

def save(fig,name):
    (F/'figures').mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(F/'figures'/(name+'.'+ext),dpi=200)
    plt.close(fig)

def equal_world(ax,x,span=None):
    mid=(x.min(0)+x.max(0))/2;span=np.ptp(x,axis=0).max()*1.2 if span is None else span
    for j,setlim in enumerate([ax.set_xlim,ax.set_ylim,ax.set_zlim]):setlim(mid[j]-span/2,mid[j]+span/2)
    ax.set_box_aspect((1,1,1));ax.view_init(16,-64);ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)');ax.tick_params(labelsize=9)

def p0_redraw():
    diagnosis=json.loads((O/'P0_TARGET_SURFACE_DIAGNOSIS.json').read_text())
    fig=plt.figure(figsize=(22,12),layout='constrained')
    for j,row in enumerate(diagnosis['targets']):
        q=load_npz(O/'p0_targets'/f'{row["case"]}_{row["component_id"]}.npz');x=q['x']*1e6;lookup={int(a):k for k,a in enumerate(q['original_atom_ids'])}
        edges=np.array([[lookup[int(a)] for a in e] for e in q['live_bond_atom_ids']]);ax=fig.add_subplot(2,4,j+1,projection='3d')
        ax.add_collection3d(Line3DCollection(x[edges],colors='#6d8797',linewidths=.45,alpha=.45));ax.scatter(*x.T,s=10,c='#215570',depthshade=False)
        for mask,color,marker in [(q['active'],GREEN,'^'),(~q['active']&~q['new_open'], '#de9d28','x'),(~q['active']&q['new_open'],'#9863b3','D')]:
            if mask.any():
                faces=np.array([[lookup[int(a)] for a in face] for face in q['face_atom_ids'][mask]])
                ax.add_collection3d(Poly3DCollection(x[faces],facecolors=color,edgecolors=color,alpha=.45,linewidths=.7))
                points=q['centroid'][mask]*1e6;ax.scatter(*points.T,c=color,s=42,marker=marker,depthshade=False)
        equal_world(ax,x);c=row['counts'];ax.set_title(f'{row["case"].title()} ID {row["component_id"]}: {row["nodes"]} nodes\n{c["supported_tetra"]} supported tets; {c["active"]} loaded facets',pad=4)
    ax=fig.add_subplot(248);ax.axis('off');ax.text(.02,.98,'Actual V2 end states\n\nBlue: real material points\nand surviving original bonds\n\nOrange crosses: unopened candidates\nPurple diamonds: open, inaccessible\nGreen triangles: loaded face centroids\n(markers are identification symbols)\n\n4 fragments: no supported tetra\n1 zero-load fragment: opening test fails\n2 others: opening/accessibility limits\n\nNo geometric shell was invented.\nNo loaded facet ≠ no physical drag.\nOriginal world coordinates throughout.',va='top',fontsize=12,linespacing=1.4)
    fig.suptitle('FIG 01 | Why connected fragments lost the old surface-load representation',fontsize=21)
    save(fig,'FIG_01_SURFACE_REPRESENTATION_FAILURE')

def quadrature_figure():
    with (O/'SURFACE_QUADRATURE_AUDIT.csv').open() as f:rows=list(csv.DictReader(f))
    coverage=json.loads((O/'SURFACE_COVERAGE_DEFINITIONS.json').read_text())['targets'];names=[('U' if r['case']=='uniform' else 'H')+' '+r['component_id'] for r in rows];idx=np.arange(len(rows))
    fig=plt.figure(figsize=(17,12),layout='constrained');ax=fig.add_subplot(221)
    ax.bar(idx-.18,[float(r['old_area_m2'])*1e6 for r in rows],.36,color=OLD,label='V2 effective area');ax.bar(idx+.18,[float(r['new_area_m2'])*1e6 for r in rows],.36,color=NEW,label='V3 effective area')
    ax.set(xticks=idx,xticklabels=names,ylabel='Effective loaded area (mm²)',title='Same saved fragments, different quadrature');ax.tick_params(axis='x',rotation=25);ax.legend()
    ax=fig.add_subplot(222);ax.bar(idx-.24,[float(r['old_node_coverage'])*100 for r in rows],.24,color=OLD,label='V2 face vertices')
    ax.bar(idx,[r['quadrature_owner_or_vertex_fraction']*100 for r in coverage],.24,color=GREEN,label='V3 site owners / vertices')
    ax.bar(idx+.24,[r['mapped_force_support_fraction']*100 for r in coverage],.24,color=NEW,label='V3 force-weight support')
    ax.set(xticks=idx,xticklabels=names,ylabel='Fraction of fragment nodes (%)',ylim=(0,105),title='Force support and geometric sites are distinct');ax.tick_params(axis='x',rotation=25);ax.legend(loc='upper left',fontsize=10)
    ref=load_npz(V2_O/'provenance/reference.npz');s=load_npz(V2_O/'cases/biphasic/keys/final.npz');reg=json.loads((V2_O/'cases/biphasic/COMPONENT_REGISTRY.json').read_text());ids=np.array(reg['3769']['atom_ids'])-1
    surf=CurrentSurface(ref,'uniform');z=surf.geometry(s['x'],s['d']);member=np.zeros(surf.N,bool);member[ids]=True;internal=member[ref['bonds']].all(1)&(s['d']<1)
    for col,new in enumerate([False,True]):
        ax=fig.add_subplot(223+col,projection='3d');x=s['x']*1e6
        ax.add_collection3d(Line3DCollection(x[ref['bonds'][internal]],colors='#8a9da7',linewidths=.4,alpha=.4));ax.scatter(*x[ids].T,s=9,c='#28536a',depthshade=False)
        mask=z['active'][:surf.L]&member[surf.base.tri].all(1)
        if mask.any():ax.add_collection3d(Poly3DCollection(x[surf.base.tri[mask]],facecolors=GREEN,edgecolors=GREEN,alpha=.9))
        if new:
            patch=np.flatnonzero(z['active'][surf.L:]&member[surf.owner]);quads=[]
            for p in patch:
                i=surf.owner[p];axis=surf.axis[p];j,k=(axis+1)%3,(axis+2)%3;F0=z['l2_deformation'][i]
                corners=np.array([surf.offset[p]+DX*(a*np.eye(3)[j]+b*np.eye(3)[k]) for a,b in [(-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5)]])
                quads.append((s['x'][i]+corners@F0.T)*1e6)
            ax.add_collection3d(Poly3DCollection(quads,facecolors='#4ba9bf',edgecolors='#14566a',linewidths=.5,alpha=.7))
        equal_world(ax,x[ids],620);ax.set_title('H ID 3769, 148 actual nodes\n'+('V3: retained L1 + partial dual-face patches' if new else 'V2: one active supported triangle'))
    fig.suptitle('FIG 03 | Surface quadrature replay — approximate, conservative, and still incomplete',fontsize=20)
    audit=json.loads((O/'SURFACE_GEOMETRY_AND_MAPPING_AUDIT.json').read_text())
    unresolved='; '.join(f'{name}: {r["unresolved_patches"]}/{r["candidate_patches"]}' for name,r in zip(names,rows))
    errors={key:max(c['mapping'][key] for c in audit['checks']) for key in ['force_N','moment_Nm','power_W']}
    fig.get_layout_engine().set(rect=(0,.075,1,.925))
    fig.text(.5,.048,'Unresolved / candidate patches — '+unresolved,ha='center',fontsize=10)
    fig.text(.5,.028,f'Max mapping errors: force {errors["force_N"]:.2e} N; moment {errors["moment_Nm"]:.2e} N m; power {errors["power_W"]:.2e} W.',ha='center',fontsize=11)
    fig.text(.5,.008,'Bounded affine weights preserve resultant, moment, and power; neither coverage metric is a continuum wetted-area fraction.',ha='center',fontsize=11)
    save(fig,'FIG_03_SURFACE_QUADRATURE_BEFORE_AFTER')

def native_figures():
    result=json.loads((O/'NATIVE_V3_RESULT.json').read_text());lineage=json.loads((O/'FRAGMENT_LINEAGE.json').read_text());new=read_history(O/'cases/biphasic/history.csv');old=read_history(V2_O/'cases/biphasic/history.csv')
    def plot(ax,rows,key,color,label,scale=1):ax.plot([r['time_s']*1e3 for r in rows],[r[key]*scale for r in rows],color=color,label=label,lw=1.6)
    def timeaxis(ax):ax.axvline(2.8,color='#666',ls='--',lw=.8);ax.set_xlabel('Physical time (ms)')
    fig,axs=plt.subplots(3,2,figsize=(16,12),layout='constrained');axs=axs.ravel()
    plot(axs[0],old,'effective_loaded_area_m2',OLD,'V2 total',1e6);plot(axs[0],new,'effective_loaded_area_m2',NEW,'V3 total',1e6);plot(axs[0],new,'L2_effective_area_m2',GREEN,'V3 supplement',1e6);axs[0].set(ylabel='Effective area (mm²)',title='Evolving load representation');axs[0].legend()
    plot(axs[1],new,'L2_active_patches',NEW,'active dual patches');plot(axs[1],new,'L2_unresolved_patches',GRAY,'unresolved dual candidates');plot(axs[1],new,'fluid_accessible_new_faces',GREEN,'accessible new faces, both methods');axs[1].set(ylabel='Patch count',title='Active and unresolved surface candidates');axs[1].legend()
    plot(axs[2],new,'shear_peak_Pa',NEW,'actual peak slip shear');plot(axs[2],new,'nominal_shear_Pa',GRAY,'prescribed nominal drive');axs[2].set(ylabel='Shear (Pa)',title='Relative velocity changes actual traction');axs[2].legend()
    plot(axs[3],new,'net_force_x_N',NEW,'total applied Fx',1e6);axs[3].set(ylabel='Resultant force (µN)',title='Whole-material surface resultant')
    target=result['final_target_components'];chosen=max(target,key=lambda q:q['nodes'])['component_id'] if target else None
    fluid=read_history(O/'FRAGMENT_FORCE_AND_POWER.csv') if target else []
    track=[r for r in fluid if r.get('component_id')==chosen]
    if track:
        t=[r['time_s']*1e3 for r in track];force=[np.linalg.norm(json.loads(r['fluid_force_N']))*1e6 for r in track];moment=[np.linalg.norm(json.loads(r['fluid_moment_COM_Nm']))*1e9 for r in track]
        axs[4].plot(t,force,color=NEW);axs[5].plot(t,moment,color=GREEN);axs[4].set(ylabel='Force magnitude (µN)',title=f'Final target ID {chosen}: exact-set lifetime');axs[5].set(ylabel='COM torque magnitude (nN m)',title='Torque uses the current mass center')
        slip=axs[4].twinx();slip.plot(t,[r['fluid_area_mean_tangential_slip_m_s'] for r in track],color=OLD,ls='--',lw=1.2);slip.set_ylabel('Area-mean relative slip (m/s)',color=OLD);slip.grid(False)
    else:
        for ax in axs[4:]:ax.text(.5,.5,'No final 20–200-node target observed',transform=ax.transAxes,ha='center')
    for ax in axs:timeaxis(ax)
    fig.suptitle('FIG 04 | V3 native dynamics and relative fluid loading',fontsize=21);save(fig,'FIG_04_NATIVE_SURFACE_LOADING')

    hist=json.loads((V2_O/'NATIVE_BIPHASIC_RESULT.json').read_text());old_sizes=np.array(hist['final_free_sizes']);new_sizes=np.array(result['final_free_sizes'])
    fig,axs=plt.subplots(2,2,figsize=(16,10),layout='constrained');names=['1–2','3–13','14–19','20–33','34–113','114–200','>200'];idx=np.arange(7)
    for sizes,delta,color,label in [(old_sizes,-.18,OLD,'V2 H'),(new_sizes,.18,NEW,'V3 H')]:
        counts=[];mass=[]
        for lo,hi in BINS:
            q=sizes[(sizes>=lo)&((sizes<=hi) if hi is not None else True)];counts.append(len(q));mass.append(q.sum()/8125*100)
        axs[0,0].bar(idx+delta,np.where(np.array(counts)>0,counts,np.nan),.36,color=color,label=label)
        axs[0,1].bar(idx+delta,mass,.36,color=color,label=label)
        axs[1,1].plot(np.arange(1,len(sizes)+1),np.sort(sizes)[::-1],'.-',ms=3,color=color,label=label)
    for ax in axs[0]:ax.set_xticks(idx,names);ax.legend()
    axs[0,0].set(yscale='log',ylabel='Free component count (zero bins omitted)',title='Unchanged original-material definition');axs[0,1].set(ylabel='Fraction of initial material mass (%)',title='Count and mass answer different questions')
    plot(axs[1,0],old,'N20_200',OLD,'V2 H');plot(axs[1,0],new,'N20_200',NEW,'V3 H');axs[1,0].axhline(5,color='#555',ls=':');axs[1,0].set(ylabel='20–200-node free components',title='Inclusive target interval');axs[1,0].legend();timeaxis(axs[1,0])
    axs[1,1].axhspan(20,200,color='#a0c8cd',alpha=.2);axs[1,1].set(xscale='log',yscale='log',xlabel='Final component rank',ylabel='Original nodes',title='All final component sizes');axs[1,1].legend()
    fig.suptitle('FIG 05 | Historical H versus new H: surface algorithms differ',fontsize=21);save(fig,'FIG_05_FRAGMENT_SIZE_DISTRIBUTION')

    fig,axs=plt.subplots(1,2,figsize=(17,9),layout='constrained');tracks=sorted(lineage['target_tracks'],key=lambda q:(q['birth_s'],q['component_id']))
    for j,q in enumerate(tracks):
        axs[0].plot([q['birth_s']*1e3,q['end_s']*1e3],[j,j],color=NEW if q['survives_exact_set_to_end'] else GRAY,lw=2)
        if q['survives_exact_set_to_end']:axs[0].plot(q['end_s']*1e3,j,'o',color=NEW,ms=4)
    axs[0].axvline(2.8,color='#555',ls='--');axs[0].set(xlabel='Physical time (ms)',ylabel='20–200-node exact member-set index',title=f'All {len(tracks)} target member-set lifetimes')
    unload=lineage['unload_targets']
    if unload:
        labels=[f'ID {q["component_id_at_drive_off"]}\n{q["nodes_at_drive_off"]} nodes' for q in unload];axs[1].bar(np.arange(len(unload)),[q['max_member_retention']*100 for q in unload],color=GREEN)
        axs[1].set(xticks=np.arange(len(unload)),xticklabels=labels,ylim=(0,105),ylabel='Largest-descendant member retention (%)',title='Drive-off targets at the last accepted state')
        for i,q in enumerate(unload):axs[1].text(i,q['max_member_retention']*100+2,f'{q["target_descendants"]} target descendants',ha='center',fontsize=9,rotation=15)
    else:axs[1].text(.5,.5,'No target at drive-off / drive-off not reached',transform=axs[1].transAxes,ha='center')
    fig.suptitle('FIG 06 | Exact-set lifetime, erosion, and descendant survival are distinct',fontsize=20);save(fig,'FIG_06_FRAGMENT_LINEAGE_AND_SURVIVAL')

    energy=read_history(O/'ENERGY_LEDGER.csv');fig,axs=plt.subplots(2,2,figsize=(16,10),layout='constrained')
    for key,label,color in [('W_drive_J','drive',NEW),('W_drag_J','passive / relative drag',OLD),('W_J','net external work','#222')]:plot(axs[0,0],energy,key,color,label,1e9)
    axs[0,0].set(ylabel='Work (nJ)',title='Signed external work');axs[0,0].legend()
    for key,label,color in [('DeltaU_J','ΔU',NEW),('DeltaK_J','ΔK',GREEN),('D_J','D',OLD),('DeltaUc_J','ΔUc','#9d78a4'),('Jc_J','Jc (subtract in balance)',GRAY)]:plot(axs[0,1],energy,key,color,label,1e9)
    axs[0,1].set(ylabel='Energy (nJ)',title='Stored, kinetic, damage, and contact channels');axs[0,1].legend()
    plot(axs[1,0],energy,'R_J',OLD,'R',1e9);axs[1,0].set(ylabel='Residual R (nJ)',title=f'Separate residual; max normalized {100*max(r["relative_R"] for r in energy):.4f}%')
    plot(axs[1,1],energy,'W_accessible_new_faces_J',NEW,'accessible new-face work',1e9);plot(axs[1,1],energy,'W_L2_surface_J',GREEN,'L2 work',1e9)
    axs[1,1].set(ylabel='Diagnostic subchannel work (nJ)',title='Overlapping subsets of total work — do not add again');axs[1,1].legend()
    for ax in axs.ravel():timeaxis(ax)
    fig.suptitle('FIG 07 | Wdrive + Wdrag = ΔU + ΔK + D + ΔUc − Jc + R',fontsize=21);save(fig,'FIG_07_FULL_ENERGY_LEDGER')

    fig,axs=plt.subplots(2,2,figsize=(17,11),layout='constrained');axs[0,0].axis('off');q=result['objective'];oldq=hist['objective']
    table=[['20–200-node targets',str(oldq['n_fragments_20_200']),str(q['n_fragments_20_200'])],['largest free component',str(max(old_sizes)),str(max(new_sizes,default=0))],['target / initial mass',f'{100*oldq["mass_fraction_target_initial"]:.3f}%',f'{100*q["mass_fraction_target_initial"]:.3f}%'],['1–2 fines / detached mass',f'{100*oldq["fine_debris_mass_fraction_detached"]:.2f}%',f'{100*(q["fine_debris_mass_fraction_detached"] or 0):.2f}%'],['detached / initial mass',f'{100*oldq["total_detached_mass_fraction"]:.2f}%',f'{100*q["total_detached_mass_fraction"]:.2f}%']]
    tab=axs[0,0].table(cellText=table,colLabels=['Measure','V2 H','V3 H'],cellLoc='center',loc='center',colWidths=[.54,.23,.23]);tab.auto_set_font_size(False);tab.set_fontsize(12);tab.scale(1,2.2);axs[0,0].set_title('Actual final outcomes')
    cost=json.loads((O/'COMPUTE_COST.json').read_text());axs[0,1].bar(['V2 H CPU','V3 H CPU','V3 H wall'],[cost['old_H_reference_cpu_s']/3600,cost['cpu_s']/3600,cost['wall_s']/3600],color=[OLD,NEW,GRAY]);axs[0,1].set(ylabel='Hours',title='Measured cost; V2 CPU excludes its user pause')
    if target:
        idx=np.arange(len(target));axs[1,0].bar(idx-.18,[100*t['quadrature_owner_or_vertex_fraction'] for t in target],.36,color=GREEN,label='site owners / vertices');axs[1,0].bar(idx+.18,[100*t['fluid_loaded_node_fraction'] for t in target],.36,color=NEW,label='force-weight support');axs[1,0].set(xticks=idx,xticklabels=[str(t['component_id']) for t in target],ylim=(0,105),ylabel='Fraction of fragment nodes (%)',xlabel='Final target ID',title='Coverage remains a discretization diagnostic');axs[1,0].legend()
    else:axs[1,0].text(.5,.5,'No final target in the specified interval',transform=axs[1,0].transAxes,ha='center')
    axs[1,1].axis('off');axs[1,1].text(0,.98,'Interpretation limits\n\n• Same H material and prescribed load; surface method changed.\n• Partial control-volume surface, bounded affine extrapolation.\n• Fines and uncertain geometry: FLUID_LOADING_UNRESOLVED.\n• Mechanical shear proxy; no validated bubble-flow coupling.\n• Material and correlation length remain uncalibrated.\n• Temporal association does not prove drag caused a split.\n• New H versus old Uniform is not a material-only control.',va='top',fontsize=12,linespacing=1.55)
    fig.suptitle('FIG 08 | Outcomes, computational cost, and remaining model limits',fontsize=21);save(fig,'FIG_08_OUTCOMES_COST_AND_LIMITATIONS')

def main():
    p=argparse.ArgumentParser();p.add_argument('--pre',action='store_true');args=p.parse_args()
    p0_redraw();quadrature_figure()
    if not args.pre:native_figures()

if __name__=='__main__':main()
