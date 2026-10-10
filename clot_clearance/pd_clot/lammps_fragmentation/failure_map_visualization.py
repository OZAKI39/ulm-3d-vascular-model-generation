"""Scientific figures and three actual-state movies; no invented outcomes."""
from pathlib import Path
import os,json,time
os.environ.setdefault('LIBGL_ALWAYS_SOFTWARE','1')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter
from .coherent_failure_card import OUTPUT as O,RESULTS as R
from .config import dump,PATTERNS

MODES=['NO_FAILURE','WHOLE_CLOT_DETACHMENT','RESOLVED_FRAGMENTATION','PULVERIZATION','HYBRID_DEBONDING_AND_FRAGMENTATION','CLASSIFICATION_UNCERTAIN']
SHORT=['No failure','Whole detachment','Resolved fragmentation','Pulverization','Hybrid','Uncertain']
COLORS=['#55a868','#4c72b0','#ddaa33','#c44e52','#8172b3','#888888']
PNAMES={'LOCAL_GAUSSIAN_SHEAR':'Local Gaussian shear','OPPOSING_SHEAR':'Opposing shear','VORTEX_TANGENTIAL':'Tangential vortex'}

def savefig(fig,name):
    (R/'figures').mkdir(exist_ok=True)
    for ext in ['png','pdf']:fig.savefig(R/'figures'/f'{name}.{ext}',dpi=180,bbox_inches='tight')
    plt.close(fig)

def failure_figures(data,c):
    b=json.loads((O/'BULK_CALIBRATION.json').read_text());i=json.loads((O/'INTERFACE_CALIBRATION.json').read_text())
    fig=plt.figure(figsize=(13,7));gs=fig.add_gridspec(2,2,height_ratios=[1.15,1]);a=fig.add_subplot(gs[0,:]);a.axis('off')
    rows=[['Whole-clot bulk','18.27 ± 1.45 J/m²','~4 kPa','4–8 kPa','6–8 J/m²'],['Clot / porcine skin','Lowest-rate peeling model','~1.7 kPa','1.7–3.5 kPa','1.3–1.6 J/m²']]
    table=a.table(cellText=rows,colLabels=['Source condition','Measured energy / test','Observed strength','Inferred CZ strength','Inferred intrinsic work'],loc='center',cellLoc='center',colWidths=[.18,.27,.17,.19,.19]);table.auto_set_font_size(False);table.set_fontsize(10);table.scale(1,2.4)
    a.text(.01,.96,'Liu et al., JMPS 193 (2024), 105858 | DOI 10.1016/j.jmps.2024.105858',fontsize=12,transform=a.transAxes)
    a.text(.01,.08,'Inner-aorta adhesion 1.35 ± 0.39 J/m² is a separate substrate measurement.\nThe inferred interface ranges above belong to the porcine-skin model.',fontsize=10,transform=a.transAxes)
    ax=fig.add_subplot(gs[1,0]);ax.set_title('Three fixed anchor layers; all clot nodes movable',fontsize=12)
    for j in range(3):
        for k in range(5):
            ax.scatter(k,j+.5,c='#4e86ae',s=65,zorder=3)
            if k<4:ax.plot([k,k+1],[j+.5,j+.5],c='#4e86ae',lw=1)
            if j<2:ax.plot([k,k],[j+.5,j+1.5],c='#4e86ae',lw=1)
    for j in range(3):
        for k in range(5):
            ax.scatter(k,-j-.5,c='#555555',s=65,marker='s',zorder=3)
            if k<4:ax.plot([k,k+1],[-j-.5,-j-.5],c='#555555',lw=1)
            if j<2:ax.plot([k,k],[-j-.5,-j-1.5],c='#555555',lw=1)
            ax.plot([k,k],[-j-.5,.5],c='#df8a23',lw=.8,alpha=.6)
    ax.text(4.4,1.6,'Clot / BULK',fontsize=10);ax.text(4.4,0,'Breakable INTERFACE',fontsize=10,color='#a65e00');ax.text(4.4,-2,'Fixed WALL\nexcluded from statistics',fontsize=10)
    ax.set(xlim=(-.5,7.5),ylim=(-3.2,3));ax.axis('off');ax.text(0,-3.1,'Topology schematic; not a solver state',fontsize=9)
    ax=fig.add_subplot(gs[1,1]);ax.axis('off');ax.text(0,.98,'Undamaged elastic mapping from SI Table S2',fontsize=13,weight='bold',va='top')
    ax.text(0,.81,f"G₀ = {c['G_Pa']:.0f} Pa; K₀ = {c['K_Pa']/1000:.3f} kPa\nR = 0.5625 mm; dx = 50 µm; horizon = 3.01dx",fontsize=11,va='top')
    ax.text(0,.53,f"Bulk: s00 = {b['s00']:.6f}; peak = {b['calibrated_macro_strength_Pa']/1000:.3f} kPa\nInterface: s00 = {i['s00']:.6f}; peak = {i['calibrated_macro_strength_Pa']/1000:.3f} kPa",fontsize=11,color='#206b2e',va='top')
    ax.text(0,.27,'Strength mappings in prescribed small cells;\nno fitted energy dissipation, Yeoh nonlinearity or Mullins effect.\nDensity and numerical drag remain engineering assumptions.',fontsize=10,va='top')
    fig.suptitle('Coherent failure source and breakable interface | EXPLORATORY ONLY',fontsize=16)
    fig.text(.5,.02,'LAMMPS critical-stretch parameters calibrated to strength ranges; fracture energies are source constraints, not directly inserted as s00.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.92,bottom=.09,hspace=.2,wspace=.16);savefig(fig,'coherent_failure_card_and_interface_model')
    fig,axs=plt.subplots(1,3,figsize=(14,5),sharex=True,sharey=True)
    for ax,pat in zip(axs,PATTERNS):
        for z in data:
            if z['pattern_id']!=pat:continue
            k=MODES.index(z['failure_mode']);marker='D' if z['map_stage']=='BOUNDARY_REFINEMENT' else 's'
            ax.scatter(z['duration_over_Tcross'],z['tau_budget'],c=COLORS[k],s=145 if marker=='s' else 90,marker=marker,edgecolors='black',linewidths=.6,zorder=3)
        ax.set(xscale='log',yscale='log',xlabel='Duration / T_cross',title=PNAMES[pat]);ax.set_xticks([.5,1,2,8],['0.5','1','2','8']);ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(which='both',alpha=.2);ax.axhline(10,color='#333333',ls='--',lw=.8,zorder=1)
    axs[0].set_ylabel('Reference mean tangential traction budget (Pa)')
    handles=[Line2D([],[],marker='s',linestyle='',color=color,label=name) for color,name in zip(COLORS,SHORT)]
    fig.legend(handles=handles,loc='lower center',ncol=3,bbox_to_anchor=(.5,.02),fontsize=9)
    fig.suptitle('LPS DURATION MAP — VISCOELASTIC RATE DEPENDENCE NOT INCLUDED',fontsize=13)
    fig.text(.5,.87,'Each symbol is an actual simulation; diamonds refine category changes. Dashed ~10 Pa line: literature scale, not a universal threshold.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.78,bottom=.25,wspace=.1);savefig(fig,'force_duration_failure_maps')
    fig,axs=plt.subplots(2,3,figsize=(13,8),sharex=True,sharey=True)
    for j,pat in enumerate(PATTERNS):
        subset=[z for z in data if z['pattern_id']==pat]
        for k,kind in enumerate(['bulk','interface']):
            ax=axs[k,j]
            for z in subset:
                if z['failure_mode']=='NO_FAILURE':ax.scatter(z['duration_over_Tcross'],z['tau_budget'],c=[z['U_'+kind+'_P99_9']],cmap='viridis',vmin=0,vmax=1,s=110,marker='s',edgecolors='black',lw=.4)
                else:ax.scatter(z['duration_over_Tcross'],z['tau_budget'],c='#cccccc',marker='x',s=35)
            ax.set(xscale='log',yscale='log',title=f'{PNAMES[pat]} | {kind}',xlabel='Duration / T_cross');ax.set_xticks([.5,1,2,8],['0.5','1','2','8']);ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(alpha=.2)
            if j==0:ax.set_ylabel('Traction budget (Pa)')
    fig.suptitle('Subcritical utilization: temporal peak of surviving-bond P99.9 | EXPLORATORY ONLY',fontsize=13)
    fig.text(.5,.025,'Colored squares: NO_FAILURE cases. Grey crosses: other categories. Saved-time maxima; no viscoelastic rate dependence.',ha='center',fontsize=10)
    fig.subplots_adjust(top=.9,bottom=.1,right=.9,hspace=.35,wspace=.18);cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='viridis'),cax=fig.add_axes([.92,.17,.015,.62]));cb.set_label('Tensile utilization P99.9 / critical stretch')
    savefig(fig,'subcritical_utilization_maps')

def select_representatives(data):
    selected={}
    for mode in MODES:
        group=[z for z in data if z['failure_mode']==mode]
        if not group:continue
        if mode=='NO_FAILURE':z=max(group,key=lambda a:max(a.get('U_bulk_P99_9') or 0,a.get('U_interface_P99_9') or 0))
        elif mode=='WHOLE_CLOT_DETACHMENT':z=max(group,key=lambda a:(a['interface_broken_fraction'],a['largest_free_fraction']))
        elif mode=='PULVERIZATION':z=max(group,key=lambda a:a['debris_fraction'])
        else:z=max(group,key=lambda a:(a['resolved_fragment_count'],-a['debris_fraction']))
        selected[mode]=z['case_id']
    chosen=[]
    transition=selected.get('NO_FAILURE') or min(data,key=lambda a:abs(1-max(a.get('U_bulk_P99_9') or 0,a.get('U_interface_P99_9') or 0)))['case_id']
    chosen.append(transition)
    for mode in ['RESOLVED_FRAGMENTATION','WHOLE_CLOT_DETACHMENT','HYBRID_DEBONDING_AND_FRAGMENTATION','PULVERIZATION','CLASSIFICATION_UNCERTAIN']:
        if len(chosen)==3:break
        if mode in selected and selected[mode] not in chosen:chosen.append(selected[mode])
    for z in sorted(data,key=lambda a:(-a['FI_free'],a['case_id'])):
        if len(chosen)==3:break
        if z['case_id'] not in chosen:chosen.append(z['case_id'])
    return selected,chosen

def load_scene(caseid,g):
    p=O/'cases'/caseid;z=np.load(p/'accepted_trajectory.npz');e=np.load(p/'broken_bonds.npz')
    return dict(case_id=caseid,trajectory={k:z[k] for k in z.files},events=e['events'],indices=e['initial_bond_indices'],request=json.loads((p/'request.json').read_text()),result=json.loads((p/'result.json').read_text()),geometry=g,metrics=np.genfromtxt(p/'fragment_metrics.csv',delimiter=',',names=True,dtype=None,encoding='utf8'))

def add_state(plot,s,index,bounds=None,overview=False):
    import pyvista as pv
    g=s['geometry'];z=s['trajectory'];X=z['position_m'][index]*1e3;n=int(z['clot_count']);dx=s['request']['dx_m']*1e3
    label=z['component_id'][index];rgb=(plt.cm.hsv(((label*0.61803398875)%1))[:,:3]*255).astype(np.uint8)
    cloud=pv.PolyData(X[:n]);cloud['rgb']=rgb
    if overview:plot.add_mesh(cloud,scalars='rgb',rgb=True,point_size=5,render_points_as_spheres=True)
    else:
        cubes=cloud.glyph(geom=pv.Cube(x_length=.9*dx,y_length=.9*dx,z_length=.9*dx),orient=False,scale=False);plot.add_mesh(cubes,scalars='rgb',rgb=True,show_scalar_bar=False)
    wall=pv.PolyData(X[n:]);plot.add_mesh(wall,color='#555555',point_size=5,render_points_as_spheres=True)
    alive=np.ones(len(g['bonds']),bool);ev=s['events'];alive[s['indices'][ev[:,1]<=z['step'][index]]]=False
    bonds=g['bonds'][alive];lines=pv.PolyData(X);lines.lines=np.c_[np.full(len(bonds),2),bonds].ravel();plot.add_mesh(lines,color='#525b66',opacity=.20,line_width=.35)
    if bounds is None:
        lo=X.min(0);hi=X.max(0);pad=np.maximum(.1*(hi-lo),dx);bounds=np.c_[lo-pad,hi+pad].ravel()
    plot.set_background('white');plot.view_isometric(render=False);plot.enable_parallel_projection();plot.reset_camera(bounds=bounds,render=False)
    # view_isometric(negative=False) clears PyVista's is_set flag. Freeze the
    # requested bounds so the first screenshot does not refit all actors.
    plot.renderer.camera.is_set=True
    plot.add_axes()
    return X

def render_all():
    import pyvista as pv
    import imageio.v2 as imageio
    from PIL import Image
    start=time.perf_counter();data=json.loads((O/'FAILURE_MAP.json').read_text());assert isinstance(data,list) and data
    c=json.loads((O/'CONFIGURATION_FROZEN.json').read_text());zg=np.load(O/'INTERFACE_GEOMETRY.npz');g={k:zg[k] for k in zg.files}
    for d in ['figures','animations','previews']: (R/d).mkdir(exist_ok=True)
    failure_figures(data,c);representatives,movies=select_representatives(data)
    dump(R/'REPRESENTATIVES.json',dict(modes=representatives,movies=movies,selection='One most-utilized NO_FAILURE transition-side case, resolved if observed, whole/hybrid if observed; otherwise fill actual distinct modes/strongest remaining cases',absent_modes=[m for m in MODES if m not in representatives]))
    fig,axs=plt.subplots(2,3,figsize=(15,9))
    for ax,mode,name in zip(axs.ravel(),MODES,SHORT):
        ax.axis('off');ax.set_title(name,fontsize=12)
        if mode not in representatives:ax.text(.5,.5,'NOT OBSERVED',ha='center',va='center',fontsize=16,color='#777777');continue
        s=load_scene(representatives[mode],g);plot=pv.Plotter(off_screen=True,window_size=(800,600));add_state(plot,s,-1);im=plot.screenshot();plot.close();Image.fromarray(im).save(R/'previews'/f'{mode}.png');ax.imshow(im)
        q=s['result'];ax.text(.5,-.02,f"{PNAMES[q['pattern_id']]} | {q['tau_budget']:.3g} Pa | {q['duration_over_Tcross']:.3g} T_cross\nFI_free={q['FI_free']:.3f}; debris={q['debris_fraction']:.3f}; resolved={q['resolved_fragment_count']}",ha='center',va='top',transform=ax.transAxes,fontsize=9)
    fig.suptitle('Observed failure modes | actual final positions and all surviving bonds | anchors grey',fontsize=14);fig.subplots_adjust(top=.91,bottom=.07,hspace=.3);savefig(fig,'representative_failure_modes')
    manifests=[]
    for caseid in movies:
        s=load_scene(caseid,g);z=s['trajectory'];q=s['result'];X=z['position_m']*1e3;lo=X.min((0,1));hi=X.max((0,1));pad=np.maximum(.12*(hi-lo),c['dx_m']*1e3);bounds=np.c_[lo-pad,hi+pad].ravel();RR=c['radius_m']*1e3
        local=np.array([-1.6*RR,1.6*RR,-1.6*RR,1.6*RR,-.4*RR,2.*RR]);mp4=R/'animations'/f'{caseid}.mp4';gifs=[];manifest=[]
        writer=imageio.get_writer(mp4,fps=12,codec='libx264',quality=7,macro_block_size=1,ffmpeg_params=['-threads','2'])
        for k,step in enumerate(z['step']):
            plot=pv.Plotter(shape=(1,2),off_screen=True,window_size=(1400,800));row=s['metrics'][k];camera_scales=[]
            for j in [0,1]:
                plot.subplot(0,j);add_state(plot,s,k,local if j==0 else bounds,overview=j==1)
                title=f"{PNAMES[q['pattern_id']]} | {q['tau_budget']:.3g} Pa | {q['duration_over_Tcross']:.3g} T_cross" if j==0 else 'Complete cloud | true connected-component colors'
                heading=plot.add_text(title+f"\nt={row['time_s']*1e3:.3f} ms; FI_free={row['FI_free']:.3f}; debris={row['debris_fraction']:.3f}",position='upper_left',font_size=10,color='#263747')
                footer=plot.add_text(('Fixed original-clot neighborhood' if j==0 else 'All clot particles and fixed anchors retained')+'\nPositions unscaled; all surviving bonds; exploratory only',position='lower_left',font_size=10,color='#263747')
                for label in [heading,footer]:
                    label.GetTextProperty().SetBackgroundColor(1,1,1);label.GetTextProperty().SetBackgroundOpacity(1)
                camera_scales.append(float(plot.renderer.camera.parallel_scale))
            im=plot.screenshot();plot.close();writer.append_data(im)
            if k%3==0 or k==len(z['step'])-1:gifs.append(Image.fromarray(im).resize((875,500)))
            if k in [0,len(z['step'])//2,len(z['step'])-1]:Image.fromarray(im).save(R/'previews'/f'{caseid}_frame_{k:04d}.png')
            manifest.append(dict(frame=k,LAMMPS_step=int(step),simulation_time_s=float(step*c['dt_s']),synthetic_state=False,camera_parallel_scale_mm=camera_scales))
            if k%60==0:print('RENDER',caseid,k,'/',len(z['step']),flush=True)
        writer.close();gifs[0].save(R/'animations'/f'{caseid}.gif',save_all=True,append_images=gifs[1:],duration=250,loop=0,optimize=False)
        info=dict(case_id=caseid,final_failure_mode=q['failure_mode'],fps=12,frames=manifest,all_surviving_bonds_drawn=True,all_debris_retained=True,positions_unscaled=True,glyph_width_dx=.9,overview_point_markers_enlarged=True,playback_not_physical_fracture_speed=True);dump(R/'animations'/f'{caseid}_MANIFEST.json',info);manifests.append(dict(case_id=caseid,MP4=str(mp4),frames=len(manifest),GIF_frames=len(gifs)))
    dump(R/'FIGURE_DATA.json',dict(source_card=json.loads((O/'COHERENT_FAILURE_CARD_LIU2024_V1.json').read_text()),elastic=json.loads((O/'ELASTIC_MAPPING_STATUS.json').read_text()),calibrations={x:json.loads((O/(x+'_CALIBRATION.json')).read_text()) for x in ['BULK','INTERFACE']},map_points=[{k:v for k,v in x.items() if not isinstance(v,(dict,list))} for x in data],representatives=representatives))
    dump(O/'VISUALIZATION_QA.json',dict(status='AWAITING_VISUAL_INSPECTION',figures=4,movies=manifests,render_wall_s=time.perf_counter()-start))
    print('RENDER_COMPLETE',time.perf_counter()-start,flush=True)
