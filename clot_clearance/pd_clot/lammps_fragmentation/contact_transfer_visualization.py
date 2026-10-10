"""Three bounded scientific figures and one movie from actual saved states."""
import os,csv,json,time
os.environ.setdefault('LIBGL_ALWAYS_SOFTWARE','1')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from .targeted_contact_transfer import O,R,geometry
from .config import dump

MAIN=['GAUSSIAN_100_MODIFIED','OPPOSING_100_D2','GAUSSIAN_200_D2']
LABELS={'GAUSSIAN_100_MODIFIED':'Gaussian 100 Pa / 1.0dx','OPPOSING_100_D2':'Opposing 100 Pa / 1.0dx','GAUSSIAN_200_D2':'Gaussian 200 Pa / 1.0dx','GAUSSIAN_100_NATIVE':'Gaussian 100 Pa / native 1.35dx','OPPOSING_100_NATIVE':'Opposing 100 Pa / native 1.35dx'}
COLORS={'GAUSSIAN_100_MODIFIED':'#2878b5','OPPOSING_100_D2':'#2b9169','GAUSSIAN_200_D2':'#b64577','GAUSSIAN_100_NATIVE':'#2878b5','OPPOSING_100_NATIVE':'#2b9169'}
MODES=['NO_FAILURE','RESOLVED_FRAGMENTATION','PULVERIZATION','WHOLE_CLOT_DETACHMENT','HYBRID_DEBONDING_AND_FRAGMENTATION','CLASSIFICATION_UNCERTAIN']
MCOLORS=['#78b47b','#e6ae41','#c84b51','#517db1','#9365a5','#cccccc']

def table(name):return np.genfromtxt(O/'analysis'/name/'timeseries.csv',delimiter=',',names=True,dtype=None,encoding='utf8')
def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(R/'figures'/f'{name}.{ext}',dpi=160)
    plt.close(fig)

def figures(summaries):
    for d in ['figures','previews','animations']:(R/d).mkdir(exist_ok=True)
    names=MAIN+['GAUSSIAN_100_NATIVE','OPPOSING_100_NATIVE'];data={n:table(n) for n in names}
    fig,axes=plt.subplots(3,3,figsize=(16,12))
    fields=[('bulk_broken_fraction','Broken bulk fraction',1),('interface_broken_fraction','Broken interface fraction',1),('debris_fraction','1-2 node unresolved debris fraction',1),('largest_attached_component_fraction','Largest attached clot fraction',1),('clot_KE_J','Clot kinetic energy (nJ)',1e9),('active_contact_pairs','Active nonbonded contact pairs',1),('min_nonbonded_distance_over_dx','Minimum nonbonded distance / dx',1),('deep_overlap_clot_fraction','Clot nodes at distance <0.5dx (%)',100)]
    for ax,(field,title,factor) in zip(axes.ravel(),fields):
        for name in names:
            a=data[name];ax.plot(a['time_s']*1e3,a[field]*factor,color=COLORS[name],ls='--' if 'NATIVE' in name else '-',lw=1.5,alpha=.75 if 'NATIVE' in name else 1,label=LABELS[name])
        ax.set(title=title,xlabel='Physical time (ms)',xlim=(0,2.33));ax.grid(alpha=.18)
        if field=='min_nonbonded_distance_over_dx':ax.axhline(.5,color='#444444',lw=.8,ls=':')
        if field=='deep_overlap_clot_fraction':ax.axhline(1,color='#444444',lw=.8,ls=':')
    ax=axes.ravel()[-1]
    for y,name in enumerate(names):
        a=data[name];tt=a['time_s']*1e3
        for k in range(len(tt)-1):ax.barh(y,tt[k+1]-tt[k],left=tt[k],height=.68,color=MCOLORS[MODES.index(a['failure_mode'][k])],edgecolor='none')
    ax.set(yticks=range(len(names)),yticklabels=['G100 new','O100 new','G200 new','G100 native','O100 native'],xlabel='Physical time (ms)',xlim=(0,2.33),title='Inherited >=3-node category timeline');ax.invert_yaxis()
    fig.suptitle('Contact transfer and stronger load | EXPLORATORY ONLY',fontsize=16,y=.99)
    handles=[plt.Line2D([],[],color=COLORS[n],ls='--' if 'NATIVE' in n else '-',label=LABELS[n]) for n in names];fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.96),ncol=3,fontsize=9)
    mode_handles=[plt.Line2D([],[],color=color,lw=6,label=mode.replace('_',' ').title()) for mode,color in zip(MODES,MCOLORS) if any(mode in data[n]['failure_mode'] for n in names)]
    fig.legend(handles=mode_handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=4,fontsize=9)
    fig.text(.5,.02,'Native 200 Pa was not run and is not inferred. Contact diagnostics are sampled every 10 us; gaps mean no nearby nonbonded pair.\nAll loads are reference-geometry engineering budgets; actual local peaks differ by pattern. No biological quantitative validation.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.88,bottom=.14,left=.07,right=.97,hspace=.38,wspace=.35);save(fig,'contact_transfer_and_stronger_load')

    fig,axes=plt.subplots(4,3,figsize=(16,13),gridspec_kw={'height_ratios':[1,1,1.15,.9]})
    lifetime_plot_selection={}
    for j,name in enumerate(MAIN):
        a=data[name];q=summaries[name];tt=a['time_s']*1e3
        for d,color in [('200','#ca7d22'),('300','#a33a96')]:
            axes[0,j].step(tt,a['N_detached_ge'+d],where='post',color=color,label='Detached >= '+d+' um')
            axes[1,j].plot(tt,a['mass_fraction_detached_ge'+d]*100,color=color)
        axes[0,j].set_title(LABELS[name],fontsize=11);axes[0,j].set_ylabel('Simultaneous detached count');axes[1,j].set_ylabel('Detached mass / initial clot (%)')
        for ax in axes[:2,j]:ax.set(xlim=(0,2.33),xlabel='Physical time (ms)');ax.grid(alpha=.2)
        with (O/'analysis'/name/'lifetimes.csv').open() as f:life=[r for r in csv.DictReader(f) if int(r['node_count'])>=34 and int(r['detached_frame_count'])>0]
        chosen=sorted(life,key=lambda r:(-int(r['node_count']),float(r['first_detached_time_s'])))[:8];chosen=sorted(chosen,key=lambda r:float(r['first_detached_time_s']));lifetime_plot_selection[name]=[int(r['track_id']) for r in chosen]
        ax=axes[2,j]
        for y,row in enumerate(chosen):
            first=float(row['first_detached_time_s'])*1e3;last=float(row['last_detached_seen_s'])*1e3;end=float(row['end_exclusive_or_censor_time_s'])*1e3;color='#a33a96' if int(row['node_count'])>=114 else '#ca7d22'
            ax.plot([first,last],[y,y],lw=4,color=color,solid_capstyle='butt');ax.scatter([first,last],[y,y],s=13,color=color);ax.plot([last,end],[y,y],color=color,ls=':',lw=1)
        ax.set(yticks=range(len(chosen)),yticklabels=[f"T{r['track_id']}: {r['node_count']} nodes" for r in chosen],xlabel='Physical time (ms)',title='Exact-membership episodes (largest 8)',xlim=(0,2.33));ax.grid(axis='x',alpha=.2)
        if not chosen:ax.text(.5,.5,'No detached component >=200 um',ha='center',transform=ax.transAxes)
        axes[3,j].axis('off');lines=[]
        for d in ['200','300']:
            z=q['thresholds'][d];lines.append(f">={d} um: ever {'YES' if z['EVER_OBSERVED'] else 'NO'} | peak N={z['MAX_SIMULTANEOUS_COUNT']} | final N={z['AT_FINAL_FRAME']}\nPeak mass {100*z['PEAK_DETACHED_MASS_FRACTION']:.2f}% | final {100*z['FINAL_MASS_FRACTION']:.2f}%\nPersistent membership episodes (>=3 frames): {z['persistent_membership_episodes']}")
        axes[3,j].text(0,1,'\n\n'.join(lines),va='top',fontsize=9)
    axes[0,0].annotate('49 nodes = 227.0 um\n0.50-0.54 ms, 5 frames',xy=(.52,1),xytext=(.8,1.6),arrowprops=dict(arrowstyle='->',color='#555555'),fontsize=9);axes[0,0].set_ylim(-.05,2.2)
    handles=[plt.Line2D([],[],color=c,lw=2,label='Detached >= '+d+' um') for d,c in [('200','#ca7d22'),('300','#a33a96')]];fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.957),ncol=2)
    fig.suptitle('Detached fragment size and lifetime | reference-volume equivalent diameter',fontsize=16,y=.99)
    fig.text(.5,.031,'200 um = literature reporting scale; 300 um = dx-based engineering display threshold\nAll attached components excluded here. Solid spans show first-to-last observation; dotted tails end at the next saved frame.\nEvery membership-changing split creates child episodes; episode counts are not independent embolus counts. All episodes remain in CSV.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.90,bottom=.105,left=.08,right=.97,hspace=.53,wspace=.35);save(fig,'detached_fragment_size_and_lifetime')
    with (O/'TARGETED_OLD_NEW_MAP_COMPARISON.csv').open() as f:points=list(csv.DictReader(f))
    fig,axs=plt.subplots(1,2,figsize=(14,5),sharex=True,sharey=True)
    for ax,model in zip(axs,['NATIVE_1P35DX','MODIFIED_1P0DX']):
        for p in points:
            if p['contact_model_id']!=model:continue
            y=0 if p['pattern_id']=='LOCAL_GAUSSIAN_SHEAR' else 1;x=float(p['tau_budget_Pa']);mode=p['terminal_mode'];ax.scatter(x,y,c=MCOLORS[MODES.index(mode)],s=180,edgecolor='#333333',linewidth=.6,zorder=3);ax.annotate(f'{x:.3g} Pa\n'+mode.replace('RESOLVED_FRAGMENTATION','Resolved').replace('PULVERIZATION','Pulverized').replace('NO_FAILURE','No failure'),(x,y),xytext=(0,14),textcoords='offset points',ha='center',fontsize=8)
        ax.set(xscale='log',xlim=(7,1600),ylim=(-.45,1.55),yticks=[0,1],yticklabels=['Gaussian','Opposing'],xlabel='Reference mean traction budget (Pa)',title=model+' | actual 2T_cross points');ax.grid(axis='x',alpha=.15);ax.set_xticks([10,31.6227766,100,200,1000],['10','31.6','100','200','1000']);ax.tick_params(axis='y',labelleft=True);ax.text(.5,.04,'Blank / unmeasured regions = UNKNOWN',ha='center',transform=ax.transAxes,fontsize=9,color='#666666')
    fig.suptitle('Targeted contact-model comparison | no boundary interpolation',fontsize=15)
    fig.text(.5,.025,'Separate contact models. Old 37-case map is unchanged. No new duration boundary, no fitted threshold and no assumed monotonicity.',ha='center',fontsize=9);fig.subplots_adjust(top=.80,bottom=.16,wspace=.22);save(fig,'targeted_contact_map_refresh')
    dump(R/'FIGURE_DATA.json',dict(summaries=summaries,map_points=points,lifetime_plot_selection=lifetime_plot_selection,all_lifetimes_retained_in_csv=True,figures=3))

def add_scene(plot,scene,index,bounds):
    import pyvista as pv
    g=scene['geometry'];z=scene['trajectory'];x=z['position_m'][index]*1e3;n=g['clot_count'];dx=scene['request']['dx_m']*1e3
    labels=scene['graph']['component_id'][index];attached=scene['graph']['attached'][index];_,inv,count=np.unique(labels,return_inverse=True,return_counts=True);sizes=count[inv]
    rgb=np.tile([46,151,106],(n,1));rgb[(~attached)&(sizes>=34)]=[223,135,29];rgb[(~attached)&(sizes>=114)]=[173,57,156];rgb[attached]=[41,112,176]
    cloud=pv.PolyData(x[:n]);cloud['rgb']=rgb.astype(np.uint8);cubes=cloud.glyph(geom=pv.Cube(x_length=.9*dx,y_length=.9*dx,z_length=.9*dx),orient=False,scale=False);plot.add_mesh(cubes,scalars='rgb',rgb=True,show_scalar_bar=False)
    wall=pv.PolyData(x[n:]).glyph(geom=pv.Cube(x_length=.9*dx,y_length=.9*dx,z_length=.9*dx),orient=False,scale=False);plot.add_mesh(wall,color='#777777')
    alive=np.ones(len(g['bonds']),bool);alive[scene['indices'][scene['events'][:,1]<=z['step'][index]]]=False;b=g['bonds'][alive];lines=pv.PolyData(x);lines.lines=np.c_[np.full(len(b),2),b].ravel();plot.add_mesh(lines,color='#4a535c',opacity=.20,line_width=.35)
    plot.set_background('white');plot.view_isometric(render=False);plot.enable_parallel_projection();plot.reset_camera(bounds=bounds,render=False);plot.renderer.camera.is_set=True;plot.add_axes()

def text_overlay(plot,text,where='upper_left',size=9):
    a=plot.add_text(text,position=where,font_size=size,color='#263747');a.GetTextProperty().SetBackgroundColor(1,1,1);a.GetTextProperty().SetBackgroundOpacity(1)

def render():
    import pyvista as pv
    import imageio.v2 as imageio
    from PIL import Image
    start=time.perf_counter();summaries=json.loads((O/'ALL_CASE_SUMMARIES.json').read_text());figures(summaries);figure_time=time.perf_counter()-start;print('THREE_FIGURES_COMPLETE',flush=True)
    scenes=[];g=geometry()
    for name in MAIN:
        p=Path(summaries[name]['source_path']);a=np.load(p/'accepted_trajectory.npz');e=np.load(p/'broken_bonds.npz');z=np.load(O/'analysis'/name/'graph_reconstruction.npz')
        scenes.append(dict(name=name,trajectory={k:a[k] for k in a.files},graph={k:z[k] for k in z.files},geometry=g,request=json.loads((p/'request.json').read_text()),events=e['events'],indices=e['initial_bond_indices'],metrics=table(name)))
    assert all(np.array_equal(s['trajectory']['step'],scenes[0]['trajectory']['step']) for s in scenes)
    lo=np.min([s['trajectory']['position_m'].min((0,1)) for s in scenes],axis=0)*1e3;hi=np.max([s['trajectory']['position_m'].max((0,1)) for s in scenes],axis=0)*1e3;pad=np.maximum(.08*(hi-lo),.05);bounds=np.c_[lo-pad,hi+pad].ravel();radius=.5625;local=np.array([-1.6*radius,1.6*radius,-1.6*radius,1.6*radius,-.4*radius,2*radius])
    ids=sorted(set(range(0,len(scenes[0]['trajectory']['step']),4))|{50,len(scenes[0]['trajectory']['step'])-1});writer=imageio.get_writer(R/'animations/contact_transfer_fragment_size.mp4',fps=10,codec='libx264',quality=7,macro_block_size=1,ffmpeg_params=['-threads','2']);gif=[];frames=[]
    for frame,k in enumerate(ids):
        plot=pv.Plotter(shape=(2,3),off_screen=True,window_size=(1500,1020));scales=[]
        for row in range(2):
            for col,s in enumerate(scenes):
                plot.subplot(row,col);add_scene(plot,s,k,local if row==0 else bounds);m=s['metrics'][k]
                text_overlay(plot,f"{LABELS[s['name']]}\nt={m['time_s']*1000:.3f} ms; detached >=200/300 um: {m['N_detached_ge200']}/{m['N_detached_ge300']}\nLargest detached d_eq={m['largest_detached_d_eq_um']:.1f} um; debris={m['debris_fraction']*100:.1f}%",size=9)
                text_overlay(plot,('Fixed original neighborhood' if row==0 else 'Complete cloud; all nodes / bonds retained')+'\nBlue: attached | grey: fixed anchors\nDetached: green <200 | orange >=200 | magenta >=300 um','lower_left',8)
                scales.append(float(plot.renderer.camera.parallel_scale))
        im=plot.screenshot();plot.close();writer.append_data(im);gif.append(Image.fromarray(im).resize((1000,680)))
        if k==50 or frame in [0,len(ids)//2,len(ids)-1]:Image.fromarray(im).save(R/'previews'/f'movie_frame_{frame:03d}.png')
        frames.append(dict(frame=frame,saved_frame_index=k,step=int(scenes[0]['trajectory']['step'][k]),time_s=float(scenes[0]['metrics']['time_s'][k]),camera_parallel_scale_mm=scales,synthetic=False))
        if frame%15==0:print('MOVIE',frame,'/',len(ids),flush=True)
    writer.close();gif[0].save(R/'animations/contact_transfer_fragment_size.gif',save_all=True,append_images=gif[1:],duration=100,loop=0,optimize=False)
    dump(R/'animations/FRAME_MANIFEST.json',dict(cases=MAIN,frames=frames,fps=10,interpolation=False,all_nodes_and_surviving_bonds=True,positions_unscaled=True,glyph_width_dx=.9,small_debris_not_enlarged=True,reference_neighborhood_bounds_mm=local.tolist(),complete_cloud_bounds_mm=bounds.tolist(),color_meaning='Actual attachment and reference-volume size, not damage or artificial cuts',playback_not_physical_experiment_time=True))
    dump(O/'VISUALIZATION_QA.json',dict(status='AWAITING_INSPECTION',figures=3,figure_wall_s=figure_time,movie_frames=len(frames),movie_render_wall_s=time.perf_counter()-start-figure_time,total_render_wall_s=time.perf_counter()-start))
    print('RENDER_COMPLETE',time.perf_counter()-start,flush=True)
