"""One synchronized comparison figure/movie, rendered from accepted native states."""
import os,json,time
os.environ.setdefault('LIBGL_ALWAYS_SOFTWARE','1')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .config import dump
from .contact_counterfactual import O,R,BASE,geometry
from .failure_map_visualization import add_state

NAMES=['CONTACT_ON','CONTACT_OFF','CONTACT_MODIFIED']
TITLES=['CONTACT_ON | native 1.35dx','CONTACT_OFF | DIAGNOSTIC_ONLY','CONTACT_MODIFIED | 1.0dx']
COLORS=['#c33f49','#4575b4','#17866d']

def scenes():
    g=geometry();result=[]
    for name in NAMES:
        p=BASE if name=='CONTACT_ON' else O/name
        z=np.load(p/'accepted_trajectory.npz');e=np.load(p/'broken_bonds.npz')
        result.append(dict(case_id=name,trajectory={k:z[k] for k in z.files},events=e['events'],indices=e['initial_bond_indices'],request=json.loads((p/'request.json').read_text()),result=json.loads((p/'result.json').read_text()),geometry=g,metrics=np.genfromtxt(p/'fragment_metrics.csv',delimiter=',',names=True,dtype=None,encoding='utf8')))
    assert all(np.array_equal(result[0]['trajectory']['step'],s['trajectory']['step']) for s in result)
    return result

def opaque_text(plot,text,where='upper_left',size=10):
    actor=plot.add_text(text,position=where,font_size=size,color='#263747')
    actor.GetTextProperty().SetBackgroundColor(1,1,1)
    actor.GetTextProperty().SetBackgroundOpacity(1)

def render():
    import pyvista as pv
    import imageio.v2 as imageio
    from PIL import Image
    start=time.perf_counter();ss=scenes();c=json.loads((O/'CAUSAL_COMPARISON.json').read_text())
    all_lo=np.min([s['trajectory']['position_m'].min((0,1)) for s in ss],axis=0)*1e3
    all_hi=np.max([s['trajectory']['position_m'].max((0,1)) for s in ss],axis=0)*1e3
    pad=np.maximum(.08*(all_hi-all_lo),.05);bounds=np.c_[all_lo-pad,all_hi+pad].ravel()
    radius=ss[0]['request']['radius_m']*1e3
    local=np.array([-1.6*radius,1.6*radius,-1.6*radius,1.6*radius,-.4*radius,2*radius])
    for d in ['figures','animations','previews']:(R/d).mkdir(exist_ok=True)
    fig=plt.figure(figsize=(15,13.8));gs=fig.add_gridspec(4,6,height_ratios=[.9,.9,1.35,1.35])
    variables=[('bulk_broken_fraction','Broken bulk bond fraction'),('interface_broken_fraction','Broken interface bond fraction'),('debris_fraction','Debris fraction (1-2 clot nodes)'),('resolved_fragment_count','Components with >=3 clot nodes')]
    first=c['cases']['CONTACT_ON']['first_any_break_time_s']*1e3
    onset=c['cases']['CONTACT_ON']['categories']['PULVERIZATION']['first_saved_time_s']*1e3
    end=ss[0]['metrics']['time_s'][-1]*1e3
    for k,(field,label) in enumerate(variables):
        ax=fig.add_subplot(gs[k//2,3*(k%2):3*(k%2)+3])
        for s,title,color in zip(ss,TITLES,COLORS):ax.plot(s['metrics']['time_s']*1e3,s['metrics'][field],label=title,color=color,lw=1.7)
        for value,style in [(first,':'),(onset,'--'),(end,'-.')]:ax.axvline(value,color='#555555',ls=style,lw=.8)
        ax.set(xlim=(0,end+.04),xlabel='Physical time (ms)',ylabel=label);ax.grid(alpha=.18)
        if field!='resolved_fragment_count':ax.set_ylim(-.02,1.03)
        if k==0:
            for v,label2 in [(first,'First bulk break'),(onset,'ON pulverization'),(end,'End')]:ax.text(v,.98,label2,rotation=90,ha='right',va='top',fontsize=8,transform=ax.get_xaxis_transform())
    shot_records=[]
    for row,target in enumerate([.0008,.00233]):
        for col,s in enumerate(ss):
            k=int(np.argmin(abs(s['metrics']['time_s']-target)));actual=float(s['metrics']['time_s'][k]);assert abs(actual-target)<1e-12
            plot=pv.Plotter(off_screen=True,window_size=(800,620));add_state(plot,s,k,bounds)
            im=plot.screenshot();scale=float(plot.renderer.camera.parallel_scale);plot.close()
            preview=R/'previews'/f'{s["case_id"]}_{int(round(target*1e6))}us.png';Image.fromarray(im).save(preview)
            ax=fig.add_subplot(gs[row+2,2*col:2*col+2]);ax.imshow(im);ax.axis('off');ax.set_title(f'{TITLES[col]}\nt = {actual*1e3:.2f} ms',fontsize=10,color=COLORS[col])
            m=s['metrics'][k];ax.text(.5,-.01,f"Debris {100*m['debris_fraction']:.2f}% | largest {100*m['largest_free_fraction']:.2f}%\n>=3-node components {m['resolved_fragment_count']}",ha='center',va='top',transform=ax.transAxes,fontsize=9)
            shot_records.append(dict(case=s['case_id'],step=int(s['trajectory']['step'][k]),time_s=actual,camera_parallel_scale_mm=scale,common_bounds_mm=bounds.tolist()))
    handles=[plt.Line2D([],[],color=color,lw=2,label=title) for title,color in zip(TITLES,COLORS)]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.947),ncol=3,fontsize=10)
    fig.suptitle('Post-break contact cascade: one frozen hemisphere case | EXPLORATORY ONLY',fontsize=16,y=.985)
    fig.text(.5,.958,'100 Pa local Gaussian shear; pulse 2 T_cross; dt 0.5 us; same 2.33 ms end; baseline reused',ha='center',fontsize=10)
    limits='; '.join(f'{d} [{bounds[2*i]:.2f}, {bounds[2*i+1]:.2f}]' for i,d in enumerate('xyz'))
    fig.text(.5,.027,'All six 3D views use the same fixed camera and bounds (mm): '+limits+'\nActual unscaled positions, all clot nodes and surviving bonds; fixed anchors grey. Component colors are graph labels.\nThe >=3-node count includes the attached main remnant; it does not mean similarly sized large fragments.',ha='center',fontsize=9)
    fig.subplots_adjust(top=.90,bottom=.085,left=.065,right=.98,hspace=.45,wspace=.50)
    for ext in ['png','pdf']:fig.savefig(R/'figures'/f'contact_cascade_causal_comparison.{ext}',dpi=160)
    plt.close(fig);print('MAIN_FIGURE_COMPLETE',flush=True)
    indices=sorted(set(range(0,len(ss[0]['trajectory']['step']),4))|{len(ss[0]['trajectory']['step'])-1})
    writer=imageio.get_writer(R/'animations/contact_cascade_synchronized.mp4',fps=10,codec='libx264',quality=7,macro_block_size=1,ffmpeg_params=['-threads','2'])
    gif=[];frames=[]
    for frame,k in enumerate(indices):
        plot=pv.Plotter(shape=(2,3),off_screen=True,window_size=(1500,960));scales=[]
        for row in range(2):
            for col,s in enumerate(ss):
                plot.subplot(row,col);add_state(plot,s,k,local if row==0 else bounds,overview=row==1)
                m=s['metrics'][k];opaque_text(plot,f"{TITLES[col]}\nt={m['time_s']*1e3:.3f} ms; debris={100*m['debris_fraction']:.1f}%",size=9)
                opaque_text(plot,('Original clot neighborhood (fixed view)' if row==0 else 'Complete cloud (fixed view); all nodes retained')+'\nAll surviving bonds | EXPLORATORY ONLY','lower_left',8)
                scales.append(float(plot.renderer.camera.parallel_scale))
        im=plot.screenshot();plot.close();writer.append_data(im);gif.append(Image.fromarray(im).resize((1000,640)))
        if frame in [0,len(indices)//2,len(indices)-1]:Image.fromarray(im).save(R/'previews'/f'movie_frame_{frame:03d}.png')
        frames.append(dict(frame=frame,saved_frame_index=k,LAMMPS_step=int(ss[0]['trajectory']['step'][k]),physical_time_s=float(ss[0]['metrics']['time_s'][k]),camera_parallel_scale_mm=scales,synthetic_state=False))
        if frame%15==0:print('MOVIE_FRAME',frame,'/',len(indices),flush=True)
    writer.close();gif[0].save(R/'animations/contact_cascade_synchronized.gif',save_all=True,append_images=gif[1:],duration=100,loop=0,optimize=False)
    dump(R/'animations/FRAME_MANIFEST.json',dict(cases=NAMES,source='accepted native LAMMPS saved states',frames=frames,fps=10,frame_stride=4,final_frame_included=True,interpolation=False,positions_unscaled=True,all_clot_nodes_retained=True,all_surviving_bonds_drawn=True,top_row_fixed_neighborhood_bounds_mm=local.tolist(),bottom_row_fixed_complete_cloud_bounds_mm=bounds.tolist(),overview_point_markers_enlarged=True,playback_is_not_physical_fracture_speed=True))
    dump(O/'VISUALIZATION_QA.json',dict(status='AWAITING_INSPECTION',main_figures=1,formats=['PNG','PDF'],movies=1,movie_formats=['MP4','GIF'],movie_frames=len(frames),snapshots=shot_records,render_wall_s=time.perf_counter()-start))
    print('VISUALIZATION_COMPLETE',time.perf_counter()-start,flush=True)
