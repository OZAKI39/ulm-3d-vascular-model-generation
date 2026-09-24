"""Shared static/animated views of P8 replay caches. SI data stays untouched."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.patches import Circle,Polygon
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from .particle8_replay import read,write,digest,snapshot,interpolate,canonical_hash,DEFAULT_OUTPUT
from .rbc_orientation import rotation_matrix

RBC='#d85b47'; MB='#168eb1'; INK='#183449'; MUTED='#8c99a3'; GREEN='#278276'; LIGHT='#eef3f6'
FIGURES=['particle8_00_overview','particle8_01_real_inlet_flux','particle8_02_rbc_volume_scheduler',
 'particle8_03_mb_scheduler','particle8_04_real_mixed_smoke','particle8_05_synthetic_lifecycle',
 'particle8_06_pending_identity','particle8_07_restart_parity','particle8_08_animation_storyboard','particle8_09_limitations']
ANIMATIONS=['particle8_anim_01_real_single_mb_entry','particle8_anim_02_real_mixed_inlet_smoke',
 'particle8_anim_03_synthetic_open_section_lifecycle','particle8_anim_04_restart_continuity','particle8_anim_05_flux_weighted_births']
plt.rcParams.update({'font.size':11,'font.family':'DejaVu Sans','axes.titleweight':'bold','axes.labelcolor':INK,
 'text.color':INK,'axes.edgecolor':'#c3cfd8','axes.spines.top':False,'axes.spines.right':False,
 'axes.grid':True,'grid.alpha':.14,'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white'})
FOOTER='H_D = 0.45 feed volume fraction; tube Hct is a separate observation  |  C_MB = 8.5e12 m^-3 nominal estimate'
BOUNDARY='MODEL ASSUMPTION: isotropic RBC orientation  |  RBC PASSAGE NOT ESTABLISHED  |  NOT PRODUCTION TIMESTEP SELECTION'


def canvas(title,subtitle,*,size=(14.4,9)):
    fig=plt.figure(figsize=size,dpi=100)
    fig.text(.045,.957,title,fontsize=20,weight='bold',va='top')
    fig.text(.045,.906,subtitle,fontsize=11,color=GREEN,va='top')
    fig.text(.045,.050,FOOTER,fontsize=9)
    fig.text(.045,.027,BOUNDARY,fontsize=9,color='#944538')
    return fig


def axes(fig,rect): return fig.add_axes(rect)


def triangle_area_2d(triangles):
    p=np.asarray(triangles);a=p[:,1]-p[:,0];b=p[:,2]-p[:,0]
    return abs(a[:,0]*b[:,1]-a[:,1]*b[:,0])/2


def inlet_basis(data):
    tri=np.array(data['triangles_m']); origin=tri.reshape(-1,3).mean(0)
    _,_,vh=np.linalg.svd(tri.reshape(-1,3)-origin,full_matrices=False)
    return origin,vh[:2],(tri-origin)@vh[:2].T*1e6


def flux_panel(fig,ax,xy,values,label):
    pc=PolyCollection(xy,array=np.asarray(values),cmap='viridis',edgecolors='white',linewidths=.18)
    ax.add_collection(pc);ax.autoscale();ax.set_aspect('equal');ax.set(xlabel='Inlet basis 1 (µm)',ylabel='Inlet basis 2 (µm)')
    fig.colorbar(pc,ax=ax,shrink=.8,pad=.03,label=label)
    return pc


def count_text(s):
    return '\n'.join(f'{sp}: scheduled {c["scheduled"]:,}   admitted {c["admitted"]:,}   pending {c["pending"]:,}\n'
                     f'        active {c["active"]:,}   exited {c["exited"]:,}   deleted {c["deleted"]:,}' for sp,c in s['counts'].items())


def actual_outline(record,state):
    """Exact xy projection boundary of the sampled sphere/ellipsoid, in metres."""
    g=record['geometry'];theta=np.linspace(0,2*np.pi,65);d=np.stack([np.cos(theta),np.sin(theta),np.zeros_like(theta)],axis=1)
    if record['species']=='MB': p=d*g['radius_m']
    else:
        rot=rotation_matrix(state['q']);a=np.asarray(g['axes_m']);cov=rot@np.diag(a*a)@rot.T
        p=(d@cov)/np.sqrt(np.einsum('ij,ij->i',d@cov,d))[:,None]
    return p[:,:2]+np.array(state['position_m'])[:2]


def wall3d(fig,rect,geo,origin,limit=7):
    ax=fig.add_axes(rect,projection='3d',computed_zorder=False)
    wall=np.asarray(geo['nearby_wall_triangles_m']);mid=wall.mean(1)
    mask=np.linalg.norm(mid-origin,axis=1)<limit*1e-6*1.5
    chosen=wall[mask][::2]  # Visible surface triangle subset only; no geometry changes.
    ax.add_collection3d(Poly3DCollection((chosen-origin)*1e6,facecolor='#adbec9',edgecolor='#8ea4b0',linewidth=.12,alpha=.10))
    ax.add_collection3d(Poly3DCollection((np.asarray(geo['inlet_triangles_m'])-origin)*1e6,
                                       facecolor='#66b6c0',edgecolor='none',alpha=.24))
    ax.set(xlim=(-limit,limit),ylim=(-limit,limit),zlim=(-limit,limit),xlabel='Δx (µm)',ylabel='Δy (µm)',zlabel='Δz (µm)')
    ax.set_box_aspect([1,1,1]);ax.view_init(elev=22,azim=-64)
    return ax


def surface_points(shape,center):
    """Render actual sphere or P3 capsule; retained RBC q is not a capsule axis."""
    u=np.linspace(0,2*np.pi,24);v=np.linspace(-np.pi/2,np.pi/2,16)
    U,V=np.meshgrid(u,v)
    if shape['mode']=='SPHERE_MB':
        r=shape['radius_m'];p=np.stack([r*np.cos(V)*np.cos(U),r*np.cos(V)*np.sin(U),r*np.sin(V)],axis=-1)
    else:
        r=shape['capsule_radius'];length=shape['capsule_length'];axis=np.asarray(shape['capsule_axis'])
        trial=np.eye(3)[np.argmin(abs(axis))];a=np.cross(axis,trial);a/=np.linalg.norm(a);b=np.cross(axis,a)
        axial=r*np.sin(V)+np.sign(V)*length/2
        p=r*np.cos(V)[...,None]*(np.cos(U)[...,None]*a+np.sin(U)[...,None]*b)+axial[...,None]*axis
    return p+np.asarray(center)


class ReplayView:
    def __init__(self,output,index):
        self.root=Path(output);self.data=self.root/'data';self.index=index
        self.dynamic=[];self.shape_artists=[]
        if index==1:self.init_mb()
        elif index==2:self.init_smoke()
        elif index==3:self.init_synthetic()
        elif index==4:self.init_restart()
        else:self.init_flux()

    def clear_dynamic(self):
        for a in self.dynamic:a.remove()
        self.dynamic=[]

    def init_mb(self):
        self.scene=read(self.data/'real_single_mb_scene.json');self.record=self.scene['records'][0]
        self.geo=read(self.data/'13_real_geometry.json');self.origin=np.array(self.record['trajectory'][0]['position_m'])
        self.fig=canvas('01  |  A microbubble enters the real Frozen inlet',
                        'REAL GEOMETRY / SINGLE MB / NO RBC IN THIS FIXTURE  •  Saved P6.5 states; display interpolation only')
        self.a=wall3d(self.fig,[.02,.23,.46,.55],self.geo,self.origin,3.5);self.a.set_title('Real inlet / nearby wall (local view)',fontsize=12)
        self.b=axes(self.fig,[.57,.27,.36,.50]);self.b.set_aspect('equal')
        end=np.array(self.record['trajectory'][-1]['position_m']);self.forward=(end-self.origin)/np.linalg.norm(end-self.origin)
        trial=np.eye(3)[np.argmin(abs(self.forward))];self.side=np.cross(self.forward,trial);self.side/=np.linalg.norm(self.side)
        wall=(np.array(self.geo['nearby_wall_triangles_m'])-self.origin)*1e6
        projected=np.stack([wall@self.forward,wall@self.side],axis=-1)
        self.b.add_collection(PolyCollection(projected,facecolors='#c3d1db',edgecolors='#9cafbb',alpha=.1,linewidths=.2))
        inlet=(np.array(self.geo['inlet_triangles_m'])-self.origin)*1e6
        self.b.add_collection(PolyCollection(np.stack([inlet@self.forward,inlet@self.side],axis=-1),facecolors=MB,alpha=.20,edgecolors='none'))
        self.b.set(xlim=(-2,3),ylim=(-2.3,2.3),xlabel='Along entry displacement (µm)',ylabel='Transverse projection (µm)',title='Equal physical scale / original radius')
        path=(np.array([s['position_m'] for s in self.record['trajectory']])-self.origin)*1e6
        self.path=path;self.b.plot(path@self.forward,path@self.side,':',color=MUTED,lw=1,label='Saved center path')
        self.b.legend(loc='upper right',fontsize=9)
        self.clock=self.fig.text(.045,.835,'',fontsize=13,weight='bold')
        self.info=self.fig.text(.52,.16,'',fontsize=12)

    def draw_mb(self,t,label):
        self.clear_dynamic();s=snapshot(self.scene,t);birth=self.record['admit_time_s']
        self.clock.set_text(f'Physical t = {t:.9f} s   |   since birth = {(t-birth)*1000:+.3f} ms')
        if s['active']:
            state=s['active'][0];p=(np.array(state['position_m'])-self.origin)*1e6
            points=(surface_points(self.record['admitted_shape'],state['position_m'])-self.origin)*1e6
            self.dynamic.append(self.a.plot_surface(points[:,:,0],points[:,:,1],points[:,:,2],color=MB,alpha=.7,linewidth=0))
            self.dynamic+=self.a.plot([0,p[0]],[0,p[1]],[0,p[2]],color=MB,lw=3)
            circle=Circle((p@self.forward,p@self.side),self.record['geometry']['radius_m']*1e6,facecolor=MB,edgecolor=INK,alpha=.45)
            self.b.add_patch(circle);self.dynamic.append(circle)
            self.dynamic+=self.b.plot([0,p@self.forward],[0,p@self.side],color=MB,lw=3)
            self.dynamic.append(self.b.text(p@self.forward,p@self.side,'ID 1',ha='center',va='center',fontsize=10))
            self.info.set_text(f'ACTIVE  •  displacement {np.linalg.norm(p)*1000:.2f} nm\nOriginal diameter {self.record["geometry"]["diameter_m"]*1e6:.4f} µm; no inlet offset')
        else:self.info.set_text('Before deterministic birth\nNo particle placed in the lumen yet')
        return s

    def init_smoke(self):
        self.scene=read(self.data/'real_mixed_scene.json');self.geo=read(self.data/'13_real_geometry.json')
        self.origin=np.array(self.geo['inlet_triangles_m']).reshape(-1,3).mean(0)
        self.fig=canvas('02  |  Real inlet: admission and an accumulating queue',
             'REAL GEOMETRY / SMOKE / RBC PASSAGE NOT ESTABLISHED  •  Admitted obstacle is STATIC; this is not transport')
        self.a=wall3d(self.fig,[.01,.25,.47,.53],self.geo,self.origin,7);self.a.set_title('Actual P3 capsule surrogate after admission',fontsize=11)
        self.b=axes(self.fig,[.57,.38,.36,.35]);self.b.set(ylim=(0,1200),ylabel='RBC count',title='All pending RBCs included in the counter')
        self.bars=self.b.bar(['scheduled','admitted','pending','active'],[0]*4,color=[RBC,RBC,MUTED,RBC],width=.6)
        self.labels=[self.b.text(k,15,'0',ha='center',fontsize=12) for k in range(4)]
        self.clock=self.fig.text(.045,.837,'',fontsize=13,weight='bold');self.info=self.fig.text(.55,.29,'',fontsize=12)
        self.queue=self.fig.text(.55,.20,'',fontsize=10,color=MUTED)
        self.fig.text(.06,.16,'Waiting is a queue, not a spatial particle cloud.\nThe original RBC geometry and orientation are retained.',fontsize=11)

    def draw_smoke(self,t,label):
        s=snapshot(self.scene,t);self.clock.set_text(f'Physical scheduler clock t = {t:.6f} s   |   active shapes held fixed')
        r=s['counts']['RBC'];m=s['counts']['MB']
        for k,bar,txt in zip(['scheduled','admitted','pending','active'],self.bars,self.labels):
            bar.set_height(r[k]);txt.set_y(r[k]+18);txt.set_text(str(r[k]))
        self.info.set_text(f'MB: scheduled {m["scheduled"]}  |  admitted {m["admitted"]}  |  pending {m["pending"]}\nRBC exits: {r["exited"]}  •  Continuous admission remains limited')
        pending=[p['particle_id'] for p in s['pending'] if p['species']=='RBC']
        self.queue.set_text('Non-spatial queue (first 6 IDs): '+', '.join(map(str,pending[:6]))+(' …' if len(pending)>6 else '')+'\nPending shapes are not drawn inside the lumen.')
        if s['active'] and not self.shape_artists:
            r=self.scene['records'][0];p=(surface_points(r['admitted_shape'],s['active'][0]['position_m'])-self.origin)*1e6
            self.shape_artists.append(self.a.plot_surface(p[:,:,0],p[:,:,1],p[:,:,2],color=RBC,alpha=.85,linewidth=0))
        if not s['active'] and self.shape_artists:
            for a in self.shape_artists:a.remove()
            self.shape_artists=[]
        return s

    def init_synthetic(self):
        self.scene=read(self.data/'synthetic_scene.json');self.records={r['particle_id']:r for r in self.scene['records']}
        self.fig=canvas('03  |  SYNTHETIC CONTROL SECTION',
            'BOOKKEEPING / LIFECYCLE DEMONSTRATION  •  NOT REAL FROZEN LUMEN  •  Original P7 common-plug fixture')
        self.a=axes(self.fig,[.075,.28,.45,.46]);self.b=axes(self.fig,[.66,.28,.29,.46]);self.b.set_aspect('equal')
        self.a.set(xlim=(-.5,10.5),ylim=(-105,105),xlabel='Center axial position z (nm) — stretched display axis',ylabel='Center x (µm)',title='INLET → continuous center motion → OUTLET / DELETE')
        self.a.axvline(0,color=GREEN,lw=2);self.a.axvline(10,color=INK,lw=2)
        self.a.axhline(-100,color=INK,lw=1);self.a.axhline(100,color=INK,lw=1)
        self.b.set(xlim=(-105,105),ylim=(-105,105),xlabel='x (µm)',ylabel='y (µm)',title='Actual xy shape projections')
        self.b.add_patch(plt.Rectangle((-100,-100),200,200,fill=False,edgecolor=INK,lw=1))
        self.clock=self.fig.text(.045,.835,'',fontsize=12,weight='bold');self.counter=self.fig.text(.075,.145,'',fontsize=11)
        self.fig.text(.075,.785,'Width 200 µm; length 10 nm; speed 25 µm/s; residence 0.4 ms. Shapes straddle open caps.',fontsize=11)
        self.fig.text(.66,.175,'Warm = RBC; cool = MB\nLeft: ID glyphs, not physical shapes\nRight: original size and orientation',fontsize=10)

    def draw_tracks(self,axis,s,records,*,footprints=None):
        for p in s['active']:
            r=records[p['particle_id']];pos=np.array(p['position_m']);color=RBC if p['species']=='RBC' else MB
            x=pos[2]*1e9;y=pos[0]*1e6
            self.dynamic+=axis.plot([0,x],[y,y],color=color,alpha=.4,lw=2)
            self.dynamic+=axis.plot([x],[y],marker='o' if p['species']=='MB' else 'D',color=color,ms=12 if p['species']=='MB' else 7)
            self.dynamic.append(axis.annotate(str(p['particle_id']),(x,y),xytext=(3,6),textcoords='offset points',fontsize=8,color=color))
            if footprints is not None:
                patch=Polygon(actual_outline(r,p)*1e6,facecolor=color,edgecolor=color,alpha=.7)
                footprints.add_patch(patch);self.dynamic.append(patch)
                self.dynamic.append(footprints.annotate(str(p['particle_id']),pos[:2]*1e6,xytext=(4,4),textcoords='offset points',fontsize=7,color=color))

    def draw_synthetic(self,t,label):
        self.clear_dynamic();s=snapshot(self.scene,t)
        self.clock.set_text(f'Physical t = {t:.9f} s  |  {label}')
        self.draw_tracks(self.a,s,self.records,footprints=self.b);self.counter.set_text(count_text(s))
        return s

    def init_restart(self):
        self.scene=read(self.data/'synthetic_scene.json');self.other=read(self.data/'restarted_scene.json')
        self.records={r['particle_id']:r for r in self.scene['records']}
        self.fig=canvas('04  |  Checkpoint → destroy → binary restart',
            'SYNTHETIC CONTROL / AUDIT  •  Two independent saved P7 runs  •  1,180 events at checkpoint; 2,358 at final time')
        self.a=axes(self.fig,[.075,.33,.38,.40]);self.b=axes(self.fig,[.565,.33,.38,.40])
        for ax,title in [(self.a,'Continuous run'),(self.b,'Destroyed / restarted run')]:
            ax.set(xlim=(-.5,10.5),ylim=(-105,105),xlabel='Center z (nm) — stretched axis',ylabel='Center x (µm)',title=title)
            ax.axvline(0,color=GREEN);ax.axvline(10,color=INK)
        self.clock=self.fig.text(.045,.835,'',fontsize=13,weight='bold');self.counter=self.fig.text(.075,.17,'',fontsize=11)
        self.match=self.fig.text(.56,.19,'',fontsize=12,color=GREEN)
        self.fig.text(.075,.775,'ID glyphs are not physical shapes. Boundary surfaces and original geometry are unchanged.',fontsize=10)

    def draw_restart(self,t,label):
        self.clear_dynamic();s=snapshot(self.scene,t);other=snapshot(self.other,t);assert s==other
        self.draw_tracks(self.a,s,self.records);self.draw_tracks(self.b,other,self.records)
        phase='CUT TO FINAL AUDIT' if t>.126 else s['restart_segment'].upper()
        self.clock.set_text(f'Physical t = {t:.9f} s  |  {phase}  |  checkpoint t = 0.125 s')
        self.counter.set_text(count_text(s));self.match.set_text('MISMATCHES: 0\nSame IDs / times / positions / queues\nSame complete RNG and scheduler state\nNo duplicate ID / no lost particle')
        return dict(**s,continuous_state_sha256=canonical_hash(s),restarted_state_sha256=canonical_hash(other),mismatch_count=0)

    def init_flux(self):
        mesh=read(self.data/'01_inlet_mesh.json');self.origin,self.basis,self.xy=inlet_basis(mesh)
        self.samples=np.genfromtxt(self.data/'02_actual_inlet_samples.csv',delimiter=',',names=True)
        self.points=(np.column_stack([self.samples[k] for k in ['x_m','y_m','z_m']])-self.origin)@self.basis.T*1e6
        self.q=np.mean(mesh['q_m_s'],axis=1)*1e6
        self.fig=canvas('05  |  Why fast inlet regions receive more samples',
            'REAL INLET GEOMETRY / ANIMATED SAMPLING AUDIT  •  P7 diagnostic proposals, NOT admitted physical births')
        self.a=axes(self.fig,[.075,.25,.36,.51]);self.b=axes(self.fig,[.57,.25,.36,.51])
        flux_panel(self.fig,self.a,self.xy,self.q,'Inward speed (µm/s)')
        area=triangle_area_2d(self.xy)
        self.area=area;self.density=flux_panel(self.fig,self.b,self.xy,np.zeros(len(area)),'Samples per µm²')
        final=np.bincount(self.samples['triangle_id'].astype(int),minlength=len(area))/area
        self.density.set_clim(0,final.max());self.a.set_title('Saved linear inlet flux');self.b.set_title('Accumulated triangle density / actual area')
        self.clock=self.fig.text(.045,.835,'',fontsize=13,weight='bold');self.info=self.fig.text(.075,.15,'',fontsize=11)

    def draw_flux(self,t,label):
        n=min(len(self.points),max(0,int(round(t/10*len(self.points)))))
        count=np.bincount(self.samples['triangle_id'][:n].astype(int),minlength=len(self.area));self.density.set_array(count/self.area)
        self.clock.set_text(f'Display clock = {t:.2f} s  |  diagnostic samples = {n:,} / {len(self.points):,}')
        self.info.set_text('The display clock is not a simulation clock. Sample proposals do not have physical birth times.\nDensity is area-normalized; this animation does not imply finite-size admission or solve RBC passage.')
        return dict(time_s=t,time_role='DISPLAY_CLOCK_NOT_PHYSICAL_BIRTH_TIME',sample_count=n,
                    source_classification='real_frozen',diagnostic_proposals_not_admitted_births=True)

    def draw(self,t,label=''):
        return [self.draw_mb,self.draw_smoke,self.draw_synthetic,self.draw_restart,self.draw_flux][self.index-1](float(t),label)


def frame_schedule(index,data):
    if index==1:
        scene=read(data/'real_single_mb_scene.json');lo,hi=scene['time_range_s']
        return [(float(t),'Saved trajectory; linear display interpolation') for t in np.linspace(lo,hi,150)]
    if index==2:
        scene=read(data/'real_mixed_scene.json');end=scene['time_range_s'][1]
        return [(float(t),'Static obstacle smoke') for t in np.r_[np.linspace(0,end,145),np.repeat(end,15)]]
    if index==3:
        # Slow physical windows separated by explicit cuts. NEVER interpolate across cuts.
        windows=[(0,.002,'WINDOW 1 / startup'),(.1170,.1184,'CUT TO WINDOW 2 / first MB'),(.2347,.2361,'CUT TO WINDOW 3 / second MB')]
        return [(float(t),label) for a,b,label in windows for t in np.linspace(a,b,70)]+[(.25,'CUT TO FINAL AUDIT / t = 0.25 s')]*15
    if index==4:
        return [(float(t),'Restart audit') for t in np.r_[np.linspace(.124,.125,70,endpoint=False),np.repeat(.125,15),np.linspace(.125,.126,70)]]+[(.25,'CUT TO FINAL AUDIT')]*15
    return [(float(t),'Diagnostic sampling') for t in np.linspace(0,10,150)]


def export_animations(output=DEFAULT_OUTPUT,indices=None):
    import imageio_ffmpeg
    from PIL import Image
    output=Path(output);folder=output/'animations';frames=output/'frames'
    folder.mkdir(exist_ok=True);frames.mkdir(exist_ok=True)
    for idx in (range(1,6) if indices is None else indices):
        view=ReplayView(output,idx);schedule=frame_schedule(idx,output/'data');name=ANIMATIONS[idx-1]
        path=folder/(name+'.mp4');encoder=imageio_ffmpeg.write_frames(str(path),(1440,900),fps=15,
            codec='libx264',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',quality=8,macro_block_size=2,
            output_params=['-movflags','+faststart'],ffmpeg_log_level='error')
        encoder.send(None);manifest=[];gif=[]
        key_indices=sorted(set([0,len(schedule)//3,len(schedule)//2,2*len(schedule)//3,len(schedule)-1]))
        try:
            for k,(t,label) in enumerate(schedule):
                state=view.draw(t,label);view.fig.canvas.draw()
                rgb=np.ascontiguousarray(np.asarray(view.fig.canvas.buffer_rgba())[:,:,:3]);encoder.send(rgb)
                receipt=dict(frame_index=k,video_time_s=k/15,physical_or_display_time_s=t,window_label=label,
                    clock_role='DISPLAY_CLOCK_NOT_PHYSICAL_TIME' if idx==5 else 'PHYSICAL_SIMULATION_CLOCK',
                    state_sha256=canonical_hash(state),state=state)
                manifest.append(receipt)
                if k in key_indices:
                    Image.fromarray(rgb).save(frames/f'{name}_frame_{k:04d}.png')
                if k%5==0:
                    gif.append(Image.fromarray(rgb).resize((864,540)))
            encoder.close()
        finally:
            encoder.close();plt.close(view.fig)
        gif[0].save(folder/(name+'.gif'),save_all=True,append_images=gif[1:],duration=round(1000*5/15),loop=0,optimize=False)
        write(output/'data'/(name+'_frames.json'),dict(animation=path.name,fps=15,size_px=[1440,900],
            frame_count=len(manifest),duration_s=len(manifest)/15,keyframe_indices=key_indices,
            gif_role='DOWNSAMPLED_PREVIEW_3FPS',time_schedule_role='EXPLICIT_TIME_WINDOWS_NO_INTERPOLATION_ACROSS_CUTS',
            frames=manifest))
        print('Exported',path.name,len(manifest),'frames',flush=True)


def export_figures(output=DEFAULT_OUTPUT):
    output=Path(output);data=output/'data';out=output/'figures';out.mkdir(exist_ok=True)
    def save(i,fig,sources):
        name=FIGURES[i];fig.savefig(out/(name+'.png'),dpi=140);plt.close(fig)
        write(data/(name+'_plot.json'),dict(figure=name+'.png',source_sha256={s:digest(output/s) for s in sources},
            renderer='particle_3d.particle8_visuals',labels=[FOOTER,BOUNDARY],units='EXPLICIT_AXIS_UNITS',
            data_policy='READ_ONLY_P7_DERIVATIVE_NO_NEW_PHYSICS'))
        print('Figure',name,flush=True)
    fig=canvas('Particle-8  |  Visualization and replay',
        'Two evidence lines, one traceable replay layer  •  Inlet population infrastructure validated within its stated scope')
    ax=axes(fig,[.045,.15,.91,.69]);ax.axis('off')
    boxes=[(.02,.55,.45,.37,'A  /  REAL FROZEN GEOMETRY',
       'Official inlet flux and sample density\nSingle MB entry: original P6.5 saved trajectory\nMixed smoke: limited admission, growing pending queue',MB),
      (.52,.55,.45,.37,'B  /  SYNTHETIC CONTROL',
       'Continuous center motion, exits and deletions\nActual LAMMPS checkpoint / destroy / restart evidence\n10 nm open section; bookkeeping, not vessel passage',RBC),
      (.02,.13,.95,.32,'SHARED REPLAY / AUDIT',
       'Stable ID + original geometry + orientation + physical clock + event ledger\nStatic figures and animations consume the same saved scene data. No resampling.\nScientific boundary: real RBC passage unresolved; no full suspension, PK, CFD or Particle-7.5.',INK)]
    for x,y,w,h,title,body,c in boxes:
        ax.add_patch(plt.Rectangle((x,y),w,h,transform=ax.transAxes,fc=LIGHT,ec=c,lw=1.5))
        ax.text(x+.02,y+h-.05,title,transform=ax.transAxes,color=c,fontsize=14,weight='bold',va='top')
        ax.text(x+.02,y+h-.13,body,transform=ax.transAxes,fontsize=11,va='top',linespacing=1.8)
    save(0,fig,['data/upstream_provenance.json'])

    mesh=read(data/'01_inlet_mesh.json');origin,basis,xy=inlet_basis(mesh);audit=read(data/'01_flux_audit.json')
    counts=np.genfromtxt(data/'02_actual_triangle_counts.csv',delimiter=',',names=True)
    samples=np.genfromtxt(data/'02_actual_inlet_samples.csv',delimiter=',',names=True)
    areas=triangle_area_2d(xy)
    density=np.bincount(samples['triangle_id'].astype(int),minlength=len(xy))/areas
    fig=canvas('01  |  Official inlet flux and sampling density',
        'REAL GEOMETRY / AUDIT  •  Sampling proposals, not finite-size admitted births')
    for x,values,title,label in [(.06,np.mean(mesh['q_m_s'],axis=1)*1e6,'Inward normal speed','µm/s'),
       (.38,np.array(mesh['weights_m3_s'])*1e18,'Per-triangle flow','fL/s'),(.70,density,'100,000 original P7 samples','samples / µm²')]:
        ax=axes(fig,[x,.32,.25,.43]);flux_panel(fig,ax,xy,values,label);ax.set_title(title,fontsize=11)
    Q=audit['INLET']['positive_Q_m3_s'];res=audit['mass_balance']['signed_in_minus_out_m3_s']
    fig.text(.075,.18,f'Q_in = {Q:.12g} m³/s   |   signed inlet − outlet sum = {res:.5g} m³/s\nTriangle density is divided by triangle area. Faster regions have greater expected density.',fontsize=13)
    save(1,fig,['data/01_inlet_mesh.json','data/01_flux_audit.json','data/02_actual_inlet_samples.csv'])

    sched=np.genfromtxt(data/'03_04_scheduler.csv',delimiter=',',names=True)
    fig=canvas('02  |  Each sampled RBC volume determines its own event',
        'REAL-Q SCHEDULER AUDIT  •  Deterministic cumulative flux  •  No mean-volume replacement')
    a=axes(fig,[.08,.30,.40,.46]);b=axes(fig,[.59,.30,.35,.46])
    a.plot(sched['time_s'],sched['rbc_target_m3']*1e18,color=INK,label='0.45 × cumulative blood volume')
    a.plot(sched['time_s'],sched['rbc_scheduled_m3']*1e18,'--',color=RBC,label='Sum of actual sampled RBC volumes')
    a.set(xlabel='Physical scheduler time (s)',ylabel='RBC volume (fL)');a.legend(fontsize=9)
    b.plot(sched['time_s'],sched['next_rbc_volume_m3']*1e18,color=MUTED,lw=1,label='Next sampled RBC volume')
    b.plot(sched['time_s'],sched['rbc_residual_m3']*1e18,color=RBC,label='Unscheduled residual')
    b.set(xlabel='Physical scheduler time (s)',ylabel='Volume (fL)');b.legend(fontsize=9)
    fig.text(.08,.17,'At each recorded checkpoint: 0 ≤ unscheduled residual < next RBC volume.\nFeed volume fraction controls arrival volume; instantaneous lumen Hct is a separate diagnostic.',fontsize=12)
    save(2,fig,['data/03_04_scheduler.csv'])

    fig=canvas('03  |  Microbubble scheduling is deterministic',
        'REAL-Q SCHEDULER AUDIT  •  C_MB = 8.5e12 m^-3  •  NOT Poisson arrivals')
    a=axes(fig,[.08,.31,.40,.45]);b=axes(fig,[.59,.31,.35,.45]);rate=8.5e12*Q;end=float(sched['time_s'][-1])
    t=np.arange(1,int(rate*end)+1)/rate;times=np.r_[0,t,end];values=np.r_[0,np.arange(1,len(t)+1),len(t)]
    a.plot(sched['time_s'],sched['mb_expected'],color=INK,label='C × cumulative blood volume')
    a.step(times,values,where='post',color=MB,label='Scheduled count at exact crossing')
    a.set(xlabel='Physical scheduler time (s)',ylabel='MB count');a.legend(fontsize=9)
    b.plot(sched['time_s'],sched['mb_error'],color=MB);b.axhline(-1,ls='--',color=RBC,label='Strict lower bound (excluded)');b.axhline(0,color=MUTED,lw=1)
    b.set(xlabel='Physical scheduler time (s)',ylabel='Scheduled − expected',ylim=(-1.07,.07));b.legend(fontsize=9)
    fig.text(.08,.17,f'Rate = {rate:.10g} /s; first MB at t = {1/rate:.9f} s.\nEvent times do not depend on the plotting frame rate or the scheduler call partition.',fontsize=12)
    save(3,fig,['data/03_04_scheduler.csv','data/01_flux_audit.json','data/08_timestep_parity.json'])

    view=ReplayView(output,2);view.draw(view.scene['time_range_s'][1]);
    view.fig.texts[0].set_text('04  |  Real inlet: admission and an accumulating queue')
    save(4,view.fig,['data/real_mixed_scene.json','data/13_real_geometry.json'])

    scene=read(data/'synthetic_scene.json');times=np.linspace(0,.25,401);snap=[snapshot(scene,float(t)) for t in times]
    fig=canvas('05  |  Continuous lifecycle and closed population accounting',
        'SYNTHETIC CONTROL / AUDIT  •  P7 actual LAMMPS fixture  •  NOT REAL FROZEN LUMEN')
    a=axes(fig,[.08,.30,.39,.47]);b=axes(fig,[.59,.30,.35,.47])
    for key,color,ls in [('scheduled',INK,'-'),('admitted',RBC,'--'),('exited',GREEN,'-'),('deleted',MB,':')]:
        a.plot(times,[s['counts']['RBC'][key] for s in snap],label=key,color=color,ls=ls,lw=2)
    a.set(xlabel='Physical time (s)',ylabel='Cumulative RBC count');a.legend(fontsize=9)
    b.plot(times,[s['counts']['RBC']['active'] for s in snap],color=RBC,label='Active RBC (sampled)')
    # Show both short MB presence intervals explicitly; coarse frames can miss them.
    for i,r in enumerate(x for x in scene['records'] if x['species']=='MB'):
        b.axvspan(r['admit_time_s'],r['exit_time_s'],color=MB,alpha=.6,label='Exact MB active interval' if i==0 else None)
    b.set(xlabel='Physical time (s)',ylabel='Active RBC count');b.legend(fontsize=9)
    fig.text(.08,.17,'2 MB + 2,356 RBC scheduled and admitted; 2 MB + 2,352 RBC exited/deleted; 4 RBC active.\nDuplicate IDs = 0; lost particles = 0. Finite shapes straddle the 10 nm open caps.',fontsize=12)
    save(5,fig,['data/synthetic_scene.json'])

    queue=read(data/'pending_identity.json');orig=queue['original'];admit=queue['admitted'];g=orig['geometry']
    fig=canvas('06  |  Pending means waiting, not resampling',
        'SYNTHETIC CONTROL / IDENTITY AUDIT  •  Same original particle across failed and successful admission')
    a=axes(fig,[.08,.32,.36,.43]);a.step(queue['times_s'],queue['queue_sizes'],where='post',color=MUTED,lw=3)
    a.set(xlabel='Physical time (s)',ylabel='Pending RBC count',ylim=(-.1,1.2));a.annotate('Same RBC admitted',(.01,0),xytext=(.002,.45),arrowprops=dict(arrowstyle='->',color=RBC))
    a=axes(fig,[.51,.24,.44,.55]);a.axis('off')
    text=f'ID {orig["particle_id"]}  →  ID {admit["particle_id"]}      EXACT IDENTITY\n\nOriginal diameter: {orig["provenance"]["D_um"]:.6f} µm\nOriginal volume: {orig["volume_m3"]*1e18:.6f} fL\nOriginal aspect ratio c/a: {g["c_m"]/g["a_m"]:.6f}\nQuaternion wxyz: {np.array2string(np.array(orig["q"]),precision=4)}\n\nScheduled: {orig["scheduled_time_s"]:.9f} s\nAdmitted: {admit["admitted_time_s"]:.9f} s\nWaiting time: {admit["admitted_time_s"]-orig["scheduled_time_s"]:.9f} s'
    a.text(0,.95,text,va='top',fontsize=12,linespacing=1.55)
    save(6,fig,['data/pending_identity.json'])

    parity=read(data/'10_restart_parity.json');fields=['events','births','exits','pending','active','scheduler']
    fig=canvas('07  |  Restart retains the future',
        'SYNTHETIC CONTROL / RESTART AUDIT  •  Actual saved binary restart; no history regeneration from seed')
    a=axes(fig,[.09,.30,.38,.46]);matrix=np.array([[int(not x[k]) for k in fields] for x in parity['comparisons']])
    im=a.imshow(matrix,aspect='auto',vmin=0,vmax=1,cmap='Reds');a.set(xticks=range(6),xticklabels=fields,ylabel='Saved comparison checkpoint',title='Mismatch matrix: every cell = 0');a.tick_params(axis='x',rotation=25)
    fig.colorbar(im,ax=a,ticks=[0,1],shrink=.8,label='Mismatch')
    a=axes(fig,[.58,.26,.35,.50]);a.axis('off');a.text(0,.97,
      'CHECKPOINT t = 0.125 s\n1,180 scheduled events\n\nDestroy actual LAMMPS instance\nRead binary state + full sidecar\n\nContinue to 2,358 events\nSame stable IDs and event times\nSame pending/admitted/exited records\nSame RNG, next RBC and scheduler\nNo duplicate ID / no particle loss',va='top',fontsize=13,linespacing=1.5)
    save(7,fig,['data/10_restart_parity.json','data/restart_exact_fields.json'])

    from PIL import Image
    fig=canvas('08  |  Storyboard from the exported animations',
        'Actual encoded-scene keyframes  •  A: real entry and congestion  •  B: synthetic lifecycle and restart  •  C: sampling audit',size=(18,13))
    chosen=[]
    for idx,name in enumerate(ANIMATIONS,1):
        manifest=read(data/(name+'_frames.json'));keys=manifest['keyframe_indices']
        key=keys[2] if idx in [1,3] else keys[-1]
        chosen.append((idx,name,key))
    for panel,(idx,name,key) in enumerate(chosen):
        col=panel%2;row=panel//2;ax=axes(fig,[.045+col*.48,.64-row*.265,.44,.235]);ax.axis('off')
        path=output/'frames'/f'{name}_frame_{key:04d}.png';ax.imshow(Image.open(path));ax.set_title(f'Animation {idx:02d} / frame {key}',fontsize=12)
    ax=axes(fig,[.53,.125,.40,.19]);ax.axis('off');ax.text(0,.90,'Read the clock and classification first.\n\nReal smoke does not show RBC transport.\nSynthetic centers move through a 10 nm section.\nFlux samples are proposals, not admitted births.\nAll display cuts and interpolation are recorded.',va='top',fontsize=14,linespacing=1.5)
    save(8,fig,[f'frames/{name}_frame_{key:04d}.png' for idx,name,key in chosen])

    fig=canvas('09  |  Scientific boundaries remain visible',
        'AUDIT / CARRY-FORWARD LIMITATIONS  •  Visualization improves inspection; it does not complete missing physics')
    a=axes(fig,[.075,.17,.85,.64]);a.axis('off')
    items=[('REAL RBC PASSAGE','Continuous admission and passage remain unresolved.'),
      ('SYNTHETIC CONTROL','Lifecycle demonstration is not physiological lumen proof.'),
      ('NUMERICAL CHOICES','Production timestep and neighbor settings are NOT FROZEN.'),
      ('NEAR-FIELD MODEL','Sphere-normal regularization only; non-spherical lubrication not frozen.'),
      ('TIME ACCURACY','P6.5 handoff event-time convergence has not been established.'),
      ('SUSPENSION / PK','Full mixed suspension physiology and full MB PK are not established.'),
      ('THIS STAGE','No CFD. No Particle-7.5. No new RBC deformation physics.'),
      ('DISPLAY vs PHYSICS','Interpolation, time cuts and stretched axes are explicitly labeled.')]
    for i,(left,right) in enumerate(items):
        y=.96-i*.12;a.text(0,y,left,color=RBC,weight='bold',fontsize=12,va='top');a.text(.28,y,right,fontsize=12,va='top')
    save(9,fig,['data/upstream_provenance.json','data/synthetic_scene.json','data/real_mixed_scene.json'])
