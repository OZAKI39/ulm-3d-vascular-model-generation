"""Full vessel camera-orbit rendering. Coordinates from P8 remain authoritative."""
from pathlib import Path
import math
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image,ImageDraw,ImageFont
from .particle8_replay import read,write,digest,canonical_hash
from .particle8_visuals import RBC,MB,INK,MUTED,GREEN,surface_points
from .rbc_orientation import rotation_matrix
from .particle8_full3d_data import OUTPUT,CASES,LABELS,safe_output,scene_from_name,schedule,camera_angles,replay_frame

FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD='/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'
SIZE=(1600,1000);FPS=15


def text(draw,xy,value,size=18,color=INK,bold=False):
    draw.text(xy,value,font=ImageFont.truetype(BOLD if bold else FONT,size),fill=color)


def surface(points,faces):
    return pv.PolyData(np.asarray(points,dtype=float),np.column_stack([np.full(len(faces),3),faces]).ravel())


def display_position(position,origin,stretch=False):
    q=(np.asarray(position,dtype=float)-origin)*1e6
    if stretch:q=q*np.array([1.,1.,20000.])
    return q


def shape_mesh(record,state,origin):
    """Render original shape, never use capsule's retained quaternion as its axis."""
    center=display_position(state['position_m'],origin);shape=record['admitted_shape']
    if shape['mode']=='SPHERE_MB':
        sphere=pv.Sphere(radius=1,theta_resolution=24,phi_resolution=18)
        sphere.points=np.asarray(sphere.points,dtype=float)*(shape['radius_m']*1e6)+center
        return sphere
    if shape['mode']=='FREE_OBLATE':
        sphere=pv.Sphere(radius=1,theta_resolution=24,phi_resolution=18)
        sphere.points=(np.asarray(sphere.points)*np.array(shape['axes_m'])*1e6)@rotation_matrix(state['q']).T+center
        return sphere
    points=(surface_points(shape,state['position_m'])-origin)*1e6
    return pv.StructuredGrid(points[:,:,0],points[:,:,1],points[:,:,2]).extract_surface(algorithm='dataset_surface')


def draw_camera(plotter,center,scale,angles):
    a=math.radians(angles['azimuth_deg']);e=math.radians(angles['elevation_deg'])
    position=np.asarray(center)+4*scale*np.array([math.cos(a)*math.cos(e),math.sin(a)*math.cos(e),math.sin(e)])
    plotter.camera_position=[position,center,[0,0,1]];plotter.camera.parallel_projection=True
    plotter.camera.parallel_scale=scale;plotter.reset_camera_clipping_range()
    return dict(position_display_units=position.tolist(),focal_point_display_units=np.asarray(center).tolist(),
                up=[0,0,1],parallel_scale_display_units=scale,**angles)


def context_note(case,state,scene):
    if case.startswith('A_'):
        if not state['particles']:return 'Before MB birth | no active particles'
        record=scene['records'][0]
        distance=np.linalg.norm(np.asarray(state['particles'][0]['position_m'])-record['trajectory'][0]['position_m'])*1e9
        return f'Displacement: {distance:.2f} nm | original MB size'
    if case.startswith('B_'):
        return 'Admitted RBC: STATIC capsule surrogate' if state['particles'] else 'No admitted RBC yet | pending count only'
    if case.startswith('D_'):return 'MISMATCH = 0 | saved P7 restart replay'
    return '10 nm open section; shapes straddle caps'


class View:
    def __init__(self,root,case,tail_mode='recent'):
        self.root=Path(root);self.case=case;self.config=CASES[case];self.tail_mode=tail_mode
        self.scene=scene_from_name(root,self.config['scene']);self.records={r['particle_id']:r for r in self.scene['records']}
        self.synthetic=self.scene['source_classification']=='synthetic_control';self.restart=case.startswith('D_')
        self.other=scene_from_name(root,'synthetic_restarted') if self.restart else None
        self.main=pv.Plotter(off_screen=True,window_size=(1000,640));self.secondary=pv.Plotter(off_screen=True,window_size=(500,365))
        for p in [self.main,self.secondary]:
            p.set_background('white');p.render_window.SetMultiSamples(0);p.enable_depth_peeling(number_of_peels=8,occlusion_ratio=0)
            p.add_axes(line_width=2,xlabel='x',ylabel='y',zlabel='z',color=INK)
        self.dynamic=[[],[]]
        if self.synthetic:self.synthetic_geometry()
        else:self.real_geometry()
        # Fixed origin and focal targets are independent of animation time.
        for p in [self.main,self.secondary]:p.show(auto_close=False,interactive=False)

    def real_geometry(self):
        geo=np.load(self.root/'data/full_frozen_geometry.npz');wall=geo['WALL_points_m']
        self.origin=(wall.min(0)+wall.max(0))/2;extent=(wall.max(0)-wall.min(0))*1e6
        self.main_center=np.zeros(3);self.main_scale=np.linalg.norm(extent)*.56
        centers=[];roles=[]
        for role in ['WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
            points=display_position(geo[role+'_points_m'],self.origin);mesh=surface(points,geo[role+'_faces'])
            for j,p in enumerate([self.main,self.secondary]):
                p.add_mesh(mesh.copy(),color='#a8bbc8' if role=='WALL' else GREEN if role=='INLET' else '#876399',
                    opacity=(.22 if role=='WALL' else .52) if j==0 else (.08 if role=='WALL' else .06),
                    smooth_shading=True,show_scalar_bar=False,render=False)
            if role!='WALL':centers.append(points.mean(0));roles.append(role)
        self.main.add_point_labels(np.array(centers),roles,font_size=14,text_color=INK,shape_color='white',
            shape_opacity=.75,show_points=False,always_visible=True,render=False)
        inlet=display_position(geo['INLET_points_m'],self.origin).mean(0)
        self.secondary_center=display_position(self.scene['records'][0]['trajectory'][0]['position_m'],self.origin) if self.case.startswith('A_') else inlet
        self.secondary_scale=1.75 if self.case.startswith('A_') else 8.3
        # Reference ring identifies where the local inset sits on the full vessel.
        self.main.add_mesh(pv.Sphere(radius=3.6,center=inlet,theta_resolution=16,phi_resolution=12),
            style='wireframe',color=GREEN,opacity=.12,line_width=.7,render=False)

    def synthetic_geometry(self):
        self.origin=np.zeros(3);self.main_center=np.array([0.,0.,100.]);self.main_scale=182.
        self.secondary_center=self.main_center if self.restart else np.array([0,0,.005])
        self.secondary_scale=182. if self.restart else 180.
        for j,p in enumerate([self.main,self.secondary]):
            stretched=(j==0 or self.restart);length=200. if stretched else .01
            box=pv.Box(bounds=(-100,100,-100,100,0,length))
            p.add_mesh(box,color='#9bb4c3',style='wireframe',opacity=.55,line_width=1.5,render=False)
            for z,c in [(0,GREEN),(length,'#876399')]:
                plane=pv.Plane(center=(0,0,z),direction=(0,0,1),i_size=200,j_size=200,i_resolution=1,j_resolution=1)
                p.add_mesh(plane,color=c,opacity=.08,render=False)
            if stretched:
                p.add_point_labels(np.array([[0,-110,0],[0,-110,length]]),['INLET','OUTLET / DELETE'],
                    font_size=12,show_points=False,always_visible=True,shape_opacity=.8,text_color=INK,render=False)

    def add_particles(self,plotter,j,state,stretched):
        for name in self.dynamic[j]:plotter.remove_actor(name,reset_camera=False,render=False)
        self.dynamic[j]=[];pos=[];labels=[]
        for particle in state['particles']:
            pid=particle['particle_id'];record=self.records[pid];color=MB if particle['species']=='MB' else RBC
            center=display_position(particle['position_m'],self.origin,stretched)
            mesh=pv.Sphere(radius=2.6 if particle['species']=='MB' else 2.1,center=center,theta_resolution=14,phi_resolution=12) if stretched else shape_mesh(record,particle,self.origin)
            wire=(not self.synthetic and j==1 and particle['species']=='MB')
            name=f'particle_{pid}';plotter.add_mesh(mesh,color=color,opacity=.06 if wire else .9 if stretched else .72,
                smooth_shading=True,name=name,render=False,reset_camera=False);self.dynamic[j].append(name)
            if wire:
                radius=record['geometry']['radius_m']*1e6;theta=np.linspace(0,2*np.pi,129)
                for k,(a,b) in enumerate([(0,1),(0,2),(1,2)]):
                    pts=np.zeros((len(theta),3));pts[:,a]=radius*np.cos(theta);pts[:,b]=radius*np.sin(theta);pts+=center
                    name=f'shell_ring_{pid}_{k}';plotter.add_mesh(pv.lines_from_points(pts),color=MB,opacity=.3,line_width=1.,name=name,render=False,reset_camera=False);self.dynamic[j].append(name)
                name=f'center_{pid}';plotter.add_mesh(pv.PolyData(np.array([center])),color='#125971',point_size=6,
                    render_points_as_spheres=True,name=name,render=False,reset_camera=False);self.dynamic[j].append(name)
            tail=particle['tail'];points=display_position([s['position_m'] for s in tail],self.origin,stretched) if tail else np.empty((0,3))
            if len(points)>1 and np.any(np.linalg.norm(np.diff(points,axis=0),axis=1)>0):
                name=f'tail_{pid}';plotter.add_mesh(pv.lines_from_points(points),color='#125971' if wire else color,line_width=3.5 if wire else 2.4 if not stretched else 1.8,
                    opacity=.95 if wire else .7,name=name,render=False,reset_camera=False);self.dynamic[j].append(name)
            if self.synthetic or j==1:
                offset=np.array([0.,0.,record['geometry'].get('radius_m',0.)*1.5e6]) if wire else np.zeros(3)
                pos.append(center+offset);labels.append(f'{particle["species"]} {pid}')
        if pos:
            name='active_ids';plotter.add_point_labels(np.array(pos),labels,font_size=12 if j else 14,show_points=False,
                text_color=INK,shape_color='white',shape_opacity=.65,always_visible=True,name=name,render=False,reset_camera=False)
            self.dynamic[j].append(name)

    def frame(self,index,total,t,window):
        state=replay_frame(self.scene,t,self.config['tail_window_s'],self.tail_mode)
        other=replay_frame(self.other,t,self.config['tail_window_s'],self.tail_mode) if self.restart else state
        if self.restart and state!=other:raise ValueError('Restart physical state mismatch')
        self.add_particles(self.main,0,state,self.synthetic)
        self.add_particles(self.secondary,1,other,self.restart)
        angles=camera_angles(self.case,index,total)
        camera_main=draw_camera(self.main,self.main_center,self.main_scale,angles)
        inset_angles=angles if self.synthetic else dict(angles,azimuth_deg=angles['azimuth_deg']+90)
        camera_second=draw_camera(self.secondary,self.secondary_center,self.secondary_scale,inset_angles)
        self.main.render();self.secondary.render()
        canvas=Image.new('RGB',SIZE,'white');canvas.paste(Image.fromarray(self.main.screenshot(return_img=True)[:,:,:3]),(25,180))
        canvas.paste(Image.fromarray(self.secondary.screenshot(return_img=True)[:,:,:3]),(1060,185))
        d=ImageDraw.Draw(canvas);color=RBC if self.synthetic else MB
        text(d,(42,24),'Particle-8 | '+self.config['title'],30,bold=True)
        d.rectangle((42,72,1558,76),fill=color)
        text(d,(42,87),self.config['classification'],19,color,bold=True)
        clock=f'Physical t = {t:.9f} s  |  frame {index+1}/{total}  |  {window}'
        if self.case.startswith('A_'):clock+=f'  |  since birth {(t-self.scene["records"][0]["admit_time_s"])*1000:+.3f} ms'
        text(d,(42,122),clock,18)
        left='FULL FROZEN VESSEL | all wall faces + inlet + 3 outlets | equal physical scale'
        right='LOCAL INLET ZOOM | equal scale'
        if self.synthetic:
            left='3D CONTROL | z ×20,000 DISPLAY ONLY | centers are ID glyphs, not finite shapes'
            right='UNSTRETCHED | original shapes'
            if self.restart:left='CONTINUOUS | center ID glyphs; z ×20,000 DISPLAY ONLY';right='RESTARTED | same center glyphs'
        text(d,(42,158),left,15,bold=True);text(d,(1060,158),right,15,bold=True)
        text(d,(1060,557),'COUNTS — pending is not drawn in-lumen',15,bold=True)
        text(d,(1065,585),'Status',15);text(d,(1360,585),'RBC',15,RBC,bold=True);text(d,(1470,585),'MB',15,MB,bold=True)
        for i,key in enumerate(['scheduled','admitted','pending','active','exited','deleted']):
            y=612+i*23;text(d,(1065,y),key,16);text(d,(1360,y),str(state['counts']['RBC'][key]),16,RBC);text(d,(1470,y),str(state['counts']['MB'][key]),16,MB)
        note=context_note(self.case,state,self.scene)
        text(d,(1060,768),note,15,color,bold=True)
        # Orthographic scale bars are derived from actual camera parallel scale.
        for x,y,pixels,scale,length in [(250,797,640,self.main_scale,50 if self.synthetic else 20),(1220,534,365,self.secondary_scale,50 if self.synthetic else 1 if self.case.startswith('A_') else 5)]:
            n=length*pixels/(2*scale);d.line((x,y,x+n,y),fill=INK,width=3)
            label=f'{length} µm'+(' transverse scale' if self.synthetic else '')
            text(d,(x,y-23),label,13)
        text(d,(42,827),f'RBC: warm  |  MB: cool  |  thin line: {self.tail_mode} pathline  |  wall: blue-grey  |  inlet: green  |  outlets: violet',16)
        text(d,(42,853),f'Camera {self.config["camera"]}: main azimuth {angles["azimuth_deg"]:.1f}°, inset {inset_angles["azimuth_deg"]:.1f}°, elevation 28°; DISPLAY ONLY. Coordinates do not rotate.',15)
        text(d,(42,878),f'{self.tail_mode.upper()} tails: active particles only. DISPLAY INTERPOLATION between saved states; no new physical solution.',14)
        for y,label in [(907,LABELS[2]),(930,LABELS[3]),(953,LABELS[4])]:text(d,(42,y),label,13,'#8a4036')
        if self.case.startswith('B_'):text(d,(42,802),'Camera motion must not be read as RBC motion. The saved admitted shape is held fixed.',14,color)
        receipt=dict(frame_index=index,video_time_s=index/FPS,physical_time_s=t,window_label=window,
            case=self.case,classification=self.config['classification'],state=state,state_sha256=canonical_hash(state),
            rendered_ids_main=[p['particle_id'] for p in state['particles']],rendered_ids_secondary=[p['particle_id'] for p in other['particles']],
            camera_main=camera_main,camera_secondary=camera_second,labels=LABELS,context_note=note,
            coordinates_rotated=False,synthetic_z_display_factor=20000 if self.synthetic else 1,
            restart_mismatches=0 if self.restart else None,restarted_state_sha256=canonical_hash(other) if self.restart else None)
        return canvas,receipt

    def close(self):
        self.main.close();self.secondary.close()


def render(output=OUTPUT,cases=None,tail_mode='recent',preview=False):
    root=safe_output(output)
    for folder in ['mp4','gif','frames','keyframes','storyboard']: (root/folder).mkdir(parents=True,exist_ok=True)
    for case in (CASES if cases is None else cases):
        v=View(root,case,tail_mode);times=schedule(case,v.scene);name='particle8_full3d_'+case+('_full_tail' if tail_mode=='full' else '')
        if preview:
            k=len(times)//2;im,receipt=v.frame(k,len(times),*times[k]);im.save(root/'keyframes'/(name+'_preview.png'));v.close();continue
        keyindices=sorted(set([0,len(times)//2,len(times)-1]));frames=[];gif=[]
        encoder=imageio_ffmpeg.write_frames(str(root/'mp4'/(name+'.mp4')),SIZE,fps=FPS,codec='libx264',quality=8,
            pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=2,output_params=['-movflags','+faststart'],ffmpeg_log_level='error')
        encoder.send(None)
        try:
            for k,(t,label) in enumerate(times):
                im,record=v.frame(k,len(times),t,label);encoder.send(np.asarray(im));frames.append(record)
                if k in keyindices:im.save(root/'keyframes'/f'{name}_frame_{k:04d}.png')
                if k%5==0:gif.append(im.resize((960,600)))
                if k%60==0:print(case,k,'/',len(times),flush=True)
        finally:encoder.close();v.close()
        gif[0].save(root/'gif'/(name+'.gif'),save_all=True,append_images=gif[1:],duration=333,loop=0,optimize=False)
        write(root/'frames'/(name+'_manifest.json'),dict(name=name,case=case,source_scene=read(root/'data/discovered_inputs.json')['scenes'][CASES[case]['scene']],
            fps=FPS,size=SIZE,frame_count=len(frames),keyframe_indices=keyindices,camera_policy=CASES[case]['camera'],tail_mode=tail_mode,frames=frames))
        print('Written',name,flush=True)
