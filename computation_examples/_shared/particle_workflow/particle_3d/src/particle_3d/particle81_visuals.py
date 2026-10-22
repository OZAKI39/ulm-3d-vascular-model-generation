"""Academic full-vessel ULM displays, exclusively from the saved P8.1 scene."""
from pathlib import Path
import hashlib,os
# Limit the software rasterizer while independent integrations use other cores.
os.environ.setdefault("LP_NUM_THREADS","2")
os.environ.setdefault("VTK_SMP_MAX_THREADS","2")
from collections import OrderedDict
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image,ImageDraw
from .particle8_replay import read,digest,canonical_hash
from .particle8_full3d_visuals import text,surface,draw_camera
from .particle81_simulation import OUTPUT,dump
from .particle81_replay import Scene,TITLE,INDEPENDENT,FOOTERS,OUTLETS,frame_schedule

COLORS={'OUTLET_01':'#c76b46','OUTLET_02':'#208ca6','OUTLET_03':'#8065a6'}
INK='#183449';MUTED='#78909c';ACTIVE='#d66a3c';SIZE=(1600,1000);FPS=15
ANIMATIONS=['01_full_vessel_replay','02_accumulation','03_single_journeys','04_outlet_ensemble']


class VesselView:
    def __init__(self,scene,size=(1090,665),wall_opacity=.10):
        self.scene=scene;self.size=size;self.cap_labels={};self.active=[];self.plotter=pv.Plotter(off_screen=True,window_size=size)
        p=self.plotter;p.set_background('white');p.render_window.SetMultiSamples(0)
        p.enable_depth_peeling(number_of_peels=8,occlusion_ratio=0)
        geo=np.load(scene.root/'data/full_frozen_geometry.npz');wall=geo['WALL_points_m']
        self.origin=(wall.max(0)+wall.min(0))/2;self.scale=np.linalg.norm((wall.max(0)-wall.min(0))*1e6)*.51
        for role in ['WALL','INLET',*OUTLETS]:
            mesh=surface((geo[role+'_points_m']-self.origin)*1e6,geo[role+'_faces'])
            p.add_mesh(mesh,color='#acbdc7' if role=='WALL' else '#278276' if role=='INLET' else COLORS[role],
                opacity=wall_opacity if role=='WALL' else .4,smooth_shading=True,render=False)
            if role!='WALL':
                center=((geo[role+'_points_m']-self.origin)*1e6).mean(0)
                self.cap_labels[role]=center
        p.add_axes(color=INK,xlabel='x',ylabel='y',zlabel='z',line_width=2)
        p.show(auto_close=False,interactive=False);self.dynamic=[];self.path_cache=OrderedDict()

    def clear(self):
        self.active=[]
        for name in self.dynamic:self.plotter.remove_actor(name,reset_camera=False,render=False)
        self.dynamic=[]

    def paths(self,ids,name,opacity=.65,width=1.0,color_by_outlet=True,partial=None):
        # Static history geometry is reused across orbit frames. No physical
        # samples are decimated; at most two history meshes stay in memory.
        key=(tuple(ids),color_by_outlet) if partial is None else None
        if key is not None and key in self.path_cache:
            mesh,kw=self.path_cache[key];self.path_cache.move_to_end(key)
        else:
            points=[];lines=[];codes=[];offset=0
            for pid in ids:
                a=self.scene.arrays[pid][:,1:4] if partial is None else partial[pid]
                if len(a)<2:continue
                points.append((a-self.origin)*1e6);lines.extend([len(a),*range(offset,offset+len(a))]);offset+=len(a)
                outlet=self.scene.entries[pid]['exit_outlet'];codes.extend([OUTLETS.index(outlet)+1 if outlet in OUTLETS else 0]*len(a))
            if not points:return
            mesh=pv.PolyData(np.vstack(points),lines=np.asarray(lines,dtype=np.int64))
            if color_by_outlet:
                mesh['outlet_code']=np.asarray(codes)
                kw=dict(scalars='outlet_code',cmap=['#a7b4ba',*COLORS.values()],clim=(0,3),show_scalar_bar=False)
            else:kw=dict(color='#478ea4')
            if key is not None:
                self.path_cache[key]=(mesh,kw)
                while len(self.path_cache)>2:self.path_cache.popitem(last=False)
        self.plotter.add_mesh(mesh,**kw,name=name,line_width=width,opacity=opacity,render=False,reset_camera=False)
        self.dynamic.append(name)

    def current(self,active):
        self.active=active
        for a in active:
            center=(np.asarray(a['position_m'])-self.origin)*1e6
            s=pv.Sphere(radius=1.,theta_resolution=16,phi_resolution=12)
            s.points=np.asarray(s.points,dtype=float)*(a['radius_m']*1e6)+center
            name=f'active_{a["particle_id"]}'
            self.plotter.add_mesh(s,color=ACTIVE,opacity=1,name=name,smooth_shading=True,render=False,reset_camera=False);self.dynamic.append(name)


    def capture(self,azimuth=65):
        camera=draw_camera(self.plotter,np.zeros(3),self.scale,dict(azimuth_deg=float(azimuth),elevation_deg=28,role='DISPLAY_ONLY'))
        self.plotter.render();pixels=self.plotter.screenshot(return_img=True)[:,:,:3]
        image=Image.fromarray(pixels);draw=ImageDraw.Draw(image)
        def project(point):
            renderer=self.plotter.renderer
            renderer.SetWorldPoint(*point,1.);renderer.WorldToDisplay()
            x,y,_=renderer.GetDisplayPoint()
            return x,self.size[1]-y
        for role,point in self.cap_labels.items():
            x,y=project(point);dx,dy=(-85,8) if role=='OUTLET_03' else (12,-20)
            label_x=float(np.clip(x+dx,5,self.size[0]-100));label_y=float(np.clip(y+dy,5,self.size[1]-22))
            text(draw,(label_x,label_y),role,14,INK,True)
        for a in self.active:
            x,y=project((np.asarray(a['position_m'])-self.origin)*1e6)
            # A screen-space ID callout identifies the original-size sphere;
            # this line is an annotation, never a saved or invented trajectory.
            draw.line((x,y,x-14,y-14),fill=ACTIVE,width=1)
            text(draw,(x-80,y-32),f'MB {a["particle_id"]}',14,ACTIVE,True)
        return image,camera

    def close(self):self.plotter.close()


def compose_frame(scene,view,kind,index,total,spec):
    state=scene.snapshot(spec['time_s'],spec['only_ids']);view.clear()
    finished=state['completed_ids'];partial={p['particle_id']:scene.partial_path(p['particle_id'],p['elapsed_time_s']) for p in state['active']}
    if kind=='03_single_journeys':
        view.paths(finished,'history',opacity=.9,width=3.)
        view.paths(list(partial),'partial',opacity=.9,width=3.,partial=partial)
    else:
        view.paths(finished,'history',opacity=.55 if kind=='02_accumulation' else .70 if kind=='04_outlet_ensemble' else .35,width=1.5 if kind in ['02_accumulation','04_outlet_ensemble'] else 1.)
        if kind!='02_accumulation':view.paths(list(partial),'partial',opacity=.95,width=3.,partial=partial)
    view.current(state['active'])
    az=30.+90.*index/max(1,total-1);picture,camera=view.capture(az)
    image=Image.new('RGB',SIZE,'white');image.paste(picture,(20,174));d=ImageDraw.Draw(image)
    names={'01_full_vessel_replay':'Full-vessel rotating replay','02_accumulation':'ULM accumulation / trajectory reconstruction',
        '03_single_journeys':'Representative inlet-to-outlet journeys','04_outlet_ensemble':'Outlet-colored ensemble replay'}
    text(d,(40,22),'Particle-8.1 | '+names[kind],29,bold=True)
    text(d,(40,66),TITLE,18,'#208ca6',True);d.line((40,100,1560,100),fill='#208ca6',width=3)
    text(d,(40,111),f'Physical acquisition t = {spec["time_s"]:.6f} s  |  frame {index+1}/{total}  |  {spec["label"]}',15)
    text(d,(40,140),'REAL FROZEN VESSEL + FLOW | complete microbubble trajectory visualization layer | original SI coordinates',16,bold=True)
    count_scope='SELECTED MB ACCOUNTING' if spec['only_ids'] is not None else 'ACQUISITION ACCOUNTING'
    text(d,(1132,187),count_scope,18,bold=True)
    for i,(label,key) in enumerate([('Scheduled','scheduled'),('Admitted','admitted'),('Active saved tracks','active'),('Completed exits','completed'),
        ('Unresolved admission','unresolved'),('Stopped incomplete','stopped')]):
        y=230+i*34;text(d,(1132,y),label,17);text(d,(1515,y),str(state['counts'][key]),18,bold=True)
    text(d,(1132,455),'PATH COLOR = FINAL OUTLET',16,bold=True)
    for i,o in enumerate(OUTLETS):
        text(d,(1132,491+i*31),o,17,COLORS[o],True)
        text(d,(1515,491+i*31),str(sum(scene.entries[pid]['exit_outlet']==o for pid in finished)),17,COLORS[o])
    for i,line in enumerate(['Active MB: original sphere + ID callout','History: saved completed polylines',
                            'Stopped tracks: no extrapolation','No resampling of MB attributes']):text(d,(1132,615+i*26),line,15)
    if state['active']:
        a=state['active'][0];e=scene.entries[a['particle_id']]
        text(d,(1132,738),f'MB {a["particle_id"]} | D = {e["diameter_um"]:.3f} µm',17,ACTIVE,True)
        text(d,(1132,768),f'Age = {a["elapsed_time_s"]*1000:.3f} ms',16)
    elif kind=='03_single_journeys':text(d,(1132,738),'Outlet reached; MB removed',17,'#208ca6',True)
    else:text(d,(1132,738),'No active MB at this physical time',15)
    covered=sum(bool(n) for n in scene.catalog['outlet_counts'].values())
    text(d,(1132,802),f'Dataset outlet coverage: {covered}/3',15,'#884b43',True)
    length=20.;px=length*665/(2*view.scale);d.line((240,816,240+px,816),fill=INK,width=3);text(d,(240,788),'20 µm | equal physical scale',14)
    text(d,(40,850),f'Camera orbit {az:.1f}° / elevation 28° — DISPLAY ONLY | display interpolation / slow motion / compression / cuts; raw data unchanged',14)
    for i,line in enumerate(FOOTERS):text(d,(40,883+25*i),line,13,'#884b43')
    receipt=dict(frame_index=index,physical_time_s=spec['time_s'],window_label=spec['label'],only_ids=spec['only_ids'],
        state=state,state_sha256=canonical_hash(state),count_scope=count_scope,rendered_completed_ids=finished,rendered_active_ids=[p['particle_id'] for p in state['active']],
        active_partial_tails_drawn=kind!='02_accumulation',camera=camera,labels=[TITLE,INDEPENDENT,*FOOTERS],
        physical_coordinates_rotated=False,physical_axis_stretching=False,source_scene_sha256=scene.sha256,
        render_window_backend=view.plotter.render_window.GetClassName(),screen_space_id_callout=True,
        polyline_decimation=False,display_interpolation_only=True,time_compression_or_cuts_explicit=True,
        time_mapping='ACQUISITION_COMPRESSION' if kind=='02_accumulation' else 'SLOW_MOTION_EVENT_WINDOWS_WITH_EXPLICIT_IDLE_CUTS')
    return image,receipt


def render_animations(root=OUTPUT,kinds=None,preview=False):
    scene=Scene(root);root=Path(root)
    for folder in ['animations','frames','keyframes','inspection']:(root/folder).mkdir(exist_ok=True)
    for kind in kinds or ANIMATIONS:
        schedule=frame_schedule(scene,kind);view=VesselView(scene,wall_opacity=.045 if kind=='02_accumulation' else .12)
        name='particle8_1_anim_'+kind
        if preview:
            i=len(schedule)//2;im,_=compose_frame(scene,view,kind,i,len(schedule),schedule[i]);im.save(root/'keyframes'/(name+'_preview.png'));view.close();continue
        encoder=imageio_ffmpeg.write_frames(str(root/'animations'/(name+'.mp4')),SIZE,fps=FPS,codec='libx264',quality=8,
            pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=2,output_params=['-movflags','+faststart'],ffmpeg_log_level='error');encoder.send(None)
        receipts=[];gifs=[];keys=[0,len(schedule)//2,len(schedule)-1]
        try:
            for i,spec in enumerate(schedule):
                im,record=compose_frame(scene,view,kind,i,len(schedule),spec);encoder.send(np.asarray(im));receipts.append(record)
                if i in keys:im.save(root/'keyframes'/f'{name}_frame_{i:04d}.png')
                if i%5==0 or i in keys:gifs.append(im.resize((960,600)))
                if i%60==0:print(kind,i,'/',len(schedule),flush=True)
        finally:encoder.close();view.close()
        gifs[0].save(root/'animations'/(name+'.gif'),save_all=True,append_images=gifs[1:],duration=333,loop=0,optimize=False)
        dump(root/'frames'/(name+'_frames.json'),dict(kind=kind,source_scene_sha256=scene.sha256,fps=FPS,size=SIZE,frame_count=len(receipts),keyframe_indices=keys,frames=receipts))
        print('WRITTEN',name,flush=True)
