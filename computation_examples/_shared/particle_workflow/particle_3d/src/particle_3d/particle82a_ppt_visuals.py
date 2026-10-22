"""Clean black-background PPT views of saved finite-size MBs and their tracks."""
from pathlib import Path
import argparse,json,math,os,socket,time
os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
import numpy as np
import pyvista as pv
import matplotlib
from matplotlib.colors import ListedColormap
from PIL import Image,ImageDraw,ImageFont
import imageio_ffmpeg
from .particle82_provenance import atomic_json,sha256
from .particle82a_media import decode_video

CMAP=ListedColormap(matplotlib.colormaps['turbo'](np.linspace(.10,.96,256)))
SIZE=(1920,1080);FPS=24;FRAMES=432
ANIMATIONS=['01_full_vessel_microbubbles','02_branch_microbubbles','03_single_microbubble_journey']


def sample_at_age(samples,age):
    if not len(samples) or age<samples[0,0] or age>samples[-1,0]:raise ValueError('Replay extrapolation forbidden')
    k=min(int(np.searchsorted(samples[:,0],age,side='right'))-1,len(samples)-1)
    if k==len(samples)-1:return samples[k,1:4].copy(),samples[k,4:7].copy(),k
    weight=(age-samples[k,0])/(samples[k+1,0]-samples[k,0])
    position=(1-weight)*samples[k,1:4]+weight*samples[k+1,1:4]
    return position,samples[k+1,4:7].copy(),k


def surface(points,faces):return pv.PolyData(points,np.column_stack([np.full(len(faces),3),faces]).ravel())


def text(draw,xy,value,size,color='#eef2f8'):
    draw.text(xy,value,font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',size),fill=color)


def annotate(pixels,title,footer,limit):
    image=Image.fromarray(np.asarray(pixels)[:,:,:3]);d=ImageDraw.Draw(image);w,h=image.size
    d.rectangle((0,0,w,int(h*.083)),fill='black');d.rectangle((0,int(h*.936),w,h),fill='black')
    text(d,(int(w*.025),int(h*.026)),title,int(h*.032))
    text(d,(int(w*.025),int(h*.958)),footer,int(h*.020),'#bec9d5')
    d.rectangle((int(w*.86),int(h*.085),w,int(h*.935)),fill='black')
    x=int(w*.885);y=int(h*.265);bw=int(w*.019);bh=int(h*.50)
    rgb=(CMAP(np.linspace(1,0,bh))[:,:3]*255).astype('uint8')
    image.paste(Image.fromarray(np.repeat(rgb[:,None,:],bw,axis=1)),(x,y))
    text(d,(int(w*.87),y-int(h*.063)),'Track speed',int(h*.022))
    text(d,(int(w*.87),y-int(h*.034)),'(mm/s)',int(h*.020))
    for value in np.linspace(0,limit,6):
        yy=y+bh*(1-value/limit);text(d,(x+bw+int(w*.007),int(yy-h*.012)),f'{value:.2f}',int(h*.021))
    return image


class Scene:
    def __init__(self,root):
        self.root=Path(root);self.meta=json.loads((self.root/'data/SCENE.json').read_text());self.sha=sha256(self.root/'data/SCENE.json')
        assert sha256(self.root/'data/geometry.npz')==self.meta['geometry_sha256']
        self.geo=np.load(self.root/'data/geometry.npz');self.records={r['event_id']:r for r in self.meta['records']};self.arrays={}
        for pid,r in self.records.items():
            file=self.root/r['samples_path'];assert sha256(file)==r['samples_sha256'];self.arrays[pid]=np.load(file)['samples']
        wall=self.geo['WALL_points_m'];self.origin=(wall.min(0)+wall.max(0))/2
        self.wall_um=(wall-self.origin)*1e6;self.branch=(self.geo['branch_focus_m']-self.origin)*1e6
        self.radius=float(np.linalg.norm(self.wall_um,axis=1).max())
        radial=np.linalg.norm(self.wall_um[:,:2],axis=1);phi=np.deg2rad(22)
        vertical=np.abs(self.wall_um[:,2])*np.cos(phi)+radial*np.sin(phi)
        self.scale=float(max(vertical.max(),radial.max()/1.70)/.88)
        self.limit=math.ceil(self.meta['summary']['speed_max_mm_s']*10)/10
        completed=[pid for pid,r in self.records.items() if r['completed']]
        self.replay_ids=completed[:48]
        # Display selection is explicit and separate from the unchanged 500-track cohort.
        for outlet in ['OUTLET_01','OUTLET_03']:
            extras=[pid for pid in completed if self.records[pid]['outlet']==outlet][:12]
            self.replay_ids.extend(pid for pid in extras if pid not in self.replay_ids)
        self.single=max(completed,key=lambda pid:self.records[pid]['path_length_um'])
        self.window=float(np.quantile([self.records[pid]['last_age_s'] for pid in completed],.9))

    def states(self,progress,single=False):
        states=[];ids=[self.single] if single else self.replay_ids
        for rank,pid in enumerate(ids):
            r=self.records[pid]
            age=progress*r['last_age_s'] if single else progress*self.window+(.8-1.6*rank/max(1,len(ids)-1))*self.window
            if not 0<=age<r['last_age_s']:continue
            position,velocity,k=sample_at_age(self.arrays[pid],age)
            states.append(dict(event_id=pid,age_s=float(age),position_m=position.tolist(),radius_m=r['radius_m'],
                saved_interval_index=k,speed_mm_s=float(np.linalg.norm(velocity)*1000)))
        return states

    def paths(self,ids,through=None):
        positions=[];speed=[];lines=[];offset=0
        for pid in ids:
            a=self.arrays[pid]
            if through is not None:
                pos,vel,k=sample_at_age(a,through);a=np.vstack([a[:k+1],a[k].copy()]);a[-1,0]=through;a[-1,1:4]=pos;a[-1,4:7]=vel
            if len(a)<2:continue
            positions.append((a[:,1:4]-self.origin)*1e6);speed.extend(np.linalg.norm(a[:,4:7],axis=1)*1000)
            lines.extend([len(a),*range(offset,offset+len(a))]);offset+=len(a)
        if not positions:return None
        mesh=pv.PolyData(np.vstack(positions),lines=np.array(lines));mesh['speed']=np.asarray(speed);return mesh


class View:
    def __init__(self,scene,detail=False,single=False,size=SIZE):
        self.scene=scene;self.detail=detail;self.single=single;self.size=size
        self.p=pv.Plotter(off_screen=True,window_size=size,shape=(1,2),col_weights=[.86,.14],border=False)
        p=self.p;p.set_background('black',all_renderers=True);p.subplot(0,0)
        p.renderer.SetViewport(0,.075,.86,.91);p.render_window.SetMultiSamples(0)
        p.enable_depth_peeling(number_of_peels=8,occlusion_ratio=0)
        wall=surface(scene.wall_um,scene.geo['WALL_faces'])
        p.add_mesh(wall,color='#abb7c5',opacity=.22,smooth_shading=True,
            show_scalar_bar=False,ambient=.7,diffuse=.3)
        for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
            cap=surface((scene.geo[role+'_points_m']-scene.origin)*1e6,scene.geo[role+'_faces'])
            p.add_mesh(cap,color='#b5c8d9',opacity=.13,show_scalar_bar=False)
        completed=[pid for pid,r in scene.records.items() if r['completed']]
        stopped=[pid for pid,r in scene.records.items() if not r['completed']]
        if single:
            mesh=scene.paths([scene.single]);p.add_mesh(mesh,color='#68819a',opacity=.20,line_width=1.,show_scalar_bar=False)
        else:
            mesh=scene.paths(completed)
            p.add_mesh(mesh,scalars='speed',cmap=CMAP,clim=(0,scene.limit),opacity=.65,line_width=1.4,show_scalar_bar=False,
                lighting=False)
            if stopped:
                mesh=scene.paths(stopped);p.add_mesh(mesh,color='#566371',opacity=.24,line_width=1.,show_scalar_bar=False)
        self.center=scene.branch if detail else np.zeros(3);self.scale=18. if detail else scene.scale
        if single:
            track=(scene.arrays[scene.single][:,1:4]-scene.origin)*1e6
            self.center=(track.min(0)+track.max(0))/2
            self.scale=max(12.,float(np.linalg.norm(np.ptp(track,axis=0)))*.65)
        self.dynamic=[];p.show(auto_close=False,interactive=False)

    def camera(self,angle):
        p=self.p;theta=np.deg2rad(angle);phi=np.deg2rad(22)
        direction=np.array([np.cos(theta)*np.cos(phi),np.sin(theta)*np.cos(phi),np.sin(phi)])
        p.camera.position=self.center+4.5*self.scene.radius*direction;p.camera.focal_point=self.center;p.camera.up=(0,0,1)
        p.camera.parallel_projection=True;p.camera.parallel_scale=self.scale;p.reset_camera_clipping_range()
        points=self.scene.wall_um-self.center;right=np.array([-np.sin(theta),np.cos(theta),0.]);up=np.cross(direction,right)
        aspect=self.size[0]*.86/(self.size[1]*.835)
        x=points@right/(self.scale*aspect);y=points@up/self.scale
        bounds=[float(x.min()),float(x.max()),float(y.min()),float(y.max())]
        if not self.detail and not self.single:assert max(map(abs,bounds))<.95,'Full-vessel clipping'
        return dict(azimuth_deg=float(angle),elevation_deg=22,parallel_scale_um=self.scale,
            position_um=list(p.camera.position),center_um=self.center.tolist(),projected_bounds=bounds,
            physical_coordinates_rotated=False,axis_stretching=False,
            view_role='SAVED_JOURNEY_CLOSEUP' if self.single else 'BRANCH_CLOSEUP' if self.detail else 'FULL_VESSEL')

    def frame(self,progress,angle,title):
        for actor in self.dynamic:self.p.remove_actor(actor,reset_camera=False,render=False)
        self.dynamic=[];states=self.scene.states(progress,self.single)
        if states:
            centers=(np.array([r['position_m'] for r in states])-self.scene.origin)*1e6
            points=pv.PolyData(centers);points['radius_um']=np.array([r['radius_m'] for r in states])*1e6
            spheres=points.glyph(orient=False,scale='radius_um',factor=1.,geom=pv.Sphere(radius=1.,theta_resolution=20,phi_resolution=16))
            self.p.add_mesh(spheres,color='#fff5d8',smooth_shading=True,opacity=1.,ambient=.3,diffuse=.7,specular=.55,
                specular_power=24,show_scalar_bar=False,name='bubbles',render=False,reset_camera=False);self.dynamic.append('bubbles')
        if self.single and progress>0:
            age=min(progress,1.)*self.scene.records[self.scene.single]['last_age_s']
            mesh=self.scene.paths([self.scene.single],through=age)
            if mesh is not None:
                self.p.add_mesh(mesh,scalars='speed',cmap=CMAP,clim=(0,self.scene.limit),line_width=3.,lighting=False,
                    show_scalar_bar=False,name='active_tail',render=False,reset_camera=False);self.dynamic.append('active_tail')
        cam=self.camera(angle);self.p.render();pixels=self.p.screenshot(return_img=True)
        if self.single:
            footer=f'Saved journey | age {progress*self.scene.records[self.scene.single]["last_age_s"]*1000:.1f} ms | true MB size | size-conditioned ensemble'
        else:footer=f'{len(self.scene.records)} saved tracks | independent phase overlay | true MB size | size-conditioned ensemble'
        image=annotate(pixels,title,footer,self.scene.limit)
        return image,dict(progress=float(progress),camera=cam,active_bubbles=states,
            physical_radius_scale=1.,source_scene_sha256=self.scene.sha,
            time_role='SINGLE_SAVED_TRAJECTORY_AGE' if self.single else 'DISPLAY_ALIGNED_INDEPENDENT_AGES_NOT_COMMON_ACQUISITION_TIME',
            no_extrapolation=True)

    def close(self):self.p.close()


def render(root,preview=False):
    root=Path(root);scene=Scene(root);start=time.time()
    for sub in ['figures','animations','frames','inspection']:(root/sub).mkdir(exist_ok=True)
    stills=[('00_full_vessel_microbubble_tracks',False,False,'Microbubbles and trajectories | Full vessel'),
        ('01_branch_microbubble_detail',True,False,'Microbubbles and trajectories | Branch detail'),
        ('03_single_microbubble_journey',False,True,'One microbubble | Saved inlet-to-outlet journey')]
    receipts={}
    for name,detail,single,title in stills:
        view=View(scene,detail,single,size=SIZE if preview else (3840,2160))
        im,receipt=view.frame(.45,45,title);im.save(root/'figures'/(name+'.png'));receipts[name]=receipt;view.close()
    views=[]
    for angle in [25,115,205,295]:
        view=View(scene,size=SIZE);im,receipt=view.frame(.45,angle,f'Microbubble trajectories | View {angle} deg')
        views.append(im);receipts['view_'+str(angle)]=receipt;view.close()
    board=Image.new('RGB',(3840,2160),'black')
    for k,im in enumerate(views):board.paste(im,((k%2)*1920,(k//2)*1080))
    board.save(root/'figures/02_multiview.png')
    atomic_json(root/'frames/STILL_RECEIPTS.json',receipts)
    atomic_json(root/'PPT_RENDER_MANIFEST.json',dict(scene_sha256=scene.sha,render_source_sha256=sha256(__file__),
        cohort_method='B',total_tracks=len(scene.records),background='BLACK',colorbar='RIGHT_OUTSIDE_3D_VIEWPORT',
        color_scale_mm_s=[0,scene.limit],colors='SAVED_MB_SPEED_ON_COMPLETED_PATHS; STOPPED_PATHS_GREY; TRUE_SIZE_BUBBLES_IVORY',
        moving_display_ids=scene.replay_ids,moving_selection='FIRST_48_COMPLETED_PLUS_UP_TO_12_OBSERVED_01_AND_03_COMPLETIONS',
        display_selection_does_not_change_500_cohort=True,single_journey_id=scene.single,
        single_selection='LONGEST_SAVED_COMPLETED_PATH_FOR_ILLUSTRATION',physical_radius_scale=1.,
        original_coordinates_preserved=True,extrapolation=False,render_host=socket.gethostname(),
        independent_phase_overlay_not_a_simultaneous_physical_population=True,preview_only=preview))
    if preview:return
    video_records=[]
    for kind,name in enumerate(ANIMATIONS):
        detail=kind==1;single=kind==2;view=View(scene,detail,single)
        title=['Microbubbles and trajectories | Full vessel','Microbubbles and trajectories | Branch detail','One microbubble | Saved inlet-to-outlet journey'][kind]
        path=root/'animations'/(name+'.mp4');temp=path.with_name(path.stem+'.tmp.mp4')
        writer=imageio_ffmpeg.write_frames(str(temp),SIZE,fps=FPS,codec='libx264',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',quality=8,
            macro_block_size=8,output_params=['-preset','medium','-movflags','+faststart']);writer.send(None);frames=[]
        try:
            for index in range(FRAMES):
                progress=index/(FRAMES-1);angle=45+(95 if single else 360)*index/FRAMES
                im,receipt=view.frame(progress,angle,title);receipt['frame']=index;frames.append(receipt)
                writer.send(np.asarray(im))
                if index%72==0:print('PPT_RENDER',name,index,'/',FRAMES,flush=True)
        finally:writer.close();view.close()
        os.replace(temp,path);atomic_json(root/'frames'/(name+'.json'),dict(frame_count=FRAMES,fps=FPS,size=SIZE,frames=frames))
        video_records.append(decode_video(path,FRAMES,root/'inspection'))
    for r in video_records:assert r['distinct_frames']>=FRAMES//2
    sheet=Image.new('RGB',(2880,1800),'black');d=ImageDraw.Draw(sheet)
    for row,r in enumerate(video_records):
        text(d,(25,row*600+12),r['file'].replace('_',' '),25)
        for col,p in enumerate(r['inspection_frames']):
            im=Image.open(p).convert('RGB');im.thumbnail((960,540));sheet.paste(im,(col*960,row*600+55))
    sheet.save(root/'figures/04_animation_storyboard.png')
    atomic_json(root/'PPT_MEDIA_VALIDATION.json',dict(all_pass=True,records=video_records,
        every_exported_frame_decoded=True,scene_sha256=scene.sha,render_seconds=time.time()-start,
        remote_hostname=socket.gethostname(),manual_visual_review='PENDING_USER_REVIEW'))
    print('PPT_VISUALS_COMPLETE',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--preview',action='store_true')
    a=p.parse_args();render(a.root,a.preview)
