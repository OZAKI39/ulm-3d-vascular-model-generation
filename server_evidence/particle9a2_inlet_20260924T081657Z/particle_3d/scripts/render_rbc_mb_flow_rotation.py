"""Shared translation and hydrodynamic rotation, with body-fixed MB markers."""
import os,json,hashlib,time,argparse
os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
import numpy as np
import pyvista as pv
import vtk
from PIL import Image,ImageDraw
import imageio_ffmpeg
from prepare_rbc_mb_flow_rotation import *
from particle_3d.rbc_orientation import rotation_matrix
from render_normal_rbc_encounter import CMAP,text,RED,CYAN,GRAY,WHITE,SIZE,FPS

class View:
    def __init__(self):
        self.scene=json.loads((OUT/'data/SCENE.json').read_text());self.data=np.load(OUT/'data/motion.npz')
        self.records=self.scene['records'];self.p=pv.Plotter(off_screen=True,window_size=SIZE,shape=(1,2),col_weights=[.86,.14],border=False)
        p=self.p;p.set_background('black',all_renderers=True);p.render_window.SetMultiSamples(0)
        p.subplot(0,0);p.renderer.SetViewport(.015,.14,.855,.86);p.enable_depth_peeling(number_of_peels=8)
        wall=pv.Cylinder(center=(0,0,0),direction=(1,0,0),radius=RADIUS*1e6,height=40,resolution=192,capping=False)
        p.add_mesh(wall,color='#b8c7d7',opacity=.10,smooth_shading=True,ambient=.65,diffuse=.3,show_scalar_bar=False)
        pts=[];cells=[];speeds=[]
        for y in np.linspace(-5.5,5.5,9):
            z=-2.5;v=PoiseuilleField(RADIUS,MEAN).sample(np.array([0,y,z])*1e-6).velocity_m_s[0]*1000
            k=len(pts);pts.extend([[-20,y,z],[20,y,z]]);cells.extend([2,k,k+1]);speeds.extend([v,v])
        lines=pv.PolyData(np.array(pts),lines=np.array(cells));lines['Speed']=np.array(speeds)
        p.add_mesh(lines,scalars='Speed',cmap=CMAP,clim=(0,4),line_width=1,opacity=.22,lighting=False,show_scalar_bar=False)
        self.planes=vtk.vtkPlaneCollection()
        for x,n in [(-20,(1,0,0)),(20,(-1,0,0))]:
            plane=vtk.vtkPlane();plane.SetOrigin(x,0,0);plane.SetNormal(*n);self.planes.AddItem(plane)
        base=pv.read(OUT/'data/normal_rbc_si.vtp');base.points=np.asarray(base.points,dtype=float)*1e6
        sphere=pv.Sphere(radius=MB_RADIUS*1e6,theta_resolution=64,phi_resolution=48)
        normal=sphere.points/(MB_RADIUS*1e6)
        rgb=np.tile(np.array([103,227,246],dtype=np.uint8),(sphere.n_points,1))
        rgb[np.abs(normal[:,1])<.15]=[242,248,255]
        rgb[normal[:,0]>.86]=[255,206,88]
        sphere['Body_fixed_rotation_marks']=rgb
        self.meshes=[];self.refs=[];self.actors=[];self.trails=[];self.trail_actors=[]
        for r in self.records:
            is_rbc=r['species']=='RBC';mesh=(base if is_rbc else sphere).copy(deep=True);self.refs.append(mesh.points.copy());self.meshes.append(mesh)
            kwargs=dict(color=RED) if is_rbc else dict(scalars='Body_fixed_rotation_marks',rgb=True)
            actor=p.add_mesh(mesh,**kwargs,smooth_shading=True,ambient=.19 if is_rbc else .48,
                diffuse=.78 if is_rbc else .5,specular=.16 if is_rbc else .2,specular_power=30,show_scalar_bar=False)
            actor.mapper.SetClippingPlanes(self.planes);self.actors.append(actor)
            line=pv.PolyData(np.array([[0.,0,0],[1e-6,0,0]]),lines=[2,0,1]);self.trails.append(line)
            la=p.add_mesh(line,color=RED if is_rbc else CYAN,line_width=1.8,opacity=.65,lighting=False,show_scalar_bar=False)
            la.mapper.SetClippingPlanes(self.planes);self.trail_actors.append(la)
        p.remove_all_lights()
        p.add_light(pv.Light(position=(-20,-30,35),focal_point=(0,0,0),intensity=.95,light_type='scene light'))
        p.add_light(pv.Light(position=(25,10,15),focal_point=(0,0,0),intensity=.22,light_type='scene light'))
        p.camera.position=(25,-55,80);p.camera.focal_point=(0,0,0);p.camera.up=(0,.8,.55)
        p.camera.parallel_projection=True;p.camera.parallel_scale=10.5;p.camera.clipping_range=(.1,400.)
        p.show(auto_close=False,interactive=False)

    def frame(self,index):
        t=float(self.data['time_s'][index]);positions=self.data['positions_m'][index];counts={'RBC':0,'MB':0};visible=[]
        for j,(r,position) in enumerate(zip(self.records,positions)):
            radius=r['radius_m'];shown=position[0]+radius>=-20e-6 and position[0]-radius<=20e-6
            self.actors[j].SetVisibility(bool(shown))
            transform=np.eye(4);transform[:3,:3]=rotation_matrix(self.data['quaternion_wxyz'][index,j])
            transform[:3,3]=position*1e6
            matrix=vtk.vtkMatrix4x4();matrix.DeepCopy(transform.ravel())
            self.actors[j].SetUserMatrix(matrix)
            if shown:counts[r['species']]+=1;visible.append(r['id'])
            begin=max(0.,t-.002)
            initial=np.array(r['initial_m']);v=np.array(r['velocity_m_s'])
            trail=np.array([initial+begin*v,position])*1e6
            self.trails[j].points=trail
            self.trail_actors[j].SetVisibility(bool(t>0 and trail[:,0].max()>=-20 and trail[:,0].min()<=20))
        self.p.render();im=Image.fromarray(self.p.screenshot());d=ImageDraw.Draw(im)
        d.rectangle((0,0,1920,150),fill='black');d.rectangle((0,922,1920,1080),fill='black')
        text(d,(48,27),'RBC and microbubble flow | Translation and rotation',35)
        text(d,(48,83),'IDEALIZED VESSEL  |  Poiseuille velocity, shear and vorticity  |  Mean flow 2.0 mm/s',25,GRAY)
        text(d,(48,173),f'Physical time  {t*1000:.2f} ms',26)
        text(d,(48,218),f'In view:  {counts["RBC"]} RBCs   |   {counts["MB"]} microbubbles',24,GRAY)
        text(d,(50,943),'Normal RBCs',25,RED);text(d,(320,943),'Microbubbles',25,CYAN)
        text(d,(655,943),'Trails: last 2 ms',24,GRAY);text(d,(1195,943),'Flow direction →',24,GRAY)
        text(d,(48,990),'Flow-driven motion | RBC: Jeffery rotation | MB: local fluid rotation',23)
        text(d,(48,1033),'30 ms physical time → 18 s playback | MB surface marks make rotation visible',22,GRAY)
        d.rectangle((1655,145,1920,921),fill='black');x=1735;y=305;bw=32;bh=440
        rgb=(CMAP(np.linspace(1,0,bh))[:,:3]*255).astype('uint8');im.paste(Image.fromarray(np.repeat(rgb[:,None,:],bw,axis=1)),(x,y))
        text(d,(1680,226),'Flow speed',24);text(d,(1680,260),'(mm/s)',22)
        for v in np.linspace(0,4,5):text(d,(x+bw+13,int(y+bh*(1-v/4)-13)),f'{v:.1f}',22)
        return im,dict(frame=index,physical_time_s=t,video_time_s=index/FPS,visible_ids=visible,visible_counts=counts,
            quaternion_source_frame=index,orientation_sha256=hashlib.sha256(self.data['quaternion_wxyz'][index].tobytes()).hexdigest())

def main():
    for directory in ['figures','animations']:(OUT/directory).mkdir(parents=True,exist_ok=True)
    parser=argparse.ArgumentParser();parser.add_argument('--still',action='store_true');a=parser.parse_args();view=View()
    im,_=view.frame(144);im.save(OUT/'figures/RBC_MB_flow_rotation_overview.png')
    if a.still:view.p.close();return
    video=OUT/'animations/RBC_MB_flow_rotation.mp4';writer=imageio_ffmpeg.write_frames(str(video),SIZE,fps=FPS,codec='libx264',quality=8,
        pix_fmt_out='yuv420p',macro_block_size=8,output_params=['-preset','medium','-movflags','+faststart']);writer.send(None)
    records=[];start=time.time()
    try:
        for k in range(FRAMES):
            im,r=view.frame(k);writer.send(np.asarray(im));records.append(r)
            if k%72==0:print('ROTATING_COFLOW_RENDER',k,'/',FRAMES,flush=True)
    finally:writer.close();view.p.close()
    write_json(OUT/'data/replay_frames.json',records)
    frames=imageio_ffmpeg.read_frames(str(video));meta=next(frames);n=0;hashes=set();keys={}
    for k,raw in enumerate(frames):
        n+=1;hashes.add(hashlib.sha256(raw).hexdigest())
        if k in [0,144,288,431]:
            arr=np.frombuffer(raw,np.uint8).reshape(1080,1920,3).copy();keys[k]=arr
            Image.fromarray(arr).save(OUT/f'figures/decoded_{k:04d}.png')
    assert n==FRAMES and len(hashes)==FRAMES and meta['size']==SIZE and meta['fps']==FPS
    assert all(min(r['visible_counts'].values())>0 for r in records)
    for arr in keys.values():
        q=arr[270:920,:1650].astype(float)
        assert ((q[:,:,0]>1.5*q[:,:,1])&(q[:,:,0]>80)).sum()>2000
        assert ((q[:,:,1]>1.5*q[:,:,0])&(q[:,:,2]>120)).sum()>500
    board=Image.new('RGB',(3840,2160),'black')
    for j,k in enumerate([0,144,288,431]):board.paste(Image.fromarray(keys[k]),((j%2)*1920,(j//2)*1080))
    board.save(OUT/'figures/RBC_MB_flow_rotation_storyboard.png')
    write_json(OUT/'MEDIA_VALIDATION.json',dict(all_pass=True,decoded_frames=n,distinct_frames=len(hashes),fps=FPS,
        size=SIZE,video_duration_s=18.,physical_duration_s=DURATION,video_sha256=sha(video),
        body_fixed_MB_markers=True,RBC_rotation_from_saved_quaternions=True,
        all_frames_show_both_species=True,minimum_visible_counts={sp:min(r['visible_counts'][sp] for r in records) for sp in ['RBC','MB']},
        maximum_visible_counts={sp:max(r['visible_counts'][sp] for r in records) for sp in ['RBC','MB']},
        source_npz_sha256=sha(OUT/'data/motion.npz'),render_seconds=time.time()-start,background='black',
        camera='FIXED',pair_contact=False,hydrodynamic_pair_coupling=False))
    print('ROTATING_COFLOW_MEDIA_COMPLETE',flush=True)

if __name__=='__main__':main()
