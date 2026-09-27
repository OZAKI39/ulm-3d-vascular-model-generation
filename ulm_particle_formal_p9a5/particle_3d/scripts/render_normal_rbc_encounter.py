"""Black scientific presentation of a normal biconcave RBC / MB encounter."""
from pathlib import Path
import argparse,json,hashlib,os,sys,time
os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
import imageio_ffmpeg
from simulate_normal_rbc_encounter import *
import matplotlib
from matplotlib.colors import ListedColormap
CMAP=ListedColormap(matplotlib.colormaps['turbo'](np.linspace(.10,.96,256)))
SIZE=(1920,1080);FRAMES=432;FPS=24;RED='#e34c50';CYAN='#67e3f6';WHITE='#eef3fa';GRAY='#b8c6d9'

def text(d,xy,value,size=24,color=WHITE):
    d.text(xy,value,font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',size),fill=color)

class Replay:
    def __init__(self):
        d=np.load(OUT/'data/fine.npz');self.times=d['times_s'];self.positions=d['positions_m']
        self.intervals=json.loads((OUT/'data/fine_intervals.json').read_text());self.end=self.times[-1]
        self.summary=json.loads((OUT/'data/fine_summary.json').read_text())

    def at(self,t):
        if not 0<=t<=self.end:raise ValueError('Replay extrapolation forbidden')
        k=min(int(np.searchsorted(self.times,t,side='right'))-1,len(self.intervals)-1)
        rec=self.intervals[k];x,y=self.positions[k];vx,vy=rec['free_velocities_m_s']
        a,b,v,partial=advance_equatorial(x,y,vx,vy,t-self.times[k]);contact=partial['contact_duration_s']>0
        phase='Approach' if t<self.summary['first_contact_time_s'] else 'Contact and sliding' if t<=self.summary['last_contact_time_s'] else 'Separation'
        return np.array([a,b]),k,partial,phase

class View:
    def __init__(self,replay):
        self.r=replay;self.p=pv.Plotter(off_screen=True,window_size=SIZE,shape=(1,2),col_weights=[.86,.14],border=False)
        p=self.p;p.set_background('black',all_renderers=True);p.render_window.SetMultiSamples(0)
        p.subplot(0,0);p.renderer.SetViewport(.015,.14,.855,.86);p.enable_depth_peeling(number_of_peels=8)
        wall=pv.Cylinder(center=(0,0,0),direction=(1,0,0),radius=TUBE_RADIUS*1e6,height=40,resolution=192,capping=False)
        self.wall=wall;self.wall_ref=wall.points.copy()
        p.add_mesh(wall,color='#b8c7d7',opacity=.10,smooth_shading=True,ambient=.65,diffuse=.3,show_scalar_bar=False)
        # Back-side reference streamlines only. These never intersect the bodies.
        pts=[];cells=[];speeds=[]
        for y in np.linspace(-5.5,5.5,9):
            z=-2.5;u=poiseuille(np.array([0,y,z])*1e-6,TUBE_RADIUS)[0]*1000
            i=len(pts);pts.extend([[-20,y,z],[20,y,z]]);cells.extend([2,i,i+1]);speeds.extend([u,u])
        lines=pv.PolyData(np.array(pts),lines=np.array(cells));lines['Speed']=np.array(speeds)
        self.lines=lines;self.lines_ref=lines.points.copy()
        p.add_mesh(lines,scalars='Speed',cmap=CMAP,clim=(0,4),line_width=1.1,opacity=.35,lighting=False,show_scalar_bar=False)
        self.rbc=pv.read(OUT/'data/normal_rbc_si.vtp');self.rbc.points=np.asarray(self.rbc.points,dtype=float)*1e6
        self.rbc_ref=self.rbc.points.copy()
        self.bubble=pv.Sphere(radius=MB_RADIUS*1e6,theta_resolution=64,phi_resolution=64)
        self.bubble_ref=self.bubble.points.copy()
        p.add_mesh(self.rbc,color=RED,smooth_shading=True,ambient=.14,diffuse=.82,specular=.18,specular_power=30,show_scalar_bar=False)
        p.add_mesh(self.bubble,color=CYAN,smooth_shading=True,ambient=.32,diffuse=.58,specular=.9,specular_power=40,show_scalar_bar=False)
        self.trails=[]
        for j,color in [(0,RED),(1,CYAN)]:
            mesh=pv.PolyData(np.repeat(replay.positions[0,j][None,:]*1e6,2,axis=0),lines=[2,0,1]);self.trails.append(mesh)
            p.add_mesh(mesh,color=color,line_width=2.6,lighting=False,opacity=.85,show_scalar_bar=False)
        p.remove_all_lights()
        self.key=pv.Light(position=(-20,-30,35),focal_point=(0,0,0),intensity=.95,light_type='scene light')
        self.fill=pv.Light(position=(25,10,15),focal_point=(0,0,0),intensity=.22,light_type='scene light')
        p.add_light(self.key);p.add_light(self.fill)
        p.camera.position=(25,-55,80);p.camera.focal_point=(0,0,0);p.camera.up=(0,.80,.55)
        p.camera.parallel_projection=True;p.camera.parallel_scale=10.5;p.camera.clipping_range=(.1,400.)
        p.show(auto_close=False,interactive=False)

    def frame(self,index):
        t=self.r.end*index/(FRAMES-1);pos,k,partial,phase=self.r.at(t)
        self.rbc.points=self.rbc_ref+pos[0]*1e6;self.bubble.points=self.bubble_ref+pos[1]*1e6
        for j,mesh in enumerate(self.trails):
            start=max(0,int(np.searchsorted(self.r.times,max(0,t-.003))))
            pts=np.vstack([self.r.positions[start:k+1,j],pos[j]])*1e6
            if len(pts)<2:pts=np.repeat(pts,2,axis=0)
            mesh.copy_from(pv.PolyData(pts,lines=np.r_[len(pts),np.arange(len(pts))]))
        center=np.array([pos[:,0].mean()*1e6,0.,0.]);self.p.camera.position=center+[25,-55,80];self.p.camera.focal_point=center
        self.wall.points=self.wall_ref+center;self.lines.points=self.lines_ref+center
        self.key.position=center+[-20,-30,35];self.key.focal_point=center
        self.fill.position=center+[25,10,15];self.fill.focal_point=center
        self.p.render();im=Image.fromarray(self.p.screenshot());d=ImageDraw.Draw(im)
        d.rectangle((0,0,1920,150),fill='black');d.rectangle((0,922,1920,1080),fill='black')
        text(d,(48,27),'Microbubble–RBC interaction | Normal biconcave RBC',35)
        text(d,(48,83),'IDEALIZED VESSEL  |  Diameter 14 µm  |  Mean flow 2.0 mm/s',25,GRAY)
        text(d,(48,173),f'Physical time  {t*1000:.2f} ms',26)
        text(d,(48,218),phase,26,'#ffd581')
        text(d,(1180,173),f'Surface gap  {max(0,partial["gap_m"])*1e6:.2f} µm',23)
        text(d,(50,943),'Normal RBC',25,RED);text(d,(305,943),'Microbubble',25,CYAN)
        text(d,(605,943),'Trails: last 3 ms',24,GRAY);text(d,(1155,943),'Flow direction  →',24,GRAY)
        text(d,(48,990),'Normal shape retained | Frictionless rim contact | Analytical background flow',24)
        text(d,(48,1033),'30 ms physical time → 18 s playback | Camera follows the pair | Original particle sizes',22,GRAY)
        d.rectangle((1655,145,1920,921),fill='black');x=1735;y=305;bw=32;bh=440
        rgb=(CMAP(np.linspace(1,0,bh))[:,:3]*255).astype('uint8');im.paste(Image.fromarray(np.repeat(rgb[:,None,:],bw,axis=1)),(x,y))
        text(d,(1680,226),'Flow speed',24);text(d,(1680,260),'(mm/s)',22)
        for v in np.linspace(0,4,5):text(d,(x+bw+13,int(y+bh*(1-v/4)-13)),f'{v:.1f}',22)
        return im,dict(frame=index,video_time_s=index/FPS,physical_time_s=t,phase=phase,positions_m=pos.tolist(),
            source_interval=k,surface_gap_m=partial['gap_m'],camera_focal_point_um=center.tolist(),rbc_deformed=False)

    def close(self):self.p.close()

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--still',action='store_true');args=parser.parse_args()
    replay=Replay();view=View(replay)
    im,_=view.frame(100);im.save(OUT/'figures/Normal_RBC_interaction_overview.png')
    if args.still:view.close();return
    path=OUT/'animations/Normal_RBC_microbubble_interaction.mp4';writer=imageio_ffmpeg.write_frames(str(path),SIZE,fps=FPS,
        codec='libx264',quality=8,pix_fmt_out='yuv420p',macro_block_size=8,output_params=['-preset','medium','-movflags','+faststart'])
    writer.send(None);records=[];start=time.time()
    try:
        for k in range(FRAMES):
            im,rec=view.frame(k);writer.send(np.asarray(im));records.append(rec)
            if k%72==0:print('NORMAL_RBC_RENDER',k,'/',FRAMES,flush=True)
    finally:writer.close();view.close()
    write_json(OUT/'data/replay_frames.json',records)
    reader=imageio_ffmpeg.read_frames(str(path));meta=next(reader);keys={};unique=set();n=0
    for k,raw in enumerate(reader):
        n+=1;unique.add(hashlib.sha256(raw).hexdigest())
        if k in [0,100,220,431]:
            arr=np.frombuffer(raw,np.uint8).reshape(1080,1920,3).copy();keys[k]=arr
            Image.fromarray(arr).save(OUT/f'figures/decoded_{k:04d}.png')
    assert n==FRAMES and len(unique)==FRAMES and meta['size']==SIZE and meta['fps']==FPS
    # Both particles must remain visible in every inspected phase.
    for arr in keys.values():
        a=arr[250:900,:1650].astype(float)
        assert ((a[:,:,0]>1.5*a[:,:,1])&(a[:,:,0]>80)).sum()>1500
        assert ((a[:,:,1]>1.5*a[:,:,0])&(a[:,:,2]>120)).sum()>250
    board=Image.new('RGB',(3840,2160),'black')
    for j,k in enumerate([0,100,220,431]):board.paste(Image.fromarray(keys[k]),((j%2)*1920,(j//2)*1080))
    board.save(OUT/'figures/Normal_RBC_interaction_storyboard.png')
    record=dict(all_pass=True,frames=n,distinct_frames=len(unique),fps=FPS,size=SIZE,video_duration_s=18.,physical_duration_s=float(replay.end),
        source_trajectory_sha256=sha(OUT/'data/fine.npz'),source_intervals_sha256=sha(OUT/'data/fine_intervals.json'),
        video_sha256=sha(path),geometry_sha256=sha(OUT/'data/normal_rbc_si.vtp'),normal_shape_retained_in_all_frames=True,
        scene_type='EXPLICITLY_LABELLED_IDEALIZED_VESSEL',camera='TRANSLATES_WITH_PAIR_WITH_FIXED_VIEW_DIRECTION',
        field='ANALYTICAL_POISEUILLE_NOT_ORIGINAL_FEM',render_seconds=time.time()-start,
        replay='EXACT_PARTIAL_CONTACT_SOLUTION_NO_LINEAR_CHORD_THROUGH_RBC',speed_range_mm_s=[0,4],background='black')
    write_json(OUT/'MEDIA_VALIDATION.json',record);print(json.dumps(record,indent=2),flush=True)

if __name__=='__main__':main()
