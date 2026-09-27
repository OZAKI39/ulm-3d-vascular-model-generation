"""Replay saved, finite-size P4 encounter; match the black FEM presentation style."""
from pathlib import Path
import argparse,hashlib,json,os,sys,time
os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
import imageio_ffmpeg
from simulate_rbc_mb_interaction import ROOT,CASE,OUT,sha,initial_shapes
from particle_3d.particle3_cases import write_json
from particle_3d.pair_geometry import pair_gap
from particle_3d.rbc_mb_encounter import advance_encounter,capsule_sphere_gap
from rbc_render_style import CMAP,camera

SIZE=(1920,1080);FPS=24;FRAMES=432
RED='#e46059';CYAN='#6ae3f2';GRAY='#aebdce';WHITE='#eef3fa'

def txt(d,xy,s,n=23,c=WHITE):
    d.text(xy,s,font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',n),fill=c)

def capsule_mesh(shape):
    # Exact capsule surface sampled for rendering only; solver uses convex support.
    axis=np.asarray(shape.axis_world);e=np.eye(3)[np.argmin(abs(axis))]
    e=np.cross(axis,e);e/=np.linalg.norm(e);f=np.cross(axis,e)
    theta=np.linspace(0,np.pi,65);phi=np.linspace(0,2*np.pi,65)
    z=shape.radius_m*np.cos(theta)+np.where(theta<np.pi/2,1,-1)*shape.cylindrical_length_m/2
    rr=shape.radius_m*np.sin(theta)
    points=shape.center_m+z[:,None,None]*axis+rr[:,None,None]*(np.cos(phi)[None,:,None]*e+np.sin(phi)[None,:,None]*f)
    return pv.StructuredGrid(points[:,:,0]*1e6,points[:,:,1]*1e6,points[:,:,2]*1e6).extract_surface(algorithm='dataset_surface')

class Replay:
    def __init__(self):
        self.init=json.loads((OUT/'data/initialization.json').read_text())
        self.states=json.loads((OUT/'data/fine_states.json').read_text())
        self.times=np.array([r['time_s'] for r in self.states]);self.end=self.times[-1]
        self.centers=np.array([[p['center_m'] for p in r['particles']] for r in self.states])
        self.velocities=np.array([[p['velocity_m_s'] for p in r['particles']] for r in self.states])
        self.contact=np.array([bool(r['projection'].get('canonical_ids')) for r in self.states])
        self.shapes=initial_shapes(self.init)

    def at(self,t):
        if not self.times[0]<=t<=self.end:raise ValueError('No replay extrapolation')
        k=min(max(0,np.searchsorted(self.times,t,side='right')-1),len(self.times)-2)
        a=self.shapes[101].moved(self.centers[k,0]);b=self.shapes[203].moved(self.centers[k,1])
        rec=self.states[k+1]['projection']
        a,b,v,partial=advance_encounter(a,b,rec['cap_free_velocity_m_s'],rec['mb_free_velocity_m_s'],t-self.times[k])
        pos=np.array([a.center_m,b.center_m]);gap=capsule_sphere_gap(a,b)
        assert gap.gap_m>=-gap.roundoff_budget_m
        return pos,np.array([v[101],v[203]]),k,gap,partial

class View:
    def __init__(self,replay):
        self.r=replay;self.p=pv.Plotter(off_screen=True,window_size=SIZE,shape=(1,3),col_weights=[.30,.56,.14],border=False)
        p=self.p;p.set_background('black',all_renderers=True);p.render_window.SetMultiSamples(0)
        wall=pv.read(CASE/'SV_MESH/mesh-surfaces/WALL.vtp');wall.points=np.asarray(wall.points,dtype=float)*1e6
        lines=pv.read(CASE/'streamlines/data/streamlines_si.vtp');lines.points=np.asarray(lines.points,dtype=float)*1e6
        path=self.r.centers*1e6;self.focus=(path.min((0,1))+path.max((0,1)))/2
        p.subplot(0,0);p.renderer.SetViewport(.015,.32,.30,.86);p.enable_depth_peeling(number_of_peels=8)
        p.add_mesh(wall,color=GRAY,opacity=.18,smooth_shading=True,ambient=.65,show_scalar_bar=False)
        p.add_mesh(lines,scalars='Speed_mm_s',cmap=CMAP,clim=(0,7.5),line_width=.7,opacity=.38,lighting=False,show_scalar_bar=False)
        center=np.array(wall.center);radius=np.linalg.norm(wall.points-center,axis=1).max()
        # The narrow context viewport needs a larger orthographic scale.
        camera(p,center,radius,45,radius*1.15)
        p.add_mesh(pv.Sphere(radius=2.4,center=self.focus),color='#ffd57b',style='wireframe',line_width=1.3,opacity=.9)
        p.subplot(0,1);p.renderer.SetViewport(.32,.17,.86,.88);p.enable_depth_peeling(number_of_peels=8)
        lo=path.min((0,1))-8.;hi=path.max((0,1))+8.
        bounds=np.column_stack((lo,hi)).ravel()
        local=wall.clip_box(bounds,invert=False)
        p.add_mesh(local,color=GRAY,opacity=.17,smooth_shading=True,ambient=.7,diffuse=.3,show_scalar_bar=False)
        axis=np.array(replay.init['rbc']['axis_world']);normal=replay.at(0)[3].normal_j_to_i
        direction=np.cross(axis,normal);direction/=np.linalg.norm(direction)
        if direction[2]<0:direction=-direction
        direction=direction+.25*normal;direction/=np.linalg.norm(direction)
        p.camera.position=self.focus+80*direction;p.camera.focal_point=self.focus
        p.camera.up=(0,0,1);p.camera.parallel_projection=True;p.camera.parallel_scale=9.
        p.camera.clipping_range=(.1,250.)
        self.rbc_mesh=capsule_mesh(replay.shapes[101]);self.rbc_original=self.rbc_mesh.points.copy()
        self.mb_mesh=pv.Sphere(radius=replay.shapes[203].radius_m*1e6,center=replay.shapes[203].center_m*1e6,
                               theta_resolution=64,phi_resolution=48)
        self.mb_original=self.mb_mesh.points.copy()
        p.add_mesh(self.rbc_mesh,color=RED,smooth_shading=True,ambient=.40,diffuse=.65,specular=.30,specular_power=22,show_scalar_bar=False)
        p.add_mesh(self.mb_mesh,color=CYAN,smooth_shading=True,ambient=.4,diffuse=.55,specular=.9,specular_power=38,opacity=.95,show_scalar_bar=False)
        # Paths grow with saved physical time. No future or invented extension.
        self.trails=[]
        for j,color in [(0,RED),(1,CYAN)]:
            mesh=pv.PolyData(np.repeat(path[0,j][None,:],2,axis=0),lines=[2,0,1]);self.trails.append(mesh)
            p.add_mesh(mesh,color=color,line_width=3.,lighting=False,show_scalar_bar=False)
        p.show(auto_close=False,interactive=False)

    def frame(self,index):
        t=self.r.end*index/(FRAMES-1);pos,vel,k,gap,partial=self.r.at(t)
        self.rbc_mesh.points=self.rbc_original+(pos[0]-self.r.centers[0,0])*1e6
        self.mb_mesh.points=self.mb_original+(pos[1]-self.r.centers[0,1])*1e6
        for j,mesh in enumerate(self.trails):
            pts=np.vstack((self.r.centers[:k+1,j],pos[j]))*1e6
            mesh.copy_from(pv.PolyData(pts,lines=np.r_[len(pts),np.arange(len(pts))]))
        self.p.render();im=Image.fromarray(self.p.screenshot());d=ImageDraw.Draw(im)
        d.rectangle((0,0,1920,125),fill='black');d.rectangle((0,958,1920,1080),fill='black')
        txt(d,(48,27),'Microbubble–RBC interaction | Inlet mean 2.0 mm/s',35)
        txt(d,(48,83),'Local two-body encounter in the actual FEM vessel',23,GRAY)
        txt(d,(48,165),'FULL-VESSEL LOCATION',22,GRAY)
        txt(d,(635,140),'INTERACTION DETAIL',22,GRAY)
        txt(d,(65,790),f'Physical time   {t*1000:.3f} ms',26)
        txt(d,(65,837),f'Surface gap     {max(0,gap.gap_m)*1e9:.1f} nm',24)
        contact=partial['contact_duration_s']>0
        phase='Contact / deflection' if contact else ('Approach' if not self.r.contact[:k+1].any() else 'After contact')
        txt(d,(65,883),phase,24,'#ffd57b')
        txt(d,(635,918),'RBC capsule',23,RED);txt(d,(905,918),'Microbubble',23,CYAN)
        txt(d,(1195,918),'Trails: saved paths',23,GRAY)
        txt(d,(48,982),'P4 frictionless contact | Fixed, volume-preserving RBC capsule | One-way frozen flow',24)
        txt(d,(48,1024),f'Slow motion: {self.r.end*1000:.2f} ms physical time shown over 18 s | Original particle sizes',22,GRAY)
        # Same color map and range as the new flow figures; applies to flow lines.
        d.rectangle((1660,125,1920,950),fill='black');x=1732;y=330;bh=465;bw=32
        rgb=(CMAP(np.linspace(1,0,bh))[:,:3]*255).astype('uint8');im.paste(Image.fromarray(np.repeat(rgb[:,None,:],bw,axis=1)),(x,y))
        txt(d,(1685,220),'Context flow',22,GRAY)
        txt(d,(1685,253),'Flow speed',23);txt(d,(1685,286),'(mm/s)',22)
        for v in np.linspace(0,7.5,6):txt(d,(x+bw+14,int(y+bh*(1-v/7.5)-14)),f'{v:.1f}',22)
        return im,dict(frame=index,video_time_s=index/FPS,physical_time_s=t,accepted_interval=k,
            positions_m=pos.tolist(),pair_gap_m=gap.gap_m,phase=phase,contact_active=contact)

    def close(self):self.p.close()

def decode(path):
    frames=imageio_ffmpeg.read_frames(str(path));meta=next(frames);unique=set();keys={};count=0
    for k,raw in enumerate(frames):
        unique.add(hashlib.sha256(raw).hexdigest());count+=1
        if k in [0,144,288,431]:keys[k]=np.frombuffer(raw,np.uint8).reshape(1080,1920,3).copy()
    assert count==FRAMES and len(unique)>FRAMES*.90 and meta['size']==SIZE and meta['fps']==FPS
    for k,arr in keys.items():
        Image.fromarray(arr).save(OUT/f'figures/decoded_frame_{k:04d}.png')
        region=arr[170:945,620:1640].astype(int)
        assert ((region[:,:,0]>region[:,:,1]*1.3)&(region[:,:,0]>90)).sum()>1000
        assert ((region[:,:,1]>region[:,:,0]*1.4)&(region[:,:,2]>100)).sum()>200
    return dict(all_pass=True,frames=count,distinct_frames=len(unique),fps=FPS,width=1920,height=1080,duration_s=18.,sha256=sha(path)),keys

def main():
    p=argparse.ArgumentParser();p.add_argument('--still',action='store_true');args=p.parse_args()
    replay=Replay();view=View(replay)
    image,_=view.frame(288);image.save(OUT/'figures/Interaction_overview.png')
    if args.still:view.close();return
    video=OUT/'animations/Microbubble_RBC_interaction.mp4';writer=imageio_ffmpeg.write_frames(str(video),SIZE,fps=FPS,
        codec='libx264',quality=8,pix_fmt_out='yuv420p',macro_block_size=8,output_params=['-preset','medium','-movflags','+faststart'])
    writer.send(None);records=[];begin=time.time()
    try:
        for k in range(FRAMES):
            im,r=view.frame(k);writer.send(np.asarray(im));records.append(r)
            if k%72==0:print('RBC_MB_RENDER',k,'/',FRAMES,flush=True)
    finally:writer.close();view.close()
    write_json(OUT/'data/replay_frames.json',records);validation,keys=decode(video)
    board=Image.new('RGB',(2880,1080),'black');d=ImageDraw.Draw(board)
    for col,k in enumerate([0,144,431]):
        # Focus on the actual decoded local view, preserving image aspect ratio.
        crop=Image.fromarray(keys[k]).crop((600,125,1645,953));crop.thumbnail((940,850))
        board.paste(crop,(col*960+10,100));txt(d,(col*960+35,35),f't = {records[k]["physical_time_s"]*1000:.3f} ms',31)
        txt(d,(col*960+35,900),records[k]['phase'],28,'#ffd57b')
    txt(d,(35,993),'Microbubble–RBC encounter | P4 kinematic contact | Fixed capsule surrogate',30)
    board.save(OUT/'figures/Interaction_storyboard.png')
    validation.update(physical_duration_s=replay.end,render_elapsed_s=time.time()-begin,source_states_sha256=sha(OUT/'data/fine_states.json'),
        original_particle_sizes=True,background='black',flow_speed_range_mm_s=[0,7.5],full_vessel_and_local_detail=True,
        replay_method='EXACT_PARTIAL_REPLAY_OF_ACCEPTED_P4_CONTACT_INTERVALS_NO_EXTRAPOLATION',fixed_camera=True)
    write_json(OUT/'MEDIA_VALIDATION.json',validation);print('RBC_MB_MEDIA_COMPLETE',flush=True)

if __name__=='__main__':main()
