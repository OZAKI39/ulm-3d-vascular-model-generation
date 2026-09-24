"""Scientific black-background views of the new, independently verified field.

Interior colors are P1 velocity evaluated on an inset sampling surface, never
velocity painted onto the no-slip wall. Camera orbits leave physical data intact.
"""
from pathlib import Path
import argparse, json, math, time, csv, hashlib, socket, platform, importlib.metadata
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from PIL import Image, ImageDraw, ImageFont
import imageio_ffmpeg
ROOT=Path(__file__).resolve().parents[2]
CASE=ROOT/'flow_cases/mean-2p0-mmps'
CMAP=ListedColormap(matplotlib.colormaps['turbo'](np.linspace(.10,.96,256)))
FPS=24; FRAMES=432; SIZE=(1920,1080)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')

def annotated(frame,title,surface_only=False,limit=None):
    """Keep annotations outside the 3D viewport, including branch close-ups."""
    im=Image.fromarray(np.asarray(frame)[:,:,:3]);draw=ImageDraw.Draw(im)
    width,height=im.size
    font_path='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
    title_font=ImageFont.truetype(font_path,max(20,int(height*.032)))
    foot_font=ImageFont.truetype(font_path,max(14,int(height*.022)))
    draw.rectangle((0,0,width,int(height*.083)),fill='black')
    draw.rectangle((0,int(height*.938),width,height),fill='black')
    draw.text((int(width*.025),int(height*.026)),title,font=title_font,fill='#f0f3fa')
    footer='Actual outer surface; no-slip wall: 0 mm/s' if surface_only else 'Interior velocity sampling surface; transparent no-slip wall'
    draw.text((int(width*.025),int(height*.958)),footer,font=foot_font,fill='#c0c8d3')
    if limit is not None:
        # Rasterize the legend outside the 3D renderer: VTK translucent passes
        # can interfere with text overlays. This uses the exact same color map.
        draw.rectangle((int(width*.86),int(height*.09),width,int(height*.935)),fill='black')
        x=int(width*.885);y=int(height*.26);bar_w=int(width*.019);bar_h=int(height*.50)
        colors=(CMAP(np.linspace(1,0,bar_h))[:,:3]*255).astype(np.uint8)
        gradient=np.repeat(colors[:,None,:],bar_w,axis=1)
        im.paste(Image.fromarray(gradient),(x,y))
        draw.text((int(width*.872),y-int(height*.058)),'Speed (mm/s)',font=foot_font,fill='#f0f3fa')
        for value in np.linspace(0,limit,6):
            yy=y+bar_h*(1-value/limit)
            draw.text((x+bar_w+10,int(yy-height*.012)),f'{value:.1f}',font=foot_font,fill='#f0f3fa')
    return np.asarray(im)

def camera(p,center,radius,angle,scale):
    theta=np.deg2rad(angle);phi=np.deg2rad(22)
    direction=np.array([np.cos(theta)*np.cos(phi),np.sin(theta)*np.cos(phi),np.sin(phi)])
    p.camera.position=center+4.5*radius*direction
    p.camera.focal_point=center
    p.camera.up=(0,0,1)
    p.camera.parallel_projection=True
    p.camera.parallel_scale=scale
    p.camera.clipping_range=(.01,1000)

def projection_bounds(points,center,angle,scale,aspect):
    theta=np.deg2rad(angle);phi=np.deg2rad(22)
    right=np.array([-np.sin(theta),np.cos(theta),0.])
    up=np.array([-np.cos(theta)*np.sin(phi),-np.sin(theta)*np.sin(phi),np.cos(phi)])
    q=points-center
    x=q@right/(scale*aspect);y=q@up/scale
    return [float(x.min()),float(x.max()),float(y.min()),float(y.max())]

class Scene:
    def __init__(self,field,limit):
        self.g=pv.read(field)
        self.original_points=np.array(self.g.points,copy=True)
        self.g.points=np.asarray(self.g.points,dtype=float)*1e6
        self.g['Speed']=np.linalg.norm(self.g['Velocity'],axis=1)*1000
        self.surface=self.g.extract_surface(algorithm='dataset_surface').compute_normals(
            auto_orient_normals=True,consistent_normals=True)
        distance=self.g.compute_implicit_distance(self.surface)
        assert distance['implicit_distance'].min()<-.35
        self.inner=distance.contour([-.35],scalars='implicit_distance')
        self.inner['Speed']=np.linalg.norm(self.inner['Velocity'],axis=1)*1000
        assert self.inner.n_points>1000 and np.isfinite(self.inner['Speed']).all()
        self.limit=limit
        self.center=np.array(self.g.center)
        self.radius=float(np.max(np.linalg.norm(self.g.points-self.center,axis=1)))
        q=self.g.points-self.center;radial=np.linalg.norm(q[:,:2],axis=1)
        vertical=np.abs(q[:,2])*np.cos(np.deg2rad(22))+radial*np.sin(np.deg2rad(22))
        self.full_scale=float(max(vertical.max(),radial.max()/1.37)/.90)
        self.branch=np.array(self.g.points[np.argmax(self.g['Speed'])])

    def plotter(self,title,detail=False,surface_only=False,size=SIZE,ports=False):
        p=pv.Plotter(off_screen=True,window_size=size,shape=(1,2),col_weights=[.86,.14],border=False)
        p.set_background('black',all_renderers=True)
        p.subplot(0,0);p.enable_depth_peeling(number_of_peels=8)
        p.renderer.SetViewport(0,.075,.86,.91)
        args=dict(cmap=CMAP,clim=(0,self.limit),smooth_shading=True,show_scalar_bar=False,
                  ambient=.55,diffuse=.45,specular=.10)
        if surface_only:
            actor=p.add_mesh(self.surface,scalars='Speed',**args)
        else:
            # The outer geometry is the actual wall/caps; their actual velocity is
            # also used here. Interior is a separate geometric sampling surface.
            p.add_mesh(self.surface,scalars='Speed',opacity=.13,**args)
            actor=p.add_mesh(self.inner,scalars='Speed',**args)
        font=max(13,int(size[1]/45))
        if ports:
            centers=[];labels=[]
            for name in ('INLET','OUTLET_01','OUTLET_02','OUTLET_03'):
                face=pv.read(CASE/'SV_MESH/mesh-surfaces'/f'{name}.vtp')
                centers.append(np.array(face.center)*1e6);labels.append(name.replace('OUTLET_','Out ').replace('INLET','Inlet'))
            p.add_point_labels(np.array(centers),labels,font_size=16,text_color='white',point_color='white',
                point_size=4,shape=None,show_points=True,always_visible=True)
        p.subplot(0,1)
        p.subplot(0,0)
        center=self.branch if detail else self.center
        scale=24. if detail else self.full_scale
        camera(p,center,self.radius,45,scale)
        return p,center,scale

def stats_figure(physics,path):
    m=physics['measurements'];comparison=physics['comparison']
    with plt.style.context('dark_background'):
        fig,axes=plt.subplots(1,2,figsize=(15,7),gridspec_kw={'width_ratios':[1.25,1]})
        fig.patch.set_facecolor('black')
        flows=[m['Q_in_m3_s'],*m['outlet_flows_m3_s'].values(),m['Q_out_total_m3_s']]
        labels=['Inlet\n(magnitude)','Outlet 01','Outlet 02','Outlet 03','Outlet sum']
        bars=axes[0].bar(labels,np.array(flows)*1e15,color=['#d9e6f6','#5389ef','#efba42','#58c9b2','#9cb2cc'],width=.62)
        axes[0].set_ylabel('Flow rate (pL/s)',fontsize=13)
        axes[0].set_ylim(0,max(flows)*1e15*1.22)
        for b,v in zip(bars,flows):axes[0].text(b.get_x()+b.get_width()/2,b.get_height()+.3,f'{v*1e15:.3f}',ha='center',fontsize=11)
        axes[0].spines[['top','right']].set_visible(False)
        axes[0].set_title('Measured boundary integrals',fontsize=16,pad=22)
        axes[1].axis('off')
        rows=[('Target inlet mean','2.000 mm/s'),('Measured inlet mean',f"{m['inlet_actual_mean_m_s']*1000:.6f} mm/s"),
              ('Speed range',f"{m['velocity_min_m_s']*1000:.3f} – {m['velocity_max_m_s']*1000:.3f} mm/s"),
              ('Volume mean speed',f"{m['velocity_volume_mean_m_s']*1000:.3f} mm/s"),
              ('Relative inflow error',f"{m['epsilon_Q']:.2e}"),('Relative mass imbalance',f"{m['epsilon_mass']:.2e}"),
              ('Maximum wall speed',f"{m['wall_velocity_max_m_s']*1000:.1f} mm/s"),
              ('Maximum-speed increase',f"{comparison['maximum_speed_ratio']:.3f} times")]
        axes[1].text(0,1,'Independent checks',transform=axes[1].transAxes,fontsize=16,va='top')
        for i,(k,v) in enumerate(rows):
            y=.86-i*.105;axes[1].text(0,y,k,color='#b9c5d5',fontsize=11);axes[1].text(.98,y,v,ha='right',fontsize=12)
        fig.suptitle('Inlet mean 2.0 mm/s | Flow verification',fontsize=20,y=.97)
        fig.subplots_adjust(left=.08,right=.96,bottom=.18,top=.81,wspace=.28)
        fig.savefig(path,dpi=160,facecolor='black');plt.close(fig)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stills-only',action='store_true');args=parser.parse_args()
    physics=json.loads((CASE/'reports/physics_validation.json').read_text());assert physics['status']=='PASS'
    field=CASE/'frozen_flow/steady_flow_mean_2p0_mmps.vtu';source_sha=sha(field)
    figures=CASE/'figures';animations=CASE/'animations';figures.mkdir(exist_ok=True);animations.mkdir(exist_ok=True)
    limit=math.ceil(physics['measurements']['velocity_max_m_s']*1000*2)/2
    scene=Scene(field,limit);start=time.time()
    capabilities=[]
    for name,title,surface in [('Figure_00_overview.png','Inlet mean 2.0 mm/s | Velocity field',False),
                               ('Figure_01_surface_speed.png','Inlet mean 2.0 mm/s | Surface speed',True)]:
        p,center,scale=scene.plotter(title,surface_only=surface,ports=surface)
        Image.fromarray(annotated(p.screenshot(),title,surface,limit)).save(figures/name)
        if not capabilities:
            capabilities=[line for line in p.ren_win.ReportCapabilities().splitlines() if any(k in line for k in ('OpenGL vendor','OpenGL renderer','OpenGL version'))]
        p.close()
    panels=[]
    for angle in (20,110,200,290):
        p,center,scale=scene.plotter(f'View {angle:03d} deg',size=(1280,960))
        camera(p,center,scene.radius,angle,scale);panels.append(Image.fromarray(annotated(p.screenshot(),f'View {angle:03d} deg',limit=limit)));p.close()
    board=Image.new('RGB',(2560,1920),'black')
    for i,im in enumerate(panels):board.paste(im,((i%2)*1280,(i//2)*960))
    board.save(figures/'Figure_02_multiview.png')
    stats_figure(physics,figures/'Figure_03_flow_validation.png')
    p,center,scale=scene.plotter('Branch detail | Inlet mean 2.0 mm/s',detail=True)
    Image.fromarray(annotated(p.screenshot(),'Branch detail | Inlet mean 2.0 mm/s',limit=limit)).save(figures/'Figure_05_branch_detail.png');p.close()
    manifest=dict(source_flow_sha256=source_sha,visual_units='mm/s',display_coordinates='micrometres; copy only',
        color_scale_mm_s=[0,limit],colormap='turbo truncated to [0.10,0.96]',background='black',colorbar='separate right viewport',
        inset_um=.35,internal_method='P1 vector interpolation on contour of nodal signed distance=-0.35 um, followed by vector norm; no wall-speed substitution',
        surface_method='norm of actual solver nodal velocity on exterior including inlet/outlet caps; wall no-slip remains zero',
        model_transforms_during_animation=False,view='orthographic camera orbit at 22 deg elevation',
        center_um=scene.center.tolist(),branch_focus_um=scene.branch.tolist(),branch_focus_rule='actual maximum-speed node, near proximal bifurcation',
        size=SIZE,fps=FPS,frames=FRAMES,duration_s=FRAMES/FPS,animations=[],
        render_host=socket.gethostname(),opengl=capabilities,renderer_script_sha256=sha(__file__))
    if not args.stills_only:
        for name,title,detail in [('Animation_01_full_vessel.mp4','Inlet mean 2.0 mm/s | Velocity field',False),
                                  ('Animation_02_branch_detail.mp4','Branch detail | Inlet mean 2.0 mm/s',True)]:
            p,center,scale=scene.plotter(title,detail=detail)
            trace=[]
            writer=imageio_ffmpeg.write_frames(str(animations/name),size=SIZE,fps=FPS,codec='libx264',quality=8,
                macro_block_size=2,output_params=['-movflags','+faststart'])
            writer.send(None)
            for i in range(FRAMES):
                angle=45+360*i/FRAMES;camera(p,center,scene.radius,angle,scale)
                # screenshot() can reuse the previous framebuffer after the first
                # show. Explicit rendering is required after moving the camera.
                p.render()
                bounds=projection_bounds(scene.g.points,center,angle,scale,SIZE[0]*.86/(SIZE[1]*.835))
                if not detail:assert max(abs(v) for v in bounds)<.94,'Full vessel would be clipped'
                writer.send(np.ascontiguousarray(annotated(p.screenshot(),title,limit=limit)))
                trace.append(dict(frame=i,time_s=i/FPS,azimuth_deg=angle,parallel_scale_um=scale,
                                  projected_bounds=bounds,camera_position_um=list(p.camera.position)))
                if i%72==0:print(name,i,'/',FRAMES,flush=True)
            writer.close();p.close()
            dump(animations/(Path(name).stem+'_camera.json'),trace)
            manifest['animations'].append(dict(file=name,sha256=sha(animations/name),detail=detail,frames=FRAMES,fps=FPS))
    manifest['render_wall_time_s']=time.time()-start
    manifest['source_unchanged']=sha(field)==source_sha;assert manifest['source_unchanged']
    manifest['figures']=[dict(file=p.name,sha256=sha(p)) for p in sorted(figures.glob('*.png'))]
    dump(CASE/'reports/render_manifest.json',manifest)
    dump(CASE/'reports/software_versions.json',dict(hostname=socket.gethostname(),python=platform.python_version(),
        packages={name:importlib.metadata.version(name) for name in ('numpy','scipy','pyvista','vtk','matplotlib','pillow','imageio-ffmpeg')},
        opengl=capabilities))
    print('Rendering complete',flush=True)

if __name__=='__main__':main()
