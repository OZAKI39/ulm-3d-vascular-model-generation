"""Fixed-camera, fixed-scale rigid turntables about one vessel longitudinal axis."""
import os
os.environ.setdefault('LP_NUM_THREADS', '2')
os.environ.setdefault('VTK_SMP_MAX_THREADS', '2')
from pathlib import Path
import argparse
import hashlib
import json
import time
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import imageio_ffmpeg
from compute_field_diagnostics import CASE, OUT, dump, sha
from render import CMAP

SIZE = (1920, 1080)
FPS, FRAMES = 24, 432
VIEWPORT = (0., .08, .86, .915)
ASPECT = SIZE[0]*(VIEWPORT[2]-VIEWPORT[0])/(SIZE[1]*(VIEWPORT[3]-VIEWPORT[1]))
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
SPECS = {
    'pressure': dict(title='Pressure field', scalar='Pressure_Pa', unit='Pressure (Pa)',
        clim=(-5., 2500.), ticks=[0, 500, 1000, 1500, 2000, 2500],
        footer='Steady FEM pressure on the vessel surface | Inlet mean velocity: 2.0 mm/s',
        file='pressure_surface_si.vtp'),
    'wss': dict(title='Wall shear stress', scalar='WSS_display_Pa', unit='WSS (Pa)',
        clim=(0., 55.), ticks=[0, 10, 20, 30, 40, 50, 55],
        footer='Derived from FEM velocity gradients | Area-weighted nodal display | Physical wall only',
        file='wall_wss_si.vtp'),
    'streamlines': dict(title='FEM streamlines', scalar='Speed_mm_s', unit='Speed (mm/s)',
        clim=(0., 7.5), ticks=[0, 1.5, 3, 4.5, 6, 7.5],
        footer='96 saved streamlines | All three outlets covered | Line counts do not represent flow fractions',
        file=None),
}


def axis_rotation(axis, angle_degrees):
    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    skew = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
    angle = np.deg2rad(angle_degrees)
    return np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(axis,axis)+np.sin(angle)*skew


class FixedAxisRotation:
    def __init__(self, points_um):
        self.points = np.asarray(points_um)
        self.center = (self.points.min(0)+self.points.max(0))/2
        # Compute the principal geometric axis once, never separately per frame.
        _, _, vt = np.linalg.svd(self.points-self.points.mean(0), full_matrices=False)
        self.axis = vt[0]
        if self.axis[2] < 0:self.axis *= -1
        heights = (self.points-self.center)@self.axis
        self.center += self.axis*(heights.min()+heights.max())/2
        self.radius = float(np.linalg.norm(self.points-self.center, axis=1).max())
        self.angles = np.arange(FRAMES)*360/FRAMES
        direction = np.array([1.,1.,0.])
        direction -= np.dot(direction,self.axis)*self.axis
        if np.linalg.norm(direction)<1e-8:
            direction=np.array([0.,0.,1.])-self.axis[2]*self.axis
        self.direction=direction/np.linalg.norm(direction)
        self.right=np.cross(self.axis,self.direction)
        self.up=self.axis.copy()
        self.camera_position=self.center+4.5*self.radius*self.direction
        q = self.points-self.center
        height=q@self.axis
        radial=np.linalg.norm(q-height[:,None]*self.axis,axis=1)
        # Continuous 360-degree bound: one constant scale fits all rotations.
        self.scale=float(max(np.max(np.abs(height)),radial.max()/ASPECT)/.945)
        self.scales=np.full(FRAMES,self.scale)
        self.centers=np.tile(self.center,(FRAMES,1))
        self.transforms=[]
        self.trace=[]
        for index in range(FRAMES):
            rotation=axis_rotation(self.axis,self.angles[index])
            matrix=np.eye(4);matrix[:3,:3]=rotation;matrix[:3,3]=self.center-rotation@self.center
            self.transforms.append(matrix)
            transformed=q@rotation.T
            x=transformed@self.right/(self.scale*ASPECT);y=transformed@self.up/self.scale
            bounds=[float(x.min()),float(x.max()),float(y.min()),float(y.max())]
            assert max(abs(v) for v in bounds)<.946
            assert np.allclose(rotation.T@rotation,np.eye(3),atol=1e-14)
            assert np.allclose(rotation@self.axis,self.axis,atol=1e-14)
            self.trace.append(dict(frame=index,video_time_s=index/FPS,rotation_angle_deg=float(self.angles[index]),
                rotation_axis_unit=self.axis.tolist(),rotation_center_um=self.center.tolist(),
                camera_position_um=self.camera_position.tolist(),camera_up_unit=self.up.tolist(),
                focal_point_um=self.center.tolist(),parallel_scale_um=self.scale,display_matrix=matrix.tolist(),
                projected_bounds=bounds,physical_data_modified=False,rigid_display_transform_only=True))

    def setup_camera(self, plotter):
        plotter.camera.position=self.camera_position
        plotter.camera.focal_point=self.center
        plotter.camera.up=self.up
        plotter.camera.parallel_projection=True
        plotter.camera.parallel_scale=self.scale
        plotter.camera.clipping_range=(.01,1000)

    def apply(self, actors, index):
        for actor in actors:actor.user_matrix=self.transforms[index]


def annotate(frame, kind):
    spec=SPECS[kind]
    im=Image.fromarray(np.asarray(frame)[:,:,:3]);w,h=im.size;d=ImageDraw.Draw(im)
    titlefont=ImageFont.truetype(FONT,round(h*.032))
    font=ImageFont.truetype(FONT,round(h*.021))
    small=ImageFont.truetype(FONT,round(h*.018))
    d.rectangle((0,0,w,h*.077),fill='black')
    d.rectangle((0,h*.935,w,h),fill='black')
    d.text((w*.025,h*.025),'Full vessel | '+spec['title'],font=titlefont,fill='#f0f3fa')
    d.text((w*.025,h*.958),spec['footer'],font=font,fill='#c0c8d3')
    d.rectangle((w*.86,h*.08,w,h*.935),fill='black')
    x,y,bw,bh=round(w*.883),round(h*.275),round(w*.018),round(h*.49)
    colors=(CMAP(np.linspace(1,0,bh))[:,:3]*255).astype(np.uint8)
    im.paste(Image.fromarray(np.repeat(colors[:,None,:],bw,axis=1)),(x,y))
    d.text((w*.872,h*.208),spec['unit'],font=font,fill='#f0f3fa')
    lo,hi=spec['clim']
    for value in spec['ticks']:
        yy=y+bh*(hi-value)/(hi-lo)
        d.line((x+bw,yy,x+bw+w*.005,yy),fill='#e9eef4',width=max(1,round(w/1920)))
        label=f'{value:g}'
        d.text((x+bw+w*.010,yy-h*.012),label,font=font,fill='#f0f3fa')
    d.text((w*.872,h*.84),'Steady field',font=small,fill='#a7b4c5')
    d.text((w*.872,h*.87),'Fixed-axis rotation',font=small,fill='#a7b4c5')
    return im


class View:
    def __init__(self,kind,orbit,size=SIZE):
        self.kind,self.orbit=kind,orbit
        self.actors=[]
        spec=SPECS[kind]
        self.p=pv.Plotter(off_screen=True,window_size=size,shape=(1,2),col_weights=[.86,.14],border=False)
        p=self.p;p.set_background('black',all_renderers=True);p.subplot(0,0)
        p.renderer.SetViewport(*VIEWPORT)
        p.enable_depth_peeling(number_of_peels=8,occlusion_ratio=0)
        p.render_window.SetMultiSamples(0)
        if kind=='streamlines':
            wall=pv.read(OUT/'data/wall_wss_si.vtp');wall.points*=1e6
            self.actors.append(p.add_mesh(wall,color='#b7c4d1',opacity=.20,smooth_shading=True,ambient=.65,diffuse=.35,show_scalar_bar=False))
            mesh=pv.read(CASE/'streamlines/data/streamlines_si.vtp');mesh.points*=1e6
            assert mesh.n_cells==96
            self.actors.append(p.add_mesh(mesh,scalars=spec['scalar'],cmap=CMAP,clim=spec['clim'],line_width=1.5*size[1]/1080,
                opacity=.9,lighting=False,show_scalar_bar=False))
        else:
            mesh=pv.read(OUT/'data'/spec['file']);mesh.points*=1e6
            assert mesh[spec['scalar']].min()>=spec['clim'][0] and mesh[spec['scalar']].max()<=spec['clim'][1]
            self.actors.append(p.add_mesh(mesh,scalars=spec['scalar'],cmap=CMAP,clim=spec['clim'],preference='point',
                smooth_shading=True,split_sharp_edges=False,ambient=.68,diffuse=.32,specular=.08,
                show_scalar_bar=False))
        orbit.setup_camera(p);orbit.apply(self.actors,0);p.show(auto_close=False,interactive=False)

    def frame(self,index):
        self.orbit.apply(self.actors,index);self.p.render()
        assert np.allclose(self.p.camera.position,self.orbit.camera_position,atol=1e-10)
        assert np.allclose(self.p.camera.focal_point,self.orbit.center,atol=1e-10)
        assert self.p.camera.parallel_scale==self.orbit.scale
        return annotate(self.p.screenshot(),self.kind)

    def close(self):
        self.p.close()


def validate_video(path,kind):
    reader=imageio_ffmpeg.read_frames(str(path));meta=next(reader)
    assert meta['size']==SIZE and meta['fps']==FPS and abs(meta['duration']-18)<.05
    distinct=set();keys={};count=0
    for i,raw in enumerate(reader):
        count+=1;distinct.add(hashlib.sha256(raw).hexdigest())
        if i in [0,108,216,324,431]:
            a=np.frombuffer(raw,np.uint8).reshape(SIZE[1],SIZE[0],3).copy();keys[i]=a
            Image.fromarray(a).save(OUT/'inspection'/f'{kind}_{i:04d}.png')
    assert count==FRAMES and len(distinct)==FRAMES
    for a in keys.values():
        body=a[100:950,:1640].astype(float)
        assert np.sum(body.max(2)-body.min(2)>25)>10000
        assert np.mean(a[5:70,1700:1900])<3
    sheet=Image.new('RGB',(3840,2160),'black')
    for k,i in enumerate([0,108,216,324]):
        sheet.paste(Image.fromarray(keys[i]),((k%2)*1920,(k//2)*1080))
    sheet.save(OUT/'figures'/f'{kind}_rotation_views.png')
    return dict(file=str(path.relative_to(OUT)),sha256=sha(path),decoded_frames=count,distinct_frames=len(distinct),
                duration_s=18.,fps=FPS,size=SIZE,all_pass=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stills-only',action='store_true')
    args=parser.parse_args()
    start=time.time()
    compute=json.loads((OUT/'COMPUTE_VALIDATION.json').read_text());assert compute['all_pass']
    for name in ['data/wall_wss_si.vtp','data/pressure_surface_si.vtp']:
        assert sha(OUT/name)==compute['outputs_sha256'][name]
    mesh=pv.read(OUT/'data/pressure_surface_si.vtp')
    orbit=FixedAxisRotation(mesh.points*1e6)
    dump(OUT/'CAMERA_TRACE.json',orbit.trace)
    videos=[];capabilities=[]
    for kind in SPECS:
        view=View(kind,orbit,size=(3840,2160))
        view.frame(0).save(OUT/'figures'/f'{kind}_overview_4k.png')
        if not capabilities:
            capabilities=[l for l in view.p.render_window.ReportCapabilities().splitlines()
                if any(k in l for k in ['OpenGL vendor','OpenGL renderer','OpenGL version'])]
        view.close()
        if args.stills_only:
            view=View(kind,orbit)
            for index in [108,216,324]:view.frame(index).save(OUT/'inspection'/f'{kind}_preview_{index}.png')
            view.close();continue
        path=OUT/'animations'/f'{kind}_full_vessel_enlarged.mp4'
        view=View(kind,orbit)
        temp=path.with_name(path.stem+'.tmp.mp4')
        writer=imageio_ffmpeg.write_frames(str(temp),size=SIZE,fps=FPS,codec='libx264',quality=8,
            pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=8,
            output_params=['-preset','medium','-movflags','+faststart'])
        writer.send(None)
        try:
            for index in range(FRAMES):
                writer.send(np.asarray(view.frame(index)))
                if index%72==0:print(kind,index,'/',FRAMES,flush=True)
        finally:
            writer.close();view.close()
        os.replace(temp,path)
        videos.append(validate_video(path,kind))
        print(kind,'verified',flush=True)
    if args.stills_only:
        print('STILLS_COMPLETE',flush=True);return
    source_lock=json.loads((OUT/'SOURCE_LOCK.json').read_text())
    assert all(sha(CASE/name)==value for name,value in source_lock.items())
    record=dict(all_pass=True,background='black',colormap='turbo [0.10,0.96]',legend_outside_3d_viewport=True,
        camera_mode='Stationary orthographic camera; rigid actor rotation about one fixed longitudinal axis',
        rotation_axis_method='First principal axis of original surface coordinates; computed once',
        rotation_axis_unit=orbit.axis.tolist(),rotation_center_um=orbit.center.tolist(),
        camera_position_um=orbit.camera_position.tolist(),fixed_parallel_scale_um=orbit.scale,
        camera_position_constant=True,focal_point_constant=True,zoom_constant=True,axis_constant=True,
        angular_speed_deg_s=20.,display_transform='Proper rigid rotation; source coordinates and scalar values unchanged',
        all_mesh_points_inside_view_all_frames=True,geometry_smoothed=False,physical_fields_modified=False,
        clipping_of_scalar_values=False,animation_is_display_rotation_of_frozen_steady_field=True,
        source_hashes_unchanged=True,pressure_color_scale_Pa=list(SPECS['pressure']['clim']),
        wss_color_scale_Pa=list(SPECS['wss']['clim']),wss_display=compute['wss_display'],
        enlarged_streamlines_reused_without_integration=True,videos=videos,opengl=capabilities,
        figures=[dict(file=str(p.relative_to(OUT)),sha256=sha(p)) for p in sorted((OUT/'figures').glob('*'))],
        elapsed_seconds=time.time()-start,script_sha256=sha(__file__))
    dump(OUT/'MEDIA_VALIDATION.json',record)
    print(json.dumps(record,indent=2),flush=True)


if __name__=='__main__':main()
