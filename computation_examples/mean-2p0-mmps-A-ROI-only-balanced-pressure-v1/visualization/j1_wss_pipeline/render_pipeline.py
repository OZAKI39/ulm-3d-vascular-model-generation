"""Five real-data J1 turntables using the existing rendering/annotation style."""
from pathlib import Path
import argparse, hashlib, json, os, sys, time
os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW','vtkEGLRenderWindow')
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image,ImageDraw,ImageFont
import render_visualization as shared
import render_surface_fields as base

STAGES={
 '01_wall_normals':dict(title='Junction 1 | 1. Wall and outward normals',
    footer='Gold: unit node normals | Traction recovery uses planar facet normals',scalar=None,glyph='normal',color='#ffc857'),
 '02_velocity_gradient':dict(title='Junction 1 | 2. Near-wall velocity gradient',
    footer='Color: |G n| = |du/dn| | Full velocity tensor from adjacent P1 tetrahedra',
    scalar='NormalVelocityDerivative_display_s_inv',unit='Rate (s⁻¹)',clim=(0,16000),ticks=[0,4000,8000,12000,16000],glyph=None),
 '03_viscous_traction':dict(title='Junction 1 | 3. Viscous traction',
    footer='t_mu = mu (G + G^T) n | Gold: viscous traction | No pressure term',scalar=None,glyph='traction',color='#ffc857'),
 '04_tangential_wss':dict(title='Junction 1 | 4. Tangential WSS vectors',
    footer='tau_w = t_mu - (t_mu · n)n | Cyan: wall-on-fluid tangential traction',scalar=None,glyph='wss',color='#56dce5'),
 '05_wss_magnitude':dict(title='Junction 1 | 5. WSS magnitude and vectors',
    footer='Color: |tau_w| (nodal display) | Black: raw facet vectors, wall on fluid',
    scalar='WSS_display_Pa',unit='WSS (Pa)',clim=(0,55),ticks=[0,10,20,30,40,50,55],glyph='wss',color='black'),
}

def decorate(pixels,kind):
    spec=STAGES[kind]
    if spec['scalar']:return base.annotate(pixels,kind)
    im=Image.fromarray(pixels[:,:,:3]);w,h=im.size;draw=ImageDraw.Draw(im)
    draw.text((w*.025,h*.025),spec['title'],font=ImageFont.truetype(shared.FONT,round(h*.032)),fill='#f0f3fa')
    draw.text((w*.025,h*.958),spec['footer'],font=ImageFont.truetype(shared.FONT,round(h*.021)),fill='#c0c8d3')
    return im

def make_glyph(points,vec,length):
    seed=pv.PolyData(points);seed['direction']=vec;seed['length_um']=length
    geom=pv.Arrow(tip_length=.28,tip_radius=.10,shaft_radius=.030,tip_resolution=10,shaft_resolution=10)
    return seed.glyph(orient='direction',scale='length_um',factor=1.,geom=geom)

def inputs():
    compute=json.loads((HERE/'COMPUTE_VALIDATION.json').read_text());assert compute['all_pass']
    for name,digest in compute['source_lock'].items():assert shared.sha(ROOT/name)==digest,name
    for name,digest in compute['files'].items():assert shared.sha(HERE/name)==digest,name
    mesh=pv.read(HERE/'data/J1_wss_pipeline_si.vtp');mesh.points=np.asarray(mesh.points,dtype=float)*1e6
    meta=json.loads((HERE/'data/region_and_annotations.json').read_text());samples=np.load(HERE/'data/glyph_samples.npz')
    ni=samples['node_indices'];fi=samples['face_indices'];n=mesh.cell_data['Outward_normal'];centers=mesh.cell_centers().points
    node_n=mesh.point_data['NodeOutwardNormal'];offset=compute['glyph_offset_um'];factor=compute['traction_glyph_scale_um_per_Pa']
    glyphs={};glyphs['normal']=make_glyph(mesh.points[ni]+offset*node_n[ni],node_n[ni],np.full(len(ni),compute['normal_glyph_length_um']))
    for name,array in [('traction','ViscousTraction_Pa'),('wss','WSSVector_Pa')]:
        v=mesh.cell_data[array][fi];glyphs[name]=make_glyph(centers[fi]+offset*n[fi],v,np.linalg.norm(v,axis=1)*factor)
    framing=np.vstack([mesh.points,*[v.points for v in glyphs.values()]])
    orbit=shared.Turntable(framing,[0,0,1],meta['J1_um'])
    return compute,mesh,meta,glyphs,orbit,pv.PolyData(framing)

class LocalAnnotations(shared.AnnotationPlan):
    def add(self,text,size,centers,valid,cost,kind,**extra):
        if kind=='axis_title':
            # The close-up has four branch labels in a small physical window.
            # Add margin candidates, retaining the original clear-space checks,
            # typography and periodic fade planner; do not relax overlap limits.
            positions=np.array([(x,y) for x in [100.,220.,1490.,1600.] for y in [150.,280.,410.,540.,670.,800.,920.]])
            more=np.broadcast_to(positions,(shared.FRAMES,*positions.shape)).copy()
            w,h=self.text_size(text,size)
            ok=self.in_view(more,w,h)&self.unoccupied(more,w,h,pad=12.)&(self.clearance(more,w,h)>9.)
            centers=np.concatenate([centers,more],axis=1);valid=np.concatenate([valid,ok],axis=1)
            cost=np.concatenate([cost,np.full(ok.shape,20.)],axis=1)
        return super().add(text,size,centers,valid,cost,kind,**extra)

class PipelineView:
    def __init__(self,mesh,meta,glyphs,orbit,occluder,kind,size=shared.SIZE,annotations=None,raw=False,face_normals=False):
        self.kind=kind;self.orbit=orbit;self.raw=raw;spec=STAGES[kind]
        p=self.p=pv.Plotter(off_screen=True,window_size=size,shape=(1,2),col_weights=[.90,.10],border=False)
        p.set_background('black',all_renderers=True);p.subplot(0,0)
        p.renderer.SetViewport(*shared.VIEWPORT);p.render_window.SetMultiSamples(0)
        kwargs=dict(smooth_shading=True,split_sharp_edges=False,ambient=.68,diffuse=.32,specular=.08,show_scalar_bar=False)
        if spec['scalar']:
            scalar='WSS_raw_Pa' if raw else spec['scalar']
            actor=p.add_mesh(mesh,scalars=scalar,preference='cell' if raw else 'point',cmap=shared.CMAP,clim=spec['clim'],**kwargs)
            mapper=actor.mapper
            if raw:
                mapper.SetScalarModeToUseCellFieldData();mapper.SelectColorArray(scalar);mapper.InterpolateScalarsBeforeMappingOff()
            else:
                mapper.SetScalarModeToUsePointFieldData();mapper.SelectColorArray(scalar);mapper.InterpolateScalarsBeforeMappingOn()
            mapper.Update();actual=pv.wrap(mapper.GetInput())
            data=mesh.cell_data if raw else mesh.point_data;out=actual.cell_data if raw else actual.point_data
            assert np.array_equal(data[scalar],out[scalar]);assert np.array_equal(mesh.points,actual.points)
            assert np.array_equal(mesh.faces,actual.faces)
            assert data[scalar].min()>=spec['clim'][0] and data[scalar].max()<=spec['clim'][1]
            self.mapping=dict(scalar=scalar,mode=mapper.GetScalarModeAsString(),values_unchanged=True,
                interpolate_scalars=bool(mapper.GetInterpolateScalarsBeforeMapping()),geometry_unchanged=True)
        else:
            actor=p.add_mesh(mesh,color='#a8b3bf',show_edges=True,edge_color='#657181',line_width=.45*size[1]/1080,**kwargs)
            self.mapping=dict(scalar=None,geometry_unchanged=True)
        if spec['glyph']:
            glyph=glyphs[spec['glyph']]
            if face_normals:
                compute=json.loads((HERE/'COMPUTE_VALIDATION.json').read_text());fi=np.load(HERE/'data/glyph_samples.npz')['face_indices']
                n=mesh.cell_data['Outward_normal'][fi];centers=mesh.cell_centers().points[fi]
                glyph=make_glyph(centers+compute['glyph_offset_um']*n,n,np.full(len(fi),.9))
            p.add_mesh(glyph,color=spec['color'],ambient=1.,diffuse=0.,specular=0.,show_scalar_bar=False)
        orbit.setup(p)
        self.axes=shared.CoordinateGrid(p,mesh.bounds,orbit.center,'branch_detail',size)
        ports=[(name,np.array(point)) for name,point in meta['annotations']]
        self.annotations=annotations or LocalAnnotations(orbit,occluder,self.axes,ports)
        self.annotations.record['port_arrow_tip']='J1 center or true wall facet toward the named branch; no new boundary is introduced'
        self.annotations.record['clearance_geometry']='Wall and all three glyph geometries, shared across all stages'
        p.show(auto_close=False,interactive=False)

    def frame(self,i,still=False):
        c=self.orbit.trace[i];self.p.camera.position=c['camera_position_um'];self.p.camera.up=c['camera_up_unit'];self.p.render()
        assert np.allclose(self.p.camera.focal_point,self.orbit.center)
        assert self.p.camera.parallel_scale==self.orbit.scale
        image=decorate(self.p.screenshot(),self.kind)
        if self.raw:
            draw=ImageDraw.Draw(image);w,h=image.size
            draw.rectangle((0,0,w*.9,h*.07),fill='black')
            draw.text((w*.025,h*.025),'Junction 1 | Raw facet WSS comparison',font=ImageFont.truetype(shared.FONT,round(h*.032)),fill='#f0f3fa')
            draw.rectangle((0,h*.945,w,h),fill='black')
            draw.text((w*.025,h*.958),'Color: raw facet WSS | Same view and 0–55 Pa scale | Black: raw facet vectors',font=ImageFont.truetype(shared.FONT,round(h*.021)),fill='#c0c8d3')
        return self.annotations.overlay(image,i,still=still)

    def close(self):self.p.close()

def validate_movie(path,out,kind):
    reader=imageio_ffmpeg.read_frames(str(path));meta=next(reader)
    assert meta['size']==shared.SIZE and meta['fps']==shared.FPS and abs(meta['duration']-18)<.05
    hashes=set();frames={};count=0
    for i,raw in enumerate(reader):
        count+=1;hashes.add(hashlib.sha256(raw).hexdigest())
        if i in shared.KEYS[:4]:frames[i]=Image.fromarray(np.frombuffer(raw,np.uint8).reshape(1080,1920,3))
    assert count==432 and len(hashes)==432
    sheet=Image.new('RGB',(3840,2160),'black')
    for j,i in enumerate(shared.KEYS[:4]):sheet.paste(frames[i],(j%2*1920,j//2*1080))
    sheet.save(out/'figures'/f'{kind}_rotation_views.png')
    return dict(file=str(path.relative_to(out)),sha256=shared.sha(path),frames=count,distinct_frames=len(hashes),fps=24,size=[1920,1080],duration_s=18,all_pass=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=HERE);parser.add_argument('--previews-only',action='store_true')
    args=parser.parse_args();out=args.output.resolve();start=time.time()
    for d in ['figures','animations']:(out/d).mkdir(parents=True,exist_ok=True)
    compute,mesh,meta,glyphs,orbit,occluder=inputs()
    for name,spec in STAGES.items():
        if spec['scalar']:base.SPECS[name]=spec
        for field,fraction in [('title',.032),('footer',.021)]:
            font=ImageFont.truetype(shared.FONT,round(1080*fraction))
            assert font.getlength(spec[field])+48<1720,(name,field)
    annotations=None;maps={};gpu=[]
    for kind in STAGES:
        view=PipelineView(mesh,meta,glyphs,orbit,occluder,kind,size=(3840,2160),annotations=annotations)
        annotations=view.annotations;maps[kind]=view.mapping
        for j,i in enumerate(shared.KEYS[:4]):
            frame=view.frame(i,still=True);name=f'{kind}_overview_4k' if j==0 else f'{kind}_angle_{j*90:03d}_4k'
            frame.save(out/'figures'/f'{name}.png',dpi=(300,300))
            if j==0:frame.save(out/'figures'/f'{name}.pdf',resolution=300)
        if not gpu:
            gpu=[l for l in view.p.render_window.ReportCapabilities().splitlines() if any(t in l for t in ['OpenGL vendor','OpenGL renderer','OpenGL version'])]
            assert any('NVIDIA GeForce RTX 4090' in l for l in gpu),gpu
        view.close();print(kind,'4K_READY',flush=True)
    view=PipelineView(mesh,meta,glyphs,orbit,occluder,'05_wss_magnitude',size=(3840,2160),annotations=annotations,raw=True)
    for j,i in enumerate([0,216]):view.frame(i,still=True).save(out/'figures'/f'05_raw_facet_comparison_{j*180:03d}_4k.png',dpi=(300,300))
    rawmap=view.mapping;view.close()
    annotations.export(out/'J1_annotations.json');shared.dump(out/'J1_camera.json',orbit.trace)
    if args.previews_only:
        shared.dump(out/'PREVIEW_VALIDATION.json',dict(all_pass=True,OpenGL=gpu,mapping=maps,annotations=annotations.record));return
    videos=[]
    for kind in STAGES:
        movie=out/'animations'/f'J1_{kind}.mp4'
        if movie.exists():raise FileExistsError('Preserve completed MP4 '+str(movie))
        view=PipelineView(mesh,meta,glyphs,orbit,occluder,kind,annotations=annotations)
        temp=movie.with_suffix('.tmp.mp4')
        writer=imageio_ffmpeg.write_frames(str(temp),size=shared.SIZE,fps=shared.FPS,codec='libx264',quality=8,
            pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=8,
            output_params=['-preset','medium','-movflags','+faststart','-threads','4'])
        writer.send(None)
        try:
            for i in range(shared.FRAMES):
                writer.send(np.ascontiguousarray(view.frame(i)))
                if i%108==0:print(kind,i,'/',shared.FRAMES,flush=True)
        finally:writer.close();view.close()
        os.replace(temp,movie);videos.append(validate_movie(movie,out,kind))
        shared.dump(out/'RENDER_PROGRESS.json',dict(completed=[v['file'] for v in videos],total=5))
    for name,digest in compute['source_lock'].items():assert shared.sha(ROOT/name)==digest,name
    for name,digest in compute['files'].items():assert shared.sha(HERE/name)==digest,name
    shared.dump(out/'MEDIA_VALIDATION.json',dict(all_pass=True,videos=videos,OpenGL=gpu,encoder='CPU libx264',
        same_camera_all_stages=True,same_glyph_samples_and_linear_scale_steps_3_4_5=True,
        camera_scale_um=orbit.scale,display_center_um=orbit.center.tolist(),z_up=True,elevation_deg=shared.ELEVATION_DEG,
        colorbar_WSS_Pa=[0,55],colorbar_normal_velocity_derivative_s_inv=[0,16000],mapping=maps,raw_comparison_mapping=rawmap,
        annotation_checks=annotations.record,source_files_unchanged=True,CFD_calls=0,
        seconds=time.time()-start,script_sha256=shared.sha(Path(__file__))))
    print('ALL_FIVE_STAGES_COMPLETE',flush=True)

if __name__=='__main__':main()
