"""Pressure and WSS using the enlarged, upright-Z full-vessel presentation.

All imports are from this self-contained directory or requirements.txt. Saved
pressure and WSS are read without recomputing or altering the physical fields.
"""
import argparse
import os
from pathlib import Path
import time

import render_visualization as shared
import imageio_ffmpeg
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE/'results/surface_fields'
SPECS = {
    'pressure': dict(title='Full vessel | Pressure field', scalar='Pressure_Pa',
        unit='Pressure (Pa)', clim=(-5.,2500.), ticks=[0,500,1000,1500,2000,2500],
        footer='Steady FEM pressure on the vessel surface', file='pressure_surface_si.vtp'),
    'wss': dict(title='Full vessel | Wall shear stress', scalar='WSS_display_Pa',
        unit='WSS (Pa)', clim=(0.,55.), ticks=[0,10,20,30,40,50,55],
        footer='Derived from FEM velocity gradient', file='wall_wss_si.vtp'),
}


def annotate(pixels, kind):
    """Transparent text: scalar-bar title/numbers and exactly one footer."""
    spec = SPECS[kind]
    image = Image.fromarray(pixels[:,:,:3])
    w,h = image.size
    draw = ImageDraw.Draw(image)
    title = ImageFont.truetype(shared.FONT,round(h*.032))
    font = ImageFont.truetype(shared.FONT,round(h*.021))
    draw.text((w*.025,h*.025),spec['title'],font=title,fill='#f0f3fa')
    draw.text((w*.025,h*.958),spec['footer'],font=font,fill='#c0c8d3')
    x,y,bw,bh = round(w*.922),round(h*.275),round(w*.014),round(h*.49)
    colors = (shared.CMAP(np.linspace(1,0,bh))[:,:3]*255).astype(np.uint8)
    image.paste(Image.fromarray(np.repeat(colors[:,None,:],bw,axis=1)),(x,y))
    draw.text((w*.904,h*.208),spec['unit'],font=font,fill='#f0f3fa')
    lo,hi = spec['clim']
    for value in spec['ticks']:
        yy = y+bh*(hi-value)/(hi-lo)
        draw.line((x+bw,yy,x+bw+w*.005,yy),fill='#e9eef4',width=max(1,round(w/1920)))
        draw.text((x+bw+w*.010,yy-h*.012),f'{value:g}',font=font,fill='#f0f3fa')
    return image


class SurfaceView:
    def __init__(self, mesh, full_surface, orbit, ports, kind, size=shared.SIZE, annotations=None):
        self.kind,self.orbit = kind,orbit
        self.p = pv.Plotter(off_screen=True,window_size=size,shape=(1,2),col_weights=[.90,.10],border=False)
        p = self.p
        p.set_background('black',all_renderers=True)
        p.subplot(0,0)
        p.renderer.SetViewport(*shared.VIEWPORT)
        p.render_window.SetMultiSamples(0)
        spec = SPECS[kind]
        self.actor = p.add_mesh(mesh,scalars=spec['scalar'],cmap=shared.CMAP,clim=spec['clim'],
            preference='point',smooth_shading=True,split_sharp_edges=False,
            ambient=.68,diffuse=.32,specular=.08,show_scalar_bar=False)
        orbit.setup(p)
        self.axes = shared.CoordinateGrid(p,full_surface.bounds,orbit.center,'full_vessel',size)
        self.annotations = annotations or shared.AnnotationPlan(orbit,full_surface,self.axes,ports)
        p.show(auto_close=False,interactive=False)

    def frame(self,index,still=False):
        camera = self.orbit.trace[index]
        self.p.camera.position = camera['camera_position_um']
        self.p.camera.up = camera['camera_up_unit']
        self.p.render()
        assert np.allclose(self.p.camera.position,camera['camera_position_um'])
        assert np.allclose(self.p.camera.focal_point,self.orbit.center)
        assert self.p.camera.parallel_scale == self.orbit.scale
        return self.annotations.overlay(annotate(self.p.screenshot(),self.kind),index,still=still)

    def close(self):
        self.p.close()


def validate_video(path,kind,output):
    reader = imageio_ffmpeg.read_frames(str(path))
    metadata = next(reader)
    assert metadata['size'] == shared.SIZE and metadata['fps'] == shared.FPS
    assert abs(metadata['duration']-shared.FRAMES/shared.FPS)<.05
    frames,distinct = {},set()
    count = 0
    font = ImageFont.truetype(shared.FONT,round(shared.SIZE[1]*.021))
    footer_end = int(shared.SIZE[0]*.025+font.getlength(SPECS[kind]['footer'])+30)
    for index,raw in enumerate(reader):
        count += 1
        distinct.add(shared.hashlib.sha256(raw).hexdigest())
        frame = np.frombuffer(raw,np.uint8).reshape(shared.SIZE[1],shared.SIZE[0],3)
        # Audit every decoded frame, including the forbidden extra legend text.
        assert frame[870:1000,1735:1910].mean()<1.
        assert frame[1020:1078,footer_end:1910].mean()<1.
        if index in shared.KEYS:
            frames[index] = frame.copy()
            Image.fromarray(frame).save(output/'inspection'/f'{kind}_{index:04d}.png')
    assert count == shared.FRAMES and len(distinct) == shared.FRAMES
    sheet = Image.new('RGB',(3840,2160),'black')
    for j,index in enumerate(shared.KEYS[:4]):
        sheet.paste(Image.fromarray(frames[index]),(j%2*1920,j//2*1080))
    sheet.save(output/'figures'/f'{kind}_rotation_views.png')
    return dict(file=str(path.relative_to(output)),sha256=shared.sha(path),
        decoded_frames=count,distinct_frames=len(distinct),size=shared.SIZE,
        fps=shared.FPS,duration_s=shared.FRAMES/shared.FPS,
        no_extra_colorbar_text_all_frames=True,footer=SPECS[kind]['footer'],all_pass=True)


def main():
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',type=Path,default=HERE/'input_data')
    parser.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument('--stills-only',action='store_true')
    args = parser.parse_args()
    case,output = args.case.resolve(),args.output.resolve()
    start = time.time()
    for name in ['figures','animations','inspection']:
        (output/name).mkdir(parents=True,exist_ok=True)
    compute = json.loads((case/'field_diagnostics/COMPUTE_VALIDATION.json').read_text())
    previous = json.loads((case/'field_diagnostics/MEDIA_VALIDATION.json').read_text())
    assert compute['all_pass']
    sources = ['field_diagnostics/data/'+spec['file'] for spec in SPECS.values()]
    sources += ['field_diagnostics/COMPUTE_VALIDATION.json','field_diagnostics/MEDIA_VALIDATION.json',
                'streamlines/data/streamlines_si.vtp']
    sources += [f'SV_MESH/mesh-surfaces/{role}.vtp' for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']]
    hashes = {name:shared.sha(case/name) for name in sources}
    meshes,original_values = {},{}
    for kind,spec in SPECS.items():
        relative = 'field_diagnostics/data/'+spec['file']
        assert hashes[relative] == compute['outputs_sha256']['data/'+spec['file']]
        mesh = pv.read(case/relative)
        original_values[kind] = mesh[spec['scalar']].copy()
        values = original_values[kind]
        assert np.isfinite(values).all()
        assert values.min()>=spec['clim'][0] and values.max()<=spec['clim'][1]
        mesh.points = np.asarray(mesh.points,dtype=float)*1e6
        meshes[kind] = mesh
    surface = meshes['pressure']
    # Match the exact full-vessel framing calculation, using the existing lines
    # only as additional bounds. No streamlines are computed or drawn here.
    lines = pv.read(case/'streamlines/data/streamlines_si.vtp')
    framing_points = np.vstack([surface.points,np.asarray(lines.points,dtype=float)*1e6])
    center = shared.fitted_z_center(framing_points)
    orbit = shared.Turntable(framing_points,[0.,0.,1.],center)
    ports = []
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        cap = pv.read(case/f'SV_MESH/mesh-surfaces/{role}.vtp')
        ports.append(('Inlet' if role=='INLET' else 'Outlet '+role[-2:],np.array(cap.center)*1e6))
    magnification = previous['fixed_parallel_scale_um']/orbit.scale*shared.VIEW_HEIGHT/.835
    assert magnification>1.3
    print(f'Fixed Z axis; magnification versus old pressure/WSS: {magnification:.6f}',flush=True)
    annotations,views,videos = None,{},[]
    # Complete both stills first so the layout can be inspected during encoding.
    for kind,spec in SPECS.items():
        view = SurfaceView(meshes[kind],surface,orbit,ports,kind,size=(3840,2160),annotations=annotations)
        annotations = view.annotations
        view.frame(0,still=True).save(output/'figures'/f'{kind}_overview_4k.png')
        annotations.export(output/f'{kind}_annotations.json')
        shared.dump(output/f'{kind}_camera.json',orbit.trace)
        views[kind] = dict(title=spec['title'],footer=spec['footer'],
            colorbar_title=spec['unit'],colorbar_ticks=spec['ticks'],colorbar_range_Pa=spec['clim'],
            colorbar_extra_text=False,scalar=spec['scalar'],
            scalar_min_Pa=float(original_values[kind].min()),scalar_max_Pa=float(original_values[kind].max()),
            source_surface='field_diagnostics/data/'+spec['file'],
            fixed_parallel_scale_um=orbit.scale,pixels_per_um=1080*shared.VIEW_HEIGHT/(2*orbit.scale),
            magnification_vs_previous_pressure_wss=magnification,
            coordinate_grid_reference_um=view.axes.origin.tolist(),coordinate_tick_step_um=view.axes.step,
            annotations=annotations.record,
            opengl=[line for line in view.p.render_window.ReportCapabilities().splitlines()
                    if any(t in line for t in ['OpenGL vendor','OpenGL renderer','OpenGL version'])])
        view.close()
    print('BOTH_4K_STILLS_READY',flush=True)
    for kind in SPECS:
        view = SurfaceView(meshes[kind],surface,orbit,ports,kind,annotations=annotations)
        if args.stills_only:
            for index in shared.KEYS[:4]:
                view.frame(index).save(output/'inspection'/f'{kind}_preview_{index:04d}.png')
            view.close()
            continue
        path = output/'animations'/f'{kind}_full_vessel_enlarged.mp4'
        temporary = path.with_name(path.stem+'.tmp.mp4')
        writer = imageio_ffmpeg.write_frames(str(temporary),size=shared.SIZE,fps=shared.FPS,
            codec='libx264',quality=8,pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=8,
            output_params=['-preset','medium','-movflags','+faststart'])
        writer.send(None)
        try:
            for index in range(shared.FRAMES):
                writer.send(np.ascontiguousarray(view.frame(index)))
                if index%72==0:
                    print(kind,index,'/',shared.FRAMES,flush=True)
        finally:
            writer.close()
            view.close()
        os.replace(temporary,path)
        videos.append(validate_video(path,kind,output))
    assert {name:shared.sha(case/name) for name in sources} == hashes
    for kind,spec in SPECS.items():
        assert np.array_equal(meshes[kind][spec['scalar']],original_values[kind])
    shared.dump(output/('PREVIEW_VALIDATION.json' if args.stills_only else 'MEDIA_VALIDATION.json'),
        dict(all_pass=True,source_sha256=hashes,source_values_unchanged=True,
            physical_geometry_unchanged=True,geometry_smoothing=False,scalar_clipping=False,
            input_coordinate_unit='m',display_coordinate_unit='µm',scalar_unit='Pa',
            z_axis_vertical_all_frames=True,rotation_axis=[0,0,1],camera_center_um=center.tolist(),
            fixed_parallel_scale_um=orbit.scale,camera_elevation_deg=shared.ELEVATION_DEG,
            frames=shared.FRAMES,fps=shared.FPS,display_rotation_duration_s=shared.FRAMES/shared.FPS,
            text_background='transparent',text_crossfade_duration_s=21/shared.FPS,
            still_policy='One dominant clear label position per text at full opacity',
            wss_display=compute['wss_display'],wss_method=compute['wss_method'],
            source_fields_recomputed=False,views=views,videos=videos,
            script_sha256=shared.sha(Path(__file__)),shared_renderer_sha256=shared.sha(HERE/'render_visualization.py'),
            elapsed_seconds=time.time()-start))
    if not args.stills_only:
        write_gallery(output)
    print('SURFACE_PREVIEWS_COMPLETE' if args.stills_only else 'SURFACE_RENDER_COMPLETE',flush=True)


def write_gallery(output):
    sections = []
    for kind,name in [('pressure','全血管压力场'),('wss','全血管壁面剪切应力（WSS）')]:
        sections.append(f'''<section><h2>{name}</h2>
<video controls loop preload="metadata" poster="figures/{kind}_overview_4k.png"
src="animations/{kind}_full_vessel_enlarged.mp4"></video>
<p><a href="animations/{kind}_full_vessel_enlarged.mp4" download>下载 MP4</a> ·
<a href="figures/{kind}_overview_4k.png" download>下载 4K 图片</a></p></section>''')
    (output/'OPEN_RESULTS.html').write_text('''<!doctype html>
<html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>FEM 压力与 WSS</title>
<style>body{max-width:1280px;margin:32px auto;padding:0 20px;background:#0b0d11;color:#edf2f7;font:17px/1.6 system-ui}
video{display:block;width:100%;background:black}a{color:#90caff}section{margin:32px 0}h1{font-size:27px}h2{font-size:21px}</style>
<h1>FEM 全血管压力场与 WSS</h1>
<p>放大视图 · Z 轴竖直旋转 · 透明文字与渐变切换 · 1920 × 1080 · 24 fps · 每段 18 秒</p>
'''+''.join(sections)+'''
<p>18 秒为稳态场的展示旋转时长。压力和 WSS 的单位均为 Pa；坐标单位为 µm。
WSS 显示沿用壁面三角形数值的面积加权节点平均。</p>
</html>''',encoding='utf-8')


if __name__=='__main__':
    main()
