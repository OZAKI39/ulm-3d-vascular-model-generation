"""Black-background, speed-colored views of saved, verified FEM streamlines."""
from pathlib import Path
import argparse,hashlib,json,os,socket,time
os.environ.setdefault('LP_NUM_THREADS','2')
os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
import numpy as np
import pyvista as pv
from PIL import Image,ImageDraw,ImageFont
import imageio_ffmpeg
from render import Scene,CMAP,camera,projection_bounds,annotated
from compute_streamlines import CASE,OUT,QUOTAS,sha,dump

SIZE=(1920,1080);FPS=24;FRAMES=432;LIMIT=7.5


def label_frame(pixels,title,count):
    image=Image.fromarray(annotated(pixels,title,limit=LIMIT));w,h=image.size;draw=ImageDraw.Draw(image)
    draw.rectangle((0,int(h*.938),w,h),fill='black')
    font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',int(h*.020))
    draw.text((int(w*.025),int(h*.958)),
        f'{count} FEM streamlines | outlet-stratified display | line counts do not represent flow fractions',font=font,fill='#c0c8d3')
    return image


class StreamView:
    def __init__(self,scene,lines,size=SIZE,outlet=None,ports=False):
        self.scene=scene;self.size=size
        self.p=pv.Plotter(off_screen=True,window_size=size,shape=(1,2),col_weights=[.86,.14],border=False)
        p=self.p;p.set_background('black',all_renderers=True);p.subplot(0,0)
        p.renderer.SetViewport(0,.075,.86,.91);p.enable_depth_peeling(number_of_peels=8,occlusion_ratio=0)
        p.render_window.SetMultiSamples(0)
        wall=pv.read(CASE/'SV_MESH/mesh-surfaces/WALL.vtp');wall.points=np.asarray(wall.points)*1e6
        p.add_mesh(wall,color='#b7c4d1',opacity=.20,smooth_shading=True,ambient=.65,diffuse=.35,show_scalar_bar=False)
        selected=lines if outlet is None else lines.extract_cells(lines.cell_data['Outlet_id']==outlet)
        self.count=selected.n_cells
        p.add_mesh(selected,scalars='Speed_mm_s',cmap=CMAP,clim=(0,LIMIT),line_width=1.5*size[1]/1080,
            opacity=.9,lighting=False,show_scalar_bar=False)
        if ports:
            centers=[];labels=[]
            for role in ['INLET',*QUOTAS]:
                face=pv.read(CASE/'SV_MESH/mesh-surfaces'/f'{role}.vtp')
                centers.append(np.array(face.center)*1e6);labels.append('Inlet' if role=='INLET' else 'Out '+role[-2:])
            p.add_point_labels(np.array(centers),labels,font_size=int(19*size[1]/1080),text_color='#eef2f8',
                point_color='#eef2f8',point_size=3,shape=None,show_points=True,always_visible=True)
        camera(p,scene.center,scene.radius,45,scene.full_scale);p.show(auto_close=False,interactive=False)

    def frame(self,angle,title):
        camera(self.p,self.scene.center,self.scene.radius,angle,self.scene.full_scale);self.p.render()
        bounds=projection_bounds(self.scene.g.points,self.scene.center,angle,self.scene.full_scale,self.size[0]*.86/(self.size[1]*.835))
        assert max(abs(v) for v in bounds)<.94
        return label_frame(self.p.screenshot(),title,self.count),dict(azimuth_deg=float(angle),
            projected_bounds=bounds,parallel_scale_um=self.scene.full_scale,
            camera_position_um=list(self.p.camera.position),physical_coordinates_rotated=False)

    def close(self):self.p.close()


def decode(path):
    reader=imageio_ffmpeg.read_frames(str(path));metadata=next(reader)
    assert metadata['size']==SIZE and metadata['fps']==FPS and abs(metadata['duration']-18)<.05
    distinct=set();keys={};count=0
    for index,raw in enumerate(reader):
        distinct.add(hashlib.sha256(raw).hexdigest());count+=1
        if index in [0,144,288,431]:keys[index]=np.frombuffer(raw,np.uint8).reshape(1080,1920,3).copy()
    assert count==FRAMES and len(distinct)>=FRAMES*.95
    for array in keys.values():
        body=array[100:950,:1600].astype(float);assert (body.max(2)-body.min(2)>20).sum()>1000
        right=array[200:900,1660:1900].astype(float);assert (right.max(2)-right.min(2)>30).sum()>3000
    delta=float(np.mean(np.abs(keys[0].astype(float)-keys[144].astype(float))));assert delta>.05
    inspection=OUT/'inspection';inspection.mkdir(exist_ok=True)
    for index,array in keys.items():Image.fromarray(array).save(inspection/f'frame_{index:04d}.png')
    return dict(all_pass=True,file=str(path.relative_to(CASE)),sha256=sha(path),decoded_frames=count,
        distinct_frames=len(distinct),fps=FPS,size=SIZE,duration_s=18.,keyframe_mean_pixel_difference=delta),keys


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stills-only',action='store_true');args=parser.parse_args()
    start=time.time();compute=json.loads((OUT/'COMPUTE_VALIDATION.json').read_text());assert compute['all_pass']
    path=OUT/'data/streamlines_si.vtp';assert sha(path)==compute['outputs_sha256']['data/streamlines_si.vtp']
    lines=pv.read(path);assert lines.n_cells==96
    lines.points=np.asarray(lines.points,dtype=float)*1e6
    scene=Scene(CASE/'frozen_flow/steady_flow_mean_2p0_mmps.vtu',LIMIT)
    figures=CASE/'figures';receipts={};capabilities=[]
    view=StreamView(scene,lines,size=(3840,2160),ports=True)
    image,receipts['overview']=view.frame(45,'Inlet mean 2.0 mm/s | FEM streamlines')
    image.save(figures/'Figure_06_streamlines_overview.png')
    capabilities=[x for x in view.p.render_window.ReportCapabilities().splitlines() if any(k in x for k in ['OpenGL vendor','OpenGL renderer','OpenGL version'])]
    view.close()
    board=Image.new('RGB',(3840,2160),'black')
    for k,outlet in enumerate([None,1,2,3]):
        view=StreamView(scene,lines,outlet=outlet)
        title='All outlets | 96 streamlines' if outlet is None else f'Outlet {outlet:02d} | {QUOTAS[f"OUTLET_{outlet:02d}"]} selected streamlines'
        im,receipts[str(outlet)]=view.frame(45,title);board.paste(im,((k%2)*1920,(k//2)*1080));view.close()
    board.save(figures/'Figure_07_streamlines_outlet_coverage.png')
    dump(OUT/'STILL_VIEW_RECEIPTS.json',receipts)
    if args.stills_only:return
    video=CASE/'animations/Animation_03_streamlines_full_vessel.mp4';temp=video.with_name(video.stem+'.tmp.mp4')
    view=StreamView(scene,lines);writer=imageio_ffmpeg.write_frames(str(temp),size=SIZE,fps=FPS,codec='libx264',quality=8,
        pix_fmt_in='rgb24',pix_fmt_out='yuv420p',macro_block_size=8,output_params=['-preset','medium','-movflags','+faststart'])
    writer.send(None);frames=[]
    try:
        for index in range(FRAMES):
            im,record=view.frame(45+360*index/FRAMES,'Inlet mean 2.0 mm/s | FEM streamlines')
            writer.send(np.asarray(im));record.update(frame=index,video_time_s=index/FPS);frames.append(record)
            if index%72==0:print('STREAMLINE_RENDER',index,'/',FRAMES,flush=True)
    finally:writer.close();view.close()
    os.replace(temp,video);dump(OUT/'CAMERA_TRACE.json',frames)
    video_record,keys=decode(video)
    sheet=Image.new('RGB',(2880,640),'black');d=ImageDraw.Draw(sheet)
    font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',26)
    for col,index in enumerate([0,144,288]):
        d.text((col*960+25,25),f'Streamlines | {index/FPS:.0f} s | camera {45+360*index/FRAMES:.0f} deg',font=font,fill='#d9e3ef')
        im=Image.fromarray(keys[index]);im.thumbnail((960,540));sheet.paste(im,(col*960,85))
    sheet.save(figures/'Figure_08_streamlines_storyboard.png')
    lock=json.loads((OUT/'SOURCE_LOCK.json').read_text());assert all(sha(CASE/p)==h for p,h in lock.items())
    validation=dict(all_pass=True,source_field_sha256=compute['source_field_sha256'],source_lines_sha256=sha(path),
        line_count=96,outlet_counts=QUOTAS,colorbar_units='mm/s',colorbar_range=[0,LIMIT],background='black',
        legend_outside_3D_viewport=True,geometry_scale_preserved=True,all_full_vessel_frames_inside_view=True,
        curve_smoothing=False,line_width_role='DISPLAY_ONLY_NOT_PARTICLE_DIAMETER',
        outlet_selection_is_for_coverage_not_flux_proportional=True,
        animation_role='CAMERA_ORBIT_OF_STEADY_FIELD_NOT_18_SECONDS_OF_TRANSIENT_FLOW',video=video_record,
        figures=[dict(file=str((figures/name).relative_to(CASE)),sha256=sha(figures/name)) for name in
            ['Figure_06_streamlines_overview.png','Figure_07_streamlines_outlet_coverage.png','Figure_08_streamlines_storyboard.png']],
        hostname=socket.gethostname(),opengl=capabilities,render_seconds=time.time()-start,
        renderer_sha256=sha(__file__),manual_visual_review='PENDING_USER_REVIEW')
    dump(OUT/'MEDIA_VALIDATION.json',validation);print('STREAMLINE_MEDIA_COMPLETE',flush=True)


if __name__=='__main__':main()
