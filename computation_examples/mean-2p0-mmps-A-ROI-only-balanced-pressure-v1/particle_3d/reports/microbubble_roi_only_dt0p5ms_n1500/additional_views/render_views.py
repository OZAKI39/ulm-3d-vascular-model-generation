"""Render three new cameras from the completed cohort; never integrate motion."""
from pathlib import Path
import hashlib,json,os,sys
os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW','vtkEGLRenderWindow')
os.environ['IMAGEIO_FFMPEG_EXE']='/usr/bin/ffmpeg'
import numpy as np
import pyvista as pv
import imageio.v2 as imageio
from PIL import Image,ImageDraw,ImageFont

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
sys.path.insert(0,str(ROOT/'scripts'))
import render_results as shared


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    lock=json.loads((HERE/'audit/SOURCE_LOCK.json').read_text())
    for relative,digest in lock['protected_files'].items():assert sha(ROOT/relative)==digest,relative
    config=json.loads((ROOT/'config.json').read_text())
    assert sha(config['flow_path'])==lock['flow_sha256']
    folders=sorted((ROOT/'tracks').glob('mb_*'))
    assert len(folders)==1500
    rows=[];tracks=[];source_files={}
    for folder in folders:
        marker=json.loads((folder/'COMPLETE.json').read_text())
        for name in ['trajectory.npz','metrics.json']:
            f=folder/name;digest=sha(f);assert digest==marker['files'][name]
            source_files[str(f.relative_to(ROOT))]=digest
        rows.append(json.loads((folder/'metrics.json').read_text()))
        tracks.append(np.load(folder/'trajectory.npz')['samples'])
    summary=json.loads((ROOT/'data/final_summary.json').read_text())
    assert summary['numerical_gate'] and summary['all_completion_hashes_verified']
    groups={o:[] for o in ['O1','O2','O3']}
    for i,row in enumerate(rows):groups.setdefault(row['outlet'] or row['status'],[]).append(i)
    count_by_group={key:len(ids) for key,ids in groups.items()}
    assert {o:count_by_group[o] for o in ['O1','O2','O3']}==summary['outlets']
    surface=pv.read(config['flow_path']).extract_surface(algorithm='dataset_surface')
    surface.points*=1e6
    center=np.asarray(surface.center);radius=np.linalg.norm(np.ptp(surface.points,axis=0))/2
    original=np.array([-.72,.69,.11]);original/=np.linalg.norm(original)
    az=float(np.arctan2(original[1],original[0]));el=float(np.arcsin(original[2]))
    def direction(a,e):return np.array([np.cos(e)*np.cos(a),np.cos(e)*np.sin(a),np.sin(e)])
    views=[('side','Side view',direction(az+np.pi/2,el)),
           ('reverse','Reverse view',-original),
           ('oblique','Elevated oblique view',direction(az+np.pi/4,np.pi/4))]
    font=ROOT/'assets/NotoSansCJKsc-Regular.otf'
    title_font=ImageFont.truetype(str(font),30);small_font=ImageFont.truetype(str(font),23)
    end=max(s[-1,0] for s in tracks);frames=288;fps=24
    ages=np.linspace(0,end,frames)
    # Prepare exactly the same interpolation/exit display as the original movie.
    positions={}
    for key,ids in groups.items():
        values=np.empty((frames,len(ids),3))
        for j,i in enumerate(ids):
            s=tracks[i]
            values[:,j,:]=np.column_stack([np.interp(ages,s[:,0],s[:,k])*1e6 for k in (1,2,3)])
            if rows[i]['status']=='COMPLETED':values[ages>s[-1,0],j,:]=np.nan
        positions[key]=values
    manifests=[]
    for name,label,view_direction in views:
        movie=HERE/'animations'/f'microbubble_age_aligned_{name}.mp4'
        if movie.exists():raise FileExistsError(movie)
        plot=pv.Plotter(off_screen=True,window_size=(1920,1080));plot.set_background(shared.BACKGROUND)
        plot.add_mesh(surface,color='#b2becf',opacity=.16,smooth_shading=True)
        plot.camera.position=center+view_direction*radius*4;plot.camera.focal_point=center
        plot.camera.up=(0,0,1);plot.camera.parallel_projection=True;plot.camera.parallel_scale=radius*.78
        clouds={}
        for key,ids in groups.items():
            if not ids:continue
            cloud=pv.PolyData(positions[key][0].copy())
            plot.add_mesh(cloud,color=shared.COLORS[key],point_size=shared.MARKER_SIZE_PX,
                          render_points_as_spheres=True,lighting=True)
            clouds[key]=cloud
        plot.show(auto_close=False)
        layout=shared.enlarge_animation(plot,surface)
        gpu=[s for s in plot.ren_win.ReportCapabilities().splitlines()
             if any(k in s for k in ['vendor string','renderer string','version string'])]
        assert any('NVIDIA GeForce RTX 4090' in s for s in gpu),gpu
        keyframes=[]
        with imageio.get_writer(movie,fps=fps,codec='libx264',quality=8,macro_block_size=1,
                                ffmpeg_params=['-threads','2']) as writer:
            for frame,age in enumerate(ages):
                for key,cloud in clouds.items():cloud.points=positions[key][frame].copy()
                plot.render();pixels=Image.fromarray(plot.screenshot(return_img=True))
                draw=ImageDraw.Draw(pixels)
                title=f'Microbubble trajectories | {len(rows)} tracks | dt = 0.5 ms | Trajectory age: {age*1000:.1f} ms'
                footer='Age-aligned independent trajectories, not simultaneous injection. Marker size is for display; physical diameters are in the data.'
                assert title_font.getlength(title)<=1840 and small_font.getlength(footer)<=1840
                draw.text((40,25),title,font=title_font,fill='white')
                draw.text((40,1018),footer,font=small_font,fill='#cad4e2')
                x=40
                for key,ids in groups.items():
                    legend=f'{shared.LABELS[key]}: {len(ids)}'
                    assert x+small_font.getlength(legend)<=1880
                    draw.text((x,70),legend,font=small_font,fill=shared.COLORS[key]);x+=320
                writer.append_data(np.asarray(pixels))
                if frame in [24,96,160,240]:keyframes.append(pixels.copy())
                if frame==96:
                    pixels.save(HERE/'figures'/f'{name}_preview.png',dpi=(300,300))
                    pixels.convert('RGB').save(HERE/'figures'/f'{name}_preview.pdf',resolution=300)
        plot.close()
        sheet=Image.new('RGB',(1920,1080))
        for i,im in enumerate(keyframes):sheet.paste(im.resize((960,540)),((i%2)*960,(i//2)*540))
        sheet.save(HERE/'figures'/f'{name}_frames.png',dpi=(300,300))
        item=dict(name=name,label=label,video=str(movie.relative_to(HERE)),sha256=sha(movie),
                  camera_direction=view_direction.tolist(),camera_azimuth_deg=float(np.degrees(np.arctan2(view_direction[1],view_direction[0]))),
                  camera_elevation_deg=float(np.degrees(np.arcsin(view_direction[2]))),
                  layout=layout,render_backend=gpu,frames=frames,fps=fps,size=[1920,1080],seconds=frames/fps,
                  marker_size_px=shared.MARKER_SIZE_PX,zoom=shared.ANIMATION_ZOOM,
                  count_by_group=count_by_group,all_tracks_used=True,text_language='English',
                  age_end_s=end,age_grid_sha256=hashlib.sha256(ages.tobytes()).hexdigest())
        manifests.append(item);print(json.dumps(item),flush=True)
    for relative,digest in source_files.items():assert sha(ROOT/relative)==digest,relative
    for relative,digest in lock['protected_files'].items():assert sha(ROOT/relative)==digest,relative
    result=dict(views=manifests,source_flow_sha256=lock['flow_sha256'],
                source_track_files_verified=len(source_files),source_tracks_unchanged=True,
                original_renderer_sha256=sha(ROOT/'scripts/render_results.py'),
                original_movie_unchanged=sha(ROOT/'animations/microbubble_age_aligned.mp4')==lock['existing_movie_sha256'],
                particle_integrations=0,CFD_calls=0,encoder='CPU libx264',GPU='NVIDIA EGL RTX 4090',
                same_time_axis_as_original=True,same_camera_scale_for_all_views=True)
    (HERE/'RENDER_MANIFEST.json').write_text(json.dumps(result,indent=2)+'\n')
    print('THREE_ADDITIONAL_VIEWS_COMPLETE',flush=True)


if __name__=='__main__':main()
