"""P8.2 academic animation layer; original-size MBs and saved polylines only."""
from pathlib import Path
import os
os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image,ImageDraw
from .particle81_visuals import VesselView,COLORS,INK,text
from .particle81_replay import TITLE,FOOTERS,INDEPENDENT,OUTLETS
from .particle82_results import Scene
from .particle82_provenance import atomic_json
from .particle8_replay import canonical_hash

ANIMATIONS=['01_full_network_rotating','02_ulm_accumulation','03_outlet_colored','04_stop_diagnostics']
SIZE=(1600,1000);FPS=15


def schedule(scene,kind):
    if kind=='02_ulm_accumulation':
        return [dict(time_s=float(t),label='ACQUISITION TIME COMPRESSION; completed paths only',only_ids=None)
                for t in np.linspace(0,scene.catalog['replay_end_time_s'],240)]
    reps=list(scene.catalog['representative_ids'])
    if kind=='04_stop_diagnostics':
        stopped=[pid for pid,e in scene.entries.items() if e['sample_count'] and not e['completed']]
        reps=[]
        for basin in OUTLETS+['UNRESOLVED_POINT_PATH']:
            ids=[pid for pid in stopped if scene.entries[pid]['point_tracer_basin']==basin]
            if ids:reps.append(ids[len(ids)//2])
        for pid in stopped:
            if len(reps)>=3:break
            if pid not in reps:reps.append(pid)
    frames=[]
    for pid in sorted(reps):
        e=scene.entries[pid]
        for t in np.linspace(e['birth_time_s'],e['last_physical_time_s'],80):
            frames.append(dict(time_s=float(t),label=f'MB {pid} saved interval; SLOW MOTION; idle time CUT',only_ids=None))
    for _ in range(20):frames.append(dict(time_s=scene.catalog['replay_end_time_s'],label='CUT TO FINAL ACCOUNTING; no extrapolation',only_ids=None))
    return frames


def compose(scene,view,kind,index,total,spec):
    state=scene.snapshot(spec['time_s']);view.clear();finished=state['completed_ids']
    view.paths(finished,'history',opacity=.8 if kind=='03_outlet_colored' else .6 if kind=='02_ulm_accumulation' else .4,
               width=1.9 if kind=='03_outlet_colored' else 1.3)
    if kind!='02_ulm_accumulation':
        partial={p['particle_id']:scene.partial_path(p['particle_id'],p['elapsed_time_s']) for p in state['active']}
        view.paths(list(partial),'active_path',partial=partial,opacity=.95,width=3.)
    stops=[]
    if kind=='04_stop_diagnostics':
        stops=[pid for pid,e in scene.entries.items() if e['sample_count'] and not e['completed'] and e['last_physical_time_s']<=spec['time_s']]
        view.paths(stops,'stopped_paths',opacity=.25,width=1.,color_by_outlet=False)
        if stops:
            actor=view.plotter.renderer.actors.get('stopped_paths')
            if actor is not None:actor.GetProperty().SetColor(.5,.5,.5)
            points=np.array([scene.arrays[pid][-1,1:4] for pid in stops])
            mesh=pv.PolyData((points-view.origin)*1e6)
            mesh['basin']=np.array([OUTLETS.index(scene.entries[pid]['point_tracer_basin'])+1 if scene.entries[pid]['point_tracer_basin'] in OUTLETS else 0 for pid in stops])
            view.plotter.add_mesh(mesh,scalars='basin',cmap=['#999999',*COLORS.values()],clim=(0,3),point_size=5,
                render_points_as_spheres=False,show_scalar_bar=False,name='stop_markers',reset_camera=False,render=False)
            view.dynamic.append('stop_markers')
    view.current(state['active']);az=25+120*index/max(1,total-1);picture,camera=view.capture(az)
    image=Image.new('RGB',SIZE,'white');image.paste(picture,(20,174));d=ImageDraw.Draw(image)
    names={'01_full_network_rotating':'Full-vessel rotating replay / observed coverage',
        '02_ulm_accumulation':'ULM accumulation / observed network',
        '03_outlet_colored':'Outlet-colored natural trajectory ensemble',
        '04_stop_diagnostics':'Computational stop diagnostics / unknown later fate'}
    text(d,(40,22),'Particle-8.2 | '+names[kind],27,bold=True)
    text(d,(40,65),TITLE,18,'#208ca6',True);d.line((40,100,1560,100),fill='#208ca6',width=3)
    text(d,(40,110),f'Physical t = {spec["time_s"]:.6f} s | frame {index+1}/{total} | {spec["label"]}',15)
    text(d,(40,140),'NATURAL FLUX-WEIGHTED MBs | verified remote-server integration | same saved scene',16,bold=True)
    text(d,(1132,186),'ACQUISITION ACCOUNTING',18,bold=True)
    for i,(label,key) in enumerate([('Scheduled','scheduled'),('Admitted','admitted'),('Active saved tracks','active'),('Completed','completed'),('Unresolved admission','unresolved'),('Stopped incomplete','stopped')]):
        text(d,(1132,226+i*34),label,17);text(d,(1515,226+i*34),str(state['counts'][key]),18,bold=True)
    text(d,(1132,448),'COMPLETE PATH: FINAL OUTLET',15,bold=True)
    for i,o in enumerate(OUTLETS):
        text(d,(1132,486+i*31),o,17,COLORS[o],True)
        text(d,(1515,486+i*31),str(sum(scene.entries[pid]['exit_outlet']==o for pid in finished)),17,COLORS[o])
    labels=['Original-size sphere = active MB','History = saved completed paths','Grey vessel = full Frozen mesh',
            'No new MB sampling during replay']
    if kind=='04_stop_diagnostics':labels=['Stop markers: point-basin colors','Marker pixels do not encode radius','Grey trails: incomplete saved paths','Stops are NOT physiological capture']
    for i,line in enumerate(labels):text(d,(1132,606+i*26),line,15)
    covered=sum(n>0 for n in scene.catalog['outlet_counts'].values())
    text(d,(1132,741),f'Observed outlet branches: {covered}/3',17,'#884b43',True)
    text(d,(1132,774),'Unobserved branches stay grey',15)
    text(d,(1132,805),'Production timestep: NOT FROZEN',14,'#884b43')
    px=20*665/(2*view.scale);d.line((240,816,240+px,816),fill=INK,width=3);text(d,(240,788),'20 µm | equal physical scale',14)
    text(d,(40,850),f'Camera {az:.1f}° / elevation 28° | DISPLAY interpolation / slow motion / cuts only; original coordinates unchanged',14)
    for i,line in enumerate(FOOTERS):text(d,(40,883+25*i),line,13,'#884b43')
    receipt=dict(frame_index=index,physical_time_s=spec['time_s'],state=state,state_sha256=canonical_hash(state),
        source_scene_sha256=scene.sha256,rendered_completed_ids=finished,rendered_active_ids=[r['particle_id'] for r in state['active']],
        rendered_stop_ids=stops,camera=camera,physical_coordinates_rotated=False,no_resampling=True,no_extrapolation=True,
        render_window_backend=view.plotter.render_window.GetClassName(),renderer_source_sha256=__import__('hashlib').sha256(Path(__file__).read_bytes()).hexdigest(),
        display_only_interpolation=True,polyline_decimation=False,time_mapping_label=spec['label'],labels=[TITLE,INDEPENDENT,*FOOTERS])
    return image,receipt


def render(root):
    root=Path(root);scene=Scene(root)
    for folder in ['animations','frames','keyframes','inspection']:(root/folder).mkdir(exist_ok=True)
    for kind in ANIMATIONS:
        frames=schedule(scene,kind);view=VesselView(scene,wall_opacity=.055 if kind=='02_ulm_accumulation' else .12)
        name='particle8_2_anim_'+kind;keys=[0,len(frames)//2,len(frames)-1]
        enc=imageio_ffmpeg.write_frames(str(root/'animations'/(name+'.mp4')),SIZE,fps=FPS,codec='libx264',quality=8,
            macro_block_size=2,pix_fmt_in='rgb24',pix_fmt_out='yuv420p',output_params=['-movflags','+faststart'],ffmpeg_log_level='error');enc.send(None)
        records=[];gifs=[]
        try:
            for i,spec in enumerate(frames):
                image,record=compose(scene,view,kind,i,len(frames),spec);enc.send(np.asarray(image));records.append(record)
                if i in keys:image.save(root/'keyframes'/f'{name}_{i:04d}.png')
                if i%5==0 or i in keys:gifs.append(image.resize((960,600)))
                if i%60==0:print('RENDER',kind,i,'/',len(frames),flush=True)
        finally:enc.close();view.close()
        gifs[0].save(root/'animations'/(name+'.gif'),save_all=True,append_images=gifs[1:],duration=333,loop=0,optimize=False)
        atomic_json(root/'frames'/(name+'.json'),dict(kind=kind,source_scene_sha256=scene.sha256,frame_count=len(records),
            fps=FPS,size=SIZE,keyframe_indices=keys,frames=records))
