"""Additional red, Z-up visualization of the sealed V3 accepted states only."""
from pathlib import Path
import os, sys, json, csv, hashlib, argparse, subprocess, datetime, time
ROOT = Path(os.environ.get('Y_CLOT_PROJECT_ROOT', '/workspace/y_clot/project'))
sys.path.insert(0, str(ROOT / 'analysis/clot_fragment_surface_coupling_v3'))
import render_videos as old
from common import R, O, F, DX, dump, sha, load_npz
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from matplotlib import colormaps

OUT = R / 'results/clot_fragment_surface_coupling_v3/red_z_vertical_20261019'
AUDIT = R / 'verification/clot_fragment_surface_coupling_v3_red_visuals/20261019'
NAME = 'VIDEO_V3_RED_Z_VERTICAL'
W, H, FPS = 1920, 1080, 24
FORCE_SCALE = 1e8  # 2 micro-newtons -> 200 micrometres in world coordinates.
RED, LIGHT_RED, DARK_RED, BLUE = '#b51e32', '#dc4050', '#741426', '#1866a0'
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'


def font(size, scale=1):
    return ImageFont.truetype(FONT, round(size * scale))


def project(pl, point):
    pl.renderer.SetWorldPoint(*point, 1)
    pl.renderer.WorldToDisplay()
    x, y, _ = pl.renderer.GetDisplayPoint()
    return np.array([x, pl.window_size[1] - y])


def setup(bounds, size=(W, H)):
    pl = pv.Plotter(off_screen=True, window_size=size)
    pl.set_background('#ffffff')
    pl.renderer.SetViewport(.055, .13, .96, .91)
    lo, hi = np.asarray(bounds)[::2], np.asarray(bounds)[1::2]
    center = (lo + hi) / 2
    direction = np.array([1.2, -3.0, .72])
    span = max(hi-lo)
    pl.camera_position = [center + span * direction, center, (0, 0, 1)]
    pl.enable_parallel_projection()
    forward = -direction / np.linalg.norm(direction)
    right = np.cross(forward, [0., 0., 1.]); right /= np.linalg.norm(right)
    up = np.cross(right, forward)
    corners = np.array([(x, y, z) for x in bounds[:2] for y in bounds[2:4] for z in bounds[4:]])
    aspect = (.905 * size[0]) / (.78 * size[1])
    pl.camera.parallel_scale = float(1.13 * max(np.ptp(corners @ up) / 2, np.ptp(corners @ right) / (2 * aspect)))
    pl.camera.clipping_range = (span * .01, span * 20)
    grid = old.box_grid_points(bounds)
    pl.add_mesh(old.line_mesh(grid, np.arange(len(grid)).reshape(-1, 2)),
                color='#61717c', opacity=.55, line_width=1.0 * size[0] / W, lighting=False)
    return pl


def choose_arrows(s, eligible, limit):
    """Spatially spread real vectors; never average or invent surface loads."""
    norms = np.linalg.norm(s['surface_force'], axis=1)
    ids = np.flatnonzero(eligible & (norms > 1e-16))
    if len(ids) <= limit:
        return ids
    points = s['surface_centroid'][ids]
    scaled = (points-points.min(0)) / np.maximum(np.ptp(points, axis=0), DX)
    chosen = [int(np.argmax(norms[ids]))]
    dist = np.sum((scaled-scaled[chosen[0]])**2, axis=1)
    while len(chosen) < limit:
        score = dist * (.15 + .85 * norms[ids]/norms[ids].max())
        score[chosen] = -1
        j = int(np.argmax(score)); chosen.append(j)
        dist = np.minimum(dist, np.sum((scaled-scaled[j])**2, axis=1))
    return ids[chosen]


def scene(pl, ref, surf, s, geometry, selected=None, field=False, arrow_limit=38):
    n = int(ref['clot_count']); xyz = s['x'] * 1e6
    mask = np.ones(n, bool) if selected is None else selected
    member = np.r_[mask, np.zeros(len(xyz)-n, bool)]
    labels, sizes, attached = old.graph.topology(ref, s['d'])
    free = ~attached[labels]
    scale = pl.window_size[0]/W
    if selected is None:
        pl.add_points(xyz[n:], color='#7e878c', point_size=3.1*scale,
                      render_points_as_spheres=True, opacity=.32)
    for detached, color in [(False, RED), (True, LIGHT_RED)]:
        ids = mask & (free == detached)
        if ids.any():
            pl.add_points(xyz[:n][ids], color=color, point_size=(6.2 if selected is None else 13)*scale,
                          render_points_as_spheres=True, opacity=.94, ambient=.25, diffuse=.75)
    b = ref['bonds']
    alive = (s['d'] < 1) & member[b].all(1)
    # The complete live graph defines components. Only its short edges are drawn
    # globally to avoid a filled mesh; the local static view draws every live edge.
    shown = alive if selected is not None else alive & (ref['r0'] <= 1.01*DX)
    if shown.any():
        pl.add_mesh(old.line_mesh(xyz, b[shown]), color=DARK_RED, opacity=.36 if selected is not None else .20,
                    line_width=(1.0 if selected is not None else .7)*scale, lighting=False)
    active = s['surface_active']; tri = surf.base.tri; length = surf.L
    l1 = active[:length] & member[tri].all(1)
    l2 = np.flatnonzero(active[length:] & member[surf.owner])
    groups = []
    if l1.any(): groups.append((old.polygon_mesh(xyz[tri[l1]]), np.flatnonzero(l1), '#ec6970'))
    if len(l2):
        quads = []
        for p in l2:
            i, a = surf.owner[p], surf.axis[p]; j, k = (a+1)%3, (a+2)%3
            corners = np.array([surf.offset[p]+DX*(u*np.eye(3)[j]+v*np.eye(3)[k])
                                for u,v in [(-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5)]])
            quads.append((s['x'][i]+corners @ geometry['l2_deformation'][i].T)*1e6)
        groups.append((old.polygon_mesh(quads), length+l2, '#f18b80'))
    tau = np.linalg.norm(s['surface_traction'], axis=1)
    clim = float(np.ceil(max(tau[np.r_[np.flatnonzero(l1), length+l2]].max(initial=0), 1)/100)*100)
    for mesh, ids, color in groups:
        if field:
            mesh.cell_data['traction_Pa'] = tau[ids]
            pl.add_mesh(mesh, scalars='traction_Pa', preference='cell', cmap='viridis', clim=(0,clim),
                        show_scalar_bar=False, opacity=.90, lighting=False, show_edges=True,
                        edge_color='#31535c', line_width=.45*scale)
        else:
            pl.add_mesh(mesh, color=color, opacity=.46, show_edges=False, lighting=False)
    eligible = np.r_[l1, active[length:] & member[surf.owner]]
    chosen = choose_arrows(s, eligible, arrow_limit)
    if len(chosen):
        vectors = s['surface_force'][chosen] * FORCE_SCALE
        pl.add_arrows(s['surface_centroid'][chosen]*1e6, vectors, mag=1, color=BLUE, lighting=False)
    stats = dict(actual_time_s=float(s['time']), selected_nodes=int(mask.sum()),
                 selected_original_atom_ids=(np.flatnonzero(mask)+1).tolist() if selected is not None else None,
                 L1_active=int(l1.sum()), L2_active=int(len(l2)), internal_live_bonds=int(alive.sum()),
                 displayed_live_bonds=int(shown.sum()), displayed_surface_indices=chosen.tolist(),
                 active_surface_indices=np.flatnonzero(eligible).tolist(),
                 resultant_force_N=s['surface_force'][eligible].sum(0).tolist(),
                 effective_area_m2=float(np.sum(s['surface_area'][eligible]*s['surface_weight'][eligible])),
                 actual_traction_max_Pa=float(tau[eligible].max(initial=0)), color_limit_Pa=clim)
    return stats


def annotate(pl, raw, bounds, s, stats, field=False):
    # VTK leaves pixels outside a reduced renderer viewport black. Use a white
    # publication canvas around the rendered viewport before adding annotations.
    height, width = raw.shape[:2]
    viewport = (round(.055*width), round(.09*height), round(.96*width), round(.87*height))
    image = Image.new('RGB', (width, height), 'white')
    image.paste(Image.fromarray(raw).convert('RGB').crop(viewport), viewport[:2])
    d = ImageDraw.Draw(image)
    width,height = image.size; scale = width/W
    def text(at, txt, size=27, anchor='lt', fill='#28333b'):
        d.text(tuple(np.asarray(at)*scale), txt, font=font(size,scale), anchor=anchor, fill=fill)
    text((78,36), f't = {float(s["time"])*1e3:.4f} ms', 30)
    if field: text((1840,36),'Surface traction',30,'rt')
    lo,hi = np.asarray(bounds)[::2],np.asarray(bounds)[1::2]
    for axis in range(3):
        start=lo.copy(); end=lo.copy()
        if axis==1: start[0]=end[0]=hi[0]
        end[axis]=hi[axis]
        anchor=['mt','lm','rm'][axis]
        shift=np.array([(0,13),(18,5),(-18,0)][axis])*scale
        for k,value in enumerate(np.linspace(lo[axis],hi[axis],4)):
            if axis in (1,2) and k==0:continue
            p=start.copy();p[axis]=value
            pos=project(pl,p)+shift
            d.text(tuple(pos),f'{value:.0f}',font=font(25,scale),anchor=anchor,fill='#34454f')
        extra=np.array([(0,48),(130,40),(-100,0)][axis])*scale
        pos=project(pl,(start+end)/2)+extra
        d.text(tuple(pos),'XYZ'[axis]+' (µm)',font=font(28,scale),anchor=anchor,fill='#24333c')
    center=(lo+hi)/2
    z0=project(pl,center);z1=project(pl,center+[0,0,200])
    assert abs(z1[0]-z0[0])<1e-7 and z1[1]<z0[1],(z0,z1)
    stats['z_axis_horizontal_error_px']=float(abs(z1[0]-z0[0]))
    # A view-plane reference arrow uses the exact same world length scale.
    pixels_per_um=(.78*height)/(2*pl.camera.parallel_scale)
    arrow_len=200*pixels_per_um; y=height-48*scale; x=width*.72
    d.line([(x,y),(x+arrow_len,y)],fill=BLUE,width=max(2,round(3*scale)))
    d.polygon([(x+arrow_len,y),(x+arrow_len-11*scale,y-5*scale),(x+arrow_len-11*scale,y+5*scale)],fill=BLUE)
    d.text((x+arrow_len+14*scale,y),'2 µN',font=font(27,scale),anchor='lm',fill=BLUE)
    d.ellipse((78*scale,y-6*scale,90*scale,y+6*scale),fill=RED)
    d.text((104*scale,y),'Clot',font=font(27,scale),anchor='lm',fill='#28333b')
    if not field:
        d.rectangle((270*scale,y-7*scale,290*scale,y+7*scale),fill='#f18b80')
        d.text((304*scale,y),'Loaded surface',font=font(27,scale),anchor='lm',fill='#28333b')
    else:
        x0,x1,y0,y1=int(width*.28),int(width*.56),int(height-60*scale),int(height-42*scale)
        grad=(colormaps['viridis'](np.linspace(0,1,x1-x0))[:,:3]*255).astype(np.uint8)
        image.paste(Image.fromarray(np.repeat(grad[None,:,:],y1-y0,axis=0)),(x0,y0))
        d=ImageDraw.Draw(image)
        for frac in [0,.5,1]:
            d.text((x0+(x1-x0)*frac,y1+5*scale),f'{frac*stats["color_limit_Pa"]:.0f}',font=font(22,scale),anchor='mt',fill='#28333b')
        d.text(((x0+x1)/2,y0-8*scale),'|τ| (Pa)',font=font(25,scale),anchor='mb',fill='#28333b')
    return np.asarray(image)


def load_state(path):
    with np.load(path) as f:
        return {k:f[k].copy() for k in old.FIELDS if k in f} | {'surface_traction':f['surface_traction'].copy()}


def bounds_for(points, pad=.08):
    lo,hi=points.min(0),points.max(0);margin=np.maximum((hi-lo)*pad,30.)
    return np.column_stack([lo-margin,hi+margin]).ravel()


def render_one(ref,surf,path,bounds,size=(W,H),field=False,selected=None):
    s=load_state(path);geom=surf.geometry(s['x'],s['d'])
    assert np.array_equal(geom['active'],s['surface_active'])
    pl=setup(bounds,size)
    stats=scene(pl,ref,surf,s,geom,selected,field,60 if field else 38)
    pl.render();raw=pl.screenshot(return_img=True)
    result=annotate(pl,raw,bounds,s,stats,field)
    stats.update(source=str(path.relative_to(R)),source_sha256=sha(path),bounds_um=np.asarray(bounds).tolist(),
                 camera_position=[list(v) for v in pl.camera_position],parallel_scale=pl.camera.parallel_scale)
    pl.close();return result,stats


def static_figures(ref,surf):
    cases=[('FIG_V3_PEAK_LOAD_FIELD',O/'cases/biphasic/keys/frame_012543.npz',None),
           ('FIG_V3_L2_FRAGMENT_LOAD_FIELD',O/'cases/biphasic/keys/frame_014415.npz',46)]
    records=[]
    for name,path,target_nodes in cases:
        s=load_state(path);selected=None
        if target_nodes:
            lab,sizes,attached=old.graph.topology(ref,s['d'])
            labels=np.flatnonzero((sizes==target_nodes)&~attached)
            assert len(labels)==1
            selected=lab==labels[0]
        points=s['x']*1e6 if selected is None else s['x'][:len(selected)][selected]*1e6
        # Include the actual load-vector endpoints in the local viewport.
        g=surf.geometry(s['x'],s['d']);n=int(ref['clot_count'])
        if selected is not None:
            member=np.r_[selected,np.zeros(len(s['x'])-n,bool)]
            eligible=np.r_[s['surface_active'][:surf.L]&member[surf.base.tri].all(1),s['surface_active'][surf.L:]&member[surf.owner]]
            points=np.vstack([points,s['surface_centroid'][eligible]*1e6,s['surface_centroid'][eligible]*1e6+s['surface_force'][eligible]*FORCE_SCALE])
        bounds=bounds_for(points,.12)
        im,record=render_one(ref,surf,path,bounds,(3200,1800),True,selected)
        Image.fromarray(im).save(OUT/(name+'.png'),dpi=(300,300))
        Image.fromarray(im).save(OUT/(name+'.pdf'),resolution=300.)
        record.update(png=name+'.png',pdf=name+'.pdf',pdf_kind='embedded 3200x1800 raster at 300 dpi',
                      load_field='Saved actual surface_traction (Pa); arrows are saved surface_force (N), not prescribed stationary nominal traction.')
        records.append(record)
        eligible=np.array(record['active_surface_indices'],int)
        with (OUT/(name+'_DATA.csv')).open('w') as stream:
            writer=csv.writer(stream);writer.writerow(['surface_index','method','x_um','y_um','z_um','traction_x_Pa','traction_y_Pa','traction_z_Pa','force_x_N','force_y_N','force_z_N','effective_area_m2','arrow_displayed'])
            chosen=set(record['displayed_surface_indices'])
            for i in eligible:
                writer.writerow([i,'L1' if i<surf.L else 'L2',*(s['surface_centroid'][i]*1e6),*s['surface_traction'][i],*s['surface_force'][i],s['surface_area'][i]*s['surface_weight'][i],i in chosen])
        print('STATIC_COMPLETE',name,flush=True)
    dump(AUDIT/'STATIC_FIGURES.json',records)
    return records


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');args=parser.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);AUDIT.mkdir(parents=True,exist_ok=True)
    source_csv=F/'videos/VIDEO_01_SURFACE_COUPLING_BEFORE_AFTER_FRAME_TIMES.csv'
    with source_csv.open() as stream:timeline=list(csv.DictReader(stream))
    paths=list(dict.fromkeys(R/q['right_source'] for q in timeline))
    data_sources={str(p.relative_to(R)):sha(p) for p in paths}
    lo=np.full(3,np.inf);hi=-lo
    for path in paths:
        with np.load(path) as f:xyz=f['x']*1e6
        lo=np.minimum(lo,xyz.min(0));hi=np.maximum(hi,xyz.max(0))
    bounds=bounds_for(np.array([lo,hi]),.045)
    ref=load_npz(O/'provenance/reference.npz');surf=old.CurrentSurface(ref,'uniform')
    start=time.monotonic();records=[]
    if args.preview:
        samples=[timeline[i] for i in [0,88,112,125,160]]
        for index,row in enumerate(samples):
            im,record=render_one(ref,surf,R/row['right_source'],bounds)
            name=f'PREVIEW_{index:02d}.png';Image.fromarray(im).save(OUT/name)
            record['image']=name;records.append(record)
        figures=static_figures(ref,surf)
        dump(AUDIT/'PREVIEW.json',dict(status='GENERATED_FOR_REVIEW',frames=records,figures=figures))
        return
    assert json.loads((AUDIT/'PREVIEW_REVIEW.json').read_text())['status']=='PASS'
    command=['ffmpeg','-y','-loglevel','warning','-f','rawvideo','-vcodec','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}',
             '-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','slow','-crf','17','-pix_fmt','yuv420p','-movflags','+faststart',str(OUT/(NAME+'.mp4'))]
    cached_path=None;im=None;last_stats=None
    with (AUDIT/'ffmpeg.log').open('wb') as log:
        enc=subprocess.Popen(command,stdin=subprocess.PIPE,stderr=log)
        try:
            for index,row in enumerate(timeline):
                path=R/row['right_source']
                if path!=cached_path:
                    im,last_stats=render_one(ref,surf,path,bounds);cached_path=path
                enc.stdin.write(im.tobytes())
                records.append(dict(frame=index,playback_s=index/FPS,actual_time_s=last_stats['actual_time_s'],source=str(path.relative_to(R)),
                                    source_sha256=data_sources[str(path.relative_to(R))],hold=row['hold'],
                                    z_axis_horizontal_error_px=last_stats['z_axis_horizontal_error_px'],
                                    displayed_surface_indices=json.dumps(last_stats['displayed_surface_indices'],separators=(',',':'))))
                if index in [0,88,112,125,160]:Image.fromarray(im).save(OUT/f'{NAME}_FRAME_{index:03d}.png')
                if index%24==0:print('FRAME',index,len(timeline),flush=True)
            enc.stdin.close();assert enc.wait()==0
        finally:
            if enc.poll() is None:enc.kill();enc.wait()
    Image.fromarray(im).save(OUT/(NAME+'_COVER.png'))
    with (OUT/(NAME+'_FRAME_TIMES.csv')).open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-show_format','-of','json',str(OUT/(NAME+'.mp4'))],text=True))
    stream=probe['streams'][0]
    assert stream['width']==W and stream['height']==H and stream['avg_frame_rate']=='24/1' and int(stream['nb_frames'])==len(timeline)
    subprocess.run(['ffmpeg','-v','error','-xerror','-i',str(OUT/(NAME+'.mp4')),'-f','null','-'],check=True,capture_output=True)
    for rel,value in data_sources.items():assert sha(R/rel)==value,rel
    dump(AUDIT/'VALIDATION.json',dict(status='PASS',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        frames=len(records),fps=FPS,resolution=[W,H],video_sha256=sha(OUT/(NAME+'.mp4')),full_decode='PASS',
        source_frame_csv=str(source_csv.relative_to(R)),source_frame_csv_sha256=sha(source_csv),source_states=data_sources,
        source_states_unchanged=True,only_new_V3=True,physical_runs_started=0,coordinate_scale=1,
        camera_up=[0,0,1],max_z_axis_horizontal_error_px=max(q['z_axis_horizontal_error_px'] for q in records),
        fixed_bounds_um=bounds,force_arrow_scale_um_per_N=FORCE_SCALE,wall_s=time.monotonic()-start,
        script_sha256=sha(Path(__file__)),notes='Red material points and actual live bonds. All material points shown. No failed bonds drawn. Active facets use saved flags; L2 corners use frozen local deformation. No convex shell or trajectory interpolation. Partial approximate fluid surface remains uncalibrated. Selected vectors only; full static loads in CSV.'))
    print('VIDEO_COMPLETE',len(records),flush=True)


if __name__=='__main__':main()
