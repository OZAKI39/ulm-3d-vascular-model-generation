"""One camera orbit of the frozen, latest Y-clot reference discretization.

The model is never rotated/deformed; only the orthographic camera moves.
Native IDs, phase assignment and reference edges are read without modification.
"""
from pathlib import Path
from itertools import product
import argparse, datetime, hashlib, json, os, subprocess, time
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont, ImageColor

ROOT = Path(__file__).resolve().parents[3]
TASK = 'y_clot_irregular_current_method_v4'
RUN = '20261020T041100Z'
A = ROOT/'analysis'/TASK/'code'
O = ROOT/'verification'/TASK/RUN/'rotation'
F = ROOT/'results'/TASK/RUN/'rotation'
SOURCE = ROOT/'verification'/TASK/RUN
W, H, FPS, NFRAMES = 3840, 2160, 30, 540
VIEW = (530, 180, 2780, 1540)
ELEVATION = -16.
AZIMUTH = -68.
FADE_HALF = 12.
TICKS = [[-6.6, -6.0, -5.4], [-.5, 0., .5], [.7, 1., 1.3]]
C = dict(weak='#d94a5b', strong='#88152b', anchor='#8499a6',
         bond='#864451', interface='#81939d', grid='#697d89', text='#293a45')
FONT = '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
FONTS = {n: ImageFont.truetype(FONT, n) for n in [48, 50, 54, 56]}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024**2), b''): h.update(block)
    return h.hexdigest()


def dump(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=False)+'\n')


def camera(angle):
    az, el = np.deg2rad([AZIMUTH+angle, ELEVATION])
    out = np.array([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])
    right = np.cross(-out, [0., 0., 1.]); right /= np.linalg.norm(right)
    up = np.cross(right, -out)
    return out, right, up


def line_mesh(points, bonds):
    return pv.PolyData(points, lines=np.c_[np.full(len(bonds), 2), bonds].ravel())


def silhouette_axes(bounds, right, up):
    """True box edges, following the previous accepted rotation's algorithm."""
    corners = np.array(list(product(*bounds)))
    xy = np.c_[corners@right, corners@up]; center = xy.mean(0); answer = []
    for axis in range(3):
        others = [k for k in range(3) if k != axis]; candidates = []
        for vals in product(*[bounds[k] for k in others]):
            edge = np.zeros((2, 3)); edge[:, axis] = bounds[axis]
            for k, val in zip(others, vals): edge[:, k] = val
            q = np.c_[edge@right, edge@up]; d = q[1]-q[0]
            side = d[0]*(xy[:, 1]-q[0, 1])-d[1]*(xy[:, 0]-q[0, 0])
            if not (side.min() >= -1e-10 or side.max() <= 1e-10): continue
            normal = np.array([-d[1], d[0]]); normal /= np.linalg.norm(normal)
            if normal@(q.mean(0)-center) < 0: normal = -normal
            score = -q.mean(0)[0] if axis == 2 else -q.mean(0)[1]
            candidates.append((score, edge, normal))
        _, edge, normal = max(candidates, key=lambda z:z[0])
        answer.append((axis, edge, normal))
    return answer


class Axes:
    def __init__(self, bounds, focal, scale):
        self.bounds, self.focal, self.scale = bounds, focal, scale
        self.keys, self.states = [], {}
        for degree in range(360):
            _, right, up = camera(degree); keys = []
            for axis, edge, outward in silhouette_axes(bounds, right, up):
                xy = np.array([self.project(q, right, up) for q in edge])
                d = xy[1]-xy[0]; normal = np.array([-d[1], d[0]])/np.linalg.norm(d)
                sign = 1 if normal@(outward*[1., -1.]) > 0 else -1
                key = (axis, *edge.ravel(), sign)
                self.states[key] = (axis, edge, sign); keys.append(key)
            self.keys.append(keys)
        self.transitions = [[i for i in range(360) if self.keys[i][a] != self.keys[i-1][a]] for a in range(3)]

    def project(self, point, right, up):
        x, y, w, h = VIEW; q = point-self.focal
        return np.array([x+w/2+(q@right)*h/(2*self.scale), y+h/2-(q@up)*h/(2*self.scale)])

    def groups(self, angle):
        result = []
        for axis in range(3):
            for degree in self.transitions[axis]:
                delta = (angle-(degree-.5)+180) % 360-180
                if abs(delta) < FADE_HALF:
                    alpha = .5-.5*np.cos(np.pi*abs(delta)/FADE_HALF)
                    key = self.keys[(degree-1) % 360 if delta < 0 else degree][axis]
                    if alpha > 1e-10: result.append((key, self.states[key], float(alpha)))
                    break
            else:
                key = self.keys[int(np.floor(angle+.5)) % 360][axis]
                result.append((key, self.states[key], 1.))
        return result


def overlays(axes, angle):
    """Only axes and a single shared legend; no fictional time/load arrows."""
    _, right, up = camera(angle)
    layer = Image.new('RGBA', (W, H)); draw = ImageDraw.Draw(layer); records = []
    for key, (axis, edge, sign), alpha in axes.groups(angle):
        xy = np.array([axes.project(q, right, up) for q in edge])
        d = xy[1]-xy[0]; normal = sign*np.array([-d[1], d[0]])/np.linalg.norm(d)
        ink = (*ImageColor.getrgb(C['text']), round(255*alpha))
        draw.line([tuple(q) for q in xy], fill=ink, width=3)
        labels = []
        for value in TICKS[axis]:
            p = edge.mean(0); p[axis] = value; pos = axes.project(p, right, up)
            draw.line([tuple(pos), tuple(pos+normal*11)], fill=ink, width=3)
            labels.append((f'{value:.1f}' if axis == 0 else f'{value:g}', pos+normal*(142 if axis == 2 else 61), 50))
        if axis == 2:
            title_position = xy[np.argmin(xy[:, 1])]+[0., -70.]
        else:
            title_position = axes.project(edge.mean(0), right, up)+normal*218
        labels.append(('XYZ'[axis]+' (mm)', title_position, 56))
        for text, pos, size in labels:
            draw.text(tuple(pos), text, font=FONTS[size], fill=ink, anchor='mm')
            records.append(dict(key=str((key, text)), text=text, pos=pos.tolist(), alpha=alpha,
                                bbox=list(draw.textbbox(tuple(pos), text, font=FONTS[size], anchor='mm'))))
    items = [('dot','Weak phase',C['weak']), ('dot','Strong phase',C['strong']),
             ('dot','Fixed anchors',C['anchor']), ('line','Reference bonds',C['bond'])]
    widths = [62+draw.textlength(label, font=FONTS[54]) for _, label, _ in items]
    gap = 100; x = (W-sum(widths)-gap*(len(items)-1))/2; y = H-102
    for (kind, label, color), width in zip(items, widths):
        if kind == 'dot': draw.ellipse((x+9,y-12,x+33,y+12),fill=color)
        else: draw.line([(x,y),(x+43,y)],fill=color,width=4)
        draw.text((x+62,y),label,font=FONTS[54],fill=C['text'],anchor='lm')
        records.append(dict(key=label,text=label,pos=[x+62,y],alpha=1.,
                            bbox=list(draw.textbbox((x+62,y),label,font=FONTS[54],anchor='lm'))))
        x += width+gap
    return layer, records


def layout_checks(axes):
    corners = np.array(list(product(*axes.bounds))); traces = []; collisions = []
    maximum_z_error = 0.
    for frame in range(NFRAMES):
        angle = frame*360/NFRAMES; _, right, up = camera(angle)
        q = np.array([axes.project(p, right, up) for p in corners]); x,y,w,h = VIEW
        assert (q[:,0] > x).all() and (q[:,0] < x+w).all() and (q[:,1] > y).all() and (q[:,1] < y+h).all()
        z0 = axes.project(axes.focal, right, up); z1 = axes.project(axes.focal+[0,0,.2],right,up)
        maximum_z_error = max(maximum_z_error, abs(float(z1[0]-z0[0])))
        assert z1[1] < z0[1]
        _, labels = overlays(axes, angle)
        for r in labels:
            a = r['bbox']; assert a[0] > 0 and a[1] > 0 and a[2] < W and a[3] < H, (frame,r)
        for i,r in enumerate(labels):
            for s in labels[i+1:]:
                if min(r['alpha'],s['alpha']) < .15: continue
                a,b = r['bbox'],s['bbox']
                if min(a[2],b[2])-max(a[0],b[0]) > 2 and min(a[3],b[3])-max(a[1],b[1]) > 2:
                    collisions.append(dict(frame=frame,labels=[r['text'],s['text']]))
        traces.append({r['key']:r for r in labels})
    assert not collisions, collisions[:12]
    max_alpha, max_shift = 0., 0.
    for i, old in enumerate(traces):
        new = traces[(i+1)%NFRAMES]
        for key in old.keys() | new.keys():
            a,b = old.get(key),new.get(key)
            max_alpha = max(max_alpha, abs((a['alpha'] if a else 0)-(b['alpha'] if b else 0)))
            if a and b and min(a['alpha'],b['alpha']) > .02:
                max_shift = max(max_shift,float(np.linalg.norm(np.array(a['pos'])-b['pos'])))
    assert maximum_z_error < 1e-8 and max_alpha < .1 and max_shift < 45
    return dict(status='PASS',frames_checked=NFRAMES,visible_text_collisions=collisions,
                clipping_count=0,max_Z_horizontal_error_px=maximum_z_error,
                max_frame_label_opacity_change=max_alpha,max_visible_label_step_px=max_shift,
                transitions_degrees=axes.transitions,periodic_loop_checked=True)


def build_scene(ref, field):
    xyz = ref['X']*1000; n = int(ref['clot_count']); b = ref['bonds']; bt = ref['bond_types']
    bounds = np.column_stack([xyz.min(0)-.10, xyz.max(0)+.10]); focal = bounds.mean(1)
    corners = np.array(list(product(*bounds)))-focal
    scale = 0.
    for angle in np.arange(0.,360.,.5):
        _,right,up = camera(angle)
        scale=max(scale,np.ptp(corners@up)/2,np.ptp(corners@right)/2*VIEW[3]/VIEW[2])
    scale *= 1.07
    p = pv.Plotter(off_screen=True,window_size=VIEW[2:],lighting='three lights'); p.set_background('white')
    p.render_window.SetMultiSamples(0);p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=8,occlusion_ratio=.01)
    grid=[]
    for plane in range(3):
        other=[k for k in range(3) if k!=plane]
        for fixed in bounds[plane]:
            for varying in other:
                constant=next(k for k in other if k!=varying)
                for value in np.unique([*bounds[constant],*TICKS[constant]]):
                    edge=np.zeros((2,3));edge[:,plane]=fixed;edge[:,constant]=value;edge[:,varying]=bounds[varying]
                    grid.append(edge)
    grid=np.array(list({tuple(e.ravel()):e for e in grid}.values()))
    p.add_mesh(line_mesh(grid.reshape(-1,3),np.arange(2*len(grid)).reshape(-1,2)),
               color=C['grid'],line_width=1.5,opacity=.68,lighting=False)
    counts={}
    for typ,color in [(1,C['bond']),(2,C['interface'])]:
        selected=bt==typ
        p.add_mesh(line_mesh(xyz,b[selected]),color=color,opacity=.045,line_width=1.1,lighting=False)
        counts[str(typ)]=dict(all_reference_edges=int(selected.sum()),short_edges_emphasized=0)
    for ids,color,size,opacity in [(np.arange(n,len(xyz)),C['anchor'],8.,.42),
                                   (np.flatnonzero(field['material_phase']==0),C['weak'],14.,.97),
                                   (np.flatnonzero(field['material_phase']==1),C['strong'],14.,.99)]:
        points=pv.PolyData(xyz[ids]);points['radius']=np.cbrt(3*ref['V'][ids]/(4*np.pi))*1000
        glyph=points.glyph(geom=pv.Sphere(radius=1.,theta_resolution=12,phi_resolution=12),scale='radius',orient=False)
        p.add_mesh(glyph,color=color,opacity=.28 if ids[0]>=n else .82,smooth_shading=True,
                   ambient=.45,diffuse=.55,specular=.10,specular_power=18)
    p.enable_parallel_projection();p.camera.parallel_scale=scale
    out,_,up=camera(0);p.camera_position=[focal+20*out,focal,up];p.camera.parallel_scale=scale
    p.camera.clipping_range=(.1,50.)
    faces=[dict(axis='XYZ'[a],side=s,segments=int(np.isclose(grid[:,:,a],value).all(1).sum()))
           for a in range(3) for s,value in enumerate(bounds[a])]
    assert all(f['segments']>=8 for f in faces)
    return p,Axes(bounds,focal,scale),dict(reference_bond_display=counts,grid_faces=faces,
        fixed_scaffold_edges_not_displayed=int((bt==3).sum()),bounds_mm=bounds.tolist(),parallel_scale_mm=scale)


def main(preview):
    start=time.perf_counter();F.mkdir(parents=True,exist_ok=True);O.mkdir(parents=True,exist_ok=True)
    inputs=[SOURCE/'provenance/reference.npz',SOURCE/'BIPHASIC_FIELD.npz',SOURCE/'BIPHASIC_MATERIAL_SPEC.json',
            SOURCE/'cases/H_Y_IRREGULAR_LPS_V4/keys/initial.npz',SOURCE/'SOURCE_FREEZE.json',
            ROOT/'analysis/clot_fragment_surface_coupling_v3_red_visuals/render.py',
            ROOT/'analysis/clot_loading_publication_v1/render_z_rotation_v2.py']
    before={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    with np.load(inputs[0]) as q:ref=dict(q)
    with np.load(inputs[1]) as q:field=dict(q)
    with np.load(inputs[3]) as q:state=dict(q)
    n=int(ref['clot_count']);assert n==1927 and len(ref['X'])-n==1561
    assert float(state['time'])==0 and np.array_equal(state['x'],ref['X']) and not state['d'].any()
    assert np.array_equal(field['atom_ids'],state['atom_ids'])
    phase=field['material_phase'];assert int((phase==1).sum())==882 and int((phase==0).sum())==1045
    p,axes,scene=build_scene(ref,field)
    qa=layout_checks(axes);dump(O/'LAYOUT_VALIDATION.json',qa);print('LAYOUT_PASS',qa,flush=True)
    if preview:frames=[0,68,135,203,270,338,405,473]
    else:frames=range(NFRAMES)
    video=F/'VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_4K.mp4';temp=video.with_name(video.stem+'.rendering.mp4')
    encoder=None;camera_path=[];stills=[]
    if not preview:
        assert not video.exists(),'Use a new directory for a new finalized version.'
        encoder=subprocess.Popen(['ffmpeg','-y','-v','error','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}',
            '-r',str(FPS),'-i','-','-an','-c:v','libx264','-crf','17','-preset','medium','-threads','4',
            '-pix_fmt','yuv420p','-movflags','+faststart',str(temp)],stdin=subprocess.PIPE)
    try:
        for frame in frames:
            angle=frame*360/NFRAMES;out,right,up=camera(angle)
            p.camera_position=[axes.focal+20*out,axes.focal,up];p.camera.parallel_scale=axes.scale;p.render()
            raw=p.screenshot(return_img=True);im=Image.new('RGBA',(W,H),'white')
            im.paste(Image.fromarray(raw).convert('RGBA'),VIEW[:2]);layer,_=overlays(axes,angle);im.alpha_composite(layer);im=im.convert('RGB')
            if encoder:encoder.stdin.write(np.asarray(im,dtype=np.uint8).tobytes())
            if preview or frame in [0,68,135,203,270,338,405,473]:
                im.save(F/f'VIEW_{frame:03d}.png');thumb=im.copy();thumb.thumbnail((960,540));stills.append(thumb)
            camera_path.append(dict(frame=frame,video_time_s=frame/FPS,orbit_angle_deg=angle,
                eye_mm=(axes.focal+20*out).tolist(),focal_mm=axes.focal.tolist(),up=up.tolist(),physical_time_s=0.))
            if preview or (frame+1)%30==0:print('RENDER',frame+1,NFRAMES,round(time.perf_counter()-start,1),flush=True)
        if encoder:
            encoder.stdin.close();assert encoder.wait()==0;temp.replace(video)
    finally:
        p.close()
        if encoder and encoder.poll() is None:encoder.kill();encoder.wait()
    sheet=Image.new('RGB',(3840,1080),'white')
    for i,im in enumerate(stills):sheet.paste(im,((i%4)*960,(i//4)*540))
    sheet.save(F/'ROTATION_CONTACT_SHEET.jpg',quality=95)
    assert before=={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    record=dict(status='PREVIEW' if preview else 'PASS',task=TASK,run=RUN,input_sha256=before,
        renderer_sha256=sha(Path(__file__)),source_case=str(SOURCE.relative_to(ROOT)),
        physical_time_s=0.,camera_orbit_only=True,geometry_and_damage_unchanged=True,new_physics_runs=0,
        material_points=n,fixed_anchors=1561,strong_points=882,weak_points=1045,strong_fraction=882/n,
        material_model='V4 variational LPS with nonuniform control volumes and original-ID heterogeneous irreversible softening; fixed boundary and jittered interior',nodal_marker='Sphere volume equals actual Vi; discrete quadrature visualization, not a continuous clot shell',
        dx_um=50.,horizon_um=150.5,PD_quadrature_volume_mm3=float(ref['V'][:n].sum()*1e9),
        no_detached_nodes_in_initial_state=True,no_applied_force_at_initial_state=True,
        resolution=[W,H],fps=FPS,duration_s=NFRAMES/FPS,frames=NFRAMES,orbit_degrees=360.,
        elevation_degrees=ELEVATION,loop='Uniform periodic camera orbit; no duplicate endpoint frame',
        colors=C,scene=scene,layout=qa,render_seconds=time.perf_counter()-start,
        overlay_text='Axes, coordinate values and one legend only; no load scale or deformation-time animation')
    if not preview:
        for target,extra in [(video,None),(F/'VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_1080P.mp4',['-vf','scale=1920:1080:flags=lanczos'])]:
            if extra:subprocess.run(['ffmpeg','-y','-v','error','-i',str(video),*extra,'-an','-c:v','libx264','-crf','18','-preset','medium','-threads','4','-pix_fmt','yuv420p','-movflags','+faststart',str(target)],check=True)
            subprocess.run(['ffmpeg','-v','error','-i',str(target),'-f','null','-'],check=True,capture_output=True)
        record['full_video_decode']='PASS_BOTH_FILES';dump(O/'CAMERA_PATH.json',camera_path)
        (F/'OPEN_RESULTS.html').write_text('''<!doctype html><html lang="zh"><meta charset="utf-8"><title>最新血块离散模型：Z轴旋转</title><body style="max-width:1300px;margin:30px auto;font-family:sans-serif;line-height:1.7"><h1>最新 Y 血块初态：360° Z 轴旋转</h1><p>相机环绕固定初态，非断裂时程。1927 血块点、1561 固定锚点；深浅红色对应真实强弱相。V4内部点云不规则、控制体积非均匀；固定外边界和锚点，按新参考坐标重建材料图并排除穿越空隙的连接。所有本体和界面键淡线呈现；固定支架内部键仅在图中省略。球形标记体积等于Vi，仅表示离散求积。</p><video controls loop style="width:100%" src="VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_4K.mp4"></video><p><a href="VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_4K.mp4">4K MP4</a> · <a href="VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_1080P.mp4">1080p PPT版</a> · <a href="VIEW_000.png">静态封面</a></p><img style="width:100%" src="ROTATION_CONTACT_SHEET.jpg"><p>Δx=50 µm；δ=150.5 µm；PD积分体积0.240875 mm³。本视频仅环绕新计算的原生初态，不是断裂时间演化；实际3.6ms计算见上级结果索引。</p></body></html>''')
    dump(O/('PREVIEW_VALIDATION.json' if preview else 'VISUALIZATION_VALIDATION.json'),record)
    print('COMPLETE',record['status'],F,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');args=parser.parse_args();main(args.preview)
