"""Restyle accepted Y-clot V4 data. No native integration or physical edits."""
from pathlib import Path
import sys, json, hashlib, subprocess, argparse, time, html
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from itertools import product
import reference_style as style

R = Path(__file__).resolve().parents[3]
TASK = 'y_clot_v4_reference_style_v1'
RUN = '20261020T084000Z'
A = R/'analysis'/TASK
O = R/'verification'/TASK/RUN
F = R/'results'/TASK/RUN
SOURCE = R/'verification/y_clot_irregular_current_method_v4/20261020T041100Z'
SOURCE_CODE = R/'analysis/y_clot_irregular_current_method_v4/code'
REFERENCE_VIDEO = R/'results/clot_discretization_current_method_rotation_v1/20261019T205403Z/VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_1080P.mp4'
sys.path.insert(0, str(SOURCE_CODE))
from analyze_result import state_paths, topology, table, CASE, OLD
from common import load_npz, dump, sha

W, H = style.W, style.H
C = style.C
FONT = style.FONT
FIG = F/'figures'
VID = F/'videos'
REV = F/'visual_review'
DEFAULT_TICKS = [q.copy() for q in style.TICKS]


def inputs():
    paths = state_paths(CASE)
    result = json.loads((SOURCE/'NATIVE_V4_RESULT.json').read_text())
    assert result['status'] == 'COMPLETED'
    paths = [(t,p) for t,p in paths if t <= result['accepted_time_s']+1e-14]
    ref = load_npz(SOURCE/'provenance/reference.npz')
    field = load_npz(SOURCE/'BIPHASIC_FIELD.npz')
    n = int(ref['clot_count'])
    assert n == 1927 and len(ref['X']) == 3488
    assert len(ref['bonds']) == 153986
    assert np.array_equal(load_npz(paths[0][1])['x'], ref['X'])
    assert np.array_equal(load_npz(paths[0][1])['atom_ids'], field['atom_ids'])
    return ref, field, paths


def freeze(paths):
    files = [p for _,p in paths]
    files += [SOURCE/'provenance/reference.npz', SOURCE/'BIPHASIC_FIELD.npz',
              SOURCE/'ACTUAL_LOAD_FIELD.npz', SOURCE/'NATIVE_V4_RESULT.json',
              SOURCE/'SCIENTIFIC_OUTCOME.json', SOURCE/'DYNAMICS_VIDEO_PROVENANCE.json',
              SOURCE/'FRAGMENT_SIZE_SPECTRA.json', CASE/'history.csv', OLD/'history.csv',
              REFERENCE_VIDEO]
    files += list(SOURCE_CODE.glob('*.py')) + list(SOURCE_CODE.glob('*.cpp'))
    frozen = {str(p.relative_to(R)):sha(p) for p in files}
    out = O/'SOURCE_HASHES.json'
    if out.exists():
        assert json.loads(out.read_text()) == frozen, 'Source changed'
    else:
        dump(out, frozen)


def ticks(bounds, initial=False):
    if initial:
        style.TICKS = [q.copy() for q in DEFAULT_TICKS]
    else:
        from matplotlib.ticker import MaxNLocator
        style.TICKS = []
        for lo,hi in bounds:
            values = MaxNLocator(nbins=3, steps=[1,2,2.5,5,10]).tick_values(lo,hi)
            values = values[(values > lo+.04*(hi-lo)) & (values < hi-.04*(hi-lo))]
            if len(values) < 2:
                values = np.round(np.linspace(lo+.2*(hi-lo),hi-.2*(hi-lo),3),1)
            style.TICKS.append(values.tolist())


def new_scene(bounds, initial=False):
    ticks(bounds, initial)
    focal = bounds.mean(1)
    corners = np.array(list(product(*bounds)))-focal
    scale = 0.
    for angle in np.arange(0.,360.,.5):
        _,right,up = style.camera(angle)
        scale = max(scale, np.ptp(corners@up)/2,
                    np.ptp(corners@right)/2*style.VIEW[3]/style.VIEW[2])
    scale *= 1.07
    p = pv.Plotter(off_screen=True, window_size=style.VIEW[2:], lighting='three lights')
    p.set_background('white')
    p.render_window.SetMultiSamples(0)
    p.render_window.SetAlphaBitPlanes(1)
    p.enable_depth_peeling(number_of_peels=8, occlusion_ratio=.01)
    grid = []
    for axis in range(3):
        others = [k for k in range(3) if k != axis]
        for side in bounds[axis]:
            for along in others:
                across = next(k for k in others if k != along)
                for value in np.unique([*bounds[across], *style.TICKS[across]]):
                    edge = np.zeros((2,3))
                    edge[:,axis] = side
                    edge[:,across] = value
                    edge[:,along] = bounds[along]
                    grid.append(edge)
    grid = np.array(list({tuple(e.ravel()):e for e in grid}.values()))
    p.add_mesh(style.line_mesh(grid.reshape(-1,3),np.arange(2*len(grid)).reshape(-1,2)),
               color=C['grid'], line_width=1.5, opacity=.68, lighting=False)
    faces = [dict(axis='XYZ'[a],side=s,segments=int(np.isclose(grid[:,:,a],v).all(1).sum()))
             for a in range(3) for s,v in enumerate(bounds[a])]
    assert all(q['segments'] >= 8 for q in faces)
    p.enable_parallel_projection()
    axes = style.Axes(bounds,focal,scale)
    return p, axes, dict(bounds_mm=bounds,parallel_scale_mm=scale,grid_faces=faces)


def draw_state(p, ref, field, state, opacity=1.):
    x = state['x']*1000
    n = int(ref['clot_count'])
    assert np.array_equal(state['atom_ids'],field['atom_ids'])
    assert np.array_equal(state['x'][n:],ref['X'][n:])
    bonds, types = ref['bonds'],ref['bond_types']
    alive = state['d'] < 1.
    counts = {}
    for typ,color in [(1,C['bond']),(2,C['interface'])]:
        selected = alive & (types == typ)
        short = selected & (ref['r0'] <= 1.01*50e-6)
        for name,mask,alpha,width in [('all',selected,.025,1.),
                                     ('short',short,.36 if typ==1 else .24,1.6)]:
            actor = f'edges_{typ}_{name}'
            if mask.any():
                p.add_mesh(style.line_mesh(x,bonds[mask]),color=color,
                           opacity=alpha*opacity,line_width=width,lighting=False,
                           name=actor,reset_camera=False)
            else:
                p.remove_actor(actor,reset_camera=False)
        counts[str(typ)] = dict(shown_surviving_edges=int(selected.sum()),
                               emphasized_short_edges=int(short.sum()))
    for name,ids,color,size,alpha in [
        ('anchors',np.arange(n,len(x)),C['anchor'],8.,.42),
        ('weak',np.flatnonzero(field['material_phase']==0),C['weak'],14.,.97),
        ('strong',np.flatnonzero(field['material_phase']==1),C['strong'],14.,.99)]:
        p.add_points(x[ids],color=color,point_size=size,opacity=alpha*opacity,
                     render_points_as_spheres=True,ambient=.45,diffuse=.55,
                     specular=.10,specular_power=18,name=name,reset_camera=False)
    lab, attached = topology(ref,state['d'])
    sizes = np.bincount(lab)
    return dict(time_s=float(state['time']),step=int(state['step']),
                free_nodes=int((~attached[lab[:n]]).sum()),
                free_components=int((~attached).sum()),
                N20_200=int((~attached&(sizes>=20)&(sizes<=200)).sum()),
                bonds=counts,
                missing_or_invented_nodes=0,broken_edges_drawn=0)


def frame(p, axes, angle=0., physical_time=None, load=None):
    out,right,up = style.camera(angle)
    p.camera_position = [axes.focal+50*out,axes.focal,up]
    p.camera.parallel_scale = axes.scale
    p.camera.clipping_range = (.1,100.)
    p.render()
    im = Image.new('RGBA',(W,H),'white')
    im.paste(Image.fromarray(p.screenshot(return_img=True)).convert('RGBA'),style.VIEW[:2])
    layer, labels = style.overlays(axes,angle)
    # Preserve the reference legend's visual positions, but describe current d<1 edges.
    draw = ImageDraw.Draw(layer)
    bond_label = next(q for q in labels if q['text']=='Reference bonds')
    draw.rectangle(bond_label['bbox'], fill=(255,255,255,255))
    draw.text(tuple(bond_label['pos']), 'Surviving bonds',font=style.FONTS[54],
              fill=C['text'],anchor='lm')
    im.alpha_composite(layer)
    draw = ImageDraw.Draw(im)
    if physical_time is not None:
        draw.text((W-270,62),f't = {physical_time*1e3:.5f} ms',
                  font=style.FONTS[50],fill=C['text'],anchor='rt')
    if load is not None:
        import matplotlib
        gradient=(matplotlib.colormaps['YlOrRd'](np.linspace(0,1,500))[:,:3]*255).astype(np.uint8)
        im.paste(Image.fromarray(np.repeat(gradient[None,:,:],24,axis=0)),(630,60))
        draw.text((610,72),'0',font=style.FONTS[48],fill=C['text'],anchor='rm')
        draw.text((1155,72),f'{load["max_Pa"]:.1f} Pa',font=style.FONTS[48],fill=C['text'],anchor='lm')
        # A legend arrow along the screen plane, with the same force/length scale as the 3D arrows.
        length_px=2e-6*load['scale_mm_per_N']*style.VIEW[3]/(2*axes.scale)
        x,y=1630,80
        draw.line([(x,y),(x+length_px,y)],fill=C['text'],width=5)
        draw.polygon([(x+length_px,y),(x+length_px-20,y-10),(x+length_px-20,y+10)],fill=C['text'])
        draw.text((x+length_px+24,y),'2 µN',font=style.FONTS[48],fill=C['text'],anchor='lm')
    z0=axes.project(axes.focal,right,up)
    z1=axes.project(axes.focal+[0,0,.2],right,up)
    assert abs(z1[0]-z0[0]) < 1e-8 and z1[1]<z0[1]
    for q in labels:
        a=q['bbox']
        assert a[0]>0 and a[1]>0 and a[2]<W and a[3]<H
    return im.convert('RGB')


def bounds_for(points):
    return np.column_stack((points.min(0)-.1,points.max(0)+.1))


def encoder(path, fps):
    assert not path.exists(),f'Finalized output already exists: {path}'
    return subprocess.Popen(['ffmpeg','-y','-v','error','-f','rawvideo','-pix_fmt','rgb24',
        '-s',f'{W}x{H}','-r',str(fps),'-i','-','-an','-c:v','libx264','-crf','17',
        '-preset','medium','-threads','4','-pix_fmt','yuv420p','-movflags','+faststart',str(path)],
        stdin=subprocess.PIPE)


def downsample(path):
    target = path.with_name(path.name.replace('_4K','_1080P'))
    assert not target.exists()
    subprocess.run(['ffmpeg','-v','error','-i',str(path),'-vf','scale=1920:1080:flags=lanczos',
                    '-c:v','libx264','-crf','17','-preset','medium','-threads','4',
                    '-pix_fmt','yuv420p','-movflags','+faststart',str(target)],check=True)


def orbit(ref,field,paths,preview):
    s=load_npz(paths[0][1])
    p,axes,info=new_scene(bounds_for(s['x']*1e3),True)
    info['state']=draw_state(p,ref,field,s)
    info['layout']=style.layout_checks(axes)
    dump(O/'ORBIT_LAYOUT.json',info)
    frames=[0,68,135,203,270,338,405,473] if preview else range(540)
    path=VID/'VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_4K.mp4'
    enc=None if preview else encoder(path,30)
    records=[];thumbs=[]
    try:
        for i in frames:
            angle=i*360/540
            im=frame(p,axes,angle)
            if enc:enc.stdin.write(np.asarray(im).tobytes())
            if i in [0,68,135,203,270,338,405,473]:
                im.save(FIG/f'ROTATION_VIEW_{i:03d}.png')
                thumb=im.copy();thumb.thumbnail((960,540));thumbs.append(thumb)
            records.append(dict(frame=i,camera_angle_deg=angle,physical_time_s=0.))
            if i%30==0:print('ORBIT',i,flush=True)
        if enc:enc.stdin.close();assert enc.wait()==0
    finally:
        p.close()
        if enc and enc.poll() is None:enc.kill()
    sheet=Image.new('RGB',(3840,1080),'white')
    for i,thumb in enumerate(thumbs):sheet.paste(thumb,((i%4)*960,(i//4)*540))
    sheet.save(REV/'ORBIT_CONTACT_SHEET.jpg',quality=95)
    dump(O/'ORBIT_FRAMES.json',records)
    if not preview:downsample(path)


def stills(ref,field,paths):
    log=json.loads((SOURCE/'DYNAMICS_VIDEO_PROVENANCE.json').read_text())
    records=[]
    targets=[('FIG_02_INITIAL_CLOT.png',paths[0][1]),
             ('FIG_02_LARGE_FRAGMENT_CLOT.png',R/log['large_fragment_picture_source']),
             ('FIG_02_FINAL_CLOT.png',paths[-1][1])]
    for name,path in targets:
        s=load_npz(path)
        p,axes,info=new_scene(bounds_for(s['x']*1e3),float(s['time'])==0)
        info['state']=draw_state(p,ref,field,s)
        im=frame(p,axes,physical_time=float(s['time']) if float(s['time']) else None)
        im.save(FIG/name);p.close()
        records.append(dict(file=name,source=str(path.relative_to(R)),source_sha256=sha(path),**info))
    # Read the already verified actual load field; no force or solver re-evaluation.
    z=load_npz(SOURCE/'ACTUAL_LOAD_FIELD.npz')
    path=R/log['load_picture_source'];s=load_npz(path)
    assert all(np.array_equal(z[k],s[k]) for k in ['time','x','d','velocity'])
    force_err=float(np.max(abs(z['nodal_force']+z['wall_force']-s['external_force'])))
    assert force_err<1e-15
    ids=np.flatnonzero(z['active']&(np.linalg.norm(z['face_force'],axis=1)>0))
    chosen=z['shown_arrow_ids'];force=z['face_force'][chosen]
    maxforce=float(np.linalg.norm(z['face_force'][ids],axis=1).max())
    scale=.28/maxforce
    magnitude=np.linalg.norm(z['traction'][ids],axis=1)
    points=np.vstack((s['x']*1e3,z['centroid'][chosen]*1e3+scale*force))
    p,axes,info=new_scene(bounds_for(points))
    info['state']=draw_state(p,ref,field,s,opacity=.55)
    p.add_points(z['centroid'][ids]*1e3,scalars=magnitude,cmap='YlOrRd',
                 clim=(0,float(magnitude.max())),point_size=14,
                 render_points_as_spheres=True,show_scalar_bar=False,lighting=False)
    p.add_arrows(z['centroid'][chosen]*1e3,force,mag=scale,color=C['text'],show_scalar_bar=False)
    frame(p,axes,physical_time=float(s['time']),
          load=dict(max_Pa=float(magnitude.max()),scale_mm_per_N=scale)).save(FIG/'FIG_03_ACTUAL_SURFACE_LOAD.png')
    p.close()
    info.update(source=str(path.relative_to(R)),actual_face_force_ids=chosen,
                force_scale_mm_per_N=scale,force_readback_max_N=force_err,
                maximum_traction_Pa=float(magnitude.max()),
                source_field_sha256=sha(SOURCE/'ACTUAL_LOAD_FIELD.npz'))
    dump(O/'LOAD_FIELD_VISUALIZATION.json',info)
    dump(O/'STATIC_SOURCES.json',records)
    overview(ref,field)


def overview(ref,field):
    from current_surface import CurrentSurface
    wall=CurrentSurface(ref,'uniform').wall.wall.copy()
    wall.points*=1e3
    p=pv.Plotter(off_screen=True,window_size=(W,H))
    p.set_background('white')
    p.add_mesh(wall,color='#b7707c',opacity=.10,smooth_shading=True)
    n=int(ref['clot_count'])
    for phase,color in [(0,C['weak']),(1,C['strong'])]:
        ids=np.flatnonzero(field['material_phase'][:n]==phase)
        p.add_points(ref['X'][ids]*1e3,color=color,point_size=6,render_points_as_spheres=True)
    center=np.array(wall.center);span=np.ptp(wall.points,axis=0).max()
    out,_,up=style.camera(0)
    p.camera_position=[center+50*out,center,up]
    p.enable_parallel_projection();p.camera.parallel_scale=span*.36
    im=Image.fromarray(p.screenshot(return_img=True)).convert('RGB')
    draw=ImageDraw.Draw(im)
    pixels=5*H/(2*p.camera.parallel_scale)
    draw.line([(240,H-220),(240+pixels,H-220)],fill=C['text'],width=5)
    draw.text((240+pixels/2,H-185),'5 mm',font=style.FONTS[50],fill=C['text'],anchor='mt')
    draw.text((W-260,H-180),'Clot in original Y vessel',font=style.FONTS[54],fill=C['text'],anchor='rt')
    im.save(FIG/'FIG_01_Y_VESSEL_LOCATION.png');p.close()


def dynamics(ref,field,paths):
    bounds=np.array([[np.inf,-np.inf]]*3)
    for _,path in paths:
        x=load_npz(path)['x']*1e3
        bounds[:,0]=np.minimum(bounds[:,0],x.min(0))
        bounds[:,1]=np.maximum(bounds[:,1],x.max(0))
    bounds[:,0]-=.1;bounds[:,1]+=.1
    p,axes,info=new_scene(bounds)
    path=VID/'VIDEO_01_Y_CLOT_V4_FRAGMENTATION_4K.mp4'
    enc=encoder(path,30);records=[]
    try:
        for i,(t,source) in enumerate(paths):
            s=load_npz(source);stats=draw_state(p,ref,field,s)
            im=frame(p,axes,physical_time=t)
            # Keep saved states exactly; holds are presentation timing, not integration.
            repeat=30 if i in [0,len(paths)-1] else 3
            for _ in range(repeat):enc.stdin.write(np.asarray(im).tobytes())
            if i in [0,len(paths)//2,len(paths)-1]:im.save(REV/f'DYNAMICS_{i:03d}.png')
            records.append(dict(source=str(source.relative_to(R)),sha256=sha(source),
                                repeat=repeat,**stats))
            if i%10==0:print('DYNAMICS',i,len(paths),flush=True)
        enc.stdin.close();assert enc.wait()==0
    finally:
        p.close()
        if enc.poll() is None:enc.kill()
    downsample(path)
    dump(O/'DYNAMICS_FRAMES.json',dict(fps=30,frames=records,**info,
         actual_saved_states=True,no_interpolation=True,displacement_scale=1,
         phase_colors_preserved=True,markers_are_visual_not_volume_spheres=True))


def plots():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':17,'axes.labelsize':19,
        'legend.fontsize':16,'axes.grid':True,'grid.color':C['grid'],'grid.alpha':.45,
        'axes.edgecolor':C['text'],'text.color':C['text'],'axes.labelcolor':C['text'],
        'xtick.color':C['text'],'ytick.color':C['text'],'savefig.facecolor':'white'})
    histories={'Y regular':table(OLD/'history.csv'),'Y V4':table(CASE/'history.csv')}
    fig,ax=plt.subplots(2,2,figsize=(18,11),layout='constrained')
    for j,(key,rows) in enumerate(histories.items()):
        t=np.array([float(q['time_s']) for q in rows])*1e3
        for a,column,label in [(ax[0,0],'N20_200','20–200-node free components'),
              (ax[0,1],'largest_nodes','Largest free component (nodes)'),
              (ax[1,0],'actual_peak_slip_shear_Pa','Peak slip traction (Pa)'),
              (ax[1,1],'effective_loaded_area_m2','Effective loaded area (mm²)')]:
            if column in rows[0]:
                y=np.array([float(q[column]) for q in rows])*(1e6 if column.endswith('m2') else 1)
                a.plot(t,y,color=[C['anchor'],C['strong']][j],lw=2,label=key)
            a.set(xlabel='Physical time (ms)',ylabel=label)
    ax[0,0].legend(frameon=False)
    fig.savefig(FIG/'FIG_04_LOAD_AND_FRAGMENT_HISTORY.png',dpi=220);plt.close(fig)
    rows=histories['Y V4'];t=np.array([float(q['time_s']) for q in rows])*1e3
    fig,ax=plt.subplots(2,1,figsize=(18,10),layout='constrained')
    palette=[C['strong'],C['weak'],C['anchor'],'#cc8a95','#535b69','#b44c68','#a6b0b8','#50434e']
    columns=['W_drive_J','W_drag_J','W_wall_J','U_J','K_J','D_J','Uc_J','Jc_J']
    for i,key in enumerate(columns):
        ax[0].plot(t,np.array([float(q[key]) for q in rows])*1e6,label=key.replace('_J',''),
                   color=palette[i],ls=['-','--','-.',':'][i%4],lw=2)
    ax[0].set(ylabel='Energy (µJ)')
    ax[0].legend(ncol=4,frameon=False)
    ax[1].plot(t,[100*float(q['relative_R']) for q in rows],color=C['strong'],lw=2)
    ax[1].axhline(1,color=C['text'],ls='--',label='1% numerical gate')
    ax[1].set(xlabel='Physical time (ms)',ylabel='Normalized energy residual (%)')
    ax[1].legend(frameon=False)
    fig.savefig(FIG/'FIG_05_ENERGY_LEDGER.png',dpi=220);plt.close(fig)
    spectra=json.loads((SOURCE/'FRAGMENT_SIZE_SPECTRA.json').read_text())
    fig,ax=plt.subplots(1,2,figsize=(20,7.5),layout='constrained')
    for j,(name,q) in enumerate(spectra.items()):
        for k,column,scale in [(0,'counts',1),(1,'mass_kg',1e6)]:
            vals=np.array(q[column])*scale
            bars=ax[k].bar(np.arange(len(vals))+(j-.5)*.37,vals,width=.37,
                           color=[C['anchor'],C['strong']][j],label=name)
            ax[k].bar_label(bars,labels=[f'{v:.3g}' for v in vals],
                           fontsize=13,padding=5,rotation=90 if k else 0)
    for a in ax:
        a.set_xticks(np.arange(len(q['node_bins'])),q['node_bins'])
        a.set_xlabel('Original nodes per free component');a.margins(y=.27)
    ax[0].legend(frameon=False)
    ax[0].set_yscale('symlog',linthresh=1);ax[0].set_ylabel('Free component count (symlog)')
    ax[1].set_ylabel('Reference mass (mg)')
    fig.savefig(FIG/'FIG_06_FRAGMENT_SIZE_DISTRIBUTION.png',dpi=220);plt.close(fig)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('mode',choices=['preview','orbit','results'])
    args=parser.parse_args()
    for d in [O,FIG,VID,REV]:d.mkdir(parents=True,exist_ok=True)
    ref,field,paths=inputs();freeze(paths)
    if args.mode=='preview':
        orbit(ref,field,paths,True);stills(ref,field,paths)
    elif args.mode=='orbit':orbit(ref,field,paths,False)
    else:stills(ref,field,paths);plots();dynamics(ref,field,paths)
    assert all(sha(R/p)==digest for p,digest in json.loads((O/'SOURCE_HASHES.json').read_text()).items())
    dump(O/f'{args.mode.upper()}_COMPLETE.json',dict(status='PASS',mode=args.mode,
         source_hashes_unchanged=True,renderer_sha256=sha(Path(__file__))))
    print('COMPLETE',args.mode,flush=True)


if __name__=='__main__':
    main()
