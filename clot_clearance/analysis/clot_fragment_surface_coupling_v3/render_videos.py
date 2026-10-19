"""True FP64 trajectories and recorded surface tractions, 1920x1080 / 24fps."""
from common import *
from current_surface import CurrentSurface
import pyvista as pv
import matplotlib.colors as mc
from PIL import Image,ImageDraw,ImageFont
from vtkmodules.vtkCommonCore import vtkStringArray
import subprocess,csv,time,resource,argparse

graph=module('movie_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
WIDTH=1920;HEIGHT=1080;FPS=24;HOLD_FRAMES=48
FORCE_SCALE_UM_PER_N=1e8 # 2 uN -> 200 um, same physical scale in every view.
NAMES={1:'VIDEO_01_SURFACE_COUPLING_BEFORE_AFTER',2:'VIDEO_02_FRAGMENT_FLUID_TRACTION_AND_SURVIVAL'}
FIELDS=['x','velocity','d','atom_ids','time','step','surface_active','surface_centroid','surface_force','surface_area','surface_weight','surface_normal','surface_velocity','surface_outer_velocity','surface_unresolved']

def catalog(case):
    paths=sorted((case/'keys').glob('*.npz'));by_step={}
    for p in paths:
        if p.stem in ['failed_unaccepted','last_accepted'] and (case/'keys/final.npz').exists():continue
        if p.stem=='failed_unaccepted':continue
        with np.load(p) as f:step=int(f['step']);t=float(f['time'])
        if step not in by_step or p.stem in ['initial','final','first_softening','first_failure','first_multinode','upper_drive_off']:
            by_step[step]=dict(path=p,time=t,step=step)
    return sorted(by_step.values(),key=lambda q:q['time'])

def state(item):
    with np.load(item['path']) as f:return {k:f[k].copy() for k in FIELDS}|dict(source_key=item['path'].stem,render_event_caption=item.get('render_event_caption'))

def stage_caption(s):
    if s.get('render_event_caption'):return s['render_event_caption']
    explicit={'initial':'INITIAL MATERIAL','first_softening':'FIRST BULK SOFTENING','first_failure':'FIRST NATURAL FAILURE',
        'first_multinode':'FIRST MULTINODE DETACHMENT','first_secondary':'FIRST RECORDED SECONDARY SPLIT',
        'upper_drive_off':'DRIVE OFF / PASSIVE DRAG','final':'LAST ACCEPTED STATE','last_accepted':'LAST VALID ACCEPTED STATE'}
    if s['source_key'] in explicit:return explicit[s['source_key']]
    t=float(s['time'])
    return 'PASSIVE DRAG WINDOW' if t>=OFF else 'DRIVE DECREASING' if t>=HOLD else 'PRESCRIBED DRIVE PLATEAU' if t>=RISE else 'PRESCRIBED DRIVE RAMP'

def registry(case):
    reg=json.loads((case/'COMPONENT_REGISTRY.json').read_text());summary=json.loads((case/'SUMMARY.json').read_text());seeds={int(k):int(k) for k in reg}
    for e in summary['secondary_events']:
        parent=int(e['parent_component_id']);children=list(map(int,e['children_component_ids']));child=max(children,key=lambda q:reg[str(q)]['nodes'])
        if reg[str(child)]['nodes']>.5*reg[str(parent)]['nodes']:seeds[child]=seeds[parent]
    lookup={tuple(q['atom_ids']):(int(cid),seeds[int(cid)]) for cid,q in reg.items()}
    return lookup,reg,summary,seeds

def topology(ref,s,lookup):
    lab,sizes,attached=graph.topology(ref,s['d']);free=~attached[lab];rgb=np.full((len(lab),3),160,np.uint8);cids={}
    for label in np.flatnonzero(~attached):
        ids=np.flatnonzero(lab==label);cid,seed=lookup[tuple((ids+1).tolist())];cids[int(label)]=cid
        rgb[ids]=(255*mc.hsv_to_rgb([(seed*.61803398875)%1,.66,.77])).astype(np.uint8)
    return dict(labels=lab,sizes=sizes,attached=attached,free=free,rgb=rgb,cids=cids)

def line_mesh(points,bonds):
    ids,inverse=np.unique(bonds,return_inverse=True);local=inverse.reshape(-1,2)
    mesh=pv.PolyData(points[ids]);mesh.verts=np.empty(0,np.int64);mesh.lines=np.c_[np.full(len(bonds),2),local].ravel();assert mesh.n_verts==0
    return mesh

def polygon_mesh(points):
    points=np.asarray(points);count,nv,_=points.shape
    mesh=pv.PolyData(points.reshape(-1,3));mesh.verts=np.empty(0,np.int64);mesh.faces=np.c_[np.full(count,nv),np.arange(count*nv).reshape(count,nv)].ravel();return mesh

def camera_for(bounds):
    center=np.array([(bounds[2*k]+bounds[2*k+1])/2 for k in range(3)]);span=max(bounds[1]-bounds[0],bounds[3]-bounds[2],bounds[5]-bounds[4])
    direction=np.array([1.4,-2.9,1.0]);forward=-direction/np.linalg.norm(direction)
    right=np.array([1.,0.,0.]);right-=np.dot(right,forward)*forward;right/=np.linalg.norm(right)
    up=np.cross(right,forward) # Image-right is projected +X; X stays horizontal.
    return [center+span*direction,center,up]

def set_camera(pl,bounds):
    camera=camera_for(bounds);pl.camera_position=camera;pl.enable_parallel_projection()
    corners=np.array([(x,y,z) for x in bounds[:2] for y in bounds[2:4] for z in bounds[4:]])
    forward=np.asarray(camera[1])-camera[0];forward/=np.linalg.norm(forward);right=np.cross(forward,camera[2]);right/=np.linalg.norm(right);up=np.cross(right,forward)
    aspect=(WIDTH/2)/(.64*HEIGHT);scale=1.32*max(np.ptp(corners@up)/2,np.ptp(corners@right)/(2*aspect))
    pl.camera.parallel_scale=float(scale);span=max(np.ptp(corners,axis=0));pl.camera.clipping_range=(span*.01,span*20)

def box_grid_points(bounds):
    lo=np.asarray(bounds)[::2];hi=np.asarray(bounds)[1::2];lines=[]
    for normal in range(3):
        others=[a for a in range(3) if a!=normal]
        for side in [lo[normal],hi[normal]]:
            for along,across in [others,others[::-1]]:
                for value in np.linspace(lo[across],hi[across],4):
                    first=lo.copy();last=hi.copy();first[normal]=last[normal]=side;first[across]=last[across]=value
                    lines.extend([first,last])
    return np.asarray(lines)

def base_view(pl,j,bounds):
    pl.subplot(0,j);pl.renderer.SetViewport(j*.5,.18,(j+1)*.5,.82);pl.set_background('#fcfcfc');set_camera(pl,bounds)
    points=box_grid_points(bounds);grid=line_mesh(points,np.arange(len(points)).reshape(-1,2))
    pl.add_mesh(grid,name=f'fixed_grid_{j}',color='#647987',opacity=.75,line_width=.85,lighting=False)
    pl.add_axes(interactive=False);return dict(grid=grid,bounds=np.asarray(bounds))

def update_axis(widget,bounds):
    widget['grid'].points=box_grid_points(bounds);widget['bounds']=np.asarray(bounds)

def axis_annotations(pl,widgets):
    # Project true axis positions through each current VTK camera; draw screen
    # text with PIL to avoid CubeAxesActor's multi-viewport text distortion.
    annotations=[]
    for j,widget in enumerate(widgets):
        pl.subplot(0,j);bounds=widget['bounds'];lo=bounds[::2];hi=bounds[1::2]
        def project(point):
            pl.renderer.SetWorldPoint(*point,1);pl.renderer.WorldToDisplay();x,y,_=pl.renderer.GetDisplayPoint()
            return np.array([x,HEIGHT-y])
        for axis in range(3):
            start=lo.copy();end=lo.copy()
            if axis==1:start[0]=end[0]=hi[0]
            end[axis]=hi[axis];delta=[np.array([0.,15.]),np.array([16.,6.]),np.array([-16.,0.])][axis]
            anchor=['mt','lm','rm'][axis]
            for k,value in enumerate(np.linspace(lo[axis],hi[axis],4)):
                if axis in (1,2) and k==0:continue
                point=start.copy();point[axis]=value
                annotations.append(dict(position=project(point)+delta,text=f'{value:.0f}',anchor=anchor,size=21))
            extra=[np.array([0.,42.]),np.array([37.,19.]),np.array([-56.,0.])][axis]
            annotations.append(dict(position=project((start+end)/2)+extra,text='XYZ'[axis]+' (um)',anchor=anchor,size=23))
    return annotations

def draw_scene(pl,ref,s,top,surf,z,new,prefix,selected=None,detail=False):
    n=int(ref['clot_count']);xyz=s['x']*1e6;mask=np.ones(n,bool) if selected is None else selected;names=[]
    def points(key,pts,**kw):
        if len(pts):pl.add_points(pts,name=prefix+key,render_points_as_spheres=True,show_scalar_bar=False,**kw);names.append(prefix+key)
    if not detail:points('anchors',xyz[n:],color='#87929c',point_size=3,opacity=.25)
    attached=~top['free']&mask
    if attached.any():points('attached',xyz[:n][attached],color='#87a4b7',point_size=3.3,opacity=.55)
    nn=top['sizes'][top['labels']]
    for key,keep,size in [('fine',nn<=2,3),('small',(nn>2)&(nn<20),4),('target',(nn>=20)&(nn<=200),7),('large',nn>200,5)]:
        ids=mask&top['free']&keep
        args=dict(color='#b5a186') if key=='fine' else dict(scalars=top['rgb'][ids],rgb=True)
        points(key,xyz[:n][ids],point_size=size+(2 if detail else 0),**args)
    b=ref['bonds'];member=np.r_[mask,np.zeros(len(ref['X'])-n,bool)]
    if detail:keep=(s['d']<1)&member[b].all(1)&(ref['bond_types']==1)
    else:keep=(s['d']<1)&(ref['bond_types']==1)&np.r_[top['free'],np.zeros(len(ref['X'])-n,bool)][b].all(1)
    if keep.any():pl.add_mesh(line_mesh(xyz,b[keep]),name=prefix+'live',color='#365a71',opacity=.38 if detail else .15,line_width=1.1 if detail else .6,lighting=False);names.append(prefix+'live')
    soft_all=(s['d']>0)&(s['d']<1)&(ref['bond_types']==1)&member[b].all(1);soft=soft_all&(ref['r0']<=1.01*DX)
    if not soft.any():soft=soft_all
    if soft.any():pl.add_mesh(line_mesh(xyz,b[soft]),name=prefix+'soft',color='#7e659e',opacity=.65,line_width=1.,lighting=False);names.append(prefix+'soft')
    failed_all=(s['d']==1)&~ref['initial_cut']&member[b].all(1);failed=failed_all&(ref['r0']<=1.01*DX)
    if not failed.any():failed=failed_all
    if failed.any():pl.add_mesh(line_mesh(xyz,b[failed]),name=prefix+'failed',color='#bc535d',opacity=.08 if detail else .035,line_width=.8,lighting=False);names.append(prefix+'failed')
    L=surf.L;tri=surf.base.tri;active=s['surface_active'];l1=active[:L]&member[tri].all(1)
    if l1.any():pl.add_mesh(polygon_mesh(xyz[tri[l1]]),name=prefix+'L1',color='#50a182',opacity=.4,show_edges=False);names.append(prefix+'L1')
    if new:
        patches=np.flatnonzero(active[L:]&member[surf.owner]);quads=[]
        for p in patches:
            i=surf.owner[p];a=surf.axis[p];j,k=(a+1)%3,(a+2)%3
            corners=np.array([surf.offset[p]+DX*(u*np.eye(3)[j]+v*np.eye(3)[k]) for u,v in [(-.5,-.5),(.5,-.5),(.5,.5),(-.5,.5)]])
            quads.append((s['x'][i]+corners@z['l2_deformation'][i].T)*1e6)
        if quads:pl.add_mesh(polygon_mesh(quads),name=prefix+'L2',color='#38a9c3',opacity=.42,show_edges=True,edge_color='#2f7188',line_width=.4);names.append(prefix+'L2')
    else:patches=[]
    owners=np.r_[tri[:,0],surf.owner] if new else tri[:,0];face_mask=active&member[owners]
    norms=np.linalg.norm(s['surface_force'],axis=1);eligible=np.flatnonzero(face_mask&(norms>1e-16))
    if len(eligible):
        chosen=eligible[np.argsort(norms[eligible])[-min(18,len(eligible)):]]
        pl.add_arrows(s['surface_centroid'][chosen]*1e6,s['surface_force'][chosen]*FORCE_SCALE_UM_PER_N,mag=1,color='#b95336',name=prefix+'traction',lighting=False)
        names.append(prefix+'traction')
    ids=np.flatnonzero(mask);stats=dict(nodes=len(ids),original_atom_ids=ids+1,internal_live_bonds=int(keep.sum()),L1_facets=int(l1.sum()),L2_patches=len(patches),active_elements=int(face_mask.sum()),
        force_N=s['surface_force'][face_mask].sum(0),surface_area_m2=float(np.sum(s['surface_area'][face_mask]*s['surface_weight'][face_mask])),
        power_W=float(np.sum(s['surface_force'][face_mask]*s['surface_velocity'][face_mask])),COM_speed_m_s=float(np.linalg.norm(s['velocity'][ids].mean(0))))
    v=s['surface_velocity'][face_mask];normal=s['surface_normal'][face_mask]
    slip=s['surface_outer_velocity'][face_mask]-(v-normal*np.sum(v*normal,axis=1)[:,None]);area=s['surface_area'][face_mask]*s['surface_weight'][face_mask]
    stats['mean_slip_m_s']=float(np.sum(np.linalg.norm(slip,axis=1)*area)/area.sum()) if area.sum()>0 else None
    if detail and len(ids) and np.all(top['free'][ids]) and np.linalg.norm(stats['force_N'])>1e-16:
        direction=stats['force_N']/np.linalg.norm(stats['force_N'])*200.
        pl.add_arrows(xyz[ids].mean(0)[None,:],direction[None,:],mag=1,color='#bf821d',name=prefix+'resultant',lighting=False)
        names.append(prefix+'resultant')
    return names,stats

def overlay(raw,which,states,tops,details,status,focus_id,hold=False,annotations=()):
    image=Image.fromarray(raw).convert('RGB');d=ImageDraw.Draw(image);font='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
    big=ImageFont.truetype(font,32);normal=ImageFont.truetype(font,25);small=ImageFont.truetype(font,21)
    d.rectangle((0,0,1920,198),fill='#fcfcfc');d.rectangle((0,886,1920,1080),fill='#fcfcfc')
    title='Old H versus new H · fragment surface coupling' if which==1 else 'New H · real fragment traction, motion and survival'
    d.text((32,12),title,font=big,fill='#243c4b')
    d.text((32,58),'Accepted coordinates ×1 | mechanical shear surrogate | same material and prescribed drive | no trajectory interpolation',font=small,fill='#4c6372')
    for j,key in enumerate(['old','new'] if which==1 else ['new','new']):
        s=states[key];top=tops[key];sizes=top['sizes'][~top['attached']];target=int(np.sum((sizes>=20)&(sizes<=200)));det=int(sizes.sum());fine=int(sizes[sizes<=2].sum())
        label=('V2 H · supported triangles' if key=='old' else 'V3 H · triangles + partial dual faces') if which==1 else ('V3 H · whole scene' if j==0 else f'V3 H · detail (camera zoom) · {focus_id}')
        x=j*960+30;d.text((x,106),label,font=normal,fill='#253f50');d.text((x,146),f'Actual t={float(s["time"])*1e3:.5f} ms'+(' · INCOMPLETE' if status[key]!='COMPLETED' else ''),font=normal,fill='#253f50')
        d.text((x,177),stage_caption(s),font=ImageFont.truetype(font,18),fill='#647784')
        if which==2 and j==1:
            q=details[j];slip=f'{q["mean_slip_m_s"]:.3g}' if q['mean_slip_m_s'] is not None else 'unresolved'
            lines=[f'Members {q["nodes"]} | live bonds {q["internal_live_bonds"]} | L1/L2 {q["L1_facets"]}/{q["L2_patches"]}',f'|F|={np.linalg.norm(q["force_N"])*1e6:.3g} uN | slip {slip} m/s | P={q["power_W"]*1e6:.3g} uW',f'Area {q["surface_area_m2"]*1e6:.4g} mm² | COM speed {q["COM_speed_m_s"]:.3g} m/s']
            if q['surface_area_m2']==0:
                lines[1]='FLUID_SURFACE_UNRESOLVED | zero numerical load, physical drag unknown'
        else:lines=[f'Target 20–200 nodes: {target} | largest free: {int(sizes.max(initial=0))}',f'Detached: {100*det/8125:.2f}% | fines / detached: '+(f'{100*fine/det:.2f}%' if det else 'undefined'),'All original material nodes remain represented.']
        for k,line in enumerate(lines):d.text((x,902+30*k),line,font=small,fill='#3b5668')
        if which==2 and j==1:
            shown=', '.join(str(i) for i in details[j]['original_atom_ids'][:8])
            if details[j]['nodes']>8:shown+=', ...'
            d.text((x,987),'Original atom IDs: '+shown+' | full set: frame CSV',font=ImageFont.truetype(font,17),fill='#4c6372')
    d.text((32,1005),'Green L1 / cyan L2 | orange face forces: 2 uN = 200 um | gold COM arrow: resultant direction only',font=small,fill='#456273')
    d.text((32,1037),'FLUID_SURFACE_UNRESOLVED remains for unsupported / uncertain geometry; a missing face is not zero physical drag.',font=small,fill='#9b4854')
    d.text((32,1059),'Purple: live softened bonds. Faint red: failed traces, no material force. Colors do not match fragments across versions.',font=ImageFont.truetype(font,18),fill='#647784')
    for label in annotations:
        d.text(tuple(label['position']),label['text'],font=ImageFont.truetype(font,label['size']),fill='#405563',anchor=label['anchor'])
    return np.asarray(image)

def render(which,timeline,ref,surf,lookups,bounds,focus_ids,focus_width,status,smoke=False):
    out=F/'videos'/('render_smoke' if smoke else '');out.mkdir(parents=True,exist_ok=True);name=NAMES[which];start=time.perf_counter();cpu=time.process_time()
    command=['ffmpeg','-hide_banner','-loglevel','error','-y','-f','rawvideo','-pix_fmt','rgb24','-s','1920x1080','-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','fast','-crf','17','-pix_fmt','yuv420p','-movflags','+faststart',str(out/(name+'.mp4'))]
    log=(out/(name+'.ffmpeg.log')).open('w');enc=subprocess.Popen(command,stdin=subprocess.PIPE,stderr=log)
    pl=pv.Plotter(shape=(1,2),off_screen=True,window_size=(WIDTH,HEIGHT),border=False);pl.ren_win.SetMultiSamples(4);axes=[base_view(pl,j,bounds) for j in range(2)];actors=[[],[]];records=[];cache={}
    previous_focus_width=None;previous_focus_time=None
    try:
        for index,row in enumerate(timeline):
            states={};tops={}
            for key in ['old','new']:
                item=row[key];cached=cache.get(key)
                if cached is None or cached[0]!=item['path']:cache[key]=(item['path'],state(item))
                states[key]=cache[key][1];tops[key]=topology(ref,states[key],lookups[key])
            s=states['new'];z=surf.geometry(s['x'],s['d']);assert np.array_equal(z['active'],s['surface_active'])
            lab=tops['new']['labels'];label=int(np.argmax(np.bincount(lab[focus_ids])));selected=lab==label
            focus_id='attached body' if tops['new']['attached'][label] else 'ID '+str(tops['new']['cids'][label]);details=[]
            detail_points=s['x'][:len(lab)][selected]*1e6;center=detail_points.mean(0)
            required_width=max(focus_width,float(abs(detail_points-center).max())*2.5)
            current_time=float(s['time'])
            if previous_focus_width is None:detail_width=required_width
            else:
                decay=np.exp(-max(0.,current_time-previous_focus_time)/60e-6)
                detail_width=max(required_width,decay*previous_focus_width+(1-decay)*required_width)
            previous_focus_width=detail_width;previous_focus_time=current_time
            for j,key in enumerate(['old','new'] if which==1 else ['new','new']):
                pl.subplot(0,j)
                for actor in actors[j]:pl.renderer.remove_actor(actor,render=False)
                detail=which==2 and j==1
                actors[j],info=draw_scene(pl,ref,states[key],tops[key],surf,z,key=='new',f'p{j}_',selected if detail else None,detail);details.append(info)
                panel_bounds=bounds
                if detail:
                    panel_bounds=np.ravel(np.column_stack([center-detail_width/2,center+detail_width/2]))
                set_camera(pl,panel_bounds);update_axis(axes[j],panel_bounds)
            # Render after both camera and grid changes; screenshot may otherwise
            # reuse an earlier buffer while labels project through the new camera.
            pl.render()
            raw=pl.screenshot(return_img=True);im=overlay(raw,which,states,tops,details,status,focus_id,annotations=axis_annotations(pl,axes));enc.stdin.write(im.tobytes())
            if smoke or index in [0,len(timeline)//2,len(timeline)-1]:Image.fromarray(im).save(out/(name+f'_frame_{index:03d}.png'))
            if index==len(timeline)-1:Image.fromarray(im).save(out/(name+'_COVER.png'))
            left='old' if which==1 else 'new'
            records.append(dict(frame=index,playback_s=index/FPS,requested_time_s=row['requested_time_s'],left_actual_time_s=float(states[left]['time']),right_actual_time_s=float(s['time']),
                left_source=str(row[left]['path'].relative_to(R)),right_source=str(row['new']['path'].relative_to(R)),focus_label=focus_id,focus_nodes=int(selected.sum()),
                focus_original_atom_ids=json.dumps((np.flatnonzero(selected)+1).tolist(),separators=(',',':')),
                focus_camera_width_um=detail_width if which==2 else '',
                focus_member_ID_sha256=hashlib.sha256(np.asarray(np.flatnonzero(selected)+1,np.int64).tobytes()).hexdigest(),state_interpolation=False,hold=False))
            if index%24==0:print('VIDEO_FRAME',which,index,flush=True)
        for k in range(HOLD_FRAMES if not smoke else 3):
            enc.stdin.write(im.tobytes());q=records[-1].copy();q.update(frame=len(records),playback_s=len(records)/FPS,hold=True);records.append(q)
        enc.stdin.close();assert enc.wait()==0
    finally:
        pl.close();log.close()
        if enc.poll() is None:enc.kill();enc.wait()
    with (out/(name+'_FRAME_TIMES.csv')).open('w') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
    (out/(name+'_FFMPEG_COMMAND.txt')).write_text(' '.join(command)+'\n')
    probe=json.loads(subprocess.check_output(['ffprobe','-v','error','-show_streams','-show_format','-of','json',str(out/(name+'.mp4'))],text=True));stream=probe['streams'][0]
    assert stream['width']==1920 and stream['height']==1080 and stream['avg_frame_rate']=='24/1' and int(stream['nb_frames'])==len(records)
    return dict(file=name+'.mp4',sha256=sha(out/(name+'.mp4')),frames=len(records),resolution=[WIDTH,HEIGHT],fps=FPS,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,ffprobe=probe)

def main():
    p=argparse.ArgumentParser();p.add_argument('--smoke',action='store_true');args=p.parse_args()
    ref=load_npz(O/'provenance/reference.npz');surf=CurrentSurface(ref,'uniform');cases=dict(old=V2_O/'cases/biphasic',new=O/'cases/biphasic')
    catalogs={k:catalog(v) for k,v in cases.items()};lookups={};registers={};status={};colors={}
    summaries={}
    for key,case in cases.items():
        lookups[key],registers[key],summary,colors[key]=registry(case);status[key]=summary['status'];summaries[key]=summary
    result=json.loads((O/'NATIVE_V3_RESULT.json').read_text());targets=result['final_target_components']
    if targets:focus=max(targets,key=lambda q:q['nodes']);focus_ids=np.array(focus['original_atom_ids'])-1
    else:
        end=state(catalogs['new'][-1]);top=topology(ref,end,lookups['new']);free=np.flatnonzero(~top['attached']);label=int(free[np.argmax(top['sizes'][free])]) if len(free) else int(np.argmax(top['sizes']));focus_ids=np.flatnonzero(top['labels']==label)
    maximum=min(catalogs['old'][-1]['time'],catalogs['new'][-1]['time'])
    event_frames=[]
    def mark_after_event(event_time,label,selection=None):
        available=[q for q in catalogs['new'] if q['time']>=event_time-1e-14 and q['time']<=maximum]
        if not available:return
        item=available[0];item['render_event_caption']=label
        event_frames.append(dict(label=label,selection=selection,event_time_s=event_time,rendered_accepted_time_s=item['time'],
            source=str(item['path'].relative_to(R)),time_offset_s=item['time']-event_time,
            scope='First archived accepted coordinates at or after the event; not a reconstructed state at the event.'))
    with (cases['new']/'history.csv').open() as f:accepted=list(csv.DictReader(f))
    first_large=next((float(r['time_s']) for r in accepted if int(r['largest_nodes'])>=20),None)
    if first_large is not None:mark_after_event(first_large,'FIRST SAVED FRAME AFTER >=20-NODE DETACHMENT')
    major_events=[e for e in summaries['new']['secondary_events'] if sum(registers['new'][str(c)]['nodes']>=20 for c in e['children_component_ids'])>=2]
    if major_events:
        mark_after_event(major_events[0]['time_s'],'FIRST SAVED FRAME AFTER MULTI-CHILD SPLIT (>=20 EACH)',
            'First accepted split with at least two children of at least 20 nodes; children may exceed the 20–200 target band.')
        largest_secondary=max(major_events,key=lambda e:sorted((registers['new'][str(c)]['nodes'] for c in e['children_component_ids']),reverse=True)[1])
        mark_after_event(largest_secondary['time_s'],'AFTER SPLIT WITH THE LARGEST SECONDARY CHILD',
            'Maximize the second-largest child node count among recorded secondary splits; display selection only, no change to the physical run or fragment statistics.')
    first_l2=next((float(r['time_s']) for r in accepted if int(r['L2_active_patches'])>0),None)
    if first_l2 is not None:mark_after_event(first_l2,'FIRST SAVED FRAME AFTER L2 ACTIVATION')
    events=[q['time'] for arr in catalogs.values() for q in arr if (not q['path'].stem.startswith('frame_') or q.get('render_event_caption')) and q['time']<=maximum]
    requested=np.unique(np.r_[np.linspace(0,maximum,145),events]);timeline=[dict(requested_time_s=t,**{k:min(arr,key=lambda q:abs(q['time']-t)) for k,arr in catalogs.items()}) for t in requested]
    lo=np.full(3,np.inf);hi=-lo;focus_width=750.
    for key,arr in catalogs.items():
        for item in arr:
            s=state(item);xyz=s['x']*1e6;lo=np.minimum(lo,xyz.min(0));hi=np.maximum(hi,xyz.max(0))
    pad=np.maximum((hi-lo)*.06,25.);bounds=np.ravel(np.column_stack([lo-pad,hi+pad]))
    if args.smoke:timeline=[timeline[i] for i in np.linspace(0,len(timeline)-1,5).astype(int)]
    videos=[render(which,timeline,ref,surf,lookups,bounds,focus_ids,focus_width,status,args.smoke) for which in [1,2]]
    dump(O/('VIDEO_SMOKE.json' if args.smoke else 'VIDEO_GENERATION.json'),dict(status='PASS',videos=videos,actual_displacement_scale=1,accepted_states_only=True,interpolation=False,
        fixed_global_bounds_um=bounds,global_camera=camera_for(bounds),global_X_axis_horizontal=True,minimum_focus_window_width_um=focus_width,focus_original_atom_ids=focus_ids+1,
        video2_detail='Camera zoom only, explicitly labeled; fixed orientation and COM tracking of the actual connected component containing most tracked original IDs. Width fits its actual extent with a 750 um minimum; zoom-in relaxes over 60 us physical time and zoom-out always fits all selected nodes. Per-frame physical width is recorded in CSV. No material-coordinate offsets or exploded views.',
        force_arrow_scale_um_per_N=FORCE_SCALE_UM_PER_N,resultant_direction_arrow_length_um=200,resultant_direction_arrow_is_not_magnitude_scaled=True,colors=colors,color_policy='same deterministic lineage rule; >50% child inherits, no cross-version identity claim',
        face_force_arrow_display='At most the 18 strongest active face forces per panel above 1e-16 N. This selects arrows for legibility only; all active elements and full recorded forces remain in the scientific data.',
        damage_trace_display='Nearest-neighbor original bonds for readable damage traces; full-length fallback if that class has no nearest-neighbor bond. Softened traces remain live; failed traces carry no material force. All original bonds still define connectivity.',
        event_frames=event_frames,per_frame_original_members='The complete actual selected-component atom IDs are stored in each frame CSV, alongside their identity hash. An onscreen excerpt is only a readable index.',
        limitations='Partial approximate surface quadrature; unsupported cells explicitly FLUID_SURFACE_UNRESOLVED. Nearest real accepted timestamps are shown per panel.'))

if __name__=='__main__':main()
