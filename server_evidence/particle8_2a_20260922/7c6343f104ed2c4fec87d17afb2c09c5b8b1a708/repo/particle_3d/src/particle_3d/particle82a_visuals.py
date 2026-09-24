"""Three-dimensional diagnostic replay from saved inputs and computed paths."""
from pathlib import Path
import argparse,gzip,json,math,os
os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
import numpy as np
import pyvista as pv
import imageio_ffmpeg
from PIL import Image,ImageDraw,ImageFont
from .particle82a_admission import context
from .particle82a_pipeline import BASINS,atomic_npz
from .particle82a_figures import BG,INK,MUTED,COLORS,METHOD_COLORS
from .particle82_provenance import atomic_json,sha256
from .particle8_replay import canonical_hash

ANIMATIONS=['particle8_2a_anim_01_inlet_geometry_rotating',
 'particle8_2a_anim_02_three_method_birth_comparison','particle8_2a_anim_03_full_vessel_outcome_comparison']
SIZE=(1920,1080);FPS=24;FRAMES=432


def text(draw,xy,value,size=24,color=INK,bold=False):
    path='/usr/share/fonts/truetype/dejavu/DejaVuSans'+('-Bold' if bold else '')+'.ttf'
    draw.text(xy,value,font=ImageFont.truetype(path,size),fill=color)


def surface(points,faces):
    return pv.PolyData(np.asarray(points),np.column_stack([np.full(len(faces),3),faces]).ravel())


def camera(plotter,center,scale,azimuth,elevation=25,up=(0,0,1)):
    a=math.radians(azimuth);e=math.radians(elevation)
    direction=np.array([math.cos(a)*math.cos(e),math.sin(a)*math.cos(e),math.sin(e)])
    plotter.camera_position=[np.asarray(center)+4*scale*direction,center,up]
    plotter.camera.parallel_projection=True;plotter.camera.parallel_scale=scale
    plotter.reset_camera_clipping_range()
    return dict(azimuth_deg=azimuth,elevation_deg=elevation,focal_point=list(center),parallel_scale=scale,
        physical_coordinates_rotated=False)


def line_mesh(paths,origin):
    points=[];lines=[];offset=0
    for path in paths:
        if len(path)<2:continue
        p=(np.asarray(path)-origin)*1e6;points.append(p);lines.extend([len(p),*range(offset,offset+len(p))]);offset+=len(p)
    return pv.PolyData(np.vstack(points),lines=np.asarray(lines)) if points else None


def prepare_scene(admission,trajectories,output,partial=False):
    admission=Path(admission);trajectories=Path(trajectories);output=Path(output)
    (output/'scene').mkdir(parents=True,exist_ok=True);c=context();g=c.geometry
    arrays={}
    for name,s in c.env.boundaries.items():
        arrays[name+'_points_m']=np.asarray(s.points);arrays[name+'_faces']=s.faces.reshape(-1,4)[:,1:]
    arrays['inlet_origin_m']=g.origin;arrays['inlet_basis']=g.basis;arrays['inward_normal']=g.normal
    arrays['map_rows']=np.load(admission/'map_n24.npz')['rows']
    atomic_npz(output/'scene'/'geometry.npz',**arrays)
    selected={};births=[]
    # First event in each resolved anchor basin; selection never uses outcomes.
    for file in sorted((admission/'events').glob('events_*.json.gz')):
        with gzip.open(file,'rt') as f:rows=json.load(f)['rows']
        needed=[(i,row) for i,row in enumerate(rows) if row['event']['point_tracer_basin'] in BASINS[:3] and row['event']['point_tracer_basin'] not in selected]
        if needed:
            raw=np.load(file.with_name(file.name.replace('.json.gz','.npz')))
            for i,row in needed:
                b=row['event']['point_tracer_basin']
                if b in selected:continue
                selected[b]=row['event']['event_id'];pid=row['event']['event_id']
                attempts=raw['attempts'];attempts=attempts[attempts[:,1]==pid]
                path=raw['paths'][raw['offsets'][i]:raw['offsets'][i+1]]
                atomic_npz(output/'scene'/f'birth_{pid}.npz',attempts=attempts,path=path)
                births.append(row)
        if len(selected)==3:break
    birth_manifest=dict(selection='FIRST_EVENT_ID_IN_EACH_ANCHOR_POINT_BASIN_NO_OUTCOME_SELECTION',
        ids=selected,events=sorted(births,key=lambda r:r['event']['point_tracer_basin']),
        source_admission_receipts={p.name:sha256(p) for p in sorted((admission/'events').glob('*.receipt.json'))})
    atomic_json(output/'scene'/'births.json',birth_manifest)
    outcomes={};source_hashes={}
    for method in 'ABC':
        folder=trajectories/'formal'/method/'trajectories'
        metas=[p for p in sorted(folder.glob('mb_*.json')) if not p.name.endswith('.receipt.json')]
        if not metas and partial:continue
        if not metas:raise ValueError('Missing formal trajectory cohort '+method)
        records=[json.loads(p.read_text()) for p in metas]
        # Representative coverage: natural first IDs plus every minor-outlet
        # completion, explicitly labelled as a display selection only.
        keep=set(r['particle_id'] for r in records[:90])
        keep.update(r['particle_id'] for r in records if r['exit_outlet'] in ['OUTLET_01','OUTLET_03'])
        chosen=[r for r in records if r['particle_id'] in keep];positions=[];offsets=[0]
        for r in chosen:
            p=folder/f'mb_{r["particle_id"]:06d}.npz';a=np.load(p)['samples'];positions.append(a[:,1:4]);offsets.append(offsets[-1]+len(a))
            source_hashes[str(p)]=sha256(p)
        atomic_npz(output/'scene'/('outcomes_'+method+'.npz'),positions=np.concatenate(positions),offsets=np.asarray(offsets))
        outcomes[method]=dict(display_records=[dict(event_id=r['particle_id'],outlet=r['exit_outlet'],end_reason=r['end_reason'],completed=r['completed']) for r in chosen],
            total_integrated=len(records),outlet_counts={b:sum(r['exit_outlet']==b for r in records) for b in BASINS[:3]},
            safety_stops=sum(r['end_reason']=='INTEGRATION_SAFETY_STOP' for r in records),
            display_selection='FIRST_90_IDS_PLUS_ALL_OUTLET_01_OR_03_COMPLETIONS; NOT_A_POPULATION_ESTIMATOR')
    atomic_json(output/'scene'/'outcomes.json',dict(methods=outcomes,source_sample_hashes=source_hashes))
    atomic_json(output/'scene'/'SCENE_MANIFEST.json',dict(files={p.name:sha256(p) for p in sorted((output/'scene').glob('*')) if p.is_file() and p.name!='SCENE_MANIFEST.json'},
        positions='SAVED_SI_COORDINATES_NO_PHYSICAL_ROTATION_NO_AXIS_STRETCH',
        radii='ORIGINAL_OR_ACTUALLY_ACCEPTED_B_DIAMETERS',
        animation_02='DIAGNOSTIC_VISUALIZATION_NOT_PHYSICAL_ENTRY_DYNAMICS'))


class Views:
    def __init__(self,output,kind):
        self.output=Path(output);self.kind=kind;self.geo=np.load(self.output/'scene'/'geometry.npz')
        self.births=json.loads((self.output/'scene'/'births.json').read_text())
        self.outcomes=json.loads((self.output/'scene'/'outcomes.json').read_text())['methods']
        self.views=[];self.dynamic=[]
        self.origin=self.geo['inlet_origin_m'] if kind<2 else (self.geo['WALL_points_m'].max(0)+self.geo['WALL_points_m'].min(0))/2
        self.columns=2 if kind==0 else 3;self.panel=(920,740) if kind==0 else (620,740)
        self.wall=surface((self.geo['WALL_points_m']-self.origin)*1e6,self.geo['WALL_faces'])
        if kind<2:
            normal=self.geo['inward_normal'];center=normal*2.5
            bounds=np.column_stack([center-5.5,center+5.5]).ravel()
            self.wall=self.wall.clip_box(bounds,invert=False).extract_surface(algorithm='dataset_surface')
            self.center=center;self.scale=5.7
        else:
            self.center=np.zeros(3);span=np.ptp((self.geo['WALL_points_m']-self.origin)*1e6,axis=0)
            self.scale=np.linalg.norm(span)*.54
        for i in range(self.columns):
            p=pv.Plotter(off_screen=True,window_size=self.panel);p.set_background(BG);p.render_window.SetMultiSamples(0)
            p.enable_depth_peeling(number_of_peels=12,occlusion_ratio=0)
            p.add_mesh(self.wall,color='#aab6c4',opacity=.12 if kind<2 else .08,smooth_shading=True,show_scalar_bar=False)
            cap=surface((self.geo['INLET_points_m']-self.origin)*1e6,self.geo['INLET_faces'])
            p.add_mesh(cap,color='#8694a8',opacity=.13,show_edges=False)
            if kind==0:
                maps=self.geo['map_rows'];points=pv.PolyData((maps[:,:3]-self.origin)*1e6)
                if i==0:
                    points['basin']=maps[:,5].astype(int)
                    p.add_mesh(points,scalars='basin',cmap=COLORS,clim=(-.5,3.5),point_size=3,show_scalar_bar=False)
                else:
                    points['Dmax']=maps[:,8]
                    p.add_mesh(points,scalars='Dmax',cmap='viridis',point_size=3,show_scalar_bar=True,
                        scalar_bar_args=dict(title='Dmax (um)',vertical=True,position_x=.86,position_y=.22,height=.52,width=.065,color=INK,title_font_size=15,label_font_size=13))
            if kind==2:
                method='ABC'[i];dataset=np.load(self.output/'scene'/('outcomes_'+method+'.npz'))
                records=self.outcomes[method]['display_records']
                for basin,color in zip([None,*BASINS[:3]],['#667487',*COLORS]):
                    paths=[dataset['positions'][dataset['offsets'][j]:dataset['offsets'][j+1]] for j,r in enumerate(records) if r['outlet']==basin]
                    mesh=line_mesh(paths,self.origin)
                    if mesh is not None:p.add_mesh(mesh,color=color,line_width=1.7 if basin else 1.,opacity=.9 if basin else .35)
                for role,color in zip(BASINS[:3],COLORS):
                    mesh=surface((self.geo[role+'_points_m']-self.origin)*1e6,self.geo[role+'_faces'])
                    p.add_mesh(mesh,color=color,opacity=.4)
            p.add_axes(color=MUTED,xlabel='x',ylabel='y',zlabel='z',line_width=2)
            p.show(auto_close=False,interactive=False);self.views.append(p);self.dynamic.append([])

    def birth_frame(self,index,total):
        events=self.births['events'];slot=min(len(events)-1,index*len(events)//total)
        e=events[slot];local=(index-slot*total/len(events))/(total/len(events));pid=e['event']['event_id']
        raw=np.load(self.output/'scene'/f'birth_{pid}.npz');specs=[]
        for mi,p in enumerate(self.views):
            for actor in self.dynamic[mi]:p.remove_actor(actor,reset_camera=False,render=False)
            self.dynamic[mi]=[];method='ABC'[mi];result=e['methods'][method]
            anchor=(np.asarray(e['event']['anchor_m'])-self.origin)*1e6
            p.add_mesh(pv.PolyData(anchor[None,:]),color=INK,point_size=9,render_points_as_spheres=True,name='anchor',reset_camera=False,render=False)
            self.dynamic[mi].append('anchor')
            if method in 'AB':
                attempts=raw['attempts'][raw['attempts'][:,0]==mi]
                j=min(len(attempts)-1,int(local*len(attempts)))
                a=attempts[j];position=a[4:7];radius=a[3]/2;state='accepted' if int(a[8])==0 else 'rejected'
                description=f'Trial {j+1}/{len(attempts)} | D={a[3]:.3f} µm | {state}'
                spec=dict(method=method,event_id=pid,draw=int(a[2]),position_m=position.tolist(),radius_um=radius,state=state)
            else:
                path=raw['path'];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path[:,1:],axis=0),axis=1))]
                end=result.get('s_birth_m',arc[-1]);s=min(end,max(0,local)*end)
                position=np.array([np.interp(s,arc,path[:,j]) for j in [1,2,3]])
                if local>.94 and result['accepted']:position=np.asarray(result['birth_center_m'])
                radius=e['event']['first_diameter_um']/2;state=result['status'] if local>.94 else 'entry representation search'
                upto=np.searchsorted(arc,s);pts=np.vstack([path[:max(1,upto),1:],position])
                mesh=line_mesh([pts],self.origin)
                if mesh is not None:p.add_mesh(mesh,color=METHOD_COLORS[2],line_width=3,name='path',reset_camera=False,render=False);self.dynamic[mi].append('path')
                description=f's={s*1e6:.3f} µm | D={2*radius:.3f} µm'
                spec=dict(method=method,event_id=pid,position_m=position.tolist(),radius_um=radius,s_m=float(s),state=state)
            sphere=pv.Sphere(radius=1.,theta_resolution=32,phi_resolution=24)
            sphere.points=np.asarray(sphere.points,float)*radius+(position-self.origin)*1e6
            p.add_mesh(sphere,color=METHOD_COLORS[mi],opacity=.85,smooth_shading=True,name='sphere',reset_camera=False,render=False)
            self.dynamic[mi].append('sphere');spec['description']=description;specs.append(spec)
        return e,specs

    def frame(self,index,total):
        image=Image.new('RGB',SIZE,BG);draw=ImageDraw.Draw(image);specs=[]
        titles=['Open inlet: outlet basins and passable diameter','Three birth strategies for the same events','Full-vessel outcomes under A / B / C']
        text(draw,(45,25),'Particle-8.2A  |  '+titles[self.kind],32,bold=True)
        subtitle='DIAGNOSTIC VISUALIZATION  |  Camera orbit only; original physical coordinates and sphere sizes'
        text(draw,(45,82),subtitle,21,MUTED)
        if self.kind==1:
            event,specs=self.birth_frame(index,total)
            text(draw,(45,119),f'Common event {event["event"]["event_id"]}  |  anchor basin {event["event"]["point_tracer_basin"]}',20)
        else:text(draw,(45,119),'Open cap is not a wall. Grey vessel is the original Frozen surface.',20)
        receipts=[]
        for i,p in enumerate(self.views):
            x=25+i*(950 if self.kind==0 else 635)
            label=['Point-tracer basins','Maximum passable diameter'][i] if self.kind==0 else ['A · fixed size','B · fixed anchor','C · fixed anchor + size'][i]
            text(draw,(x+12,161),label,23,INK if self.kind==0 else METHOD_COLORS[i],True)
            cam=camera(p,self.center,self.scale,20+110*index/max(1,total-1),elevation=27)
            p.render();image.paste(Image.fromarray(p.screenshot(return_img=True)[:,:,:3]),(x,205));receipts.append(cam)
            if self.kind==1:
                text(draw,(x+12,968),specs[i]['description'],17)
                text(draw,(x+12,995),specs[i]['state'],14,MUTED)
            elif self.kind==2:
                m=self.outcomes['ABC'[i]];counts=[m['outlet_counts'][b] for b in BASINS[:3]]
                text(draw,(x+12,973),f'All {m["total_integrated"]}: 01 / 02 / 03 = '+ ' / '.join(map(str,counts)),16)
        if self.kind==0:
            for j,b in enumerate(BASINS):text(draw,(50+j*460,993),b.replace('_POINT_PATH',''),20,COLORS[j],True)
        else:
            note='A/B trial jumps are resampling, not motion. C line is an entry representation, not finite-size entry dynamics.' if self.kind==1 else 'Coloured lines: completed saved paths. Grey lines: stopped partial paths. Counts refer to entire cohorts.'
            text(draw,(45,1021),note,18,MUTED)
        return image,dict(frame=index,cameras=receipts,diagnostic_states=specs,coordinates_rotated=False,
            time_role='DISPLAY_ONLY_ORBIT_AND_DIAGNOSTIC_SEARCH_PROGRESS_NO_SHARED_PHYSICAL_TIME')

    def close(self):
        for p in self.views:p.close()


def render(output,kinds=(0,1,2),frames=FRAMES):
    output=Path(output)
    for sub in ['animations','frames','inspection']:(output/sub).mkdir(parents=True,exist_ok=True)
    scene_sha=sha256(output/'scene'/'SCENE_MANIFEST.json')
    for kind in kinds:
        name=ANIMATIONS[kind];view=Views(output,kind);records=[]
        path=output/'animations'/(name+'.mp4');temp=path.with_name(path.stem+'.tmp.mp4')
        writer=imageio_ffmpeg.write_frames(str(temp),SIZE,fps=FPS,codec='libx264',pix_fmt_in='rgb24',pix_fmt_out='yuv420p',quality=8,
            macro_block_size=8,output_params=['-preset','medium','-movflags','+faststart'])
        writer.send(None);hashes={}
        try:
            for i in range(frames):
                im,receipt=view.frame(i,frames);writer.send(np.asarray(im));records.append(receipt)
                if i in [0,frames//2,frames-1]:
                    im.save(output/'inspection'/f'{name}_raw_{i:04d}.png');hashes[i]=canonical_hash(receipt)
                if i%48==0:print(name,i,'/',frames,flush=True)
        finally:writer.close();view.close()
        os.replace(temp,path)
        atomic_json(output/'frames'/(name+'.json'),dict(frame_count=frames,fps=FPS,size=SIZE,frames=records,
            source_scene_sha256=scene_sha,renderer_source_sha256=sha256(Path(__file__)),
            keyframe_receipt_hashes=hashes,physical_coordinates_rotated=False,
            ORIGINAL_SPHERE_SIZE=True,DIAGNOSTIC_VISUALIZATION=True))


def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','render','audit']);p.add_argument('--output',required=True)
    p.add_argument('--admission');p.add_argument('--trajectories');p.add_argument('--partial',action='store_true')
    p.add_argument('--kind',type=int);p.add_argument('--frames',type=int,default=FRAMES);a=p.parse_args()
    if a.action=='prepare':prepare_scene(a.admission,a.trajectories,a.output,a.partial)
    elif a.action=='render':render(a.output,(a.kind,) if a.kind is not None else (0,1,2),a.frames)
    else:
        from .particle82a_media import audit_directory
        audit_directory(a.output)


if __name__=='__main__':main()
