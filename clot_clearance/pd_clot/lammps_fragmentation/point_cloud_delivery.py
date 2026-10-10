"""Default, raw-node scientific movie. Historical shell renderers stay optional."""
import argparse, json, time, subprocess
from pathlib import Path
from pd_clot.portable_paths import relocate
from itertools import product
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import cKDTree
import imageio_ffmpeg
from .config import ROOT, dump
from .y_clot_morphology_v3 import OUT as PRIOR, load_geometry
from .y_bifurcation_sources import YSources, digest
from .detached_fragment_metrics import graph_components
from analysis.clot_discretization_three_views_v1.render_three_views import silhouette_axes

TASK='y_clot_point_cloud_and_early_flow_coupling_v1'
OUT=ROOT/'verification'/TASK
RESULTS=ROOT/'results'/TASK
COLORS=dict(clot='#b85c64',anchor='#8aa7b9',free='#586d7b',bond='#718593',grid='#93a3b0',axis='#647887',text='#263746')
FONT=str(relocate('/mnt/c/Windows/Fonts/arial.ttf'))
FONTS={s:ImageFont.truetype(FONT,s) for s in [20,22,24,26,28,32,36]}
FPS=20

def lines(X,b):
    return pv.PolyData(X,lines=np.c_[np.full(len(b),2),b].ravel())

class Data:
    def __init__(self,source):
        self.source=Path(source);self.case=self.source/'cases/C1_GAUSSIAN_A200'
        self.g=load_geometry();self.n=self.g['clot_count']
        # Cache each compressed array once. NPZ indexing otherwise decompresses
        # the entire trajectory again for every component in every frame.
        with np.load(self.case/'accepted_trajectory.npz') as z:self.tr={k:z[k] for k in z.files}
        with np.load(self.case/'broken_bonds.npz') as z:self.ev={k:z[k] for k in z.files}
        self.hand=json.loads((self.source/'FRAGMENT_HANDOFF.json').read_text())['components']
        with np.load(self.source/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz') as z:self.trans={k:z[k] for k in z.files}
        self.sources=YSources()
        self.files=[self.case/'accepted_trajectory.npz',self.case/'broken_bonds.npz',self.source/'FRAGMENT_HANDOFF.json',self.source/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz']
    def schedule(self):
        records=[dict(phase='LAMMPS',index=k,time_s=float(s*.5e-6),step=int(s)) for k,s in enumerate(self.tr['step'])]
        for event in ['first_break','first_detachment']:
            p=self.case/(event+'_state.npz')
            if p.exists():
                self.files.append(p);z=np.load(p);step=int(z['step'])
                if step not in self.tr['step']:records.append(dict(phase='LAMMPS_EVENT',index=None,time_s=step*.5e-6,step=step,event=event))
                else:
                    next(r for r in records if r['step']==step)['event']=event
        records.sort(key=lambda r:r['time_s'])
        times=self.trans['relative_time_s']
        # Exact stored states only. Dense early new-transport samples, then 40ms.
        chosen=[k for k,t in enumerate(times) if t>0 and (t<.02 or k==len(times)-1 or abs(t/.04-round(t/.04))<1e-7)]
        records += [dict(phase='TRANSPORT',index=k,time_s=float(.00233+times[k]),step=4660) for k in chosen]
        for k,r in enumerate(records):
            r['frame']=k;r['video_time_s']=k/FPS
            r['physical_seconds_per_video_second']=(r['time_s']-records[k-1]['time_s'])*FPS if k else (records[1]['time_s']-r['time_s'])*FPS
        assert np.all(np.diff([r['time_s'] for r in records])>0)
        return records
    def state(self,r):
        if r['phase']=='LAMMPS_EVENT':
            z=np.load(self.case/(r['event']+'_state.npz'));X=z['position_m'].copy();attached=z['attached'].copy();alive=z['alive'].copy();cid=z['component_id']
            l,c,a=graph_components(self.n,self.g['bonds'],self.g['bond_types'],alive)
            assert np.array_equal(a[l],attached)
            return X,attached,alive,cid
        k=r['index'] if r['phase']=='LAMMPS' else -1
        X=self.tr['position_m'][k].copy();attached=self.tr['attached'][k].copy();cid=self.tr['component_id'][k]
        alive=np.ones(len(self.g['bonds']),bool);alive[self.ev['initial_bond_indices'][self.ev['events'][:,1]<=r['step']]]=False
        lab,count,a=graph_components(self.n,self.g['bonds'],self.g['bond_types'],alive)
        stable=np.full(len(count),self.n+1);np.minimum.at(stable,lab,np.arange(1,self.n+1))
        assert np.array_equal(stable[lab],cid) and np.array_equal(a[lab],attached)
        if r['phase']=='TRANSPORT':
            for j,item in enumerate(self.hand):
                ids=np.array(item['member_atom_ids'])-1
                assert self.trans['component_id'][j]==item['component_id'] and not attached[ids].any()
                X[ids]+=self.trans['position_m'][r['index'],j]-item['COM_xyz_m']
        return X,attached,alive,cid

class View:
    def __init__(self,d,full=False):
        self.full=full;self.size=(640,790) if full else (1190,790)
        self.offset=np.array([1240,150] if full else [25,150])
        self.bounds=np.array([[-14.5,11.5],[-7.3,7.3],[-1.5,1.5]]) if full else np.array([[-6.8,-5.2],[-.75,.85],[.55,1.5]])
        self.focal=self.bounds.mean(1)
        az,el=np.deg2rad([-88.,65.] if full else [-72.3498758,-20.])
        self.out=np.array([np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),np.sin(el)])
        self.right=np.cross(-self.out,[0.,0.,1.]);self.right/=np.linalg.norm(self.right);self.up=np.cross(self.right,-self.out)
        corners=np.array(list(product(*self.bounds)))-self.focal
        self.scale=max(np.ptp(corners@self.up)/2,np.ptp(corners@self.right)/2*self.size[1]/self.size[0])*1.27
        p=self.p=pv.Plotter(off_screen=True,window_size=self.size,lighting='three lights');p.set_background('white')
        p.render_window.SetMultiSamples(0);p.render_window.SetAlphaBitPlanes(1);p.enable_depth_peeling(number_of_peels=10,occlusion_ratio=.01)
        grid=[]
        for plane in range(3):
            others=[k for k in range(3) if k!=plane]
            for fixed in self.bounds[plane]:
                for varying in others:
                    constant=next(k for k in others if k!=varying)
                    for value in np.linspace(*self.bounds[constant],4):
                        edge=np.zeros((2,3));edge[:,plane]=fixed;edge[:,constant]=value;edge[:,varying]=self.bounds[varying];grid.append(edge)
        grid=list({tuple(e.ravel()):e for e in grid}.values())
        p.add_mesh(lines(np.array(grid).reshape(-1,3),np.arange(2*len(grid)).reshape(-1,2)),color=COLORS['grid'],opacity=.65 if full else .9,line_width=.8,lighting=False)
        wall=d.sources.faces['WALL'].copy();wall.points*=1000
        if not full:wall=wall.clip_box(self.bounds.ravel(),invert=False).extract_surface(algorithm='dataset_surface')
        p.add_mesh(wall,color='#bac4cc',opacity=.035,lighting=False)
        if full:
            p.add_silhouette(wall,color='#a9b4bd',line_width=1,opacity=.55)
            for name,surface in d.sources.faces.items():
                if name=='WALL':continue
                edge=surface.extract_feature_edges(boundary_edges=True,feature_edges=False,manifold_edges=False,non_manifold_edges=False);edge.points*=1000
                p.add_mesh(edge,color='#93a3b0',line_width=1,lighting=False)
            # Sparse true native CFD vectors, one fixed scale, no invented routes.
            ids=np.arange(0,len(d.sources.u),1500);u=d.sources.u[ids];keep=np.linalg.norm(u,axis=1)>.003;ids=ids[keep]
            arr=pv.PolyData(d.sources.flow.points[ids]*1000);arr['velocity']=d.sources.u[ids]*1000
            p.add_mesh(arr.glyph(orient='velocity',scale='velocity',factor=.07),color='#7996a5',opacity=.65,lighting=False)
            roi=pv.Box(bounds=(-6.8,-5.2,-.75,.85,.55,1.5)).outline()
            p.add_mesh(roi,color=COLORS['clot'],line_width=1.8,lighting=False)
        self.actors={}
        for group,size,alpha in [('clot',4.9,.94),('anchor',3.5,.57),('free',6.,.98)]:
            if full:size={'clot':3.,'anchor':2.,'free':4.5}[group]
            self.actors[group]=p.add_points(pv.PolyData([[0.,0.,0.]]),color=COLORS[group],point_size=size,opacity=alpha,render_points_as_spheres=True,ambient=.55,diffuse=.45,specular=.1,specular_power=20)
        if not full:
            for typ,alpha in [(1,.021),(2,.019),(3,.010),('emphasis',.40)]:
                self.actors[typ]=p.add_mesh(pv.Line([0,0,0],[1,0,0]),color=COLORS['bond'] if typ=='emphasis' else '#9cabb5',opacity=alpha,line_width=1. if typ=='emphasis' else .7,lighting=False)
            X=d.g['X']*1000;proj=np.c_[X@self.right,X@self.up];depth=X@self.out;near=cKDTree(proj).query_ball_point(proj,.055)
            front=np.array([depth[js].max()-depth[i] for i,js in enumerate(near)])<.045
            b=d.g['bonds'];self.emphasis=front[b].all(1)&(d.g['bond_types']<=2)&(np.linalg.norm(X[b[:,1]]-X[b[:,0]],axis=1)<=.085)
        p.camera_position=[self.focal+100*self.out,self.focal,self.up];p.enable_parallel_projection();p.camera.parallel_scale=self.scale
        p.show(auto_close=False)
    def project(self,x):
        q=np.asarray(x)-self.focal
        return self.offset+np.array([self.size[0]/2+(q@self.right)*self.size[1]/(2*self.scale),self.size[1]/2-(q@self.up)*self.size[1]/(2*self.scale)])
    def update(self,d,X,attached,alive):
        X=X*1000;n=d.n;inside=np.ones(len(X),bool) if self.full else np.all((X>=self.bounds[:,0])&(X<=self.bounds[:,1]),axis=1)
        groups={'anchor':inside&(np.arange(len(X))>=n),'clot':inside&np.r_[attached,np.zeros(len(X)-n,bool)],'free':inside&np.r_[~attached,np.zeros(len(X)-n,bool)]}
        for name,mask in groups.items():
            a=self.actors[name];a.visibility=bool(mask.any())
            if mask.any():a.mapper.dataset=pv.PolyData(X[mask])
        if not self.full:
            b=d.g['bonds'];keep=alive&inside[b].all(1)
            for typ in [1,2,3,'emphasis']:
                mask=keep&(self.emphasis if typ=='emphasis' else d.g['bond_types']==typ);a=self.actors[typ];a.visibility=bool(mask.any())
                if mask.any():a.mapper.dataset=lines(X,b[mask])
        self.p.render()
        return self.p.screenshot(return_img=True),int((inside[:n]&~attached).sum()),int(inside[:n].sum())
    def axes(self,draw):
        for item in silhouette_axes(self.bounds,self.right,self.up):
            k=item['axis'];edge=np.array(item['edge']);base=edge.mean(0);normal=np.array(item['outward_screen'])*[1,-1];screen=np.array([self.project(p) for p in edge])
            draw.line(list(map(tuple,screen)),fill=COLORS['axis'],width=1)
            ticks=([[-6.8,-6,-5.2],[-.5,0,.5],[.6,1,1.4]][k]) if not self.full else ([-10,0,10] if k==0 else [-5,0,5] if k==1 else [0])
            for v in ticks:
                x=base.copy();x[k]=v;q=self.project(x)
                draw.line([tuple(q),tuple(q+normal*5)],fill=COLORS['axis'],width=1)
                draw.text(tuple(q+normal*(29 if k!=2 else 40)),f'{v:g}',font=FONTS[22],fill=COLORS['text'],anchor='mm')
            if self.full and k==1:pos=screen[np.argmin(screen[:,1])]+[22,-33]
            elif k==2:pos=screen[np.argmin(screen[:,1])]+[-25,-70] if self.full else screen[np.argmin(screen[:,1])]+[12,-30]
            else:pos=self.project(base)+normal*(69 if not self.full else 62)
            draw.text(tuple(pos),'XYZ'[k]+' (mm)',font=FONTS[24],fill=COLORS['text'],anchor='mm')

def contract():
    (OUT/'POINT_CLOUD_RENDERING_CONTRACT.md').write_text('''# Default point-cloud rendering

DEFAULT_RENDERER = pd_clot.lammps_fragmentation.point_cloud_delivery
SMOOTH_SHELL_DEFAULT = DISABLED. Historical medical outputs are preserved, optional only.
The raw 1927 clot atom IDs and 1561 fixed anchor IDs are retained. Current surviving-bond topology determines attached carmine nodes and detached grey nodes. Anchor glyphs are pale cyan-grey. Every node uses native xyz converted from m to mm, with displacement scale 1. No shell, synthetic crack, merged fragments or orange balls.
Local glyph diameters: clot 4.9 px, anchor 3.5 px, detached 6 px. Global view: 3, 2, 4.5 px. These fixed screen symbols are NOT physical hydrodynamic diameters. All three surviving bond types are drawn faintly inside the local box; a fixed subset of true initial short front bonds receives emphasis while it survives. Global-view bonds are hidden and explicitly labelled. No edge is invented.
Both cameras, six-face coordinate grids and zoom remain fixed. Local cropping is explicit: each manifest records omitted local nodes; the whole-Y view includes all original IDs. Context is the original native open Y WALL at opacity .035, never an opaque medical vessel. Sparse arrows show the actual native velocity field at a documented constant display scale (0.07 s); the field is UNOBSTRUCTED BACKGROUND CFD.
Frames use exact saved native states, or exact accepted common-time transport states plus the model's retained member offsets. No display interpolation, motion rescaling or false timing. Physical absolute time and instantaneous physical-time/playback-time rate appear in every frame. At phase 2 the attached core is explicitly FROZEN REFERENCE SNAPSHOT AFTER 2.33 MS. Source hashes, frame indices, topology/coordinates hashes and ID counts accompany every frame.
''')

def render(source=PRIOR,name='y_clot_point_cloud_only_preview',stills=False):
    started=time.perf_counter();d=Data(source);schedule=d.schedule();contract()
    for folder in ['animations','figures','frame_manifests']:(RESULTS/folder).mkdir(exist_ok=True)
    views=[View(d),View(d,True)];records=[];gif=[];selected={0:'initial',next(i for i,r in enumerate(schedule) if r['phase']=='LAMMPS' and r['step']==4660):'handoff',len(schedule)-1:'transport_5s'}
    if source!=PRIOR:
        for i,r in enumerate(schedule):
            if 'event' in r:selected[i]=r['event']+'_exact_native_step'
    movie=RESULTS/'animations'/f'{name}.mp4';temp=movie.with_suffix('.rendering.mp4');proc=None
    if not stills:
        proc=subprocess.Popen([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-hide_banner','-loglevel','error','-f','rawvideo','-pix_fmt','rgb24','-s:v','1920x1080','-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','medium','-crf','17','-pix_fmt','yuv420p','-threads','4','-movflags','+faststart',str(temp)],stdin=subprocess.PIPE,stderr=subprocess.PIPE)
    for r in schedule:
        if stills and r['frame'] not in selected:continue
        X,attached,alive,cid=d.state(r);im=Image.new('RGB',(1920,1080),'white');dr=ImageDraw.Draw(im)
        for view in views:
            frame,free,visible=view.update(d,X,attached,alive);im.paste(Image.fromarray(frame).convert('RGB'),tuple(view.offset));view.axes(dr)
            if not view.full:local_free,local_visible=free,visible
        def text(pos,s,size=26,color='text',anchor='lm'):dr.text(pos,s,font=FONTS[size],fill=COLORS[color],anchor=anchor)
        text((50,42),'Clot discretization · raw material points',36)
        text((50,88),f'Physical time  {r["time_s"]*1000:.4f} ms     |     Playback: {r["physical_seconds_per_video_second"]:.4g} s physical / s video',26,'axis')
        text((55,140),'Local view',28);text((1260,140),'Native Y vessel',28)
        text((1260,915),'UNOBSTRUCTED BACKGROUND CFD',22,'axis')
        text((55,950),f'{int(attached.sum())} attached + {int((~attached).sum())} detached = 1927 clot nodes',24)
        text((55,987),f'Local box: {local_visible}/1927 nodes; all IDs in whole-Y view',22,'axis')
        text((1260,954),'Bonds: all types, local view only',22,'axis')
        text((1260,986),'Point sizes are screen symbols',22,'axis')
        if r['phase']=='TRANSPORT':text((960,119),'FROZEN REFERENCE SNAPSHOT AFTER 2.33 MS',22,'axis','mm')
        x=150
        for label,color in [('Attached clot','clot'),('Fixed anchors','anchor'),('Detached nodes','free'),('Surviving bonds','bond')]:
            if color=='bond':dr.line([(x,1038),(x+28,1038)],fill=COLORS[color],width=2)
            else:dr.ellipse((x+7,1031,x+21,1045),fill=COLORS[color])
            text((x+40,1038),label,26);x+=390
        import hashlib
        rec={**r,'source':str(d.case/(r['event']+'_state.npz') if r['phase']=='LAMMPS_EVENT' else d.case if r['phase']=='LAMMPS' else d.source/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz'),'original_clot_ids':[1,1927],'anchor_ids':[1928,3488],'attached':int(attached.sum()),'detached':int((~attached).sum()),'local_visible_clot':local_visible,'local_visible_detached':local_free,'surviving_bonds':int(alive.sum()),'xyz_sha256':hashlib.sha256(X.astype('<f8').tobytes()).hexdigest(),'component_id_sha256':hashlib.sha256(cid.astype('<i8').tobytes()).hexdigest(),'no_display_interpolation':True}
        records.append(rec)
        if proc:proc.stdin.write(np.asarray(im).tobytes())
        if r['frame'] in selected:
            path=RESULTS/'figures'/f'{name}_{selected[r["frame"]]}.png';im.save(path);im.save(path.with_suffix('.pdf'),resolution=150)
        # Short new-case GIF shows late fracture through the dense handoff window,
        # at the SAME playback rate as the MP4 (no hidden speedup).
        if not stills and r['frame']%4==0 and (source==PRIOR or 3.5<=r['video_time_s']<=13.5):gif.append(im.resize((960,540),Image.Resampling.LANCZOS))
        if r['frame']%40==0:print('RENDER',name,r['frame'],len(schedule),round(time.perf_counter()-started,1),flush=True)
    for view in views:view.p.close()
    if proc:
        proc.stdin.close();err=proc.stderr.read().decode();assert proc.wait()==0,err;temp.replace(movie)
        gif[0].save(movie.with_suffix('.gif'),save_all=True,append_images=gif[1:],duration=200,loop=0,optimize=False)
    suffix='_stills' if stills else ''
    dump(RESULTS/'frame_manifests'/f'{name}{suffix}.json',dict(frames=records,fps=FPS,source_files=[dict(path=str(p),sha256=digest(p)) for p in d.files],camera=[dict(bounds_mm=v.bounds,parallel_scale_mm=v.scale,outward=v.out) for v in views],render_seconds=time.perf_counter()-started))
    print('COMPLETE',name,'seconds',time.perf_counter()-started,flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--new',action='store_true');a.add_argument('--stills',action='store_true');args=a.parse_args()
    render(OUT if args.new else PRIOR,'y_clot_point_cloud_early_flow_coupling' if args.new else 'y_clot_point_cloud_only_preview',args.stills)
