"""Standalone native-resolution whole-Y movie, level X axis and real force directions.

Read-only reuse of the accepted early-flow experiment; no solver launch.
"""
import argparse, csv, hashlib, json, subprocess, sys, time
from pathlib import Path
from itertools import product
import numpy as np
import pyvista as pv
from PIL import Image, ImageDraw, ImageFont
import imageio_ffmpeg

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pd_clot.lammps_fragmentation.point_cloud_delivery import Data, OUT as SOURCE, RESULTS as PREVIOUS, COLORS, FONT, lines
from pd_clot.lammps_fragmentation.early_flow_coupling import DTYPE, drag_law
from pd_clot.lammps_fragmentation.background_flow_sampler import BackgroundFlowSampler
from pd_clot.lammps_fragmentation.y_bifurcation_sources import digest
from pd_clot.lammps_fragmentation.config import dump
from analysis.clot_discretization_three_views_v1.render_three_views import silhouette_axes

RESULTS=PREVIOUS/'whole_vessel_4k_v1'
VERIFY=SOURCE/'whole_vessel_4k_v1'
SIZE=(3840,2160)
VIEWPORT=(160,280,3520,1540)
FPS=20
NAME='y_clot_whole_vessel_force_directions_4k'
C={**COLORS,'load':'#b57920','drag':'#357f83','flow':'#829dac'}
FONTS={s:ImageFont.truetype(FONT,s) for s in [36,38,40,42,44,48,52,64]}

class Forces:
    def __init__(self,data):
        self.d=data;self.path=data.case/'saved_frame_force_audit.bin'
        self.raw=np.memmap(self.path,dtype=DTYPE,mode='r');n=data.n
        self.offset={int(self.raw[i]['step']):i for i in range(0,len(self.raw),n)}
        # Select one actual, early detached singleton that remains in the lumen
        # through the displayed interval. The selection changes display only.
        rows=list(csv.DictReader((SOURCE/'FRAGMENT_TRANSPORT_SUMMARY.csv').open()))
        onset=np.load(data.case/'ONLINE_TOPOLOGY_EVENTS.npz')['detach_step']
        candidates=[int(r['component_id']) for r in rows if int(r['node_count'])==1 and not r['final_status'].startswith('EXIT')]
        self.atom=min(candidates,key=lambda i:(onset[i-1],i));self.j=next(j for j,r in enumerate(data.hand) if r['member_atom_ids']==[self.atom])
        self.flow=BackgroundFlowSampler(data.sources);self.history=np.loadtxt(data.case/'force_history.dat');self.max_sum_error=0.
    def state(self,r,X,attached):
        center=np.asarray(json.loads((self.d.case/'request.json').read_text())['load_centers_m'][0])
        if r['phase']!='TRANSPORT':
            i=self.offset[r['step']];a=self.raw[i:i+self.d.n]
            assert np.array_equal(a['atom_id'],np.arange(1,self.d.n+1))
            assert np.max(abs(a['x_m']-X[:self.d.n]))<1e-16
            load=a['load_N'];net=load.sum(0);weight=np.linalg.norm(load,axis=1)
            origin=np.average(X[:self.d.n],axis=0,weights=weight) if weight.sum()>0 else center
            np.testing.assert_allclose(net,self.history[r['step'],4:7],rtol=6e-6,atol=1e-14)
            self.max_sum_error=max(self.max_sum_error,float(abs(net-self.history[r['step'],4:7]).max()))
            item=a[self.atom-1];drag=item['drag_N'].copy()
            active=not attached[self.atom-1] and bool(item['valid_CFD'])
            mode='RECORDED_NATIVE_POST_FORCE';valid=bool(item['valid_CFD'])
        else:
            net=np.zeros(3);origin=center;active=bool(self.d.trans['active'][r['index'],self.j]);drag=np.zeros(3);valid=False
            if active:
                pos=self.d.trans['position_m'][r['index'],self.j];v=self.d.trans['velocity_m_s'][r['index'],self.j]
                assert np.array_equal(pos,X[self.atom-1]) or np.max(abs(pos-X[self.atom-1]))<1e-17
                q=self.flow.sample([pos]);valid=bool(q['valid'][0]);assert valid
                drag=drag_law(np.array([1]),v[None,:],q['velocity_m_s'])[0][0]
            mode='SAME_DRAG_LAW_AT_EXACT_ACCEPTED_TRANSPORT_STATE'
        return dict(load_N=net,load_origin_m=origin,drag_N=drag,drag_origin_m=X[self.atom-1],drag_active=active,valid_CFD=valid,drag_mode=mode,representative_atom_id=self.atom)

class Scene:
    def __init__(self,data):
        self.bounds=np.array([[-14.5,11.5],[-7.3,7.3],[-1.5,1.5]])
        self.focal=self.bounds.mean(1);self.out=np.array([.22,-.43,.875]);self.out/=np.linalg.norm(self.out)
        ex=np.array([1.,0,0]);self.right=ex-self.out[0]*self.out;self.right/=np.linalg.norm(self.right)
        self.up=np.cross(self.out,self.right)
        # Screen-right is the projection of world +X; its vertical component
        # is zero exactly, while Y and Z retain distinct oblique projections.
        assert abs(self.up[0])<1e-15
        corners=np.array(list(product(*self.bounds)))-self.focal
        self.scale=max(np.ptp(corners@self.up)/2,np.ptp(corners@self.right)/2*VIEWPORT[3]/VIEWPORT[2])*1.12
        p=self.p=pv.Plotter(off_screen=True,window_size=VIEWPORT[2:],lighting='three lights');p.set_background('white')
        p.render_window.SetMultiSamples(0);p.render_window.SetAlphaBitPlanes(1);p.enable_depth_peeling(number_of_peels=6,occlusion_ratio=.01)
        grid=[]
        for plane in range(3):
            others=[k for k in range(3) if k!=plane]
            for fixed in self.bounds[plane]:
                for varying in others:
                    constant=next(k for k in others if k!=varying)
                    for value in np.linspace(*self.bounds[constant],5 if constant==0 else 4):
                        e=np.zeros((2,3));e[:,plane]=fixed;e[:,constant]=value;e[:,varying]=self.bounds[varying];grid.append(e)
        grid=list({tuple(e.ravel()):e for e in grid}.values())
        p.add_mesh(lines(np.array(grid).reshape(-1,3),np.arange(2*len(grid)).reshape(-1,2)),color=C['grid'],opacity=.72,line_width=1.35,lighting=False)
        wall=data.sources.faces['WALL'].copy();wall.points*=1000
        p.add_mesh(wall,color='#bac4cc',opacity=.045,lighting=False)
        p.add_silhouette(wall,color='#8d9fab',opacity=.8,line_width=2.4)
        for name,surface in data.sources.faces.items():
            if name=='WALL':continue
            edge=surface.extract_feature_edges(boundary_edges=True,feature_edges=False,manifold_edges=False,non_manifold_edges=False);edge.points*=1000
            p.add_mesh(edge,color='#8d9fab',line_width=1.8,lighting=False)
        ids=np.arange(0,len(data.sources.u),2200);ids=ids[np.linalg.norm(data.sources.u[ids],axis=1)>.003]
        self.flow_ids=ids;arr=pv.PolyData(data.sources.flow.points[ids]*1000);arr['velocity']=data.sources.u[ids]*1000
        p.add_mesh(arr.glyph(orient='velocity',scale='velocity',factor=.12),color=C['flow'],opacity=.5,lighting=False)
        self.actors={}
        for group,size,alpha in [('clot',4.8,.95),('anchor',3.2,.6),('free',9.,1.)]:
            self.actors[group]=p.add_points(pv.PolyData([[0.,0.,0.]]),color=C[group],point_size=size,opacity=alpha,render_points_as_spheres=True,ambient=.55,diffuse=.45,specular=.1,specular_power=20)
        for group in ['load','drag']:
            self.actors[group]=p.add_mesh(pv.Arrow(),color=C[group],ambient=.55,diffuse=.45,smooth_shading=True)
        p.camera_position=[self.focal+100*self.out,self.focal,self.up];p.enable_parallel_projection();p.camera.parallel_scale=self.scale;p.show(auto_close=False)
        # Verify against actual VTK projection, rather than just the label code.
        points=[]
        for x in [-14.5,11.5]:
            p.renderer.SetWorldPoint(x,-7.3,-1.5,1);p.renderer.WorldToDisplay();points.append(p.renderer.GetDisplayPoint())
        self.x_vertical_error_px=abs(points[1][1]-points[0][1]);assert self.x_vertical_error_px<1e-8
        self.x_angle_deg=float(np.degrees(np.arctan2(points[1][1]-points[0][1],points[1][0]-points[0][0])))
    def project(self,x):
        q=np.asarray(x)-self.focal
        return np.array(VIEWPORT[:2])+[VIEWPORT[2]/2+(q@self.right)*VIEWPORT[3]/(2*self.scale),VIEWPORT[3]/2-(q@self.up)*VIEWPORT[3]/(2*self.scale)]
    def update(self,X,attached,force):
        n=len(attached);Xmm=X*1000
        for group,ids in [('clot',np.flatnonzero(attached)),('anchor',np.arange(n,len(X))),('free',np.flatnonzero(~attached))]:
            actor=self.actors[group];actor.visibility=len(ids)>0
            if len(ids):actor.mapper.dataset=pv.PolyData(Xmm[ids])
        for name,length in [('load',2.9),('drag',2.1)]:
            v=np.asarray(force[name+'_N']);mag=np.linalg.norm(v);visible=bool(mag>1e-17 and (name!='drag' or force['drag_active']))
            actor=self.actors[name];actor.visibility=visible
            if visible:
                actor.mapper.dataset=pv.Arrow(start=force[name+'_origin_m']*1000,direction=v/mag,scale=length,tip_length=.23,tip_radius=.075,shaft_radius=.022,shaft_resolution=20,tip_resolution=24)
            force[name+'_arrow_visible']=visible
            force[name+'_direction']=v/mag if visible else np.zeros(3)
        self.p.render();return self.p.screenshot(return_img=True)
    def axes(self,draw,boxes):
        def text(pos,s,size=44,anchor='mm'):
            draw.text(tuple(pos),s,font=FONTS[size],fill=C['text'],anchor=anchor);boxes.append((s,draw.textbbox(tuple(pos),s,font=FONTS[size],anchor=anchor)))
        for item in silhouette_axes(self.bounds,self.right,self.up):
            k=item['axis'];edge=np.array(item['edge']);base=edge.mean(0);normal=np.array(item['outward_screen'])*[1,-1];screen=np.array([self.project(x) for x in edge])
            draw.line(list(map(tuple,screen)),fill=C['axis'],width=2)
            for v in [[-10,-5,0,5,10],[-5,0,5],[-1,1]][k]:
                x=base.copy();x[k]=v;xy=self.project(x);draw.line([tuple(xy),tuple(xy+normal*10)],fill=C['axis'],width=2);text(xy+normal*(62 if k==0 else 78),f'{v:g}')
            pos=self.project(base)+normal*(140 if k==0 else 220)
            if k==2:pos=screen[np.argmin(screen[:,1])]+[70,-60]
            text(pos,'XYZ'[k]+' (mm)',52)

def decorate(scene,image,r,force):
    im=Image.new('RGB',SIZE,'white');im.paste(Image.fromarray(image).convert('RGB'),VIEWPORT[:2]);d=ImageDraw.Draw(im);boxes=[]
    def text(pos,s,size=44,color='text',anchor='lm'):
        d.text(tuple(pos),s,font=FONTS[size],fill=C[color],anchor=anchor);boxes.append((s,d.textbbox(tuple(pos),s,font=FONTS[size],anchor=anchor)))
    text((155,87),'Native Y vessel · force directions',64)
    text((155,176),f'Physical time  {r["time_s"]*1000:.4f} ms    |    Playback  {r["physical_seconds_per_video_second"]:.4g} s physical / s video',44,'axis')
    load=np.linalg.norm(force['load_N'])*1e6
    text((3650,91),f'Applied load  {load:.3f} µN' if force['load_arrow_visible'] else 'Applied load OFF',48,'load','rm')
    text((3650,171),f'Drag direction: node {force["representative_atom_id"]}' if force['drag_arrow_visible'] else f'Node {force["representative_atom_id"]}: drag inactive / below display threshold',38,'drag','rm')
    if r['phase']=='TRANSPORT':text((1920,243),'FROZEN REFERENCE SNAPSHOT AFTER 2.33 MS',36,'axis','mm')
    scene.axes(d,boxes)
    items=[('Attached clot','clot'),('Fixed anchors','anchor'),('Detached nodes','free'),('Applied load','load'),(f'Drag: node {force["representative_atom_id"]}','drag')]
    widths=[70+d.textlength(label,font=FONTS[42]) for label,_ in items];x=(SIZE[0]-sum(widths)-95*(len(items)-1))/2;y=2010
    for (label,color),width in zip(items,widths):
        if color in ['load','drag']:
            d.line([(x,y),(x+43,y)],fill=C[color],width=7);d.polygon([(x+53,y),(x+31,y-13),(x+31,y+13)],fill=C[color])
        else:d.ellipse((x+11,y-10,x+31,y+10),fill=C[color])
        text((x+68,y),label,42,color);x+=width+95
    text((155,2100),'UNOBSTRUCTED BACKGROUND CFD',36,'axis')
    text((3650,2100),'Force arrows: direction only · point glyphs · bonds hidden',36,'axis','rm')
    assert all(x0>=0 and y0>=0 and x1<SIZE[0] and y1<SIZE[1] for _,(x0,y0,x1,y1) in boxes),boxes
    # Labels and numerical ticks must not overlap one another.
    collisions=[]
    for i,(s,a) in enumerate(boxes):
        for t,b in boxes[i+1:]:
            if min(a[2],b[2])>max(a[0],b[0])+1 and min(a[3],b[3])>max(a[1],b[1])+1:collisions.append([s,t])
    assert not collisions,collisions
    return im

def run(stills=False):
    RESULTS.mkdir(parents=True,exist_ok=True);VERIFY.mkdir(parents=True,exist_ok=True);start=time.perf_counter()
    d=Data(SOURCE);manifest_path=PREVIOUS/'frame_manifests/y_clot_point_cloud_early_flow_coupling.json';old=json.loads(manifest_path.read_text());records=old['frames'];force=Forces(d);scene=Scene(d)
    wanted={next(r['frame'] for r in records if r['phase']=='LAMMPS' and r['step']==540):'peak_loading',next(r['frame'] for r in records if r.get('event')=='first_detachment'):'first_detachment',next(r['frame'] for r in records if r['phase']=='LAMMPS' and r['step']==1000):'load_and_drag',next(r['frame'] for r in records if r['phase']=='LAMMPS' and r['step']==4660):'handoff',records[-1]['frame']:'transport_5s'}
    process=None;movie=RESULTS/(NAME+'.mp4');temp=movie.with_suffix('.rendering.mp4');rows=[]
    if not stills:process=subprocess.Popen([imageio_ffmpeg.get_ffmpeg_exe(),'-y','-hide_banner','-loglevel','error','-f','rawvideo','-pix_fmt','rgb24','-s:v','3840x2160','-r',str(FPS),'-i','-','-an','-c:v','libx264','-preset','fast','-crf','15','-pix_fmt','yuv420p','-threads','4','-movflags','+faststart',str(temp)],stdin=subprocess.PIPE,stderr=subprocess.PIPE)
    for r in records:
        if stills and r['frame'] not in wanted:continue
        X,attached,alive,cid=d.state(r)
        assert hashlib.sha256(X.astype('<f8').tobytes()).hexdigest()==r['xyz_sha256']
        f=force.state(r,X,attached);raw=scene.update(X,attached,f);im=decorate(scene,raw,r,f)
        row={**r,**f,'x_axis_vertical_error_px':scene.x_vertical_error_px};rows.append(row)
        if process:process.stdin.write(np.asarray(im).tobytes())
        if r['frame'] in wanted:im.save(RESULTS/(NAME+'_'+wanted[r['frame']]+'.png'))
        if r['frame']%40==0:print('RENDER_4K',r['frame'],len(records),round(time.perf_counter()-start,2),flush=True)
    scene.p.close()
    if process:
        process.stdin.close();error=process.stderr.read().decode();assert process.wait()==0,error;temp.replace(movie)
    info=dict(mode='STILLS' if stills else 'MOVIE',resolution_px=SIZE,fps=FPS,frames=rows,render_seconds=time.perf_counter()-start,source_frame_manifest=str(manifest_path),source_frame_manifest_sha256=digest(manifest_path),source_force_audit=str(force.path),source_force_audit_sha256=digest(force.path),source_force_history_sha256=digest(d.case/'force_history.dat'),max_net_load_vs_native_history_error_N=force.max_sum_error,camera=dict(focal_mm=scene.focal,outward=scene.out,screen_right=scene.right,screen_up=scene.up,parallel_scale_mm=scene.scale,x_axis_angle_deg=scene.x_angle_deg,x_axis_vertical_error_px=scene.x_vertical_error_px),render_contract=dict(scene_native_resolution_px=VIEWPORT[2:],no_upscaling=True,all_six_grid_faces=True,geometry_displacement_scale=1,clot_glyph_px=4.8,anchor_glyph_px=3.2,detached_glyph_px=9.,bonds='HIDDEN_IN_WHOLE_VESSEL_VIEW',load='Actual sum of stored Gaussian nodal loads. Arrow at force-weighted node position; hidden when zero.',drag=f'Actual background force on representative material node {force.atom}, recorded during LAMMPS and evaluated using the unchanged law at exact transport states. Only shown while detached and valid.',force_arrow_lengths_mm=dict(load=2.9,drag=2.1),force_arrow_length_semantics='DIRECTION_ONLY_GRAPHIC_NOT_FORCE_MAGNITUDE',force_display_floor_N=1e-17,CFD_velocity_arrows='Sparse native velocity vectors, separate muted color, scale0.12s; not external-force arrows',native_geometry_and_flow_unchanged=True,new_physics_solves=0),source_files=old['source_files'])
    dump(VERIFY/('STILLS_MANIFEST.json' if stills else 'FRAME_MANIFEST.json'),info)
    print('COMPLETE_4K',info['render_seconds'],flush=True)

def validate():
    m=json.loads((VERIFY/'FRAME_MANIFEST.json').read_text());file=RESULTS/(NAME+'.mp4');gen=imageio_ffmpeg.read_frames(str(file),pix_fmt='rgb24');meta=next(gen);count=0;last=None;duplicate=0;sample_frames={r['frame']:r for r in m['frames'] if r['phase']=='LAMMPS' and r['step'] in [540,1000,4660]};errors=[]
    for raw in gen:
        h=hashlib.sha256(raw).hexdigest();duplicate+=int(last==h);last=h
        if count in sample_frames:
            name={540:'peak_loading',1000:'load_and_drag',4660:'handoff'}[sample_frames[count]['step']];png=np.asarray(Image.open(RESULTS/(NAME+'_'+name+'.png')).convert('RGB'),np.float32);frame=np.frombuffer(raw,np.uint8).reshape(SIZE[1],SIZE[0],3).astype(np.float32);mse=np.mean((frame-png)**2);errors.append(float(10*np.log10(255**2/max(mse,1e-30))))
        count+=1
    assert count==len(m['frames'])==501 and meta['size']==SIZE and meta['fps']==FPS and duplicate==0
    assert min(errors)>35 and abs(m['camera']['x_axis_angle_deg'])<1e-10
    assert np.all(np.diff([r['time_s'] for r in m['frames']])>0)
    for r in m['frames']:
        if r['time_s']>=.0011605:assert not r['load_arrow_visible'] and np.linalg.norm(r['load_N'])==0
        for key in ['load','drag']:
            if r[key+'_arrow_visible']:assert np.linalg.norm(np.asarray(r[key+'_direction'])-np.asarray(r[key+'_N'])/np.linalg.norm(r[key+'_N']))<1e-14
    for p in m['source_files']:assert digest(p['path'])==p['sha256']
    assert digest(m['source_force_audit'])==m['source_force_audit_sha256']
    result=dict(status='PASS',resolution=meta['size'],fps=meta['fps'],frames=count,duration_s=meta['duration'],duplicate_frames=duplicate,x_axis_angle_deg=m['camera']['x_axis_angle_deg'],x_axis_vertical_error_px=m['camera']['x_axis_vertical_error_px'],sample_PSNR_dB=errors,original_positions_unchanged=True,load_off_after_pulse=True,force_arrows_use_actual_vector_directions=True,new_solves=0,video_sha256=digest(file),video_bytes=file.stat().st_size)
    dump(VERIFY/'VIDEO_VERIFICATION.json',result);print(result,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--stills',action='store_true');p.add_argument('--validate',action='store_true');args=p.parse_args()
    if args.validate:validate()
    else:run(args.stills)
