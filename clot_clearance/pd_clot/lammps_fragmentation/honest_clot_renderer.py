"""Offline component-faithful cached scenes with one absolute physical clock."""
import csv,json,time,resource,hashlib
from pathlib import Path
import numpy as np
import pyvista as pv
from scipy.sparse import coo_matrix,diags
from scipy.optimize import brentq
from PIL import Image,ImageDraw,ImageFont
import imageio.v2 as imageio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .continuous_handoff_audit import OUT,RESULTS,PRIOR,PRIOR_RESULTS,CASE,load_transport_steps
from .y_clot_morphology_v3 import load_geometry
from .y_bifurcation_sources import YSources,digest
from .medical_fragment_renderer import component_shell
from .native_y_wall_query import NativeYWallQuery
from .config import dump
from .failure_map import csv_out

FPS=25
def material_only(mesh):
    """Smooth normals along mesh edges only. All vertices, faces, volumes stay fixed."""
    before=mesh.points.copy();faces=mesh.faces.copy();volume=mesh.volume;tri=mesh.faces.reshape(-1,4)[:,1:];edges=np.vstack([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]);row=np.r_[edges[:,0],edges[:,1]];col=np.r_[edges[:,1],edges[:,0]];A=coo_matrix((np.ones(len(row)),(row,col)),shape=(mesh.n_points,mesh.n_points)).tocsr();A=diags(1/np.maximum(np.asarray(A.sum(1)).ravel(),1))@A
    normal=np.asarray(mesh.point_data['Normals'],float).copy()
    for _ in range(28):normal=.45*normal+.55*(A@normal);normal/=np.maximum(np.linalg.norm(normal,axis=1),1e-30)[:,None]
    mesh.point_data['Normals']=normal.astype(np.float32);mesh.point_data.active_normals_name='Normals'
    x,y,z=mesh.points.T*1000;variation=.055*np.sin(4.3*x+1.7*y+2.1*z)+.035*np.sin(2.1*x-3.2*y+.5*z)
    mesh['medical_rgb']=np.clip(np.array([112.,31.,43.])[None,:]*(1+variation[:,None]),0,255).astype(np.uint8)
    assert np.array_equal(before,mesh.points) and np.array_equal(faces,mesh.faces) and volume==mesh.volume
    return mesh

class Dataset:
    def __init__(self):
        self.g=load_geometry();self.n=self.g['clot_count'];self.tr=np.load(CASE/'accepted_trajectory.npz');self.ev=np.load(CASE/'broken_bonds.npz');self.hand=json.loads((PRIOR/'FRAGMENT_HANDOFF.json').read_text())['components'];self.steps=load_transport_steps();self.s=YSources();self.query=NativeYWallQuery(self.s);self.shells={};self.qa=[];self.transport_brackets=[]
        self.handoff=self.tr['position_m'][-1,:self.n];self.final_items=self.items(233);self.by_cid={r['component_id']:r for r in self.hand}
        self.states={}
        for r in self.hand:
            cid=r['component_id'];a=self.steps[cid];self.states[cid]=dict(time=np.r_[0,a[:,0]],xyz=np.vstack([r['COM_xyz_m'],a[:,2:5]]))
        (RESULTS/'figures').mkdir(exist_ok=True);(RESULTS/'animations').mkdir(exist_ok=True);(RESULTS/'paraview').mkdir(exist_ok=True);(OUT/'shell_cache').mkdir(exist_ok=True)
    def items(self,k):
        if k in self.shells:return self.shells[k]
        X=self.tr['position_m'][k,:self.n];labels=self.tr['component_id'][k];attached=self.tr['attached'][k];items=[]
        for cid in np.unique(labels):
            ids=np.flatnonzero(labels==cid);mesh=None
            if len(ids)>=3:
                # Reuse the prior component reconstruction exactly; this task changes only shading.
                cache=OUT/'shell_cache'/f'frame_{k:03d}_component_{cid}.vtp'
                if cache.exists():mesh=pv.read(cache)
                else:
                    mesh=material_only(component_shell(X[ids]));cache.parent.mkdir(exist_ok=True);mesh['component_id']=np.full(mesh.n_points,cid);mesh.save(cache)
                mesh['medical_rgb']=np.asarray(mesh['medical_rgb'],np.uint8)
                self.qa.append(dict(source_frame=k,component_id=int(cid),member_count=len(ids),attached=bool(attached[ids[0]]),closed=mesh.n_open_edges==0,display_volume_m3=float(mesh.volume),reference_volume_m3=len(ids)*(50e-6)**3,bounding_extent_m=np.ptp(mesh.points,axis=0).tolist(),geometry_changed_by_new_material=False))
            items.append(dict(cid=int(cid),ids=ids,attached=bool(attached[ids[0]]),mesh=mesh))
        self.shells[k]=items;return items
    def transport_pose(self,t,frame_index=None):
        X=self.handoff.copy();shifts={};brackets=[];outside=[]
        for r in self.hand:
            cid=r['component_id'];s=self.states[cid];ids=np.array(r['member_atom_ids'])-1
            if t>=s['time'][-1]:lo=hi=len(s['time'])-1;alpha=0.;pos=s['xyz'][-1];terminal=bool(s['time'][-1]<5-1e-10)
            else:
                hi=int(np.searchsorted(s['time'],t,side='left'));lo=max(0,hi-1)
                if abs(s['time'][hi]-t)<1e-12:lo=hi
                alpha=0. if hi==lo else (t-s['time'][lo])/(s['time'][hi]-s['time'][lo]);pos=(1-alpha)*s['xyz'][lo]+alpha*s['xyz'][hi];terminal=False
            shift=pos-r['COM_xyz_m'];shifts[cid]=shift;X[ids]+=shift
            brackets.append(dict(render_frame=frame_index,component_id=cid,elapsed_s=t,source_step_lo=lo,source_step_hi=hi,source_time_lo_s=float(s['time'][lo]),source_time_hi_s=float(s['time'][hi]),alpha=float(alpha),terminal_track=terminal,same_exact_members=True,positions_only_display_interpolation=lo!=hi))
        # Interpolation is accepted only if the displayed member centers remain in the native domain or legal port exterior.
        free=np.flatnonzero(~self.tr['attached'][-1]);q=self.query.query(X[free]);legal=0
        for j in np.flatnonzero(~q['inside_lumen']):
            p=X[free[j]];valid=False
            for name,port in self.query.ports.items():
                d=np.dot(p-port['center'],port['outward'])
                if d>=-2e-9 and self.query.port_crossing(p-port['outward']*(max(d,0)+1e-6),p+port['outward']*2e-9):valid=True;break
            if valid:legal+=1
            else:outside.append(int(free[j]+1))
        assert not outside,dict(t=t,illegal_interpolated_nodes=outside)
        if frame_index is not None:self.transport_brackets.extend(brackets)
        return X,shifts,dict(illegal_outside_nodes=0,legal_port_exterior_nodes=legal,member_count=len(free))

class CachedScene:
    def __init__(self,data,comparison=False):
        self.d=data;self.comparison=comparison;self.p=pv.Plotter(off_screen=True,shape=(1,2) if comparison else '1|2',splitting_position=.60,window_size=(1600,780 if comparison else 800),border=False)
        self.views=['raw','local'] if comparison else ['local','full','branch'];self.actors=[{} for _ in self.views];self.markers=[];self.flow_actors=[]
        wall=data.s.faces['WALL'].copy();wall.points*=1000;patch=wall.clip_box(bounds=(-7.5,-4.3,-1.6,1.6,.4,1.6),invert=False).extract_surface(algorithm='dataset_surface')
        ids=np.arange(0,len(data.s.u),max(1,len(data.s.u)//200));arrows=pv.PolyData(data.s.flow.points[ids]*1000);arrows['u']=data.s.u[ids]*1000;arrows['speed']=np.linalg.norm(arrows['u'],axis=1);arrows=arrows.glyph(orient='u',scale='speed',factor=.018)
        for j,view in enumerate(self.views):
            self.p.subplot(j) if not comparison else self.p.subplot(0,j);self.p.set_background('#f4f0eb')
            self.p.add_mesh(wall if view!='local' else patch,color='#b4c3c0',opacity=.13 if view!='local' else .075,lighting=False,show_scalar_bar=False)
            if view in ['raw','local']:
                self.p.camera_position=[(-5.2,-2.2,-2.3),(-5.9,0,1.05),(0,0,1)];self.p.camera.parallel_scale=1.1
                for position,intensity,color in [((-7,-3,-4),.45,'#fff0e5'),((-4,3,0),.26,'#e7b3a6')]:self.p.add_light(pv.Light(position=position,focal_point=(-6,0,1),color=color,intensity=intensity,light_type='scene light'))
            elif view=='full':self.p.camera_position=[(-1.5,0,40),(-1.5,0,0),(0,1,0)];self.p.camera.parallel_scale=8.6
            else:self.p.camera_position=[(1,0,25),(1,0,0),(0,1,0)];self.p.camera.parallel_scale=3.8
            self.p.enable_parallel_projection()
            arr_actor=self.p.add_mesh(arrows,color='#437d84',opacity=.38,show_scalar_bar=False,lighting=False);arr_actor.SetVisibility(False);self.flow_actors.append(arr_actor)
            marker=pv.PolyData(np.zeros((1,3)));actor=self.p.add_points(marker,color='#683d43',point_size=2.,opacity=.32,render_points_as_spheres=False,lighting=False,show_scalar_bar=False);actor.SetVisibility(False);self.markers.append(actor)
            self.p.add_text('RAW: all nodes, surviving bonds, native break events' if view=='raw' else 'Native wall cutaway | geometry unchanged' if view=='local' else 'Full native Y | 1x coordinates' if view=='full' else 'Bifurcation | original CFD in Phase 2',position='upper_left',font_size=10,color='#344443')
            if view in ['full','branch']:
                positions=[[-13.7,0,0],[10.4,6,0],[10.4,-6,0]] if view=='full' else [[3,2,0],[3,-2,0]];labels=['inlet','A','B'] if view=='full' else ['A','B'];self.p.add_point_labels(positions,labels,show_points=False,always_visible=True,shape=None,text_color='#3c6060',font_size=9)
            # Fixed ruler in the camera plane; never across the clot.
            camera=self.p.camera;target=np.asarray(camera.focal_point);forward=target-np.asarray(camera.position);forward/=np.linalg.norm(forward);right=np.cross(forward,camera.up);right/=np.linalg.norm(right);up=np.cross(right,forward);scale=camera.parallel_scale;viewport=self.p.renderer.viewport;width=(viewport[2]-viewport[0])*self.p.window_size[0];height=(viewport[3]-viewport[1])*self.p.window_size[1];aspect=width/height;length=.5 if view in ['local','raw'] else 5 if view=='full' else 2;start=target-.72*scale*aspect*right-.78*scale*up
            self.p.add_mesh(pv.Line(start,start+length*right),color='#52635e',line_width=3);self.p.add_point_labels([start+length/2*right+.08*scale*up],[f'{length:g} mm'],font_size=9,text_color='#52635e',show_points=False,shape=None,always_visible=True)
        self.p.render();self.last_raw=None
        self.fixed_cameras=[(tuple(r.camera.position),tuple(r.camera.focal_point),tuple(r.camera.up),r.camera.parallel_scale) for r in self.p.renderers]
    def raw_update(self,j,k):
        X=self.d.tr['position_m'][k];n=self.d.n;labels=self.d.tr['component_id'][k];att=self.d.tr['attached'][k];p=self.p;p.subplot(0,j)
        cloud=pv.PolyData(X*1000);rgb=np.tile([148,39,58],(len(X),1));rgb[:n][~att]=[85,102,106];rgb[n:]=[50,118,147];cloud['rgb']=rgb.astype(np.uint8);cloud['atom_id']=np.arange(1,len(X)+1);cloud['component_id']=np.r_[labels,np.zeros(len(X)-n,int)];cloud['attached']=np.r_[att,np.ones(len(X)-n,bool)].astype(np.uint8)
        alive=np.ones(len(self.d.g['bonds']),bool);sel=self.d.ev['events'][:,1]<=self.d.tr['step'][k];alive[self.d.ev['initial_bond_indices'][sel]]=False;b=self.d.g['bonds'][alive];lines=pv.PolyData(X*1000,lines=np.c_[np.full(len(b),2),b].ravel());lines['bond_type']=self.d.g['bond_types'][alive]
        p.add_mesh(lines,color='#766a68',opacity=.075,line_width=.35,name='raw_bonds',show_scalar_bar=False,render=False);p.add_points(cloud,scalars='rgb',rgb=True,point_size=3,render_points_as_spheres=False,name='raw_nodes',show_scalar_bar=False,render=False)
        if sel.any():p.add_points(self.d.ev['events'][sel,4:7]*1000,color='#b36880',point_size=1.6,opacity=.4,name='raw_breaks',render_points_as_spheres=False,render=False)
        self.last_raw=(cloud,lines)
    def update(self,k,transport_t=None,index=None):
        d=self.d;items=d.items(k);X=d.tr['position_m'][k,:d.n];shifts={};validity=None
        if transport_t is not None:X,shifts,validity=d.transport_pose(transport_t,index)
        debris=np.concatenate([r['ids'] for r in items if r['mesh'] is None]) if any(r['mesh'] is None for r in items) else np.empty(0,int)
        for j,view in enumerate(self.views):
            self.p.subplot(0,j) if self.comparison else self.p.subplot(j)
            self.flow_actors[j].SetVisibility(transport_t is not None)
            if view=='raw':self.raw_update(j,k);continue
            actors=self.actors[j];needed=set()
            for r in items:
                if r['mesh'] is None:continue
                cid=r['cid'];needed.add(cid)
                if cid not in actors:
                    mesh=r['mesh'].copy();mesh.points*=1000;actor=self.p.add_mesh(mesh,scalars='medical_rgb',rgb=True,smooth_shading=False,ambient=.36,diffuse=.72,specular=.12,specular_power=14,show_scalar_bar=False,render=False);actor.GetProperty().SetInterpolationToPhong();actors[cid]=dict(actor=actor,frame=k)
                elif actors[cid]['frame']!=k:
                    mesh=r['mesh'].copy();mesh.points*=1000;actors[cid]['actor'].mapper.dataset=mesh;actors[cid]['frame']=k
                actors[cid]['actor'].SetPosition(*(shifts.get(cid,np.zeros(3))*1000));actors[cid]['actor'].SetVisibility(True)
            for cid in list(actors):
                if cid not in needed:self.p.remove_actor(actors[cid]['actor'],render=False);del actors[cid]
            if len(debris):self.markers[j].mapper.dataset=pv.PolyData(X[debris]*1000);self.markers[j].SetVisibility(True)
            else:self.markers[j].SetVisibility(False)
        # Updating cached actors does not invalidate the off-screen window buffer.
        # Render explicitly before every capture, including phase changes.
        self.p.render()
        assert self.fixed_cameras==[(tuple(r.camera.position),tuple(r.camera.focal_point),tuple(r.camera.up),r.camera.parallel_scale) for r in self.p.renderers]
        image=self.p.screenshot(return_img=True);return image,validity
    def close(self):self.p.close()

def caption(array,row):
    h,w=array.shape[:2];out=Image.new('RGB',(w,h+128),'#eeeae4');out.paste(Image.fromarray(array),(0,84));draw=ImageDraw.Draw(out);fp='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf';font=ImageFont.truetype(fp,22);small=ImageFont.truetype(fp,17)
    fracture=row['phase']=='LAMMPS_FRACTURE';t=row['physical_time_s'];rate=row['physical_seconds_per_playback_second'];factor=f'{1/rate:,.0f}x slow motion' if rate<1 else f'{rate:.2f}x time compression'
    title=f'LAMMPS fracture  |  physical t = {t*1000:.3f} ms' if fracture else f'ONE-WAY BACKGROUND FLOW  |  physical t = 0.002330 s + {row["transport_elapsed_s"]:.6f} s = {t:.6f} s'
    draw.text((22,12),title,font=font,fill='#273d3a');draw.text((22,47),'Playback time is not physical time  |  '+factor+'  |  fixed cameras, length scale 1x',font=small,fill='#50635e')
    footer='Surface shedding, not large-block fracture  |  Fine translucent points = numerically unresolved mass markers (not particle sizes)' if fracture else 'Background drag ON only in Phase 2  |  Attached body = saved 2.330 ms snapshot  |  Exit markers = stopped tracks, not deposits'
    draw.text((22,h+94),footer,font=small,fill='#364d49');return np.asarray(out)

def timeline(data):
    result=[];times=data.tr['step']*.5e-6
    for k,t in enumerate(times):
        dt=times[k]-times[k-1] if k else times[1]-times[0];result.append(dict(movie_frame=len(result),phase='LAMMPS_FRACTURE',source_frame=k,source_step=int(data.tr['step'][k]),physical_time_s=float(t),transport_elapsed_s=None,physical_seconds_per_playback_second=float(dt*FPS),positions_interpolated=False,topology_source='Native bonds at this actual saved step'))
    # Display clock only: pace by the saved mass-weighted COM arc length so the
    # rapid100us drag relaxation does not sit in a poorly sampled20ms gap.
    # Never change any source position, velocity, material parameter or force.
    N=250;first=times[-1]-times[-2]
    grid=np.unique(np.r_[0.,first,np.linspace(.02,5,1500),*[s['time'][s['time']<=.02] for s in data.states.values()]])
    travel=np.zeros(len(grid));mass=sum(r['mass_kg'] for r in data.hand)
    for r in data.hand:
        s=data.states[r['component_id']];arc=np.r_[0.,np.cumsum(np.linalg.norm(np.diff(s['xyz'],axis=0),axis=1))];travel+=r['mass_kg']/mass*np.interp(grid,s['time'],arc)
    assert np.all(np.diff(travel)>0),'A true collective stationary interval requires explicit depiction'
    first_travel=float(np.interp(first,grid,travel));total=float(travel[-1]);growth=brentq(lambda r:first_travel*(r**N-1)/(r-1)-total,1.000001,1.1)
    targets=np.cumsum(first_travel*growth**np.arange(N));targets[-1]=total;tt=np.interp(targets,travel,grid);tt[-1]=5.
    dump(OUT/'PLAYBACK_CLOCK_DESIGN.json',dict(method='Mass-weighted cumulative COM arc length of saved accepted transport positions. Geometric display-distance ramp starts at the native10us handoff increment; only the playback clock changes.',phase1_saved_frames=len(times),phase2_frames=N,first_transport_display_elapsed_s=float(tt[0]),distance_ramp_ratio=float(growth),initial_mass_weighted_display_displacement_m=first_travel,total_mass_weighted_transport_path_m=total,physical_time_monotonic=True,collective_stationary_intervals_skipped=False,force_velocity_or_position_rescaling=False))
    for j,t in enumerate(tt):
        dt=t-(tt[j-1] if j else 0);result.append(dict(movie_frame=len(result),phase='ONE_WAY_BACKGROUND_FLOW',source_frame=233,source_step=4660,physical_time_s=float(.00233+t),transport_elapsed_s=float(t),physical_seconds_per_playback_second=float(dt*FPS),positions_interpolated=True,topology_source='Fixed exact native handoff membership; same-component interpolation only between adjacent real accepted transport states'))
    assert np.all(np.diff([r['physical_time_s'] for r in result])>0)
    return result

def preview():
    data=Dataset();scene=CachedScene(data,comparison=True);im,_=scene.update(233);Image.fromarray(im).save(OUT/'HONEST_MATERIAL_PREVIEW.png');scene.close()

def render():
    started=time.perf_counter();data=Dataset();rows=timeline(data);scene=CachedScene(data);video=RESULTS/'animations/y_clot_fracture_to_flow_seamless_v2.mp4';assert not video.exists();writer=imageio.get_writer(video,fps=FPS,codec='libx264',quality=8,macro_block_size=16);preview_images=[];camera_state=[];pose_rows=[];frame_records=[]
    previous_pixels=None
    for i,row in enumerate(rows):
        a,validity=scene.update(row['source_frame'],row['transport_elapsed_s'],i if row['transport_elapsed_s'] is not None else None)
        row['scene_pixels_sha256']=hashlib.sha256(a.tobytes()).hexdigest()
        row['changed_scene_pixel_fraction']=None if previous_pixels is None else float(np.mean(np.any(a!=previous_pixels,axis=2)))
        previous_pixels=a.copy()
        items=data.items(row['source_frame']);row['medical_surface_member_nodes']=sum(len(r['ids']) for r in items if r['mesh'] is not None);row['medical_marker_member_nodes']=sum(len(r['ids']) for r in items if r['mesh'] is None)
        assert row['medical_surface_member_nodes']+row['medical_marker_member_nodes']==data.n
        row['fixed_cameras_verified']=True
        a=caption(a,row);writer.append_data(a)
        if i%4==0:preview_images.append(Image.fromarray(a).resize((960,557)))
        if i in [0,59,116,233,234,250,300,len(rows)-1]:Image.fromarray(a).save(RESULTS/'figures'/f'movie_key_{i:03d}.png')
        row.update(display_validation=validity,all_source_nodes_retained=True,rigid_transport_translation_only=True,new_solver_states_generated=False)
        if i%50==0:print('RENDER',i+1,'/',len(rows),'physical_t',row['physical_time_s'],flush=True)
        if i<233 and i not in [0,59,116]:data.shells.pop(row['source_frame'],None)
    writer.close();scene.close();preview_images[0].save(RESULTS/'animations/y_clot_fracture_to_flow_seamless_v2.gif',save_all=True,append_images=preview_images[1:],duration=160,loop=0,optimize=False)
    movie_seconds=time.perf_counter()-started
    # Two synchronized raw/medical views, same data, same cameras, all native graph points.
    compare=CachedScene(data,comparison=True)
    for k,name in [(0,'raw_vs_honest_medical_initial'),(233,'raw_vs_honest_medical_handoff')]:
        im,_=compare.update(k);fig,ax=plt.subplots(figsize=(15,8));ax.imshow(im);ax.axis('off');fig.suptitle(f'Raw solver and medical display | same native state at {data.tr["step"][k]*.0005:.3f} ms',fontsize=15);fig.text(.12,.06,'No cross-component merging. Surface shading only; unchanged vertices and volume.\nFine translucent points retain all unresolved members; marker size is a display choice.',fontsize=10);fig.savefig(RESULTS/'figures'/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(RESULTS/'figures'/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
        cloud,lines=compare.last_raw;cloud=cloud.copy();lines=lines.copy();cloud.points/=1000;lines.points/=1000;blocks=pv.MultiBlock({'native_WALL_uncut':data.s.faces['WALL'],'raw_all_nodes':cloud,'raw_surviving_typed_bonds':lines})
        for item in data.items(k):
            if item['mesh'] is not None:blocks[f'medical_component_{item["cid"]}']=item['mesh']
        unresolved=np.concatenate([item['ids'] for item in data.items(k) if item['mesh'] is None]) if any(item['mesh'] is None for item in data.items(k)) else np.empty(0,int)
        if len(unresolved):
            markers=pv.PolyData(data.tr['position_m'][k,unresolved]);markers['atom_id']=unresolved+1;markers['component_id']=data.tr['component_id'][k,unresolved];markers['mass_kg']=np.full(len(unresolved),1.25e-10);markers['attached']=data.tr['attached'][k,unresolved].astype(np.uint8);blocks['medical_unresolved_mass_markers']=markers
        blocks.save(RESULTS/'paraview'/f'raw_medical_step_{int(data.tr["step"][k]):04d}.vtm')
    compare.close()
    csv_out(OUT/'RENDER_TRANSPORT_SOURCE_BRACKETS.csv',data.transport_brackets);dump(OUT/'SHELL_MATERIAL_QA.json',dict(status='PASS',meshes=data.qa,material_changes_vertices=False,material_changes_faces=False,material_changes_mass=False,normal_smoothing_only_within_mesh_edges=True,unresolved_marker_point_size_pixels=2,unresolved_marker_opacity=.32,marker_spheres=False,source_reconstruction='Exact existing per-component reconstruction; newly smoothed normals/material only',no_new_mechanical_geometry=True))
    dump(OUT/'RENDER_MANIFEST.json',dict(fps=FPS,frames=rows,frame_count=len(rows),duration_s=len(rows)/FPS,source_native_trajectory_sha256=digest(CASE/'accepted_trajectory.npz'),source_native_events_sha256=digest(CASE/'broken_bonds.npz'),source_transport_csv_sha256=digest(PRIOR/'FRAGMENT_TRANSPORT_TRAJECTORIES.csv'),source_handoff_sha256=digest(PRIOR/'FRAGMENT_HANDOFF.json'),marker_membership_rule='All1-2node components retained as small display mass markers, including the anchored singleton; no mass deleted',interpolation='Only fixed-member rigid transport COM positions, linear between adjacent actual accepted steps. No LAMMPS position/topology interpolation and no velocity replacement.',static_attached_body='Saved native handoff geometry, not a continued elastic dynamics prediction',camera_motion=False,artificial_handoff_hold=False,full_Y_scene_original_native_WALL=True,local_wall_clip='Display-only original WALL cutaway; raw view and VTM retain uncut WALL',length_scale=1.))
    dump(OUT/'PERFORMANCE_RENDER_ONLY.json',dict(movie_wall_seconds=movie_seconds,total_render_and_static_wall_seconds=time.perf_counter()-started,peak_renderer_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,memory_scope='Renderer process high-water RSS; encoder child excluded',source_steps_cache_used=True,static_background_scene_instances=3,static_background_rebuilt_per_frame=False,transport_meshes_reconstructed_per_frame=False,mesh_cache_directory=str(OUT/'shell_cache'),new_LAMMPS_solves=0,new_transport_solves=0,new_CFD_solves=0))
    print('RENDER_COMPLETE',len(rows),time.perf_counter()-started,flush=True)

if __name__=='__main__':
    import sys
    preview() if len(sys.argv)>1 and sys.argv[1]=='preview' else render()
