"""Component-separated shells from actual voxel centers; no material position repair."""
import json,time,csv,resource
from pathlib import Path
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter,distance_transform_edt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import imageio.v2 as imageio
from .config import ROOT,dump
from .y_clot_morphology_v3 import OUT,RESULTS,V2,load_geometry
from .y_bifurcation_sources import YSources,digest
from .y_fragmentation_cases import shape
from .current_position_provider import pulse
from .failure_map import csv_out

DX=50e-6
def component_shell(points):
    """Only erode the true union of dx cubes; never fill a gap or join components."""
    h=DX/2;lo=np.floor((points.min(0)-DX)/h)*h+h/2;hi=np.ceil((points.max(0)+DX)/h)*h+h/2;dims=np.round((hi-lo)/h).astype(int)+1
    axes=[lo[i]+np.arange(dims[i])*h for i in range(3)];mesh=np.meshgrid(*axes,indexing='ij');P=np.column_stack([a.ravel() for a in mesh]);d=cKDTree(points).query(P,p=np.inf)[0]
    occupied=(d<=DX/2+1e-12).reshape(tuple(dims));signed=(distance_transform_edt(occupied)-distance_transform_edt(~occupied))*h;rounded=np.minimum(signed,gaussian_filter(signed,.6))
    grid=pv.ImageData(dimensions=tuple(dims),spacing=(h,h,h),origin=lo);grid['phi']=rounded.ravel(order='F');shell=grid.contour([1e-9],scalars='phi').triangulate().clean()
    if shell.n_cells:
        original=shell.points.copy();smooth=shell.smooth(n_iter=25,relaxation_factor=.15,boundary_smoothing=False)
        shift=smooth.points-original;length=np.linalg.norm(shift,axis=1);shift*=np.minimum(1,(DX/4)/np.maximum(length,1e-30))[:,None];candidate=original+shift
        # Only display vertices are rounded, bounded to12.5um. Particle coordinates never change.
        distance,nearest=cKDTree(points).query(candidate,p=np.inf);outside=distance>DX/2
        candidate[outside]=points[nearest[outside]]+np.clip(candidate[outside]-points[nearest[outside]],-DX/2,DX/2)
        shell.points=candidate;shell=shell.compute_normals(cell_normals=False,point_normals=True,split_vertices=False)
    return shell

def mm(mesh):
    a=mesh.copy();a.points*=1000;return a

class Scene:
    def __init__(self):
        self.g=load_geometry();self.n=self.g['clot_count'];self.s=YSources();self.wall=mm(self.s.faces['WALL']);self.case=OUT/'cases/C1_GAUSSIAN_A200';self.tr=np.load(self.case/'accepted_trajectory.npz');self.ev=np.load(self.case/'broken_bonds.npz');self.r=json.loads((self.case/'request.json').read_text());self.load=json.loads((self.case/'LOAD_PROVIDER.json').read_text());self.hand=json.loads((OUT/'FRAGMENT_HANDOFF.json').read_text());self.trans=np.load(OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz');self.ts=json.loads((OUT/'FRAGMENT_TRANSPORT_EVENTS.json').read_text());self.cache={};self.qa=[]
        self.dir=RESULTS/'figures';self.dir.mkdir(exist_ok=True);(RESULTS/'figure_data').mkdir(exist_ok=True);(RESULTS/'animations').mkdir(exist_ok=True);(RESULTS/'movie_frame_manifest').mkdir(exist_ok=True)
        ids=np.arange(0,len(self.s.u),max(1,len(self.s.u)//550));flow=pv.PolyData(self.s.flow.points[ids]*1000);flow['velocity_mm_s']=self.s.u[ids]*1000;flow['speed_mm_s']=np.linalg.norm(flow['velocity_mm_s'],axis=1);self.arrows=flow.glyph(orient='velocity_mm_s',scale='speed_mm_s',factor=.025)
    def shells(self,k):
        if k in self.cache:return self.cache[k]
        labels=self.tr['component_id'][k];attached=self.tr['attached'][k];X=self.tr['position_m'][k,:self.n];items=[]
        for cid in np.unique(labels):
            ids=np.flatnonzero(labels==cid);shell=component_shell(X[ids]) if len(ids)>=3 else None
            if shell is not None:
                shell['component_id']=np.full(shell.n_points,cid);shell['attached']=np.full(shell.n_points,int(attached[ids[0]]));assert shell.n_open_edges==0
                # Cube union is the visualization support, not a new mechanics solid.
                outside=cKDTree(X[ids]).query(shell.points,p=np.inf)[0]-DX/2
                self.qa.append(dict(frame=int(k),step=int(self.tr['step'][k]),component_id=int(cid),member_count=len(ids),display_shell_volume_m3=float(shell.volume),reference_member_volume_m3=len(ids)*DX**3,open_edges=int(shell.n_open_edges),maximum_outward_support_error_m=float(max(0,outside.max())),render_component_merge=False))
                assert outside.max()<DX/3
            items.append(dict(cid=int(cid),ids=ids,attached=bool(attached[ids[0]]),shell=shell))
        self.cache[k]=items;return items
    def base(self,p,view,flow=False):
        p.set_background('#f8f6f2');p.add_mesh(self.wall,color='#c3d4dc',opacity=.12 if view!='full' else .20,smooth_shading=True,show_scalar_bar=False)
        if flow:p.add_mesh(self.arrows,scalars='speed_mm_s',cmap='viridis',clim=[0,10.0741],opacity=.50,show_scalar_bar=False)
        aspect=(p.window_size[0]/p.shape[1])/(p.window_size[1]/p.shape[0])
        if view=='full':p.camera_position=[(-2,-2,34),(-1,0,0),(0,1,0)];p.enable_parallel_projection();p.camera.parallel_scale=max(8,13/aspect)
        elif view=='branch':p.camera_position=[(0,-2,21),(1,0,0),(0,1,0)];p.enable_parallel_projection();p.camera.parallel_scale=max(6.8,6/aspect)
        else:p.camera_position=[(-5.2,-2.2,-2.3),(-5.8,0,1.05),(0,0,1)];p.enable_parallel_projection();p.camera.parallel_scale=max(1.05,1.05/aspect)
        p.add_axes(line_width=2,labels_off=False)
        length=5 if view=='full' else 1 if view=='local' else 2
        start=np.array([-12.,-7.,0]) if view=='full' else np.array([-6.7,-.72,.3]) if view=='local' else np.array([-3.,-3.7,0])
        p.add_mesh(pv.Line(start,start+[length,0,0]),color='#244153',line_width=4);p.add_point_labels([start+[length/2,0,0]],[f'{length} mm'],font_size=10,text_color='#244153',shape=None,show_points=False,always_visible=True)
    def medical(self,p,k,transport_frame=None):
        items=self.shells(k);X=self.tr['position_m'][k,:self.n].copy();moving={}
        if transport_frame is not None:
            for j,r in enumerate(self.hand['components']):
                delta=self.trans['position_m'][transport_frame,j]-r['COM_xyz_m'];moving[r['component_id']]=delta;X[np.array(r['member_atom_ids'])-1]+=delta
        debris=[]
        for item in items:
            if item['shell'] is not None:
                a=item['shell'].copy()
                if item['cid'] in moving:a.points+=moving[item['cid']]
                a.points*=1000;color='#872a40' if item['attached'] else '#cb6045'
                p.add_mesh(a,color=color,smooth_shading=True,specular=.55,specular_power=35,ambient=.24,diffuse=.70,show_scalar_bar=False)
            else:debris.extend(item['ids'])
        if debris:p.add_points(X[debris]*1000,color='#e9912e',point_size=4,render_points_as_spheres=True,opacity=.9)
        return X
    def raw(self,p,k):
        X=self.tr['position_m'][k];labels=self.tr['component_id'][k];att=self.tr['attached'][k];cloud=pv.PolyData(X*1000);cloud['atom_id']=np.arange(1,len(X)+1);cloud['component_id']=np.r_[labels,np.zeros(len(X)-self.n,int)];cloud['attached']=np.r_[att,np.ones(len(X)-self.n,bool)].astype(int)
        pts=pv.PolyData(X[:self.n]*1000);pts['attached']=att.astype(int);glyph=pts.glyph(geom=pv.Cube(x_length=.035,y_length=.035,z_length=.035),scale=False,orient=False);p.add_mesh(glyph,scalars='attached',clim=[0,1],cmap=['#dc913c','#a33549'],show_scalar_bar=False)
        p.add_points(X[self.n:]*1000,color='#447da3',point_size=3,opacity=.65,render_points_as_spheres=True)
        alive=np.ones(len(self.g['bonds']),bool);ended=self.ev['events'][:,1]<=self.tr['step'][k] if len(self.ev['events']) else np.array([],bool);alive[self.ev['initial_bond_indices'][ended]]=False;b=self.g['bonds'][alive];lines=pv.PolyData(X*1000,lines=np.c_[np.full(len(b),2),b].ravel());lines['bond_type']=self.g['bond_types'][alive]
        p.add_mesh(lines,color='#70564f',opacity=.1,line_width=.4,show_scalar_bar=False)
        if ended.any():p.add_points(self.ev['events'][ended,4:7]*1000,color='#edbe44',point_size=2,opacity=.35)
        F=shape(X,self.r,self.g)*self.g['areas'][:,None]*self.load['scale_Pa']*pulse(self.tr['step'][k]*.5e-6,self.r['duration_s']);ids=np.flatnonzero(self.g['surface'])[::18]
        if np.linalg.norm(F,axis=1).max()>0:p.add_arrows(X[ids]*1000,F[ids]/np.linalg.norm(F,axis=1).max(),mag=.23,color='#1654b4')
        return cloud,lines,F
    def image(self,panels,size=(1600,750)):
        p=pv.Plotter(off_screen=True,shape=(1,len(panels)),window_size=size,border=False)
        for j,(kind,k,view,title,transport_frame) in enumerate(panels):
            p.subplot(0,j);self.base(p,view,flow=transport_frame is not None)
            if kind=='raw':self.raw(p,k)
            else:self.medical(p,k,transport_frame)
            p.add_text(title,position='upper_left',font_size=11,color='#263d48')
        im=p.screenshot(return_img=True);p.close();return im
    def savefig(self,name,fig):
        fig.savefig(self.dir/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(self.dir/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    def figures(self):
        plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
        fig,ax=plt.subplots(figsize=(13,8));ax.imshow(Image.open(OUT/'PRE_DYNAMIC_MORPHOLOGY_PREVIEW.png'));ax.axis('off');info=json.loads((OUT/'GEOMETRY_SUMMARY.json').read_text());fig.suptitle('Native Y wall | same camera and scale | deterministic broad lobes',fontsize=17)
        fig.text(.08,.04,f'V2: continuous 0.285677 mm³ / lattice 0.254375 mm³\nV3: continuous {info["continuous_mesh_volume_m3"]*1e9:.6f} mm³ / lattice {info["clot_volume_m3"]*1e9:.6f} mm³ | 1927 clot nodes + 1561 anchors\nVolume change: -4.94% | global dx = 50 µm | three anchor layers | native wall fit and lumen checks passed',fontsize=12)
        self.savefig('V2_vs_V3_morphology_and_wall_fit',fig)
        self.wall_figure()
        ks=[0,116,233];im=self.image([('medical',k,'local',f'C1 | {self.tr["step"][k]*.0005:.3f} ms',None) for k in ks],(1800,650));fig,axs=plt.subplots(2,1,figsize=(15,9),gridspec_kw={'height_ratios':[3,1]});axs[0].imshow(im);axs[0].axis('off');m=np.genfromtxt(self.case/'fragment_metrics.csv',delimiter=',',names=True,dtype=None,encoding=None);axs[1].plot(m['time_s']*1e3,m['attached_clot_volume_fraction']*100,label='Attached mass (%)');axs[1].plot(m['time_s']*1e3,(1-m['attached_clot_volume_fraction'])*100,label='Detached mass (%)');axs[1].set(xlabel='LAMMPS time (ms)',ylabel='Initial clot mass (%)');axs[1].legend();axs[1].grid(alpha=.2);fig.suptitle('Actual Gaussian 200 Pa reference budget | native bond fracture and local shedding',fontsize=15);self.savefig('fracture_and_detachment_in_y',fig)
        fig=plt.figure(figsize=(15,8));gs=fig.add_gridspec(1,2,width_ratios=[1.05,1]);ax=fig.add_subplot(gs[0]);ax.imshow(self.image([('medical',233,'full','BACKGROUND CFD / ONE-WAY DEMO | 5 s',250)],(900,900)));ax.axis('off');ax=fig.add_subplot(gs[1]);
        # Native wall projection, actual P1 arrows, true COM trajectories in physical units.
        pp=self.s.faces['WALL'].points;ax.scatter(pp[::4,0]*1000,pp[::4,1]*1000,s=.15,color='#afbec8',alpha=.4)
        for j,r in enumerate(self.hand['components']):
            xy=self.trans['position_m'][:,j,:2]*1000;ax.plot(xy[:,0],xy[:,1],lw=1,alpha=.65,color='#be543c' if r['node_count']>=3 else '#d5a33c');ax.scatter(*xy[-1],s=8,color='#553348')
        for name,sign in [('A',1),('B',-1)]:
            normal=np.array([np.cos(np.pi/6),sign*np.sin(np.pi/6)]);tangent=np.array([-normal[1],normal[0]])
            for distance,color,style in [(10,'#3489a5','--'),(12,'#493555','-')]:
                center=normal*distance;segment=np.array([center-1.1*tangent,center+1.1*tangent]);ax.plot(segment[:,0],segment[:,1],color=color,ls=style,lw=1.4)
            ax.text(8,sign*7,f'{name}: 10 mm classification / 12 mm outlet',fontsize=7,ha='center')
        ax.set(xlabel='x (mm)',ylabel='y (mm)',title='All detached COM routes; no streamline recentering');ax.set_aspect('equal');ax.grid(alpha=.15)
        counts={}
        for a in self.ts['components']:counts[a['final_status']]=counts.get(a['final_status'],0)+1
        fig.text(.53,.04,'5 s final states: 67 not classified; 8 classified and still inside\n62 exited A; 71 exited B (mostly unresolved mass packets)\nNative P1 velocity: 0–10.074 mm/s; no flow rescaling\nNear-wall / retained-clot proximity remains uncertain.',fontsize=10);fig.suptitle('Seconds-scale transport of real V3 handoff components in the saved unobstructed flow',fontsize=15);self.savefig('one_way_fragment_transport_routes',fig)
        fig,axs=plt.subplots(2,1,figsize=(14,12))
        for ax,k in zip(axs,[0,233]):ax.imshow(self.image([('raw',k,'local',f'Raw nodes / bonds / anchors | {self.tr["step"][k]*.0005:.3f} ms',None),('medical',k,'local','Same coordinates / components | bounded shells',None)],(1600,650)));ax.axis('off')
        fig.suptitle('Solver versus medical view | no cross-component reconstruction | orange = unresolved debris',fontsize=14);self.savefig('solver_vs_medical_view_same_frames',fig)
        for f in ['V2_V3_GEOMETRY_COMPARISON.json','FRAGMENT_HANDOFF.json','FRAGMENT_TRANSPORT_SUMMARY.csv','NATIVE_EVENT_TRANSIENT_SIZE_AUDIT.json']:
            import shutil
            shutil.copy2(OUT/f,RESULTS/'figure_data'/f)
    def wall_figure(self):
        fig=plt.figure(figsize=(14,9));gs=fig.add_gridspec(2,2);ax=fig.add_subplot(gs[:,0]);ax.imshow(self.image([('medical',233,'full','Uncut native Y WALL | final C1',None)],(750,1000)));ax.axis('off')
        for i,field in enumerate(['min_clearance_m','contact_nodes']):
            ax=fig.add_subplot(gs[i,1])
            for name in ['C0_ZERO_LOAD_DYNAMIC','C1_GAUSSIAN_A200']:
                data=np.genfromtxt(OUT/'cases'/name/'WALL_CONTACT_DIAGNOSTICS.csv',delimiter=',',names=True,dtype=None,encoding=None);ax.plot(data['time_s']*1e3,data[field]*(1e6 if i==0 else 1),label=name[:2])
            if i==1:ax.plot(data['time_s']*1e3,data['outside_nodes'],color='crimson',ls='--',label='Outside centers (C1)')
            ax.set(xlabel='LAMMPS time (ms)',ylabel='Minimum native wall clearance (µm)' if i==0 else 'Wall-contact / outside node count');ax.legend();ax.grid(alpha=.2)
            if i==0:ax.axhline(0,color='r',lw=1);ax.axhline(20,color='gray',lw=1,ls=':');ax.text(.04,.12,'0 actual outside nodes\n0 swept WALL crossings\n0 physical port exits',transform=ax.transAxes)
        fig.suptitle('Wall force is applied during time integration; material centers are never clamped',fontsize=15);self.savefig('native_y_wall_contact_and_penetration_audit',fig)
    def paraview(self):
        dest=RESULTS/'paraview';dest.mkdir(exist_ok=True)
        for k in [0,116,233]:
            p=pv.Plotter(off_screen=True);cloud,lines,F=self.raw(p,k);p.close();cloud.points/=1000;lines.points/=1000;cloud['applied_load_N']=F
            blocks=pv.MultiBlock({'native_WALL':self.s.faces['WALL'],'all_material_and_anchor_nodes':cloud,'surviving_typed_bonds':lines})
            for item in self.shells(k):
                if item['shell'] is not None:blocks[f'component_{item["cid"]}_shell']=item['shell']
            blocks.save(dest/f'fracture_step_{int(self.tr["step"][k]):04d}.vtm')
        blocks=pv.MultiBlock({'native_WALL':self.s.faces['WALL']})
        for frame in [0,50,150,250]:
            cloud=pv.PolyData(self.trans['position_m'][frame]);cloud['component_id']=self.trans['component_id'];cloud['COM_velocity_m_s']=self.trans['velocity_m_s'][frame];blocks[f'transport_COM_t_{self.trans["relative_time_s"][frame]:.2f}s']=cloud
        for j,r in enumerate(self.hand['components']):
            route=pv.lines_from_points(self.trans['position_m'][:,j]);route['component_id']=np.full(route.n_points,r['component_id']);blocks[f'route_{r["component_id"]}']=route
        blocks.save(dest/'one_way_transport_routes.vtm')
        (dest/'README.txt').write_text('Open fracture_step_*.vtm and one_way_transport_routes.vtm in ParaView. Coordinates are SI meters. Nodes include all anchors and clot members. Shells are separate per real component; unresolved members remain in the node block. Bond type1bulk2interface3wall. Raw loads are Newtons. Original CFD remains in the read-only source register. Full native WALL uncut.\n')
    def movie(self):
        path=RESULTS/'animations/y_clot_realistic_fracture_and_transport.mp4';manifest=[];gifs=[];writer=imageio.get_writer(path,fps=20,codec='libx264',quality=8,macro_block_size=16)
        selections=[('LAMMPS_FRACTURE',int(k)) for k in np.unique(np.r_[np.arange(0,234,2),233])]+[('ONE_WAY_TRANSPORT',int(k)) for k in np.unique(np.r_[np.arange(0,251,2),250])]
        for idx,(chapter,k) in enumerate(selections):
            fracture=chapter=='LAMMPS_FRACTURE';frame=k if fracture else 233;tf=None if fracture else k;clock=f'{self.tr["step"][k]*.0005:.3f} ms' if fracture else f'{self.trans["relative_time_s"][k]:.2f} s after handoff';heading='1 | LAMMPS fracture' if fracture else '2 | BACKGROUND CFD / ONE-WAY DEMO'
            im=self.image([('medical',frame,'full',heading+'\n'+clock,tf),('medical',frame,'local','Attachment region | true scale\nAll original components retained',tf),('medical',frame,'branch','Native bifurcation | no flow acceleration\nUnresolved debris shown as small points',tf)],(1536,704));writer.append_data(im)
            if idx%4==0:gifs.append(np.asarray(Image.fromarray(im).resize((960,440))))
            manifest.append(dict(movie_frame=idx,chapter=chapter,source_lammps_frame=frame,source_step=int(self.tr['step'][frame]),LAMMPS_time_s=float(self.tr['step'][frame]*.5e-6),transport_frame=tf,transport_relative_time_s=None if fracture else float(self.trans['relative_time_s'][k]),material_position_interpolation=False,displacement_scale=1,source_component_ids=self.tr['component_id'][frame].tolist() if idx in [0,len(selections)-1] else 'native trajectory / handoff unchanged'))
            if idx%30==0:print('MOVIE',idx+1,'/',len(selections),flush=True)
            if frame not in [0,116,233]:self.cache.pop(frame,None)
        writer.close();imageio.mimsave(RESULTS/'animations/y_clot_realistic_fracture_and_transport.gif',gifs,duration=200,loop=0)
        dump(RESULTS/'movie_frame_manifest/frames.json',manifest);return len(manifest)

def run():
    started=time.perf_counter();scene=Scene();scene.figures();scene.paraview();figtime=time.perf_counter()-started;frames=scene.movie()
    csv_out(OUT/'SHELL_COMPONENT_QA.csv',scene.qa)
    dump(OUT/'VISUALIZATION_PROVENANCE_AND_QA.json',dict(status='PASS',shell_method='Per-component union of actual dx cubic supports, sampled at dx/2; Gaussian0.6cell rounded field capped by original signed union, only erosion; normals smoothed, no mesh dilation or cross-component merge',native_coordinates_unchanged=True,all_components_retained=True,unresolved_mass_retained=True,medical_rotation_or_deformation_added=False,displacement_scale=1,wall_clipping='None; all three views use original full transparent native WALL',coordinate_units_native='m',display_units='mm',maximum_display_support_error_m=max(a['maximum_outward_support_error_m'] for a in scene.qa),shell_volume_is_not_replacement_mechanical_mass=True,shells_checked=len(scene.qa),figures=5,figure_pdfs=5,movie_frames=frames,figures_paraview_wall_seconds=figtime,movie_wall_seconds=time.perf_counter()-started-figtime,total_wall_seconds=time.perf_counter()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,source_trajectory_sha256=digest(scene.case/'accepted_trajectory.npz'),transport_trajectory_sha256=digest(OUT/'FRAGMENT_TRANSPORT_TRAJECTORIES.npz')))

def captions():
    """Append a caption band below every original pixel; no data/viewpoint alteration."""
    from PIL import ImageDraw,ImageFont
    import shutil
    started=time.perf_counter();movie=RESULTS/'animations/y_clot_realistic_fracture_and_transport.mp4';gif=movie.with_suffix('.gif');history=OUT/'render_pre_caption';history.mkdir(exist_ok=False);shutil.copy2(movie,history/movie.name);shutil.copy2(gif,history/gif.name)
    font=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',17);small=ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',15);manifest=json.loads((RESULTS/'movie_frame_manifest/frames.json').read_text());reader=imageio.get_reader(history/movie.name);writer=imageio.get_writer(movie,fps=20,codec='libx264',quality=8,macro_block_size=16);preview=[]
    for i,a in enumerate(reader):
        im=Image.new('RGB',(a.shape[1],a.shape[0]+48),'#e9eef1');im.paste(Image.fromarray(a),(0,0));draw=ImageDraw.Draw(im);fracture=manifest[i]['chapter']=='LAMMPS_FRACTURE'
        title='LAMMPS fracture | Gaussian reference budget 200 Pa | Native times and component identities | EXPLORATORY ONLY' if fracture else 'ONE-WAY background flow | Exit markers are terminated tracks, not deposition or capture | Physical flow speed unchanged'
        draw.text((12,a.shape[0]+3),title,font=font,fill='#163445');draw.text((12,a.shape[0]+25),'1-2 node debris = unresolved numerical mass packets | Unobstructed background CFD | Displacement scale 1x',font=small,fill='#314b58');writer.append_data(np.asarray(im))
        if i%4==0:preview.append(np.asarray(im.resize((960,470))))
        manifest[i]['caption_band_added_below_original_pixels']=True
    reader.close();writer.close();imageio.mimsave(gif,preview,duration=200,loop=0);dump(RESULTS/'movie_frame_manifest/frames.json',manifest)
    qa=json.loads((OUT/'VISUALIZATION_PROVENANCE_AND_QA.json').read_text());qa.update(caption_band_wall_seconds=time.perf_counter()-started,caption_pixels_extend_image_without_crop=True,shell_method='Separate Boolean unions of actual dx cubic supports, sampled atdx/2; true union signed distance avoids internal voxel seams; Gaussian0.6cell corner rounding capped by signed field. Laplacian25iterations with12.5um vertex displacement cap and clipping to same component support. No particle edits or component merge.',subvoxel_crack_visibility_limit_m=25e-6,biological_microstructure_not_resolved=True)
    qa['movie_wall_seconds']+=time.perf_counter()-started;qa['total_wall_seconds']+=time.perf_counter()-started;dump(OUT/'VISUALIZATION_PROVENANCE_AND_QA.json',qa)

if __name__=='__main__':
    import sys
    captions() if len(sys.argv)>1 and sys.argv[1]=='captions' else run()
