"""Offline views: native saved trajectories only, cube shells without graph bridging."""
import json,csv,time,resource,xml.etree.ElementTree as ET
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image
import imageio.v2 as imageio
from .config import ROOT,dump
from .y_irregular_morphology import OUT,RESULTS,load_geometry
from .y_bifurcation_sources import digest
from .y_fragmentation_cases import CASES
from .y_scene_visualization import mm

FIG=RESULTS/'figures';DATA=RESULTS/'FIGURE_DATA';ANIM=RESULTS/'animations';PV=RESULTS/'paraview'
RED='#8d263c';BLUE='#83aabd';INK='#233746';GOLD='#dd952f';GREY='#9ba4ae'
NAMES=['Zero load','Gaussian 100 Pa','Gaussian 200 Pa','Opposing 100 Pa']
COLORS=['#708694','#3979ae','#bd3f49','#c58b23']
LOCAL=[(-3.7,-3.6,-1.6),(-6,0,1.02),(0,0,1)]
FULL=[(21,-30,32),(-1,0,0),(0,0,1)]

def saved():
    g=load_geometry();z=[]
    for name,_,_ in CASES:
        p=OUT/'cases'/name
        z.append(dict(name=name,path=p,result=json.loads((p/'result.json').read_text()),traj=dict(np.load(p/'accepted_trajectory.npz')),graph=dict(np.load(p/'lineage/graph_reconstruction.npz')),events=dict(np.load(p/'broken_bonds.npz')),rows=list(csv.DictReader((p/'lineage/timeseries.csv').open()))))
    return g,z

def plotter(size=(1000,800)):
    p=pv.Plotter(off_screen=True,window_size=size);p.set_background('#fafcfe');p.enable_anti_aliasing('ssaa');return p

def cropwall(wall):
    c=wall.cell_centers().points
    return wall.extract_cells((abs(c[:,0]+.006)<.00115)&(abs(c[:,1])<.0011)&(c[:,2]>.0005)).extract_surface()

def image(p,name):
    a=p.screenshot(return_img=True);p.close();Image.fromarray(a).save(DATA/(name+'.png'));return a

def scale(p):
    p.add_lines(np.array([[-6.70,-.70,.47],[-6.20,-.70,.47]]),color=INK,width=3)
    p.add_point_labels(np.array([[-6.45,-.70,.47]]),['0.5 mm'],font_size=11,text_color=INK,show_points=False,shape=None,always_visible=True)

def section(ax,mesh,color,label):
    cut=mesh.slice(normal=(0,1,0),origin=(0,0,0));lines=cut.lines;i=0;segs=[]
    while i<len(lines):
        n=lines[i];pts=cut.points[lines[i+1:i+n+1]][:,[0,2]]*1000;i+=n+1
        segs.extend(np.stack([pts[:-1],pts[1:]],axis=1))
    ax.add_collection(LineCollection(segs,colors=color,linewidths=1.8,label=label));ax.set(xlim=(-6.95,-5.1),ylim=(.55,1.5),xlabel='X (mm)',ylabel='Z (mm)');ax.set_aspect('equal');ax.grid(alpha=.15)

def savefig(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=240,facecolor='white');fig.savefig(FIG/(name+'.pdf'),dpi=240,facecolor='white');plt.close(fig)

def geometry_figures(g):
    wall=pv.read(PV/'vessel_wall.vtp');local=cropwall(wall)
    old=ROOT/'results/y_bifurcation_clot_placement_visual_integration_v1/paraview'
    comparison=json.loads((OUT/'V1_V2_GEOMETRY_COMPARISON.json').read_text())
    fig,ax=plt.subplots(2,2,figsize=(14,10),gridspec_kw={'height_ratios':[2.2,1]})
    fig.suptitle('A larger, irregular mural clot on the same native Y wall',fontsize=19,x=.06,ha='left')
    for j,(folder,label) in enumerate([(old,'V1'),(PV,'V2')]):
        top=pv.read(folder/'exposed_surface.vtp');solid=pv.read(folder/'clot_solid.vtp');attach=pv.read(folder/'attachment_surface.vtp');m=comparison[label]
        p=plotter();p.add_mesh(mm(local),color=BLUE,opacity=.20);p.add_mesh(mm(top),color=RED,smooth_shading=True,specular=.2)
        edge=attach.extract_feature_edges(boundary_edges=True,feature_edges=False,manifold_edges=False,non_manifold_edges=False);p.add_mesh(mm(edge),color=GOLD,line_width=3)
        p.camera_position=LOCAL;p.camera.parallel_projection=True;p.camera.parallel_scale=1.00;scale(p)
        ax[0,j].imshow(image(p,'geometry_'+label));ax[0,j].axis('off')
        inset=ax[0,j].inset_axes([.66,.67,.30,.27]);ep=edge.points
        inset.scatter((ep[:,0]+.006)*1000,1.3*np.arctan2(ep[:,1],ep[:,2]),s=.7,color=RED)
        ang=np.linspace(0,2*np.pi,150);aa=.5625 if j==0 else .68;bb=.5625 if j==0 else .62
        inset.plot(aa*np.cos(ang),bb*np.sin(ang),color='#aebac2',lw=.8,ls='--')
        inset.set(xlim=(-.82,.82),ylim=(-.82,.82));inset.set_aspect('equal');inset.set_xticks([-.6,0,.6]);inset.set_yticks([-.6,0,.6]);inset.tick_params(labelsize=6)
        inset.set_title('Unwrapped footprint (mm)',fontsize=8);inset.set_xlabel('x - x0; dashed: nominal ellipse',fontsize=6)

        ax[0,j].set_title(f"{label}: {m['continuous_mesh_volume_m3']*1e9:.4f} mm³ reference | {m['clot_volume_m3']*1e9:.4f} mm³ PD\nFootprint {m['footprint_axial_extent_m']*1e3:.3f} × {m['footprint_circumferential_arc_extent_m']*1e3:.3f} mm; height {m['maximum_mesh_height_m']*1e3:.3f} mm",fontsize=11)
        section(ax[1,j],wall,BLUE,'Native WALL');section(ax[1,j],solid,RED,'Clot boundary');ax[1,j].set_title('Same y = 0 plane; same axes');ax[1,j].legend(fontsize=9)
    fig.text(.06,.025,'Synthetic morphology, not patient anatomy. One frozen V2; wall patch and all displayed coordinates are native mesh data.',fontsize=10,color=INK)
    fig.subplots_adjust(left=.06,right=.97,top=.89,bottom=.10,wspace=.12,hspace=.16);savefig(fig,'y_irregular_clot_v1_vs_v2')

    X=g['X']*1000;b=g['bonds'];n=g['clot_count'];p=plotter((1400,1050));p.add_mesh(mm(local),color=BLUE,opacity=.12)
    p.add_points(X[:n],color=RED,point_size=4,render_points_as_spheres=True);p.add_points(X[n:],color='#39978e',point_size=4,render_points_as_spheres=True)
    for typ,col in [(1,RED),(2,GOLD),(3,'#39978e')]:
        lines=pv.PolyData();lines.points=X;bb=b[g['bond_types']==typ];lines.lines=np.c_[np.full(len(bb),2),bb];p.add_mesh(lines,color=col,opacity=.1 if typ!=2 else .26,line_width=.7)
    sm=np.load(OUT/'Y_CLOT_V2_SURFACE_NODES_AND_NORMALS.npz');ids=np.arange(0,len(sm['surface_node_id']),12)
    p.add_arrows(sm['position_xyz_m'][ids]*1000,sm['axial_tangent_xyz'][ids]*.12,mag=1,color='#3150a6')
    p.camera_position=LOCAL;p.camera.parallel_projection=True;p.camera.parallel_scale=1.1;scale(p)
    fig,ax=plt.subplots(figsize=(13,10));ax.imshow(image(p,'material_interface'));ax.axis('off');ax.set_title('Actual cubic nodes, three anchor layers, typed bonds and V2 tangents',fontsize=17,pad=12)
    fig.text(.07,.075,f"Clot: {n:,}   Anchors: {len(X)-n:,}   Bulk: {(g['bond_types']==1).sum():,}   Interface: {(g['bond_types']==2).sum():,}   Wall: {(g['bond_types']==3).sum():,}",fontsize=11)
    fig.text(.07,.045,'SI input; display in mm. Arrows show direction only. Orange: interface; teal: fixed wall proxy; red: clot.',fontsize=10)
    fig.text(.07,.022,'Native exposed triangle areas assigned once. Reference tangents frozen; spatial load envelope follows current coordinates.',fontsize=10)
    fig.subplots_adjust(top=.93,bottom=.10,left=.01,right=.99);savefig(fig,'y_irregular_clot_material_and_interface')

def timeline_figures(cases):
    fig,ax=plt.subplots(2,3,figsize=(15,9));fields=[('bulk_broken_fraction','Bulk broken (%)',100),('interface_broken_fraction','Interface broken (%)',100),('total_attached_fraction','Attached mass (%)',100),('detached_1_2node_debris_fraction','Detached 1–2 node debris (%)',100),('detached_component_count','Detached components',1),('actual_mean_traction_budget_Pa','Current mean traction budget (Pa)',1)]
    for k,case in enumerate(cases):
        rows=case['rows'];t=np.array([float(a['time_s']) for a in rows])*1e3
        for aa,(field,title,factor) in zip(ax.flat,fields):aa.plot(t,[float(a[field])*factor for a in rows],color=COLORS[k],label=NAMES[k],lw=1.6);aa.set(xlabel='Physical time (ms)',ylabel=title);aa.grid(alpha=.18)
    ax[0,0].legend(fontsize=9);fig.suptitle('Native bond breakage, attachment and released material',fontsize=19,y=.985)
    msg=' | '.join(f"{NAMES[k]}: reference / actual peak {a['result']['peak_local_traction']:.1f} / {a['result']['actual_peak_local_traction_Pa']:.1f} Pa" for k,a in enumerate(cases[1:],1))
    fig.text(.04,.035,msg,fontsize=9);fig.text(.04,.013,'Artificial analytic shear; 1.1605 ms pulse; modified 1.0dx contact. Same mean budget does not mean same peak, net force or moment.',fontsize=10)
    fig.tight_layout(rect=[.015,.075,.99,.95]);savefig(fig,'y_clot_fracture_response_comparison')
    fig,ax=plt.subplots(2,2,figsize=(14,9));summary=[]
    for k,case in enumerate(cases[1:],1):
        rows=case['rows'];t=np.array([float(a['time_s']) for a in rows])*1e3;ss=case['result']['lineage_summary']
        for j,d in enumerate(['200','300']):
            ax[0,j].step(t,[float(a['N_detached_ge'+d]) for a in rows],where='post',color=COLORS[k],label=NAMES[k]);ax[1,j].plot(t,[float(a['mass_fraction_detached_ge'+d])*100 for a in rows],color=COLORS[k])
            q=ss['thresholds'][d];summary.append(f"{NAMES[k]}, ≥{d} µm: "+(f"peak {q['MAX_SIMULTANEOUS_COUNT']}; longest exact-set span {q['longest_exact_membership_detached_span_s']*1e6:.0f} µs; final {q['AT_FINAL_FRAME']}" if q['EVER_OBSERVED'] else 'NOT OBSERVED'))
    for j,d in enumerate(['200','300']):
        ax[0,j].set_title(f'Detached ≥{d} µm (≥{34 if j==0 else 114} nodes)');ax[0,j].set_ylabel('Simultaneous count');ax[1,j].set_ylabel('Fraction of initial clot mass (%)')
        for aa in ax[:,j]:aa.set_xlabel('Physical time (ms)');aa.grid(alpha=.2)
        if not any(c['result']['lineage_summary']['thresholds'][d]['EVER_OBSERVED'] for c in cases):
            ax[0,j].text(.5,.5,'NOT OBSERVED',ha='center',transform=ax[0,j].transAxes,color=RED,fontsize=16)
            ax[0,j].set_ylim(0,1);ax[0,j].set_yticks([0,1]);ax[1,j].set_ylim(0,.1)

    ax[0,0].legend(fontsize=9);fig.suptitle('Detached fragment size and exact-membership lifetimes',fontsize=18,y=.985)
    fig.text(.055,.018,'\n'.join(summary),fontsize=9);fig.tight_layout(rect=[.02,.18,.98,.95]);savefig(fig,'y_clot_detached_fragment_sizes')

def cubes(X,ids,attached):
    """One isolated cube per saved material node; no polygons connect components."""
    nodes=pv.PolyData(X*1000);unique,inv,counts=np.unique(ids,return_inverse=True,return_counts=True)
    rgb=np.empty((len(X),3),np.uint8);rgb[:]=[218,145,40];rgb[attached]=[139,37,58];rgb[counts[inv]<=2]=[161,168,178]
    nodes['RGB']=rgb;nodes['component_id']=ids;nodes['attached']=attached.astype(np.uint8)
    return nodes.glyph(geom=pv.Cube(x_length=.05,y_length=.05,z_length=.05),orient=False,scale=False)

def movie(g,cases):
    vessel=pv.read(PV/'vessel_surface.vtp');wall=cropwall(pv.read(PV/'vessel_wall.vtp'));n=g['clot_count']
    p=pv.Plotter(shape=(2,3),off_screen=True,window_size=(1920,1088),border=False);p.enable_anti_aliasing('ssaa')
    for j,case in enumerate(cases[1:]):
        for row in range(2):
            p.subplot(row,j);p.set_background('#fafcfe');p.add_mesh(mm(vessel if row==0 else wall),color=BLUE,opacity=.17 if row==0 else .10,smooth_shading=True)
            p.add_text(NAMES[j+1]+(' | native Y context' if row==0 else ' | actual fracture close-up'),position='upper_left',font_size=12,color=INK)
            p.camera_position=FULL if row==0 else LOCAL;p.camera.parallel_projection=True;p.camera.parallel_scale=13.9 if row==0 else 1.15
            if row==0:
                p.add_point_labels(np.array([[-14,0,0],[10.4,6,0],[10.4,-6,0]]),['INLET','A','B'],font_size=11,show_points=False,shape=None,always_visible=True)
                p.add_text('Vessel = visual context only\nNo wall collision / fluid transport',position='lower_left',font_size=10,color=INK)
            else:scale(p)
    steps=cases[1]['traj']['step'];frames=[];manifest=[];keys={0,len(steps)//4,len(steps)//2,len(steps)-1}
    for c in cases:assert np.array_equal(c['traj']['step'],steps)
    with imageio.get_writer(ANIM/'y_clot_initial_fracture_synchronized.mp4',fps=20,codec='libx264',quality=8,macro_block_size=16) as writer:
        for k,step in enumerate(steps):
            for j,case in enumerate(cases[1:]):
                tr=case['traj'];mesh=cubes(tr['position_m'][k,:n],tr['component_id'][k],tr['attached'][k])
                stat=case['rows'][k]
                for row in range(2):
                    p.subplot(row,j);p.add_mesh(mesh,scalars='RGB',rgb=True,name='clot_state',show_scalar_bar=False,smooth_shading=False)
                    if row==1:p.add_text(f"t = {step*.0005:.3f} ms | attached {float(stat['total_attached_fraction'])*100:.1f}%\nRed: attached | amber: detached | grey: unresolved\nENGINEERING ANALYTIC SHEAR | displacement ×1",position='lower_left',font_size=10,color=INK,name='time_label')
            p.render();im=p.screenshot(return_img=True);writer.append_data(im)
            if k%4==0 or k==len(steps)-1:frames.append(Image.fromarray(im).resize((960,544)))
            if k in keys:Image.fromarray(im).save(DATA/f'movie_key_{k:03d}.png')
            manifest.append(dict(frame=k,step=int(step),physical_time_s=float(step*.5e-6)))
    p.close();frames[0].save(ANIM/'y_clot_initial_fracture_preview.gif',save_all=True,append_images=frames[1:],duration=200,loop=0)
    dump(DATA/'MOVIE_FRAME_MAP.json',dict(frames=manifest,fps=20,playback_duration_s=len(steps)/20,solver_frames_only=True,no_interpolation=True,display_transform='1000 * saved SI positions; same isotropic scale; displacement amplification=1',shell='Separate cube per actual dx-volume node. No edge/face generated between nodes, therefore disconnected components cannot be bridged. Component RGB is from native graph; grey 1-2 node unresolved debris retained at dx scale.',full_vessel_collision='NOT_MODELED',fluid_transport='NOT_MODELED'))

def export_paraview_and_audit(g,cases):
    n=g['clot_count'];b=g['bonds'];typ=g['bond_types'];auditimgs=[]
    for case in cases:
        folder=PV/case['name'];folder.mkdir(exist_ok=True);tr=case['traj'];ev=case['events'];stepbreak=ev['events'][:,1] if len(ev['events']) else []
        root=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');col=ET.SubElement(root,'Collection')
        for k,step in enumerate(tr['step']):
            X=tr['position_m'][k];alive=np.ones(len(b),bool);alive[ev['initial_bond_indices'][np.asarray(stepbreak)<=step]]=False
            nodes=pv.PolyData(X);nodes['node_id']=np.arange(1,len(X)+1);nodes['component_id']=np.r_[tr['component_id'][k],np.zeros(len(X)-n,int)];nodes['track_id']=np.r_[case['graph']['track_id'][k],np.zeros(len(X)-n,int)];nodes['attached']=np.r_[tr['attached'][k],np.ones(len(X)-n,bool)].astype(np.uint8);nodes['is_anchor']=g['anchors'].astype(np.uint8);nodes['velocity_m_s']=tr['velocity_m_s'][k]
            nodes.save(folder/f'points_{k:03d}.vtp')
            links=pv.PolyData();links.points=X;links.lines=np.c_[np.full(alive.sum(),2),b[alive]];links.cell_data['bond_type']=typ[alive];links.cell_data['original_bond_id']=np.flatnonzero(alive)+1;links.save(folder/f'bonds_{k:03d}.vtp')
            for part,kind in enumerate(['points','bonds']):ET.SubElement(col,'DataSet',timestep=f'{step*.5e-6:.17g}',group=kind,part=str(part),file=f'{kind}_{k:03d}.vtp')
        ET.ElementTree(root).write(folder/'trajectory.pvd',encoding='unicode',xml_declaration=True)
        # Raw native broken events preserve position at actual break time, not inferred cracks.
        events=ev['events'];q=pv.PolyData(events[:,4:7] if len(events) else np.empty((0,3)));q['time_s']=events[:,1]*.5e-6 if len(events) else np.empty(0);q['bond_type']=typ[ev['initial_bond_indices']];q.save(folder/'native_break_events.vtp')
        if case['name']=='C2_GAUSSIAN_A200':
            for k in [0,len(tr['step'])//4,len(tr['step'])//2,len(tr['step'])-1]:
                X=tr['position_m'][k]*1000;alive=np.ones(len(b),bool);alive[ev['initial_bond_indices'][np.asarray(stepbreak)<=tr['step'][k]]]=False
                p=plotter((1152,864));p.add_points(X[:n],scalars=tr['attached'][k].astype(int),cmap=[GOLD,RED],point_size=5,render_points_as_spheres=True,show_scalar_bar=False);p.add_points(X[n:],color='#39978e',point_size=3)
                lines=pv.PolyData();lines.points=X;lines.lines=np.c_[np.full(alive.sum(),2),b[alive]];p.add_mesh(lines,color=RED,opacity=.12,line_width=.6)
                if len(events):
                    sel=events[:,1]<=tr['step'][k]
                    if sel.any():p.add_points(events[sel,4:7]*1000,color='#2870bc',point_size=2,opacity=.35)
                p.camera_position=LOCAL;p.camera.parallel_projection=True;p.camera.parallel_scale=1.2
                p.add_text(f"Gaussian 200 Pa | raw solver frame {int(tr['step'][k])}\nNodes + all surviving bonds; blue = native break locations",font_size=13,color=INK);scale(p)
                auditimgs.append(image(p,f'solver_audit_{k:03d}'))
    canvas=Image.new('RGB',(2304,1728),'white')
    for i,im in enumerate(auditimgs):canvas.paste(Image.fromarray(im),((i%2)*1152,(i//2)*864))
    canvas.save(FIG/'solver_view_keyframes.png')
    (PV/'README.md').write_text('Coordinates SI metres. Open each case trajectory.pvd in ParaView: all 234 native times; actual material+anchor points, surviving typed bonds, native component/lineage IDs and velocities. vessel_surface.vtp is visual background only, never collision. native_break_events.vtp contains actual time/position of breaks. Reference clot_solid.vtp is INITIAL geometry only. Do not overlay it as deformed solid. No CFD traction or branch transport is represented.\n')

def main(stage):
    for d in [FIG,DATA,ANIM]:d.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter();g,cases=saved()
    if stage=='figures':geometry_figures(g);timeline_figures(cases)
    elif stage=='movie':movie(g,cases)
    elif stage=='paraview':export_paraview_and_audit(g,cases)
    dump(OUT/f'VISUALIZATION_PERFORMANCE_{stage.upper()}.json',dict(wall_seconds=time.perf_counter()-started,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,solver_invocations=0,stage=stage))
    print('VISUAL_DONE',stage,time.perf_counter()-started,flush=True)

if __name__=='__main__':
    import sys
    main(sys.argv[1])
