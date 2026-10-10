"""Offline medical-style and raw-audit views of the same saved SI scene."""
import json
import math
import resource
import time
from pathlib import Path
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import imageio.v2 as imageio
from .config import dump
from .y_bifurcation_sources import OUT,RESULTS,digest
from .y_scene_task import load_geometry
from .y_background_flow_view import prepare_background

RED='#7e1428';BLUE='#628cac';INK='#243347';TEAL='#168b8c'
CAMERA=[(30,-38,35),(-1,0,0),(0,0,1)]


def mm(mesh):
    x=mesh.copy();x.points=np.asarray(x.points,float)*1000;return x


def plotter(size=(1500,1000)):
    p=pv.Plotter(off_screen=True,window_size=size);p.set_background('#fbfcfe')
    p.enable_anti_aliasing('ssaa');return p


def axes(p,bounds=None):
    p.show_bounds(bounds=bounds,grid='back',location='outer',all_edges=False,
                  xtitle='X (mm)',ytitle='Y (mm)',ztitle='Z (mm)',font_size=15,
                  color=INK,n_xlabels=5,n_ylabels=5,n_zlabels=3,fmt='%.1f',ticks='outside')
    p.add_axes(line_width=2,labels_off=False)


def scene(p,vessel,solid,flow=None,labels=True):
    p.add_mesh(mm(vessel),color=BLUE,opacity=.13,smooth_shading=True,show_scalar_bar=False)
    p.add_mesh(mm(solid),color=RED,smooth_shading=True,specular=.20,specular_power=20,show_scalar_bar=False)
    if flow is not None:
        f=mm(flow);f['speed_mm_s']=np.asarray(flow['speed_m_s'])*1000
        p.add_mesh(f,scalars='speed_mm_s',cmap='viridis',clim=(0,10.0740172),line_width=3,
                   scalar_bar_args=dict(title='Original velocity (mm/s)',vertical=False,position_x=.27,position_y=.025,width=.5,height=.06,title_font_size=16,label_font_size=14))
        arrows=np.load(RESULTS/'FIGURE_DATA/flow_direction_arrows.npz')
        vec=arrows['velocity_m_s'];vec=vec/np.linalg.norm(vec,axis=1)[:,None]*1.4
        p.add_arrows(arrows['position_xyz_m']*1000,vec,mag=1.,color=INK)
    if labels:
        p.add_point_labels(np.array([[-14,0,0],[10.392,6,0],[10.392,-6,0]]),['INLET','OUTLET A','OUTLET B'],
                           point_size=5,font_size=18,text_color=INK,shape_color='white',shape_opacity=.85,always_visible=True)
        p.add_point_labels(np.array([[-6,0,1.3]]),['Mural clot'],point_size=6,font_size=17,text_color=RED,shape_color='white',shape_opacity=.85,always_visible=True)
    p.camera_position=CAMERA;p.camera.zoom(1.05);axes(p)
    p.add_lines(np.array([[-12.5,-2.4,0.],[-10.5,-2.4,0.]]),color=INK,width=3)
    p.add_point_labels(np.array([[-11.5,-2.4,0.]]),['2 mm'],show_points=False,font_size=14,text_color=INK,shape=None,always_visible=True)


def capture(p,name):
    p.screenshot(RESULTS/'FIGURE_DATA'/name);im=p.screenshot(return_img=True);p.close();return im


def savefig(fig,name):
    fig.savefig(RESULTS/(name+'.png'),dpi=300,facecolor=fig.get_facecolor())
    fig.savefig(RESULTS/(name+'.pdf'),dpi=300,facecolor=fig.get_facecolor())
    plt.close(fig)


def basefig(title,subtitle,shape=(1,1),size=(13,9)):
    fig,ax=plt.subplots(*shape,figsize=size,squeeze=False)
    fig.patch.set_facecolor('#fbfcfe')
    fig.suptitle(title,x=.045,y=.975,ha='left',fontsize=20,fontweight='bold',color=INK)
    fig.text(.045,.925,subtitle,ha='left',fontsize=10.5,color=INK)
    fig.subplots_adjust(left=.035,right=.97,bottom=.07,top=.885,wspace=.13)
    return fig,ax


def figures():
    pdir=RESULTS/'paraview';solid=pv.read(pdir/'clot_solid.vtp');vessel=pv.read(pdir/'vessel_surface.vtp');wall=pv.read(pdir/'vessel_wall.vtp')
    top=pv.read(pdir/'exposed_surface.vtp');attach=pv.read(pdir/'attachment_surface.vtp');g=load_geometry()
    info=g['info'];fit=json.loads((OUT/'GEOMETRY_AND_WALL_FIT_CHECK.json').read_text());clear=json.loads((OUT/'LUMEN_INTERSECTION_AND_CLEARANCE.json').read_text())
    p=plotter();scene(p,vessel,solid);im=capture(p,'overview_native_scene.png')
    fig,aa=basefig('An attached clot in the solved Y-vessel geometry',
                   'Idealized Y-vessel  |  Curved mural clot (synthetic demo morphology)  |  Static geometry; true relative scale')
    aa[0,0].imshow(im);aa[0,0].axis('off')
    fig.text(.05,.045,'Mother: L = 14 mm, D = 2.6 mm     Daughters: L = 12 mm, D = 2.0 mm     Branch angle: ±30°',fontsize=11,color=INK)
    fig.text(.05,.020,'Clot centerline projection: (−6, 0, 0) mm     Attachment: +Z inner wall     Nominal a = b = H = 0.5625 mm',fontsize=10,color=RED)
    savefig(fig,'y_bifurcation_clot_overview')

    centers=wall.cell_centers().points
    local=wall.extract_cells((abs(centers[:,0]+.006)<.0010)&(abs(centers[:,1])<.0010)&(centers[:,2]>.0005)).extract_surface().triangulate()
    p=plotter((1050,1000));p.add_mesh(mm(local),color=BLUE,opacity=.2,show_edges=True,line_width=.4)
    p.add_mesh(mm(top),color=RED,smooth_shading=True,opacity=1)
    p.add_mesh(mm(attach.extract_feature_edges(boundary_edges=True,feature_edges=False,non_manifold_edges=False,manifold_edges=False)),color='#ffb44c',line_width=4)
    p.camera_position=[(-3.7,-3.2,-1.5),(-6,0,1.1),(0,0,1)];p.camera.parallel_projection=True;p.camera.parallel_scale=1.55
    axes(p,bounds=(-7,-5,-.9,.9,.5,1.5));im=capture(p,'wall_fit_native_scene.png')
    fig,aa=basefig('Curved attachment, checked against the actual wall',
                   'Boundary-fitted reference solid and its cubic material discretization are reported separately',shape=(1,2),size=(14,8))
    aa[0,0].imshow(im);aa[0,0].axis('off')
    ax=aa[0,1]
    def plot_lines(mesh,color,lw,label):
        cut=mesh.slice(normal=(0,1,0),origin=(0,0,0));lines=cut.lines;i=0;first=True
        while i<len(lines):
            n=lines[i];q=cut.points[lines[i+1:i+1+n]]*1000;i+=n+1
            if len(q)>1:ax.plot(q[:,0],q[:,2],color=color,lw=lw,label=label if first else None);first=False
    plot_lines(wall,BLUE,2,'Native WALL section');plot_lines(solid,RED,2,'Fitted clot boundary')
    X=g['X'][:g['clot_count']];A=g['X'][g['clot_count']:]
    bx=X[abs(X[:,1])<1e-10]*1000;ba=A[abs(A[:,1])<1e-10]*1000
    ax.scatter(bx[:,0],bx[:,2],s=11,color=RED,label='Clot lattice points');ax.scatter(ba[:,0],ba[:,2],s=14,color=TEAL,marker='s',label='Fixed wall anchors')
    ax.set(xlim=(-6.85,-5.15),ylim=(.55,1.55),xlabel='X (mm)',ylabel='Z (mm)',title='Actual y = 0 section');ax.set_aspect('equal');ax.grid(alpha=.2);ax.legend(fontsize=9,loc='lower left')
    fig.text(.05,.06,f"Footpoint correction: {fit['nominal_footpoint_correction_m']*1e6:.3f} µm   |   Max attachment gap: {fit['attachment_vertex_gap_max_m']*1e6:.6f} µm",fontsize=11,color=INK)
    fig.text(.05,.028,f"Reference volume: {info['continuous_mesh_volume_m3']*1e9:.4f} mm³   |   PD quadrature volume: {info['clot_volume_m3']*1e9:.4f} mm³   |   Center-section occupancy: {clear['center_section_area_occupancy']*100:.2f}%",fontsize=11,color=INK)
    savefig(fig,'curved_mural_clot_wall_fit')

    flow=pv.read(pdir/'background_streamlines_masked.vtp');p=plotter();scene(p,vessel,solid,flow);im=capture(p,'background_native_scene.png')
    fig,aa=basefig('Original background flow with the new clot overlay',
                   'UNOBSTRUCTED BACKGROUND CFD — NOT CLOTTED FLOW FIELD')
    aa[0,0].imshow(im);aa[0,0].axis('off')
    overlay=json.loads((OUT/'BACKGROUND_FLOW_OVERLAY_CHECK.json').read_text())
    fig.text(.05,.048,'Original unobstructed background CFD; not recomputed with clot. No clot-surface traction is inferred.',fontsize=11,fontweight='bold',color=RED)
    fig.text(.05,.023,f"Original Qin = 1.58124 mL/min; A/B = 49.9964/50.0036%   |   {overlay['streamlines_terminated_at_clot']} displayed lines stop at the solid; no deflection is simulated.",fontsize=10,color=INK)
    savefig(fig,'original_background_flow_with_clot_overlay')

    p=plotter((1400,1050));p.add_mesh(mm(local),color=BLUE,opacity=.12,show_edges=True,line_width=.3)
    Xall=g['X']*1000;n=g['clot_count'];types=g['bond_types'];b=g['bonds']
    p.add_points(Xall[:n],color=RED,point_size=5,render_points_as_spheres=True)
    p.add_points(Xall[n:],color=TEAL,point_size=5,render_points_as_spheres=True)
    colors={1:'#a94f60',2:'#ff9800',3:'#358d99'}
    for typ in [1,2,3]:
        links=b[types==typ];line=pv.PolyData();line.points=Xall;line.lines=np.c_[np.full(len(links),2),links]
        p.add_mesh(line,color=colors[typ],line_width=1,opacity=.075 if typ!=2 else .15)
    p.add_points(Xall[g['surface']],color='#ffd37b',point_size=7,render_points_as_spheres=True)
    p.camera_position=[(-3.8,-2.8,3.2),(-6,0,1.15),(0,0,1)];p.camera.parallel_projection=True;p.camera.parallel_scale=1.0
    axes(p,bounds=(-6.7,-5.3,-.65,.65,.7,1.5));im=capture(p,'lammps_native_scene.png')
    fig,aa=basefig('LAMMPS assembly: material points, anchors and typed bonds',
                   'All original cubic points and all bonds shown  |  Equal spacing = 50 µm  |  No dynamics',size=(14,9))
    aa[0,0].imshow(im);aa[0,0].axis('off')
    from matplotlib.lines import Line2D
    legend=[Line2D([0],[0],marker='o',ls='',color=RED,label=f"Clot: {n:,}"),Line2D([0],[0],marker='o',ls='',color=TEAL,label=f"Anchors: {len(Xall)-n:,} (3 layers)"),Line2D([0],[0],marker='o',ls='',color='#e5af46',label='Exposed area nodes'),Line2D([0],[0],color='#ff9800',label=f"Interface bonds: {info['interface_bonds']:,}")]
    fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.5,.039),ncol=4,frameon=False,fontsize=10)
    fig.text(.05,.023,f"Bulk bonds: {info['bulk_bonds']:,}   |   Wall bonds: {info['wall_bonds']:,}   |   Total: {info['bond_count']:,}   |   LAMMPS run 0: PASS; time advanced = 0",fontsize=10,color=INK)
    savefig(fig,'lammps_curved_interface_audit')


def movie():
    solid=pv.read(RESULTS/'paraview/clot_solid.vtp');vessel=pv.read(RESULTS/'paraview/vessel_surface.vtp')
    p=plotter((960,720));scene(p,vessel,solid,labels=False)
    p.add_text('Idealized Y-vessel + curved mural clot',position='upper_left',font_size=15,color=INK)
    p.add_text('STATIC GEOMETRY | Camera orbit only\nNo CFD recomputation or fracture dynamics',position='lower_left',font_size=12,color=RED)
    frames=[];cameras=[]
    with imageio.get_writer(RESULTS/'y_clot_static_camera_orbit.mp4',fps=12,codec='libx264',quality=8,macro_block_size=16) as writer:
        for i,angle in enumerate(np.linspace(0,2*np.pi,72,endpoint=False)):
            focal=np.array([-1.,0.,0.]);pos=focal+[38*np.cos(angle),38*np.sin(angle),30]
            p.camera_position=[pos,focal,[0,0,1]];p.camera.parallel_projection=True;p.camera.parallel_scale=14
            p.render();im=p.screenshot(return_img=True);writer.append_data(im)
            frames.append(Image.fromarray(im).resize((640,480)));cameras.append(dict(frame=i,position_mm=pos.tolist(),focal_mm=focal.tolist(),up=[0,0,1]))
    p.close();frames[0].save(RESULTS/'y_clot_static_camera_orbit.gif',save_all=True,append_images=frames[1:],duration=1000/12,loop=0)
    dump(RESULTS/'FIGURE_DATA/CAMERA_ORBIT.json',dict(frames=cameras,frame_count=72,fps=12,duration_s=6,meaning='Static geometry; only camera moves'))
    return dict(frames=72,fps=12,duration_s=6)


def export_interactive():
    pdir=RESULTS/'paraview'
    datasets={name:pv.read(pdir/file) for name,file in dict(vessel='vessel_surface.vtp',clot='clot_solid.vtp',clot_points='clot_points.vtu',anchors='anchors.vtu',bonds='all_bonds.vtp',masked_background_streamlines='background_streamlines_masked.vtp').items()}
    pv.MultiBlock(datasets).save(pdir/'Y_CLOT_SCENE.vtm')
    p=plotter();scene(p,datasets['vessel'],datasets['clot'],datasets['masked_background_streamlines'])
    p.add_text('UNOBSTRUCTED BACKGROUND CFD\nNOT CLOTTED FLOW FIELD',position='upper_left',font_size=14,color=RED)
    p.trame.export_html(RESULTS/'interactive_y_clot_scene.html');p.close()


def automatic_qa():
    flow=json.loads((OUT/'BACKGROUND_FLOW_OVERLAY_CHECK.json').read_text())
    animation=json.loads((RESULTS/'FIGURE_DATA/CAMERA_ORBIT.json').read_text())
    stats=[]
    for p in RESULTS.glob('*.png'):
        im=np.asarray(Image.open(p).convert('RGB'));nonwhite=np.any(im<240,axis=2)
        stats.append(dict(path=str(p),width=im.shape[1],height=im.shape[0],nonwhite_fraction=float(nonwhite.mean()),std=float(im.std()),sha256=digest(p)))
        assert nonwhite.mean()>.035 and im.std()>10
    assert len(stats)==4
    qa=dict(status='AUTOMATED_PASS_PENDING_VISUAL_REVIEW',figures=stats,movie=animation,
            coordinate_display_transform='1000 * saved SI coordinates; isotropic mm display only',same_geometry_for_all_views=True,
            smoothing='Surface normal shading only; no point displacement, no independent cosmetic clot contour',
            streamline_mask=flow,visual_review='PENDING',CLOT_PRESENT_CFD_SOLVED='NOT_RUN')
    dump(OUT/'VISUALIZATION_QA.json',qa)


def main():
    started=time.perf_counter();(RESULTS/'FIGURE_DATA').mkdir(exist_ok=True)
    flow=prepare_background();flowtime=time.perf_counter()-started;figures();figtime=time.perf_counter()-started-flowtime
    movie();export_interactive();automatic_qa()
    dump(OUT/'VISUALIZATION_PERFORMANCE.json',dict(wall_seconds=time.perf_counter()-started,flow_display_seconds=flowtime,figures_seconds=figtime,
         peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,CFD_steps=0,fracture_steps=0))
    print('VISUALS_READY',time.perf_counter()-started,flush=True)


if __name__=='__main__':main()
