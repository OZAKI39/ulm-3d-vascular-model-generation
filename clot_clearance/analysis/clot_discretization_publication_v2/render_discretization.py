"""Publication discretization figure from the requested raw handoff point cloud.

No medical surface reconstruction and no new physical states.  The style follows
the supplied PZT figure; the viewing direction follows the supplied raw panel.
"""
import argparse
import hashlib
import json
from pathlib import Path
from pd_clot.portable_paths import relocate
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image
import pyvista as pv
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/clot_discretization_publication_v2'
PRIOR=ROOT/'results/y_clot_continuous_handoff_and_honest_render_v1'
STYLE=Path(str(relocate('/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/publication_figures/PZT_control_reference_v6')))
VTM=PRIOR/'paraview/raw_medical_step_4660.vtm'
FONT=Path(str(relocate('/mnt/c/Windows/Fonts/arial.ttf')))
COLORS=dict(clot='#b85c64',anchor='#8aa7b9',debris='#586d7b',bond='#718593',grid='#93a3b0',axis='#647887',text='#263746')

def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def load():
    blocks=pv.read(VTM);nodes=blocks['raw_all_nodes'];lines=blocks['raw_surviving_typed_bonds']
    native=np.load(ROOT/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/cases/C1_GAUSSIAN_A200/accepted_trajectory.npz')
    assert native['step'][-1]==4660 and np.max(np.abs(nodes.points-native['position_m'][-1]))<1e-17
    assert nodes.n_points==3488 and lines.n_lines==144827
    return nodes,lines

def camera_geometry(size):
    # Exact direction of the requested left raw panel; only framing is adapted.
    outward=np.array([.7,-2.2,-3.35]);outward/=np.linalg.norm(outward);forward=-outward
    right=np.cross(forward,[0.,0.,1.]);right/=np.linalg.norm(right);up=np.cross(right,forward)
    bounds=np.array([[-6.8,-5.2],[-.75,.85],[.55,1.5]])
    from itertools import product
    corners=np.array(list(product(*bounds)))
    focal=bounds.mean(axis=1)
    projected=corners-focal;u=projected@right;v=projected@up
    focal+=(u.max()+u.min())/2*right+(v.max()+v.min())/2*up
    scale=max(np.ptp(v)/1.85,np.ptp(u)/(1.85*size[0]/size[1]))
    return bounds,focal,right,up,outward,scale

def render(dpi,bonds=True):
    nodes,lines=load();X=nodes.points.astype(float)*1000
    # Local framing matches the requested raw-panel subject. Remote detached
    # material outside the local viewing box is recorded, never moved or merged.
    # More page height enlarges the unchanged 3D view at the same output width,
    # while keeping the full v1 local bounds and all of its detached nodes.
    figsize=np.array([180,168])/25.4;box=np.array([.025,.11,.94,.82]);size=np.rint(figsize*dpi*box[2:]).astype(int)
    bounds,focal,right,up,outward,scale=camera_geometry(size)
    in_box=np.all((X>=bounds[:,0])&(X<=bounds[:,1]),axis=1)
    selected=np.flatnonzero(in_box);ids=np.asarray(nodes['atom_id']);attached=np.asarray(nodes['attached']).astype(bool)
    segments=lines.lines.reshape(-1,3)[:,1:];keep=in_box[segments].all(axis=1);types=np.asarray(lines['bond_type'])
    p=pv.Plotter(off_screen=True,window_size=size.tolist(),lighting='three lights');p.set_background('white')
    p.render_window.SetMultiSamples(0);p.render_window.SetAlphaBitPlanes(1);p.enable_depth_peeling(number_of_peels=100,occlusion_ratio=0.)
    # The grid and label palette match the supplied PZT style. Rear planes are
    # adapted to the below-clot raw camera rather than rotating the specimen.
    xmin,xmax=bounds[0];ymin,ymax=bounds[1];zmin,zmax=bounds[2];grid=[]
    for x in np.linspace(xmin,xmax,5):grid.extend([[[x,ymin,zmax],[x,ymax,zmax]],[[x,ymax,zmin],[x,ymax,zmax]]])
    for y in np.linspace(ymin,ymax,5):grid.extend([[[xmin,y,zmax],[xmax,y,zmax]],[[xmin,y,zmin],[xmin,y,zmax]]])
    for z in np.linspace(zmin,zmax,4):grid.extend([[[xmin,ymin,z],[xmin,ymax,z]],[[xmin,ymax,z],[xmax,ymax,z]]])
    gg=pv.PolyData(np.asarray(grid).reshape(-1,3),lines=np.c_[np.full(len(grid),2),np.arange(2*len(grid)).reshape(-1,2)].ravel())
    p.add_mesh(gg,color=COLORS['grid'],lighting=False,line_width=max(.75,dpi/180*.62),opacity=.90)
    # Three silhouette edges keep numeric ticks outside the point cloud for
    # this below-specimen camera (unlike the above-vessel PZT camera).
    for edge in [[[xmin,ymax,zmin],[xmax,ymax,zmin]],[[xmin,ymin,zmin],[xmin,ymax,zmin]],[[xmax,ymax,zmin],[xmax,ymax,zmax]]]:p.add_mesh(pv.Line(*edge),color=COLORS['axis'],lighting=False,line_width=max(.8,dpi/180*.55))
    highlighted_bonds=0
    if bonds:
        for typ,opacity in [(1,.021),(2,.019),(3,.010)]:
            bb=segments[keep&(types==typ)];mesh=pv.PolyData(X,lines=np.c_[np.full(len(bb),2),bb].ravel())
            p.add_mesh(mesh,color='#9cabb5',opacity=opacity,line_width=max(.5,dpi/180*.42),lighting=False)
        # Keep the complete graph pale, and emphasize genuine short bonds on
        # the front layer. This gives legible links without an opaque grey mass.
        projected=np.c_[(X-focal)@right,(X-focal)@up];depth=(X-focal)@outward
        tree=cKDTree(projected[selected]);neighbors=tree.query_ball_point(projected[selected],.055)
        front=np.zeros(len(X),bool)
        front[selected]=np.array([depth[selected[js]].max()-depth[i] for i,js in zip(selected,neighbors)])<.045
        length=np.linalg.norm(X[segments[:,1]]-X[segments[:,0]],axis=1)
        emph=keep&front[segments].all(axis=1)&(types<=2)&(length<=.085)
        bb=segments[emph];highlighted_bonds=len(bb)
        mesh=pv.PolyData(X,lines=np.c_[np.full(len(bb),2),bb].ravel())
        p.add_mesh(mesh,color=COLORS['bond'],opacity=.42,line_width=max(.85,dpi/180*.70),lighting=False)
    groups=[('anchor',in_box&(ids>1927),1.45,.57),('clot',in_box&(ids<=1927)&attached,1.7,.94),('debris',in_box&(ids<=1927)&~attached,2.6,.98)]
    for name,mask,point_pt,alpha in groups:
        cloud=pv.PolyData(X[mask]);p.add_points(cloud,color=COLORS[name],point_size=dpi*point_pt/72,opacity=alpha,render_points_as_spheres=True,ambient=.55,diffuse=.45,specular=.10,specular_power=20,show_scalar_bar=False)
    p.camera_position=[(focal+outward*20).tolist(),focal.tolist(),up.tolist()];p.enable_parallel_projection();p.camera.parallel_scale=scale
    p.show(auto_close=False);p.render();image=p.screenshot(return_img=True);peeling=bool(p.renderer.GetLastRenderingUsedDepthPeeling());p.close()
    audit=dict(box=box.tolist(),figsize=figsize.tolist(),size=size.tolist(),camera_focal_mm=focal.tolist(),camera_right=right.tolist(),camera_up=up.tolist(),camera_outward=outward.tolist(),parallel_scale_mm=float(scale),local_bounds_mm=bounds.tolist(),source_points=3488,source_clot_nodes=1927,source_anchor_nodes=1561,rendered_clot_nodes=int((in_box&(ids<=1927)).sum()),rendered_attached_clot_nodes=int((in_box&(ids<=1927)&attached).sum()),rendered_detached_clot_nodes=int((in_box&(ids<=1927)&~attached).sum()),rendered_anchor_nodes=int((in_box&(ids>1927)).sum()),remote_detached_nodes_outside_local_view=int((~in_box&(ids<=1927)).sum()),rendered_bonds=int(keep.sum()) if bonds else 0,surviving_bond_counts={str(t):int((keep&(types==t)).sum()) for t in [1,2,3]},depth_peeling_used=peeling,coordinate_transformation='metres to millimetres only; no position, topology or shape editing',representation='Original discrete nodes and optional surviving bonds only; no reconstructed red surface',marker_semantics='Display symbols at real node centers; marker diameters are not particle diameters',source_state='C1 native step4660, physical time2.330ms, as in supplied raw panel')
    previous=json.loads((ROOT/'results/clot_discretization_publication_v1/PROVENANCE.json').read_text())
    for key in ['local_bounds_mm','rendered_clot_nodes','rendered_anchor_nodes','rendered_detached_clot_nodes','rendered_bonds','surviving_bond_counts','camera_outward']:
        assert audit[key]==previous[key],key
    audit['view_magnification_relative_to_v1_at_equal_page_width']=(168*.82/(132*.80))*previous['parallel_scale_mm']/scale
    audit['display_revision']=dict(clot_and_anchor_colors_unchanged=True,coordinate_bounds_and_members_unchanged=True,detached_marker_points=2.6,detached_marker_opacity=.98,background_bond_opacity_by_type={'1':.021,'2':.019,'3':.010},highlighted_front_bonds=highlighted_bonds,highlight_rule='Actual surviving bulk/interface links, both ends within45um of the front-most local layer in a55um screen-plane radius and deformed length<=85um. Entire graph retained in pale background.',highlight_bond_opacity=.42,highlight_bond_linewidth_at180dpi=.70,grid_opacity=.90,grid_linewidth_at180dpi=.62)
    return image,audit

def project(point,a):
    delta=np.asarray(point)-a['camera_focal_mm'];aspect=a['size'][0]/a['size'][1];u=.5+delta@a['camera_right']/(2*a['parallel_scale_mm']*aspect);v=.5+delta@a['camera_up']/(2*a['parallel_scale_mm']);x,y,w,h=a['box'];return np.array([x+w*u,y+h*v])

def compose(image,a,dpi):
    if FONT.exists():font_manager.fontManager.addfont(str(FONT))
    plt.rcParams.update({'font.family':'Arial' if FONT.exists() else 'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','text.color':'#20262d'})
    fig=plt.figure(figsize=a['figsize'],dpi=dpi,facecolor='white');ax=fig.add_axes(a['box']);ax.imshow(image,interpolation='none');ax.axis('off')
    fig.text(.045,.965,'Clot discretization',fontsize=12,va='top');fig.text(.955,.965,'Δx = 50 µm',fontsize=10,ha='right',va='top',color='#425665')
    xmin,xmax=a['local_bounds_mm'][0];ymin,ymax=a['local_bounds_mm'][1];zmin,zmax=a['local_bounds_mm'][2]
    entries=[([x,ymax,zmin],f'{x:.1f}',(0,-11),'center') for x in [-6.8,-6.0,-5.2]]
    entries.extend(([xmin,y,zmin],f'{y:g}',(-10,0),'right') for y in [-.5,0,.5])
    entries.extend(([xmax,ymax,z],f'{z:g}',(9,0),'left') for z in [.6,1.,1.4])
    for point,label,offset,ha in entries:
        ax.annotate(label,xy=project(point,a),xycoords=fig.transFigure,xytext=offset,textcoords='offset points',fontsize=8,color=COLORS['text'],ha=ha,va='center',annotation_clip=False)
    for point,label,offset in [([-6.,ymax,zmin],'X (mm)',(0,-25)),([xmin,.05,zmin],'Y (mm)',(-35,0)),([xmax,ymax,1.],'Z (mm)',(36,0))]:
        ax.annotate(label,xy=project(point,a),xycoords=fig.transFigure,xytext=offset,textcoords='offset points',fontsize=9,color=COLORS['text'],ha='center',va='center',annotation_clip=False)
    handles=[Line2D([],[],marker='o',linestyle='none',markersize=4,color=COLORS['clot'],label='Clot nodes'),Line2D([],[],marker='o',linestyle='none',markersize=4,color=COLORS['anchor'],label='Fixed anchors'),Line2D([],[],marker='o',linestyle='none',markersize=4.5,color=COLORS['debris'],label='Detached nodes')]
    if a['rendered_bonds']:handles.append(Line2D([],[],color=COLORS['bond'],lw=.85,label='Surviving bonds'))
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.035),ncol=len(handles),frameon=False,fontsize=8,handlelength=1.1,columnspacing=1.5)
    fig.canvas.draw();w,h=fig.canvas.get_width_height();renderer=fig.canvas.get_renderer()
    for art in fig.findobj(matplotlib.text.Text):
        if art.get_visible() and art.get_text():
            b=art.get_window_extent(renderer);assert b.x0>=0 and b.y0>=0 and b.x1<=w+1 and b.y1<=h+1,(art.get_text(),b)
    return fig

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--preview',action='store_true');args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True);dpi=180 if args.preview else 600
    sources=[STYLE/'Fig1a_PZT_control.png',STYLE/'paper_pzt_reference.py',PRIOR/'figures/raw_vs_honest_medical_handoff.png',VTM,*sorted(VTM.with_suffix('').glob('*.vtp'))];before={str(p):digest(p) for p in sources}
    image,a=render(dpi);fig=compose(image,a,dpi);stem='Fig_clot_discretization';suffix='_preview' if args.preview else '';fig.savefig(OUT/(stem+suffix+'.png'),dpi=dpi)
    if not args.preview:
        fig.savefig(OUT/(stem+'.pdf'),dpi=dpi,metadata={'Title':'Clot discretization from native handoff nodes'});fig.savefig(OUT/(stem+'.svg'),dpi=dpi);im=Image.open(OUT/(stem+'.png')).convert('RGB');im.save(OUT/(stem+'.tif'),compression='tiff_lzw',dpi=(600,600));im.thumbnail((1600,1600));im.save(OUT/(stem+'_preview.png'))
    plt.close(fig)
    assert before=={p:digest(p) for p in before}
    a.update(status='PASS',protected_source_sha256=before,style_reference=str(STYLE/'Fig1a_PZT_control.png'),shape_reference=str(PRIOR/'figures/raw_vs_honest_medical_handoff.png'),shape_reference_panel='LEFT RAW ONLY',font='Arial' if FONT.exists() else 'DejaVu Sans',colors=COLORS,dpi=dpi,new_solver_runs=0,export='Rendered scientific scene; editable vector text in PDF/SVG')
    (OUT/('PREVIEW_PROVENANCE.json' if args.preview else 'PROVENANCE.json')).write_text(json.dumps(a,indent=2))
    if not args.preview:
        (OUT/'CAPTION.txt').write_text(f'''Clot discretization. Original discrete clot and anchor nodes and their surviving connections at the2.330ms native handoff state, corresponding to the LEFT raw panel specified by the user. White background, grey coordinate cage, Arial typography and muted red/blue palette follow the supplied PZT reference. Camera direction follows the raw panel. No reconstructed medical surface is used. Coordinates are original solver coordinates in millimetres; Δx=50µm is the reference lattice spacing, not a physical particle diameter or every deformed edge length.\n\nThe local view contains {a['rendered_clot_nodes']} clot nodes ({a['rendered_attached_clot_nodes']} attached and {a['rendered_detached_clot_nodes']} detached), {a['rendered_anchor_nodes']} fixed anchors and {a['rendered_bonds']} surviving bonds. {a['remote_detached_nodes_outside_local_view']} remote detached nodes lie outside the local viewing box. They are retained in source data; this is a local discretization view, not a full-domain fragment census. Marker sizes represent graphical symbols, not physical radii. No solver data, mass or coordinates are changed.\n\nRevision2 enlarges the scene by {a['view_magnification_relative_to_v1_at_equal_page_width']:.3f} at equal page width through page layout, retaining the same camera direction and local coordinate bounds. Detached markers are darker and larger; coordinate grid lines are darker. All surviving links remain represented; {a['display_revision']['highlighted_front_bonds']} genuine short links in the front layer are emphasized for legibility. This view-dependent line emphasis is purely graphical and is not a distinct mechanical bond class.\n''')
        (OUT/'README.txt').write_text('''血栓离散化图 v2：保持原白底、Arial字体、红蓝配色与观察方向；同宽显示时主体约放大30%。游离节点加大、加深，网格线加深。全部原存活连接保留淡色显示，并强调前侧真实短连接，避免密集线条叠成灰块。节点位置、取景范围与拓扑未改变。\n\nPNG/TIF：600dpi。PDF/SVG：三维场景为高清嵌入图，文字可编辑。preview.png用于快速浏览。局部范围和计数见PROVENANCE.json，图外远端碎屑保留在原数据中。旧v1保留；没有新增求解。\n\n复现：从clot_clearance根目录运行analysis/clot_discretization_publication_v2/render_discretization.py。\n''')
        files=[p for p in OUT.iterdir() if p.suffix in ['.png','.pdf','.svg','.tif']];(OUT/'DELIVERY_SHA256.json').write_text(json.dumps({p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in files},indent=2))
    print(json.dumps({k:a[k] for k in ['status','rendered_clot_nodes','rendered_anchor_nodes','remote_detached_nodes_outside_local_view','rendered_bonds','depth_peeling_used']}),flush=True)

if __name__=='__main__':main()
