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
OUT=ROOT/'results/clot_discretization_three_views_v1'
PRIOR=ROOT/'results/y_clot_continuous_handoff_and_honest_render_v1'
STYLE=Path(str(relocate('/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles/results/stage4/publication_figures/PZT_control_reference_v6')))
VTM=PRIOR/'paraview/raw_medical_step_4660.vtm'
FONT=Path(str(relocate('/mnt/c/Windows/Fonts/arial.ttf')))
COLORS=dict(clot='#b85c64',anchor='#8aa7b9',debris='#586d7b',bond='#718593',grid='#93a3b0',axis='#647887',text='#263746')
VIEWS=[dict(id=1,name='Original oblique',outward=[.7,-2.2,-3.35]),dict(id=2,name='Side oblique',azimuth=-145.,elevation=-18.),dict(id=3,name='Opposite oblique',azimuth=50.,elevation=-35.)]

def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def load():
    blocks=pv.read(VTM);nodes=blocks['raw_all_nodes'];lines=blocks['raw_surviving_typed_bonds']
    native=np.load(ROOT/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/cases/C1_GAUSSIAN_A200/accepted_trajectory.npz')
    assert native['step'][-1]==4660 and np.max(np.abs(nodes.points-native['position_m'][-1]))<1e-17
    assert nodes.n_points==3488 and lines.n_lines==144827
    return nodes,lines

def camera_geometry(size,view):
    # Exact direction of the requested left raw panel; only framing is adapted.
    if 'outward' in view:outward=np.array(view['outward'],float)
    else:
        az,el=np.deg2rad([view['azimuth'],view['elevation']]);outward=np.array([np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),np.sin(el)])
    outward/=np.linalg.norm(outward);forward=-outward
    right=np.cross(forward,[0.,0.,1.]);right/=np.linalg.norm(right);up=np.cross(right,forward)
    bounds=np.array([[-6.8,-5.2],[-.75,.85],[.55,1.5]])
    from itertools import product
    corners=np.array(list(product(*bounds)))
    focal=bounds.mean(axis=1)
    projected=corners-focal;u=projected@right;v=projected@up
    focal+=(u.max()+u.min())/2*right+(v.max()+v.min())/2*up
    # A common scale makes the three views directly comparable with revised v2.
    scale=json.loads((ROOT/'results/clot_discretization_publication_v2/PROVENANCE.json').read_text())['parallel_scale_mm']
    return bounds,focal,right,up,outward,scale

def silhouette_axes(bounds,right,up):
    from itertools import product
    corners=np.array(list(product(*bounds)));P=np.c_[corners@right,corners@up];center=P.mean(0);result=[]
    for axis in range(3):
        other=[k for k in range(3) if k!=axis];candidates=[]
        for ab in product(*[bounds[k] for k in other]):
            endpoints=np.zeros((2,3));endpoints[:,axis]=bounds[axis]
            for k,value in zip(other,ab):endpoints[:,k]=value
            Q=np.c_[endpoints@right,endpoints@up];d=Q[1]-Q[0];side=d[0]*(P[:,1]-Q[0,1])-d[1]*(P[:,0]-Q[0,0])
            if not (np.min(side)>=-1e-10 or np.max(side)<=1e-10):continue
            normal=np.array([-d[1],d[0]]);normal/=np.linalg.norm(normal)
            if normal@(Q.mean(0)-center)<0:normal=-normal
            score=Q.mean(0)[0] if axis==2 else -Q.mean(0)[1]
            candidates.append((score,endpoints,normal))
        _,edge,normal=max(candidates,key=lambda r:r[0]);result.append(dict(axis=axis,edge=edge.tolist(),outward_screen=normal.tolist()))
    return result

def render(dpi,view,bonds=True):
    nodes,lines=load();X=nodes.points.astype(float)*1000
    # Local framing matches the requested raw-panel subject. Remote detached
    # material outside the local viewing box is recorded, never moved or merged.
    # More page height enlarges the unchanged 3D view at the same output width,
    # while keeping the full v1 local bounds and all of its detached nodes.
    figsize=np.array([180,168])/25.4;box=np.array([.025,.11,.94,.82]);size=np.rint(figsize*dpi*box[2:]).astype(int)
    bounds,focal,right,up,outward,scale=camera_geometry(size,view)
    in_box=np.all((X>=bounds[:,0])&(X<=bounds[:,1]),axis=1)
    selected=np.flatnonzero(in_box);ids=np.asarray(nodes['atom_id']);attached=np.asarray(nodes['attached']).astype(bool)
    segments=lines.lines.reshape(-1,3)[:,1:];keep=in_box[segments].all(axis=1);types=np.asarray(lines['bond_type'])
    p=pv.Plotter(off_screen=True,window_size=size.tolist(),lighting='three lights');p.set_background('white')
    p.render_window.SetMultiSamples(0);p.render_window.SetAlphaBitPlanes(1);p.enable_depth_peeling(number_of_peels=100,occlusion_ratio=0.)
    # The grid and label palette match the supplied PZT style. Rear planes are
    # adapted to the below-clot raw camera rather than rotating the specimen.
    grid=[]
    for plane in range(3):
        fixed=bounds[plane,0 if outward[plane]>0 else 1];other=[k for k in range(3) if k!=plane]
        for varying in other:
            constant=next(k for k in other if k!=varying)
            for value in np.linspace(*bounds[constant],4 if constant==2 else 5):
                edge=np.zeros((2,3));edge[:,plane]=fixed;edge[:,constant]=value;edge[:,varying]=bounds[varying];grid.append(edge)
    gg=pv.PolyData(np.asarray(grid).reshape(-1,3),lines=np.c_[np.full(len(grid),2),np.arange(2*len(grid)).reshape(-1,2)].ravel())
    p.add_mesh(gg,color=COLORS['grid'],lighting=False,line_width=max(.75,dpi/180*.62),opacity=.90)
    # Three silhouette edges keep numeric ticks outside the point cloud for
    # this below-specimen camera (unlike the above-vessel PZT camera).
    axes=silhouette_axes(bounds,right,up)
    for axis in axes:p.add_mesh(pv.Line(*axis['edge']),color=COLORS['axis'],lighting=False,line_width=max(.8,dpi/180*.55))
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
    for key in ['local_bounds_mm','rendered_clot_nodes','rendered_anchor_nodes','rendered_detached_clot_nodes','rendered_bonds','surviving_bond_counts']:
        assert audit[key]==previous[key],key
    audit['view_magnification_relative_to_v1_at_equal_page_width']=(168*.82/(132*.80))*previous['parallel_scale_mm']/scale
    audit['display_revision']=dict(clot_and_anchor_colors_unchanged=True,coordinate_bounds_and_members_unchanged=True,detached_marker_points=2.6,detached_marker_opacity=.98,background_bond_opacity_by_type={'1':.021,'2':.019,'3':.010},highlighted_front_bonds=highlighted_bonds,highlight_rule='Actual surviving bulk/interface links, both ends within45um of the front-most local layer in a55um screen-plane radius and deformed length<=85um. Entire graph retained in pale background.',highlight_bond_opacity=.42,highlight_bond_linewidth_at180dpi=.70,grid_opacity=.90,grid_linewidth_at180dpi=.62)
    audit.update(view_id=view['id'],view_name=view['name'],axis_edges=axes,camera_azimuth_deg=float(np.rad2deg(np.arctan2(outward[1],outward[0]))),camera_elevation_deg=float(np.rad2deg(np.arcsin(outward[2]))),same_scale_as_other_views=True)
    return image,audit

def project(point,a):
    delta=np.asarray(point)-a['camera_focal_mm'];aspect=a['size'][0]/a['size'][1];u=.5+delta@a['camera_right']/(2*a['parallel_scale_mm']*aspect);v=.5+delta@a['camera_up']/(2*a['parallel_scale_mm']);x,y,w,h=a['box'];return np.array([x+w*u,y+h*v])

def compose(image,a,dpi):
    if FONT.exists():font_manager.fontManager.addfont(str(FONT))
    plt.rcParams.update({'font.family':'Arial' if FONT.exists() else 'DejaVu Sans','font.size':10,'pdf.fonttype':42,'svg.fonttype':'none','text.color':'#20262d'})
    fig=plt.figure(figsize=a['figsize'],dpi=dpi,facecolor='white');ax=fig.add_axes(a['box']);ax.imshow(image,interpolation='none');ax.axis('off')
    fig.text(.045,.965,f'Clot discretization · View {a["view_id"]}',fontsize=12,va='top');fig.text(.955,.965,'Δx = 50 µm',fontsize=10,ha='right',va='top',color='#425665')
    ticks=[[-6.8,-6.0,-5.2],[-.5,0,.5],[.6,1.,1.4]]
    for item in a['axis_edges']:
        k=item['axis'];normal=np.asarray(item['outward_screen']);base=np.mean(item['edge'],axis=0)
        for value in ticks[k]:
            point=base.copy();point[k]=value;label=f'{value:.1f}' if k==0 else f'{value:g}'
            ax.annotate(label,xy=project(point,a),xycoords=fig.transFigure,xytext=normal*12,textcoords='offset points',fontsize=8,color=COLORS['text'],ha='center',va='center',annotation_clip=False)
        ax.annotate('XYZ'[k]+' (mm)',xy=project(base,a),xycoords=fig.transFigure,xytext=normal*(36 if k==2 else 30),textcoords='offset points',fontsize=9,color=COLORS['text'],ha='center',va='center',annotation_clip=False)
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
    sources=[STYLE/'Fig1a_PZT_control.png',PRIOR/'figures/raw_vs_honest_medical_handoff.png',VTM,*sorted(VTM.with_suffix('').glob('*.vtp')),ROOT/'results/clot_discretization_publication_v2/PROVENANCE.json'];before={str(p):digest(p) for p in sources};audits=[];preview_images=[]
    for view in VIEWS:
        image,a=render(dpi,view);fig=compose(image,a,dpi);stem=f'Fig_clot_discretization_view_{view["id"]}';suffix='_preview' if args.preview else '';png=OUT/(stem+suffix+'.png');fig.savefig(png,dpi=dpi)
        if not args.preview:
            fig.savefig(OUT/(stem+'.pdf'),dpi=dpi,metadata={'Title':f'Clot discretization: view {view["id"]}'});fig.savefig(OUT/(stem+'.svg'),dpi=dpi)
            im=Image.open(png).convert('RGB');im.save(OUT/(stem+'.tif'),compression='tiff_lzw',dpi=(600,600));im.thumbnail((1600,1600));im.save(OUT/(stem+'_preview.png'))
        plt.close(fig)
        a.update(status='PASS',protected_source_sha256=before,font='Arial' if FONT.exists() else 'DejaVu Sans',colors=COLORS,dpi=dpi,new_solver_runs=0,export='Rendered scientific scene; editable vector labels in individual PDF/SVG')
        audits.append(a);im=Image.open(OUT/(stem+'_preview.png')).convert('RGB');im=im.resize((960,896),Image.Resampling.LANCZOS);preview_images.append(im)
        print('VIEW_COMPLETE',view['id'],a['camera_azimuth_deg'],a['camera_elevation_deg'],a['display_revision']['highlighted_front_bonds'],flush=True)
    assert before=={p:digest(p) for p in before}
    assert len(set(a['parallel_scale_mm'] for a in audits))==1
    sheet=Image.new('RGB',(2880,896),'white')
    for j,im in enumerate(preview_images):sheet.paste(im,(j*960,0))
    sheet.save(OUT/'Fig_clot_discretization_three_views_preview.png')
    if not args.preview:
        sheet=Image.new('RGB',(6378,1984),'white')
        for j,view in enumerate(VIEWS):
            im=Image.open(OUT/f'Fig_clot_discretization_view_{view["id"]}.png').convert('RGB').resize((2126,1984),Image.Resampling.LANCZOS);sheet.paste(im,(j*2126,0))
        sheet.save(OUT/'Fig_clot_discretization_three_views.png',dpi=(300,300))
    name='PREVIEW_PROVENANCE.json' if args.preview else 'PROVENANCE.json'
    (OUT/name).write_text(json.dumps(dict(status='PASS',views=audits,same_physical_state=True,same_local_bounds=True,same_nodes_and_surviving_bonds=True,same_world_scale=True,original_positions_unchanged=True,source_sha256=before,new_solver_runs=0),indent=2))
    if not args.preview:
        (OUT/'README.txt').write_text("三视角血栓离散化图，沿用调整后的v2样式、放大比例、深色网格、较大游离节点和前侧连接强调。三个视角使用完全相同的局部节点/存活键、物理坐标、取景边界及世界长度比例；只改变相机方向。\n\nView1：原斜视角。View2：侧向斜视角。View3：另一侧斜视角。各图PNG/TIF为600dpi，PDF/SVG保留可编辑文字；three_views为并排对照。前侧强调仍只选真实存活连接，随视角重新判定，不是新力学类别。\n\n每图包含1718个clot节点（1702附着、16游离）与1561个锚节点，144811条局部存活连接。209个远端游离节点在相同局部范围之外，原数据保留。所有源文件与旧v1/v2图只读，无新增求解。\n")
        files=[p for p in OUT.iterdir() if p.suffix in ['.png','.pdf','.svg','.tif']];(OUT/'DELIVERY_SHA256.json').write_text(json.dumps({p.name:dict(bytes=p.stat().st_size,sha256=digest(p)) for p in files},indent=2))

if __name__=='__main__':main()
