"""Four audit figures; existing presentation renderer is only used read-only."""
from pathlib import Path
import sys,os,json,csv
os.environ.setdefault('LP_NUM_THREADS','2');os.environ.setdefault('VTK_SMP_MAX_THREADS','2')
sys.dont_write_bytecode=True
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
A=Path(__file__).resolve().parents[1];C=A.parent;V=C/'rotate_visualization'
sys.path.insert(0,str(V))
import render_surface_fields as renderer
shared=renderer.shared

def readcsv(name):return list(csv.DictReader((A/'data'/name).open()))
def render_audit(mesh,center=None,scale=None,raw=True,edges=False,labels=None,quality=False):
    p=pv.Plotter(off_screen=True,window_size=(1300,1000));p.set_background('white')
    p.add_mesh(mesh,scalars='minSICN' if quality else ('WSS_raw_Pa' if raw else 'WSS_display_Pa'),preference='cell' if raw else 'point',
        cmap='magma' if quality else 'viridis',clim=(0,1) if quality else (0,55),lighting=False,show_edges=edges,edge_color='#33404a',line_width=.5,show_scalar_bar=False)
    if center is None:center=np.array(mesh.center)
    direction=np.array([1.,-2.2,1.1]);direction/=np.linalg.norm(direction)
    p.camera.position=center+direction*250;p.camera.focal_point=center;p.camera.up=(0,0,1);p.enable_parallel_projection()
    p.camera.parallel_scale=scale or 65
    if labels:
        p.add_point_labels(np.array([r[1] for r in labels]),[r[0] for r in labels],point_color='red',point_size=9,font_size=18,text_color='black',shape_color='white',shape_opacity=.85,always_visible=True)
    im=p.screenshot();p.close();return im

def main():
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    # Figure 1: exact existing renderer, fields, fixed-Z camera and scalar scale.
    case=V/'input_data';wall=pv.read(case/'field_diagnostics/data/wall_wss_si.vtp');surface=pv.read(case/'field_diagnostics/data/pressure_surface_si.vtp')
    wall.points*=1e6;surface.points*=1e6;lines=pv.read(case/'streamlines/data/streamlines_si.vtp')
    framing=np.vstack([surface.points,lines.points*1e6]);center=shared.fitted_z_center(framing);orbit=shared.Turntable(framing,[0.,0.,1.],center)
    ports=[]
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        cap=pv.read(case/f'solver_mesh/mesh-surfaces/{role}.vtp');ports.append(('Inlet' if role=='INLET' else 'Outlet '+role[-2:],np.array(cap.center)*1e6))
    view=renderer.SurfaceView(wall,surface,orbit,ports,'wss',size=(3840,2160))
    image=view.frame(0,still=True);image.save(A/'figures/Figure_01_current_style_reproduced.png');view.close()
    original=Image.open(V/'results/surface_fields/figures/wss_overview_4k.png').convert('RGB');repro=image.convert('RGB')
    comparison=dict(original_size=original.size,reproduced_size=repro.size,original_path=str(V/'results/surface_fields/figures/wss_overview_4k.png'))
    if original.size==repro.size:
        delta=np.asarray(original,dtype=float)-np.asarray(repro,dtype=float);comparison['pixel_mean_abs_0_255']=float(abs(delta).mean());comparison['fraction_pixels_with_any_difference']=float(np.any(delta!=0,axis=2).mean())
    comparison['note']='Pixel equality is environment dependent; scalar equality verified separately. Meeting slide itself not available.'
    (A/'data/render_reproduction.json').write_text(json.dumps(comparison,indent=2)+'\n')
    reg=json.loads((A/'data/regions.json').read_text());critical=json.loads((A/'data/critical_locations.json').read_text());pa=np.load(A/'data/plot_arrays.npz')
    mesh=pv.read(A/'data/wall_audit.vtp');mesh.points*=1e6
    labels=[(name,np.array(reg[name]['cap_center_um'])) for name in ['INLET','O1','O2','O3']]
    labels+=[('J1 / node 3238',np.array(reg['first_junction_r5um']['center_um'])),('J2 / L: min WSS nearby',np.array(reg['junctions'][1]['xyz_um']))]
    globalim=render_audit(mesh,labels=labels)
    fig,ax=plt.subplots(figsize=(13,10));ax.imshow(globalim);ax.axis('off');ax.set_title('Region map | current H0 raw facet WSS | coordinates in micrometres')
    sm=plt.cm.ScalarMappable(norm=plt.Normalize(0,55),cmap='viridis');fig.colorbar(sm,ax=ax,fraction=.025,pad=.01,label='Raw WSS [Pa]')
    fig.savefig(A/'figures/Figure_02_region_map.png',dpi=180,bbox_inches='tight');plt.close(fig)
    # Same camera, same scale, same scalar range within each raw/display pair.
    rows=[('J1: first bifurcation',np.array(reg['first_junction_r5um']['center_um']),6.),
          ('O2: artificial extension',.5*(np.array(reg['O2']['cap_center_um'])+np.array(reg['O2']['real_cut_um'])),7.),
          ('L: global minimum',np.array(critical[0]['xyz_um']),4.)]
    md=np.load(A/'inputs/mesh_arrays.npz');tet=md['tetra'];xyz=md['points_m']*1e6
    volume=pv.UnstructuredGrid(np.c_[np.full(len(tet),4),tet],np.full(len(tet),10,np.uint8),xyz);volume.cell_data['minSICN']=md['min_sicn']
    tetcent=xyz[tet].mean(axis=1)
    fig,axs=plt.subplots(3,3,figsize=(21,17))
    for rr,(name,ctr,scale) in enumerate(rows):
        # Render only nearby facets to prevent another branch occluding the selected patch.
        select=np.linalg.norm(pa['centers_um']-ctr,axis=1)<scale*1.7
        local=mesh.extract_cells(select).extract_surface(algorithm='dataset_surface')
        for cc in range(2):
            im=render_audit(local,ctr,scale,raw=cc==0,edges=cc==0);axs[rr,cc].imshow(im);axs[rr,cc].axis('off')
            axs[rr,cc].set_title(name+(' | raw facets + actual mesh' if cc==0 else ' | area-weighted nodal display'))
        selected_volume=volume.extract_cells(np.linalg.norm(tetcent-ctr,axis=1)<scale*1.7)
        cut=selected_volume.slice(normal=(1.,-2.2,1.1),origin=ctr)
        axs[rr,2].imshow(render_audit(cut,ctr,scale,raw=True,edges=True,quality=True));axs[rr,2].axis('off');axs[rr,2].set_title(name+' | actual tetrahedral section / minSICN')
    fig.colorbar(sm,ax=axs[:,:2].ravel().tolist(),fraction=.015,pad=.01,label='WSS [Pa], shared linear scale 0-55')
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='magma'),ax=axs[:,2].tolist(),fraction=.03,pad=.01,label='minSICN [dimensionless], 0-1')
    fig.savefig(A/'figures/Figure_03_local_mesh_and_mapping.png',dpi=160,bbox_inches='tight');plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained')
    r=readcsv('poiseuille_validation.csv');xs=[float(a['radial_step_um']) for a in r];ys=[float(a['area_relative_L2_error_pct']) for a in r]
    axs[0,0].plot(xs,ys,'o-');axs[0,0].set(xlabel='Radial interval [um]',ylabel='WSS relative L2 error [%]',title='Analytic Poiseuille nodal field; NOT a CFD solve')
    for xx,yy in zip(xs,ys):axs[0,0].annotate(f'{yy:.2f}%',(xx,yy),xytext=(5,5),textcoords='offset points')
    r=readcsv('outlet2_path_profile.csv')
    for case,label in [('H0','Current H0'),('zero_outlets','Old: zero outlet pressures')]:
        rr=[a for a in r if a['case']==case];xx=[float(a['s_mid_um'])for a in rr];yy=[float(a['mean_Pa'])for a in rr]
        axs[0,1].plot(xx,yy,'o-',ms=3,label=label)
        if case=='H0':axs[0,1].fill_between(xx,[float(a['p05_Pa']) for a in rr],[float(a['p95_Pa'])for a in rr],alpha=.2)
    path=json.loads((A/'data/outlet2_path.json').read_text());axs[0,1].axvline(path['real_cut_arc_um'],c='gray',ls='--',label='Real cut / extension begins')
    axs[0,1].set(xlabel='Distance from J1 toward O2 cap [um]',ylabel='Raw WSS [Pa]',title='O2: area-weighted 1-um bins (H0 shade: P5-P95)',ylim=(0,55));axs[0,1].legend(fontsize=9)
    r=readcsv('region_summary.csv');names=['O1_extension','O2_extension','O3_extension'];xx=np.arange(3)
    for off,case,label in [(-.18,'zero_outlets','Old: zero outlets'),(.18,'H0','Current H0')]:
        yy=[float(next(a['mean_Pa'] for a in r if a['case']==case and a['region']==name)) for name in names]
        axs[1,0].bar(xx+off,yy,.36,label=label)
    axs[1,0].set(xticks=xx,xticklabels=['O1','O2','O3'],ylabel='Area-weighted mean WSS [Pa]',title='Same mesh / inlet / material; outlet-pressure pattern changed');axs[1,0].legend()
    r=readcsv('time_convergence.csv')[:-1];axs[1,1].semilogy([int(a['step'])for a in r],[float(a['wss_area_L2_relative_to71']) for a in r],'o-',label='WSS')
    axs[1,1].semilogy([int(a['step'])for a in r],[float(a['velocity_L2_relative_to71']) for a in r],'s-',label='Velocity')
    axs[1,1].set(xlabel='Saved solver step (dt = 1.50406e-7 s)',ylabel='Relative L2 difference from step 71',title='Time-iteration convergence; NOT spatial convergence');axs[1,1].legend()
    for ax in axs.flat:ax.grid(alpha=.2)
    fig.savefig(A/'figures/Figure_04_validation_and_sensitivity.png',dpi=200);plt.close(fig)
    print('Four figures generated inside wss_audit/figures.');print(comparison)

if __name__=='__main__':main()
