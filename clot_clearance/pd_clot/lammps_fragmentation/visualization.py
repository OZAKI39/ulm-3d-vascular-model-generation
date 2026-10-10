"""Actual particles and every surviving bond; no synthetic fragmentation."""
from pathlib import Path
import json,time,colorsys
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
import imageio.v2 as imageio
from .config import OUTPUT,RESULTS,PATTERNS,dump

SHORT=dict(LOCAL_GAUSSIAN_SHEAR='Local Gaussian shear',OPPOSING_SHEAR='Opposing shear',VORTEX_TANGENTIAL='Tangential vortex')

def _label(plot,text,**kwargs):
    actor=plot.add_text(text,**kwargs)
    prop=actor.GetTextProperty();prop.SetBackgroundColor(1,1,1);prop.SetBackgroundOpacity(.96)
    return actor

def colors(ids,base):
    palette={int(i):np.array(colorsys.hsv_to_rgb((int(i)*.61803398875)%1,.67,.83))*255 for i in np.unique(ids)}
    rgb=np.array([palette[int(i)] for i in ids],np.uint8);rgb[base]=[80,80,85];return rgb

def scene(case_path,index=-1,load=False,window=(1000,800),bounds=None):
    p=Path(case_path);r=json.loads((p/'request.json').read_text());z=np.load(p/'accepted_trajectory.npz');g=np.load(OUTPUT/'HEMISPHERE_POINT_CLOUD.npz');data=np.load(p/'broken_bonds.npz')
    X=(g['X'] if load else z['position_m'][index])*1e3;ids=np.ones(len(X),int) if load else z['component_id'][index];base=z['base'];bonds=data['initial_bonds'];alive=np.ones(len(bonds),bool)
    if not load:alive[data['initial_bond_indices'][data['events'][:,1]<=z['step'][index]]]=False
    plot=pv.Plotter(off_screen=True,window_size=window);plot.set_background('white');plot.enable_anti_aliasing('ssaa')
    pts=pv.PolyData(X);pts['rgb']=colors(ids,base);dx=r['dx_m']*1e3
    cube=pv.Cube(x_length=.90*dx,y_length=.90*dx,z_length=.90*dx);glyph=pts.glyph(orient=False,scale=False,geom=cube)
    
    if load:plot.add_mesh(glyph,scalars='rgb',rgb=True,show_scalar_bar=False,opacity=.22,ambient=.45,diffuse=.65)
    else:plot.add_mesh(pts,style='points',scalars='rgb',rgb=True,show_scalar_bar=False,point_size=4)
    e=bonds[alive]
    if not load and len(e):
        lines=pv.PolyData(X);lines.lines=np.c_[np.full(len(e),2),e].ravel();plot.add_mesh(lines,color='#334155',opacity=.13,line_width=.6)
    R=r['radius_m']*1e3;plot.add_mesh(pv.Plane(center=(0,0,0),direction=(0,0,1),i_size=2.3*R,j_size=2.3*R),color='#d5d7db',opacity=.32)
    if load:
        table=np.load(p/'surface_forces.npz');ii=table['surface_node_id']-1;F=table['force_N'];take=np.arange(0,len(ii),7);length=np.linalg.norm(F[take],axis=1);take=take[length>1e-18]
        common=max(np.linalg.norm(np.load(OUTPUT/'cases'/name/'surface_forces.npz')['force_N'],axis=1).max() for name in PATTERNS)
        plot.add_arrows(X[ii[take]],F[take]/common*.6*R,mag=1,color='#bd3d22')
        title=SHORT[r['pattern_id']]+'\nReference tangential load; equal total |F|'
    else:
        nc=len(np.unique(ids));vmax=np.bincount(ids,weights=z['volume_m3']).max()/z['volume_m3'].sum();title=f'{SHORT[r["pattern_id"]]}\nt = {z["step"][index]*r["dt_s"]*1e3:.3f} ms | FI = {1-vmax:.3f} | N = {nc}'
    _label(plot,title,font_size=12,color='#182532',position='upper_left');_label(plot,'Component colors; grey base\nPositions unscaled; axes in mm',font_size=10,color='#334155',position='lower_left')
    plot.add_axes();plot.show_bounds(grid='back',location='outer',all_edges=False,xtitle='x (mm)',ytitle='y (mm)',ztitle='z (mm)',font_size=9)
    plot.camera_position='iso'
    if bounds is None:bounds=(X[:,0].min()-.15*R,X[:,0].max()+.15*R,X[:,1].min()-.15*R,X[:,1].max()+.15*R,min(0.,X[:,2].min())-.15*R,X[:,2].max()+.15*R)
    plot.reset_camera(bounds=bounds);plot.camera.zoom(.85);plot.show(auto_close=False);arr=plot.screenshot(return_img=True);plot.close();return arr

def render_case(case):
    return scene(case['path'])

def render_all():
    start=time.perf_counter();F=RESULTS/'figures';A=RESULTS/'animations';P=RESULTS/'previews'
    for d in [F,A,P]:d.mkdir(parents=True,exist_ok=True)
    results=[json.loads((OUTPUT/'cases'/name/'result.json').read_text()) for name in PATTERNS]
    images=[]
    for name in PATTERNS:
        load=scene(OUTPUT/'cases'/name,load=True);final=scene(OUTPUT/'cases'/name)
        Image.fromarray(load).save(P/(name+'_load.png'));Image.fromarray(final).save(P/(name+'_final.png'));images.append((load,final))
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig,axs=plt.subplots(3,2,figsize=(12,13),layout='constrained')
    for i,(r,pair) in enumerate(zip(results,images)):
        for ax,img in zip(axs[i],pair):ax.imshow(img);ax.axis('off')
        axs[i,0].set_title(f'{SHORT[r["pattern_id"]]} | peak traction {r["peak_traction"]:.0f} Pa',fontsize=11)
        axs[i,1].set_title(f'Final: FI={r["FI"]:.3f}; N={r["fragment_count"]}; Vmax/V0={r["largest_fragment_fraction"]:.3f}',fontsize=11)
    fig.suptitle(f'Equal-budget streaming patterns | clot-like LPS | exploratory only\nTotal tangential force budget = {results[0]["load_budget"]*1e3:.3f} mN; component colors from surviving bonds\nEndpoint ranking uncertain: FI near anchored-base limit',fontsize=14)
    fig.savefig(F/'streaming_pattern_fragmentation_comparison.png',dpi=300);fig.savefig(F/'streaming_pattern_fragmentation_comparison.pdf',dpi=300);fig.savefig(P/'comparison_preview.png',dpi=90);plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(12,4.2),layout='constrained')
    for name in PATTERNS:
        data=np.genfromtxt(OUTPUT/'cases'/name/'fragment_metrics.csv',delimiter=',',names=True)
        for ax,key in zip(axs,['broken_bond_fraction','FI','N_fragments']):ax.plot(data['time_s']*1e3,data[key],label=SHORT[name],lw=2)
    for ax,label in zip(axs,['Broken bond fraction','FI = 1 − Vmax/V0','Raw connected components']):ax.set_xlabel('Simulation time (ms)');ax.set_ylabel(label);ax.grid(alpha=.2)
    axs[0].legend(fontsize=8);fig.suptitle('Actual bond deletion and connected-component evolution | no debris removed')
    fig.savefig(F/'fragmentation_metrics_vs_time.png',dpi=300);fig.savefig(F/'fragmentation_metrics_vs_time.pdf',dpi=300);fig.savefig(P/'metrics_preview.png',dpi=120);plt.close(fig)
    strongest=max(results,key=lambda r:r['FI']);p=OUTPUT/'cases'/strongest['pattern_id'];r=json.loads((p/'request.json').read_text());z=np.load(p/'accepted_trajectory.npz');data=np.load(p/'broken_bonds.npz');Xall=z['position_m']*1e3;lo=Xall.min((0,1));hi=Xall.max((0,1));pad=np.maximum((hi-lo)*.07,r['dx_m']*1e3*2);bounds=(lo[0]-pad[0],hi[0]+pad[0],lo[1]-pad[1],hi[1]+pad[1],min(0.,lo[2])-pad[2],hi[2]+pad[2])
    # Two simultaneous views retain all particles: local formation and complete cloud.
    plot=pv.Plotter(off_screen=True,shape=(1,2),window_size=(1400,800));R=r['radius_m']*1e3;dx=r['dx_m']*1e3
    local_bounds=(-2*R,3*R,-2*R,2*R,-.3*R,3*R)
    for col in [0,1]:
        plot.subplot(0,col);plot.set_background('white');plot.add_mesh(pv.Plane(center=(0,0,0),direction=(0,0,1),i_size=2.3*R,j_size=2.3*R),color='#b8bec7',opacity=.4);plot.add_axes()
        _label(plot,'Original clot neighborhood (fixed view)\nUnscaled positions; all surviving bonds' if col==0 else 'Complete particle cloud (all debris retained)\nPositions unscaled; points enlarged for visibility',position='lower_left',font_size=10,color='#334155')
        plot.camera_position='iso';plot.reset_camera(bounds=local_bounds if col==0 else bounds);plot.camera.zoom(.9)
    plot.show(auto_close=False);cube=pv.Cube(x_length=.90*dx,y_length=.90*dx,z_length=.90*dx);gif=[];manifest=[];bonds=data['initial_bonds']
    with imageio.get_writer(A/'hemisphere_fragmentation.mp4',fps=12,codec='libx264',quality=8,macro_block_size=2) as writer:
        for i,step in enumerate(z['step']):
            X=Xall[i];ids=z['component_id'][i];pts=pv.PolyData(X);pts['rgb']=colors(ids,z['base']);glyph=pts.glyph(orient=False,scale=False,geom=cube)
            alive=np.ones(len(bonds),bool);alive[data['initial_bond_indices'][data['events'][:,1]<=step]]=False;e=bonds[alive];line=pv.PolyData(X)
            if len(e):line.lines=np.c_[np.full(len(e),2),e].ravel()
            counts=np.unique(ids,return_counts=True)[1];FI=1-counts.max()/len(ids)
            for col in [0,1]:
                plot.subplot(0,col)
                if col==0:plot.add_mesh(glyph,scalars='rgb',rgb=True,show_scalar_bar=False,name='particles',ambient=.45,diffuse=.65,reset_camera=False)
                else:plot.add_mesh(pts,style='points',scalars='rgb',rgb=True,show_scalar_bar=False,name='particles',point_size=3,reset_camera=False)
                if len(e):plot.add_mesh(line,name='surviving_bonds',color='#334155',opacity=.13,line_width=.6,reset_camera=False)
                else:plot.remove_actor('surviving_bonds')
                title=(SHORT[r['pattern_id']]+' | exploratory only\n' if col==0 else 'Color = actual connected component\n')+f't = {step*r["dt_s"]*1e3:.3f} ms | FI = {FI:.3f} | N = {len(counts)}'
                _label(plot,title,name='state',position='upper_left',font_size=11,color='#182532')
            arr=plot.screenshot(return_img=True);writer.append_data(arr)
            if i%3==0 or i==len(z['step'])-1:gif.append(Image.fromarray(arr).resize((875,500)))
            if i in [0,len(z['step'])//4,len(z['step'])//2,len(z['step'])-1]:Image.fromarray(arr).save(P/f'animation_frame_{i:04d}.png')
            manifest.append(dict(frame=i,LAMMPS_step=int(step),simulation_time_s=float(step*r['dt_s']),synthetic_state=False))
            if i%50==0:print('RENDER',i,'/',len(z['step']),flush=True)
    plot.close();gif[0].save(A/'hemisphere_fragmentation_preview.gif',save_all=True,append_images=gif[1:],duration=250,loop=0)
    dump(A/'FRAME_MANIFEST.json',dict(pattern_id=r['pattern_id'],fps=12,physical_fracture_speed_claim=False,all_surviving_bonds_drawn=True,no_removed_debris=True,no_interpolated_solver_states=True,frames=manifest))
    qa=dict(status='AWAITING_VISUAL_INSPECTION',render_wall_seconds=time.perf_counter()-start,animation_pattern=r['pattern_id'],frames=len(z['step']),figures=2,all_surviving_bonds_drawn=True,actual_1x=True)
    dump(OUTPUT/'VISUALIZATION_QA.json',qa);return qa
