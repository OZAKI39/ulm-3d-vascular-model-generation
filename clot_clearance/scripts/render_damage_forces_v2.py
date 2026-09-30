"""V2 black layout with enlarged damage view and compact graphical legends.

Adapted from render_fragmentation.py; neither that file nor physics results change.
All force vectors are macro-state samples, not invented subcycle oscillations.
"""
import argparse,csv,hashlib,json,os,shutil
from pathlib import Path
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import proj3d
import imageio.v2 as imageio
import imageio_ffmpeg


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def spaced_ids(X,ids,count):
    """Fixed material IDs, selected by deterministic reference-space coverage."""
    q=X[ids];chosen=[int(np.argmin(q[:,0]+q[:,1]+q[:,2]))]
    d=np.full(len(q),np.inf)
    for _ in range(min(count,len(ids))-1):
        d=np.minimum(d,np.sum((q-q[chosen[-1]])**2,axis=1));d[chosen]=-1
        chosen.append(int(np.argmax(d)))
    return ids[chosen]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();root=Path(__file__).resolve().parents[1];run=a.run.resolve();out=a.output.resolve()
    if out.exists():raise FileExistsError(f'Use a new output directory: {out}')
    # Protect all existing output and core files; do not overwrite the previous movies.
    protect=[root/'scripts/render_fragmentation.py',root/'scripts/render_damage_forces_v1.py',*root.glob('pd_clot/*.py'),*root.glob('configs/*.json')]
    protect.extend(p for p in run.rglob('*') if p.is_file())
    protect.extend(p for p in (root/'visualization/streaming_fragmentation_demo').rglob('*') if p.is_file())
    before={str(p):sha(p) for p in sorted(set(protect))}
    out.mkdir(parents=True);(out/'frames').mkdir();(out/'provenance').mkdir()
    (out/'provenance/INPUT_FILES_SHA256.json').write_text(json.dumps(before,indent=2)+'\n')
    shutil.copy2(root/'scripts/render_fragmentation.py',out/'provenance/render_fragmentation_original.py')
    shutil.copy2(root/'scripts/render_damage_forces_v1.py',out/'provenance/render_damage_forces_v1_original.py')
    shutil.copy2(Path(__file__),out/'provenance/render_damage_forces_v2.py')
    z=np.load(run/'states.npz');c=json.loads((run/'CONFIG.json').read_text());history=json.loads((run/'history.json').read_text())
    X=z['X'];x=z['x'];v=z['v'];fixed=z['fixed'];vol=z['volume'];t=np.array([row['mechanical_time_s'] for row in history])
    area=c['clot']['particle_spacing_m']**2*c['surface']['area_factor'];mass=vol*c['clot']['density_kg_m3']
    surface=[];relax=[];saved_uf=[];vtk_matches=True
    for k in range(len(t)):
        mesh=pv.read(run/f'vtk/particles_{k:04d}.vtp')
        vtk_matches &= np.array_equal(mesh.points,x[k]) and np.array_equal(mesh['velocity'],v[k]) and np.array_equal(mesh['damage'],z['damage'][k])
        fs=np.array(mesh['surface_traction_Pa'])*area;fh=np.zeros_like(fs)
        free=(~z['attached'][k])&~fixed
        uf=np.array(mesh['fluid_velocity_m_s']);saved_uf.append(uf)
        if c['transport']['enabled']:fh[free]=mass[free,None]*(uf[free]-v[k,free])/c['transport']['tau_h_s']
        surface.append(fs);relax.append(fh)
    surface=np.array(surface);relax=np.array(relax);force=surface+relax
    assert vtk_matches and np.all(force[:,fixed]==0)
    # Read-only provider evaluation checks that saved traction and velocity match
    # the exact source model. No trajectory integration or solver is run.
    from pd_clot.geometry import make_cloud
    from pd_clot.fragment_fluid import FragmentFluid
    identity=json.loads((run/'IDENTITY.json').read_text())
    source_matches=all(sha(root/p)==s for p,s in identity['source_sha256'].items())
    assert source_matches
    cloud=make_cloud(c['clot']);fluid=FragmentFluid(c);max_traction_error=0.;max_velocity_error=0.
    for k in range(len(t)):
        calc,_,_=fluid.surface_force(cloud,x[k],z['integrity'][k],z['surface'][k],z['attached'][k],t[k])
        max_traction_error=max(max_traction_error,float(np.max(np.abs(calc-surface[k]))))
        max_velocity_error=max(max_velocity_error,float(np.max(np.abs(fluid.velocity(x[k],t[k])-saved_uf[k]))))
    assert max_traction_error<1e-18 and max_velocity_error<1e-14
    ids=spaced_ids(X,np.flatnonzero(~fixed),54)
    # One GLOBAL linear scale, no framewise normalization, clipping or saturation.
    arrow_mm_per_uN=.35
    vectors=force*1e6*arrow_mm_per_uN;positions=x*1e3
    tips=positions[:,ids]+vectors[:,ids]
    lower=np.minimum(positions.min(axis=(0,1)),tips.min(axis=(0,1)))-.12
    upper=np.maximum(positions.max(axis=(0,1)),tips.max(axis=(0,1)))+.12
    lower[0]=min(lower[0],-.8);lower[1]=min(lower[1],-.55);lower[2]=min(lower[2],-1.2)
    upper[1]=max(upper[1],.55);upper[2]=max(upper[2],.02)
    net=force.sum(axis=1)*1e6;peak=np.linalg.norm(force,axis=2).max(axis=1)*1e6
    np.savez_compressed(out/'force_samples.npz',cycles=z['cycles'],mechanical_time_s=t,positions_m=x,particle_damage=z['damage'],
        surface_force_N=surface,relaxation_force_N=relax,external_fluid_force_N=force,net_external_force_N=net*1e-6,
        displayed_arrow_particle_ids=ids,displayed_arrow_vectors_mm=vectors[:,ids],fixed=fixed)
    rows=[]
    for k,row in enumerate(history):
        rows.append(dict(state=k,represented_cycles=int(z['cycles'][k]),mechanical_time_s=t[k],mean_damage=row['mean_damage'],maximum_damage=row['maximum_damage'],
            net_force_x_uN=net[k,0],net_force_y_uN=net[k,1],net_force_z_uN=net[k,2],maximum_particle_force_uN=peak[k],
            visible_force_arrows=int(np.sum(np.linalg.norm(force[k,ids],axis=1)>0))))
    with (out/'force_history.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    dark='#000000';white='#ffffff';muted='#b0bfd2';cyan='#49d5ff'
    colors=['#49d5ff','#bc8aff','#ffad6b']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14})
    frame_paths=[];layout_checks=[];main_rect=[.00,.17,.765,.75];main_zoom=1.50
    for k,row in enumerate(rows):
        fig=plt.figure(figsize=(19.2,10.8),facecolor=dark)
        ax=fig.add_axes(main_rect,projection='3d',facecolor=dark,computed_zorder=False)
        ax.set_proj_type('ortho')
        point_colors=plt.colormaps['inferno'](z['damage'][k])
        ax.scatter(*positions[k,~fixed].T,c=point_colors[~fixed],s=36,edgecolors='#acb8c7',linewidths=.35,depthshade=False,zorder=5)
        ax.scatter(*positions[k,fixed].T,c=point_colors[fixed],s=36,marker='s',edgecolors='#e6eaf0',linewidths=.65,depthshade=False,zorder=5)
        nz=np.linalg.norm(force[k,ids],axis=1)>0;shown=ids[nz]
        ax.quiver(*positions[k,shown].T,*vectors[k,shown].T,color=cyan,length=1,normalize=False,arrow_length_ratio=.22,linewidth=1.35,alpha=.88,zorder=8)
        ax.set(xlim=(lower[0],upper[0]),ylim=(lower[1],upper[1]),zlim=(lower[2],upper[2]),xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)')
        ax.set_box_aspect(upper-lower,zoom=main_zoom);ax.view_init(elev=24,azim=-68)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.label.set_color(white);axis.label.set_size(16);axis.line.set_color('#42556c');axis.set_pane_color((0,0,0,1))
            axis._axinfo['grid']['color']='#24364c'
        ax.tick_params(colors=white,labelsize=13,pad=2)
        ax.xaxis.labelpad=16;ax.yaxis.labelpad=12;ax.zaxis.labelpad=12
        fig.text(.042,.946,'Particle damage & external forces',color=white,fontsize=29,weight='bold')
        fig.text(.042,.892,f'N = {row["represented_cycles"]}',color=white,fontsize=19)
        fig.text(.79,.823,'Net external force',color=white,fontsize=17,weight='bold')
        fig.text(.79,.786,'All particles · saved states',color=muted,fontsize=13)
        chart=fig.add_axes([.795,.48,.18,.255],facecolor=dark)
        for j,label in enumerate(['Fx','Fy','Fz']):
            chart.plot(z['cycles']/1000,net[:,j],color=colors[j],lw=1.4,label=label,alpha=.8)
            chart.scatter(z['cycles'][k]/1000,net[k,j],color=colors[j],s=35,zorder=8)
        chart.axvline(z['cycles'][k]/1000,color='white',lw=1,ls='--',alpha=.7);chart.axhline(0,color=muted,lw=.5,alpha=.4)
        chart.set(xlim=(-.6,25.6),xlabel='N (×1000)',ylabel='Force (μN)');chart.tick_params(colors=muted,labelsize=11)
        chart.xaxis.label.set_color(muted);chart.yaxis.label.set_color(muted)
        for spine in chart.spines.values():spine.set_color('#42556c')
        chart.grid(alpha=.12);chart.legend(loc='upper right',frameon=False,labelcolor=white,fontsize=11,ncol=3,handlelength=1,columnspacing=.6)
        for j,label in enumerate(['Fx','Fy','Fz']):fig.text(.80,.406-j*.037,f'{label}   {net[k,j]:+.3f} μN',color=colors[j],fontsize=16)
        fig.text(.80,.249,f'Peak particle force\n{peak[k]:.3f} μN',color=white,fontsize=15,linespacing=1.6)
        cax=fig.add_axes([.155,.065,.47,.018]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='inferno'),cax=cax,orientation='horizontal')
        cb.ax.tick_params(colors=white,labelsize=12,pad=3);cb.set_label('Particle damage',color=white,fontsize=15,labelpad=4)
        # A graphic legend with an EXACT projected +X reference-force length.
        # The orthographic camera makes this reference independent of position.
        # Matplotlib's 3-D axes retain a square clip box even in this wide panel.
        # Use the verified figure-space bounds, preserving ALL enlarged particles.
        for collection in ax.collections:collection.set_clip_on(False)
        fig.canvas.draw()
        def pixel_coordinates(points):
            u,w,_=proj3d.proj_transform(*np.asarray(points).T,ax.get_proj())
            return ax.transData.transform(np.column_stack((u,w)))
        anchor=(lower+upper)/2
        ends=pixel_coordinates(np.array([anchor,anchor+[arrow_mm_per_uN,0,0]]))
        delta=(ends[1]-ends[0])/np.array([1920,1080]);start=np.array([.18,.144])
        fig.add_artist(FancyArrowPatch(start,start+delta,transform=fig.transFigure,color=cyan,arrowstyle='-|>',mutation_scale=15,linewidth=2,shrinkA=0,shrinkB=0))
        fig.text(float(start[0]+delta[0]+.012),.134,'External force · 1 μN',color=white,fontsize=15)
        fig.add_artist(Line2D([.475],[.143],transform=fig.transFigure,marker='s',markersize=9,markerfacecolor='none',markeredgecolor=white,linestyle='none'))
        fig.text(.49,.134,'Fixed base',color=white,fontsize=15)
        projected=pixel_coordinates(np.concatenate((positions[k],tips[k]),axis=0))
        bounds=np.array([projected.min(axis=0),projected.max(axis=0)])
        # Main data stay within the left plot region and away from header/legend.
        safe=bool(np.all(bounds[0]>[20,205]) and np.all(bounds[1]<[1420,925]))
        layout_checks.append(dict(state=k,data_bounds_pixels_bottom_origin=bounds.tolist(),data_inside_safe_region=safe,
            collection_clipping_disabled=all(not q.get_clip_on() for q in ax.collections)))
        assert safe,layout_checks[-1]
        path=out/f'frames/frame_{k:04d}.png';fig.savefig(path,dpi=100,facecolor=dark);plt.close(fig);frame_paths.append(path)
    os.environ['IMAGEIO_FFMPEG_EXE']=imageio_ffmpeg.get_ffmpeg_exe()
    with imageio.get_writer(out/'particle_damage_external_forces.mp4',fps=12,codec='libx264',quality=8,macro_block_size=1) as w:
        for p in frame_paths:
            im=imageio.imread(p)[:,:,:3]
            for _ in range(6):w.append_data(im)
    from PIL import Image
    images=[Image.open(p).convert('RGB') for p in frame_paths]
    images[0].save(out/'particle_damage_external_forces.gif',save_all=True,append_images=images[1:],duration=500,loop=0)
    for im in images:im.close()
    shutil.copy2(frame_paths[-1],out/'preview.png')
    changed=[p for p,s in before.items() if not Path(p).is_file() or sha(Path(p))!=s]
    manifest=dict(run=str(run),renderer=str(Path(__file__).resolve()),renderer_sha256=sha(Path(__file__)),
        input_state_sha256=sha(run/'states.npz'),saved_state_count=len(t),fps=12,repeats_per_state=6,duration_s=13,
        particle_coloring='damage only, including fixed-base particles',bubble_marker=False,bond_lines=False,fragment_ID_panel=False,clearance_plane=False,
        all_particles_retained=True,displacement_scale=1,equal_axis_units=True,
        camera=dict(elevation_deg=24,azimuth_deg=-68,projection='orthographic',limits_mm=np.column_stack((lower,upper)).tolist()),
        external_force_definition='saved surface traction * h^2 * area_factor, plus rho*V*(saved_uf-v)/tau on detached mobile particles; anchors zero',
        force_scope='Applied fluid force only; internal PD forces, numerical damping and constraint reactions excluded',
        sampling='26 saved macro endpoints after damage/topology update, all at the same 500 Hz load phase; no subcycle interpolation or invented oscillation',
        arrow_selection='54 fixed non-base material IDs by reference-space coverage; exactly zero forces omitted',
        arrow_particle_ids=ids.tolist(),arrow_mm_per_uN=arrow_mm_per_uN,arrows_normalized_per_frame=False,arrows_clipped=False,
        curve='Sum of external fluid force over ALL particles at each saved state',
        layout=dict(background=dark,main_title='Particle damage & external forces',right_title='Net external force',cycle_text_color=white,
            main_axes_rectangle=main_rect,projection_zoom=main_zoom,previous_projection_zoom=1.18,
            projected_linear_enlargement=(.75/.61)*(1.50/1.18),footer_explanatory_text=False,footer_statistics=False,
            legend='Projected 1 microNewton force arrow along +X and fixed-base square',legend_force_arrow_pixels=(delta*np.array([1920,1080])).tolist()),
        layout_checks=layout_checks,
        checks=dict(VTK_matches_NPZ=bool(vtk_matches),simulation_source_identity_matches=bool(source_matches),anchors_zero_force=bool(np.all(force[:,fixed]==0)),
            maximum_surface_force_reconstruction_error_N=max_traction_error,maximum_saved_fluid_velocity_error_m_s=max_velocity_error,
            protected_preexisting_files=len(before),changed_preexisting_files=changed))
    (out/'RENDER_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    assert not changed
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
