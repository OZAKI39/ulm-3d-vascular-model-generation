"""V3 rendering: Arial regular, source glyph, 60 fps display interpolation.

Derived from a copy of render_damage_forces_v2.py. Positions are interpolated
only for display. Damage remains the preceding original saved value.
"""
import argparse,csv,hashlib,json,os,shutil,time
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import proj3d
import imageio.v2 as imageio
import imageio_ffmpeg
from PIL import Image


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--font',type=Path,default=Path('/mnt/c/Windows/Fonts/arial.ttf'))
    a=ap.parse_args();root=Path(__file__).resolve().parents[1];run=a.run.resolve();out=a.output.resolve()
    if out.exists():raise FileExistsError(f'Use a new directory: {out}')
    font_manager.fontManager.addfont(str(a.font));prop=font_manager.FontProperties(family='Arial',weight='normal')
    font=Path(font_manager.findfont(prop,fallback_to_default=False));assert sha(font)==sha(a.font)
    protect=[root/'scripts/render_damage_forces_v2.py',*root.glob('pd_clot/*.py'),*root.glob('configs/*.json')]
    protect.extend(p for p in run.rglob('*') if p.is_file())
    protect.extend(p for p in (root/'visualization/streaming_fragmentation_demo').rglob('*') if p.is_file())
    before={str(p):sha(p) for p in sorted(set(protect))}
    out.mkdir(parents=True);(out/'keyframes').mkdir();(out/'provenance').mkdir()
    (out/'provenance/INPUT_FILES_SHA256.json').write_text(json.dumps(before,indent=2)+'\n')
    shutil.copy2(root/'scripts/render_damage_forces_v2.py',out/'provenance/render_damage_forces_v2_original.py')
    shutil.copy2(Path(__file__),out/'provenance/render_damage_source_v3.py')
    z=np.load(run/'states.npz');c=json.loads((run/'CONFIG.json').read_text());history=json.loads((run/'history.json').read_text())
    x=z['x'];D=z['damage'];fixed=z['fixed'];n=len(x);source=np.array(c['streaming']['bubble_center_m'])
    prior=json.loads((root/'visualization/streaming_fragmentation_demo/damage_forces_v2/RENDER_MANIFEST.json').read_text())
    bounds=np.array(prior['camera']['limits_mm']);lower=bounds[:,0];upper=bounds[:,1]
    fps=60;stride=30;nf=(n-1)*stride+1;index=np.arange(nf)
    left=np.minimum(index//stride,n-1);right=np.minimum(left+1,n-1);alpha=(index%stride)/stride;alpha[-1]=0
    display_x=(1-alpha[:,None,None])*x[left]+alpha[:,None,None]*x[right];display_D=D[left]
    display_x[::stride]=x
    display_x[:,fixed]=z['X'][fixed]  # Preserve exact fixed coordinates, including roundoff.
    assert np.array_equal(display_x[::stride],x) and np.array_equal(display_D[::stride],D)
    assert np.all(display_x[:,fixed]==z['X'][fixed])
    times=np.array([h['mechanical_time_s'] for h in history])
    np.savez_compressed(out/'display_samples.npz',display_positions_m=display_x,display_damage=display_D,source_state_left=left,
        source_state_right=right,interpolation_alpha=alpha,video_time_s=index/fps,display_proxy_time_s=(1-alpha)*times[left]+alpha*times[right],
        original_cycles=z['cycles'],fixed=fixed,source_position_m=source)
    with (out/'frame_map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['frame','video_time_s','left_state','right_state','position_alpha','damage_state','is_original_state'])
        for i in index:w.writerow([int(i),i/fps,int(left[i]),int(right[i]),float(alpha[i]),int(left[i]),bool(i%stride==0)])
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','font.size':14,'axes.labelweight':'normal','axes.titleweight':'normal'})
    fig=plt.figure(figsize=(19.2,10.8),facecolor='#000000')
    ax=fig.add_axes([.015,.16,.95,.77],projection='3d',facecolor='#000000',computed_zorder=False)
    ax.set_proj_type('ortho');ax.view_init(elev=24,azim=-68)
    ax.set(xlim=tuple(bounds[0]),ylim=tuple(bounds[1]),zlim=tuple(bounds[2]),xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)')
    ax.set_box_aspect(upper-lower,zoom=1.70)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
        axis.label.set_color('white');axis.label.set_size(17);axis.line.set_color('#42556c');axis.set_pane_color((0,0,0,1))
        axis._axinfo['grid']['color']='#24364c'
    ax.tick_params(colors='white',labelsize=14,pad=2)
    ax.xaxis.labelpad=16;ax.yaxis.labelpad=12;ax.zaxis.labelpad=12
    colors=plt.colormaps['inferno'](D[0])
    moving=ax.scatter(*(x[0,~fixed]*1000).T,c=colors[~fixed],s=38,edgecolors='#acb8c7',linewidths=.35,depthshade=False,zorder=5)
    base=ax.scatter(*(x[0,fixed]*1000).T,c=colors[fixed],s=38,marker='s',edgecolors='#e6eaf0',linewidths=.65,depthshade=False,zorder=5)
    # Only the existing localized loading centre; no sphere or force vector.
    marker=ax.scatter(*(source*1000),marker='D',s=74,facecolors='none',edgecolors='white',linewidths=1.6,depthshade=False,zorder=9)
    for item in [moving,base,marker]:item.set_clip_on(False)
    title=fig.text(.042,.944,'Particle damage & load source',color='white',fontsize=29,fontproperties=prop,weight='normal')
    cax=fig.add_axes([.245,.06,.51,.018]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap='inferno'),cax=cax,orientation='horizontal')
    cb.ax.tick_params(colors='white',labelsize=13,pad=3);cb.set_label('Particle damage',color='white',fontsize=16,labelpad=4)
    handles=[Line2D([],[],marker='D',markersize=9,markerfacecolor='none',markeredgecolor='white',linestyle='none',label='Load source'),
        Line2D([],[],marker='s',markersize=9,markerfacecolor='none',markeredgecolor='white',linestyle='none',label='Fixed base')]
    fig.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.136),ncol=2,frameon=False,labelcolor='white',fontsize=16,handletextpad=.5,columnspacing=3)
    fig.canvas.draw();assert title.get_fontproperties().get_weight()=='normal' and title.get_fontproperties().get_name()=='Arial'
    allpoints=np.concatenate((x.reshape(-1,3)*1000,(source*1000)[None,:]))
    u,vp,_=proj3d.proj_transform(*allpoints.T,ax.get_proj());px=ax.transData.transform(np.column_stack((u,vp)))
    screen=np.array([px.min(axis=0),px.max(axis=0)])
    assert np.all(screen[0]>[20,200]) and np.all(screen[1]<[1870,945]),screen
    os.environ['IMAGEIO_FFMPEG_EXE']=imageio_ffmpeg.get_ffmpeg_exe();started=time.perf_counter();gif_images=[]
    with imageio.get_writer(out/'particle_damage_load_source_60fps.mp4',fps=fps,codec='libx264',quality=8,macro_block_size=1,
            ffmpeg_params=['-pix_fmt','yuv420p','-movflags','+faststart']) as writer:
        for i in index:
            pos=display_x[i]*1000;col=plt.colormaps['inferno'](display_D[i])
            moving._offsets3d=tuple(pos[~fixed].T);moving.set_facecolor(col[~fixed])
            base._offsets3d=tuple(pos[fixed].T);base.set_facecolor(col[fixed])
            fig.canvas.draw();im=np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy();writer.append_data(im)
            if i%stride==0:imageio.imwrite(out/f'keyframes/state_{i//stride:04d}.png',im)
            if i%3==0:gif_images.append(Image.fromarray(im).resize((1280,720),Image.Resampling.LANCZOS))
            if i%60==0 or i==nf-1:print(json.dumps(dict(frame=int(i),total=nf,elapsed_s=time.perf_counter()-started)),flush=True)
    plt.close(fig)
    gif_images[0].save(out/'particle_damage_load_source_20fps.gif',save_all=True,append_images=gif_images[1:],duration=50,loop=0)
    for im in gif_images:im.close()
    shutil.copy2(out/f'keyframes/state_{n-1:04d}.png',out/'preview.png')
    changed=[p for p,h in before.items() if not Path(p).is_file() or sha(Path(p))!=h];assert not changed
    manifest=dict(run=str(run),renderer=str(Path(__file__).resolve()),renderer_sha256=sha(Path(__file__)),input_states_sha256=sha(run/'states.npz'),
        original_state_count=n,particle_count=x.shape[1],fps=fps,frame_count=nf,duration_s=nf/fps,repeated_hold_frames_inserted=0,
        positions='Piecewise linear display interpolation; exact saved positions every 30 frames; no extrapolation',
        damage='Original recorded value held until the next macro endpoint; no damage interpolation',
        interpolation_is_new_simulation=False,original_endpoint_indices=(np.arange(n)*stride).tolist(),
        glyph='Small hollow white diamond at configured streaming.bubble_center_m; localized load centre, not a bubble surface',
        source_position_m=source.tolist(),source_modulation_visualized=False,
        background_flow='Existing distributed pipe field remains in the simulation; no invented point emitter drawn',
        force_arrows=False,force_chart=False,force_numbers=False,cycle_label=False,bubble_sphere=False,
        title='Particle damage & load source',font=dict(family='Arial',weight='normal',path=str(font),sha256=sha(font)),background='#000000',
        camera=dict(elevation_deg=24,azimuth_deg=-68,projection='orthographic',limits_mm=bounds.tolist(),zoom=1.70),
        displacement_scale=1,all_particles_retained=True,screen_bounds_pixels_bottom_origin=screen.tolist(),collection_clipping=False,
        elapsed_render_s=time.perf_counter()-started,protected_files=len(before),changed_preexisting_files=changed)
    (out/'RENDER_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
