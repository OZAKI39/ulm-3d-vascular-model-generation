"""Scientific trajectory plots and RTX/EGL age-aligned motion animation."""
from pathlib import Path
import argparse,json,os
os.environ.setdefault('VTK_DEFAULT_OPENGL_WINDOW','vtkEGLRenderWindow')
import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from PIL import Image,ImageDraw,ImageFont

HERE=Path(__file__).resolve().parents[1]
COLORS={'O1':'#b994ff','O2':'#39d8e6','O3':'#ffb153','SUPPORTED_STATIONARY':'#bac2cf',
        'LONG_RESIDENCE_CENSORED':'#f672a7','SOLVER_FAILURE':'#ff4b4b'}
BACKGROUND='#101824'
ANIMATION_ZOOM=1.30
MARKER_SIZE_PX=7*ANIMATION_ZOOM
LABELS={'O1':'O1','O2':'O2','O3':'O3',
        'SUPPORTED_STATIONARY':'Contact-supported stationary',
        'LONG_RESIDENCE_CENSORED':'Long-residence censored','SOLVER_FAILURE':'Solver failure'}


def projected_bounds(plot,points):
    matrix=plot.camera.GetCompositeProjectionTransformMatrix(plot.renderer.GetTiledAspectRatio(),-1,1)
    matrix=np.array([[matrix.GetElement(i,j) for j in range(4)] for i in range(4)])
    clip=np.column_stack((points,np.ones(len(points))))@matrix.T
    ndc=clip[:,:2]/clip[:,3,None]
    width,height=plot.window_size
    pixels=np.column_stack(((ndc[:,0]+1)*width/2,(1-ndc[:,1])*height/2))
    return np.array([pixels.min(axis=0),pixels.max(axis=0)])


def enlarge_animation(plot,surface):
    """Exact 1.3x orthographic display zoom; retain direction and physical data."""
    old_scale=plot.camera.parallel_scale
    before=projected_bounds(plot,surface.points)
    plot.camera.parallel_scale=old_scale/ANIMATION_ZOOM
    enlarged=projected_bounds(plot,surface.points)
    # Center within the existing header/footer whitespace, without changing view direction.
    target=np.array([960.,558.])
    delta=target-enlarged.mean(axis=0)
    forward=np.asarray(plot.camera.direction);forward/=np.linalg.norm(forward)
    right=np.cross(forward,np.asarray(plot.camera.up));right/=np.linalg.norm(right)
    up=np.cross(right,forward)
    shift=(-delta[0]*right+delta[1]*up)*(2*plot.camera.parallel_scale/plot.window_size[1])
    plot.camera.position=np.asarray(plot.camera.position)+shift
    plot.camera.focal_point=np.asarray(plot.camera.focal_point)+shift
    after=projected_bounds(plot,surface.points)
    ratios=(after[1]-after[0])/(before[1]-before[0])
    assert np.allclose(ratios,ANIMATION_ZOOM,rtol=0,atol=1e-10),ratios
    assert np.all(after[0]>=[40,108]) and np.all(after[1]<=[1880,1008]),after
    return dict(original_bounds_px=before.tolist(),enlarged_bounds_px=after.tolist(),
                projected_size_ratios=ratios.tolist(),original_parallel_scale=old_scale,
                parallel_scale=plot.camera.parallel_scale,display_translation_px=delta.tolist(),
                camera_position=list(plot.camera.position),camera_focal_point=list(plot.camera.focal_point))


def main(preview=False,animation_smoke=False):
    config=json.loads((HERE/'config.json').read_text())
    source=Path(config['source_root'])
    font=source/'particle_3d/reports/particle9a5_formal_trajectories/assets/NotoSansCJKsc-Regular.otf'
    if not font.exists():font=HERE/'assets/NotoSansCJKsc-Regular.otf'
    font_manager.fontManager.addfont(str(font));plt.rcParams['font.family']='sans-serif'
    plt.rcParams['font.sans-serif']=[font_manager.FontProperties(fname=font).get_name()]
    plt.rcParams['axes.unicode_minus']=False
    # Preview may use completed rows only; it is never published as full output.
    folders=sorted((HERE/'tracks').glob('mb_*'))
    folders=[p for p in folders if (p/'COMPLETE.json').exists()]
    rows=[json.loads((p/'metrics.json').read_text()) for p in folders]
    tracks=[np.load(p/'trajectory.npz')['samples'] for p in folders]
    if not preview:assert len(rows)==1500
    if not rows:raise ValueError('No completed track records')
    out=HERE/('preview' if preview else 'figures');out.mkdir(exist_ok=True)
    mesh=pv.read(config['flow_path']);surface=mesh.extract_surface(algorithm='dataset_surface');surface.points*=1e6
    center=np.asarray(surface.center);radius=np.linalg.norm(np.ptp(surface.points,axis=0))/2
    direction=np.array([-.72,.69,.11]);direction/=np.linalg.norm(direction)
    plot=pv.Plotter(off_screen=True,window_size=(2400,1500));plot.set_background(BACKGROUND)
    plot.add_mesh(surface,color='#b2becf',opacity=.11,smooth_shading=True)
    plot.camera.position=center+direction*radius*4
    plot.camera.focal_point=center;plot.camera.up=(0,0,1)
    plot.camera.parallel_projection=True;plot.camera.parallel_scale=radius*.82
    groups={o:[] for o in ['O1','O2','O3']}
    for row,track in zip(rows,tracks):
        key=row['outlet'] or row['status'];groups.setdefault(key,[]).append(track[:,1:4]*1e6)
    for key,paths in groups.items():
        if not paths:continue
        points=np.concatenate(paths);off=0;lines=[]
        for p in paths:
            lines.extend([len(p),*range(off,off+len(p))]);off+=len(p)
        poly=pv.PolyData(points);poly.lines=np.array(lines)
        plot.add_mesh(poly,color=COLORS[key],line_width=1.35,opacity=.68,lighting=False)
    plot.show(auto_close=False)
    capabilities=plot.ren_win.ReportCapabilities()
    caps=Path(config['flow_path']).parents[1]/'SV_MESH/mesh-surfaces'
    import vtk
    projected={}
    for name in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        path=caps/(name+'.vtp')
        if path.exists():
            xyz=np.asarray(pv.read(path).center)*1e6
            point=vtk.vtkCoordinate();point.SetCoordinateSystemToWorld();point.SetValue(*xyz)
            xy=point.GetComputedDoubleDisplayValue(plot.renderer)
            projected[name]=[xy[0]/2400,1-xy[1]/1500]
    image=plot.screenshot(return_img=True)
    fig,ax=plt.subplots(figsize=(12,8),facecolor=BACKGROUND)
    ax.imshow(image);ax.set_axis_off()
    for role,(x,y) in projected.items():
        label='INLET' if role=='INLET' else 'O'+str(int(role[-2:]))
        xx=.05 if x<.5 else .95;yy=np.clip(y,.08,.92)
        ax.annotate(label,xy=(x,y),xycoords='axes fraction',xytext=(xx,1-yy),
            textcoords='axes fraction',color='white',fontsize=12,
            ha='left' if xx<.5 else 'right',arrowprops=dict(arrowstyle='->',color='white'))
    # image coordinate has inverted y, while axes-fraction y increases upward.
    for annotation in list(ax.texts):
        if hasattr(annotation,'xy'):
            annotation.xy=(annotation.xy[0],1-annotation.xy[1])
    fig.suptitle(f'Finite-size microbubble trajectories | {len(rows)} tracks | dt = 0.5 ms',color='white',fontsize=17)
    labels=[Line2D([0],[0],color=COLORS[k],lw=2,label=f'{LABELS[k]}: {len(v)}') for k,v in groups.items()]
    fig.subplots_adjust(bottom=.14)
    fig.legend(handles=labels,loc='lower center',bbox_to_anchor=(.5,.055),ncol=3,labelcolor='white',frameon=False,fontsize=10)
    fig.text(.5,.025,'Color indicates final outlet or terminal state; original 3D coordinates and near-wall contact trajectories',ha='center',color='#c0cad8',fontsize=10)
    fig.savefig(out/'01_trajectories.png',dpi=300,facecolor=fig.get_facecolor())
    fig.savefig(out/'01_trajectories.pdf',facecolor=fig.get_facecolor());plt.close(fig)
    plot.close()
    # Outcome and diameter statistics use every trajectory, not display subsamples.
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    keys=list(groups);counts=[len(groups[k]) for k in keys]
    bars=axes[0].bar([LABELS[k] for k in keys],counts,color=[COLORS[k] for k in keys]);axes[0].bar_label(bars)
    axes[0].set(ylabel='Number of tracks',title='Final outlets and terminal states (all tracks)')
    axes[0].tick_params(axis='x',labelrotation=20)
    bins=np.linspace(.75,4,27)
    for key in keys:
        values=[r['diameter_um'] for r in rows if (r['outlet'] or r['status'])==key]
        axes[1].hist(values,bins=bins,histtype='step',lw=1.8,color=COLORS[key],label=LABELS[key])
    axes[1].set(xlabel='Microbubble diameter (μm)',ylabel='Number of tracks',title='Diameter and outcome of the fixed inlet cohort');axes[1].legend(fontsize=8)
    fig.savefig(out/'02_outlet_and_size.png',dpi=300);fig.savefig(out/'02_outlet_and_size.pdf');plt.close(fig)
    gpu_strings=[s for s in capabilities.splitlines() if any(k in s for k in ['vendor string','renderer string','version string'])]
    if preview and not animation_smoke:
        (out/'preview.json').write_text(json.dumps(dict(count=len(rows),render_backend=gpu_strings,projected_ports=projected),indent=2))
        return
    # Explicitly age-aligned ensemble playback: not a physical infusion movie.
    os.environ['IMAGEIO_FFMPEG_EXE']='/usr/bin/ffmpeg'
    import imageio.v2 as imageio
    video=(out/'animations' if preview else HERE/'animations');video.mkdir(exist_ok=True)
    plot=pv.Plotter(off_screen=True,window_size=(1920,1080));plot.set_background(BACKGROUND)
    plot.add_mesh(surface,color='#b2becf',opacity=.16,smooth_shading=True)
    plot.camera.position=center+direction*radius*4;plot.camera.focal_point=center;plot.camera.up=(0,0,1)
    plot.camera.parallel_projection=True;plot.camera.parallel_scale=radius*.78
    clouds={};selection={}
    for key in groups:
        ids=[i for i,r in enumerate(rows) if (r['outlet'] or r['status'])==key]
        if not ids:continue
        cloud=pv.PolyData(np.array([tracks[i][0,1:4]*1e6 for i in ids]))
        actor=plot.add_mesh(cloud,color=COLORS[key],point_size=MARKER_SIZE_PX,render_points_as_spheres=True,lighting=True)
        clouds[key]=cloud;selection[key]=ids
    plot.show(auto_close=False)
    animation_layout=enlarge_animation(plot,surface)
    end=max(s[-1,0] for s in tracks);frames=24 if animation_smoke else 288;fps=24
    title_font=ImageFont.truetype(str(font),30);small_font=ImageFont.truetype(str(font),23)
    # This instance's NVENC advertises an encoder but actual initialization
    # returns "unsupported device (2)". Keep real NVIDIA EGL rendering and
    # use the available CPU encoder; do not claim unavailable hardware encoding.
    with imageio.get_writer(video/'microbubble_age_aligned.mp4',fps=fps,codec='libx264',
                           quality=8,macro_block_size=1,ffmpeg_params=['-threads','2']) as writer:
        for age in np.linspace(0,end,frames):
            for key,ids in selection.items():
                xyz=[]
                for i in ids:
                    s=tracks[i]
                    if age>s[-1,0] and rows[i]['status']=='COMPLETED':xyz.append([np.nan]*3)
                    else:xyz.append([np.interp(age,s[:,0],s[:,j])*1e6 for j in (1,2,3)])
                clouds[key].points=np.array(xyz)
            plot.render();pixels=Image.fromarray(plot.screenshot(return_img=True))
            draw=ImageDraw.Draw(pixels)
            title=f'Microbubble trajectories | {len(rows)} tracks | dt = 0.5 ms | Trajectory age: {age*1000:.1f} ms'
            footer='Age-aligned independent trajectories, not simultaneous injection. Marker size is for display; physical diameters are in the data.'
            assert title_font.getlength(title)<=1840 and small_font.getlength(footer)<=1840
            draw.text((40,25),title,font=title_font,fill='white')
            draw.text((40,1018),footer,font=small_font,fill='#cad4e2')
            x=40
            for key in groups:
                legend=f'{LABELS[key]}: {len(groups[key])}'
                assert x+small_font.getlength(legend)<=1880
                draw.text((x,70),legend,font=small_font,fill=COLORS[key]);x+=320
            writer.append_data(np.asarray(pixels))
    plot.close()
    manifest=dict(count=len(rows),render_backend=gpu_strings,projected_ports=projected,
        video_frames=frames,video_fps=fps,video_seconds=frames/fps,age_end_s=end,
        video_semantics='INDEPENDENT_TRAJECTORY_AGE_ALIGNED_NOT_SIMULTANEOUS_PHYSICAL_INFUSION',
        marker_size=f'{MARKER_SIZE_PX:g} pixels display only; physical diameter in cohort',all_tracks_used=True,
        text_language='English',animation_zoom_factor=ANIMATION_ZOOM,animation_layout=animation_layout,
        encoder='libx264_CPU',ffmpeg_executable='/usr/bin/ffmpeg',smoke_only=animation_smoke,
        nvenc_available=False,nvenc_policy='Reuse previously verified libx264; no new NVENC test')
    if animation_smoke:
        (out/'animation_smoke_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(json.dumps(manifest),flush=True);return
    (HERE/'data/render_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    html='''<!doctype html><html lang="en"><meta charset="utf-8"><title>Microbubble trajectories</title>
<style>body{max-width:1300px;margin:35px auto;background:#101824;color:#edf2fa;font:18px sans-serif}img,video{width:100%}a{color:#50d9ef}</style>
<h1>Microbubble trajectories: 1500 tracks, dt = 0.5 ms</h1>
<p>Independent finite-size trajectories, aligned by trajectory age for playback; not simultaneous injection.</p>
<p><a href="MICROBUBBLE_RESULTS_ZH.md">Results report (Chinese)</a> · <a href="data/trajectory_catalog.csv">Track catalog CSV</a> · <a href="data/outlet_summary.csv">Outlet statistics</a> · <a href="data/trajectories_dt0p5ms.npz">0.5 ms position data</a></p>
<video controls src="animations/microbubble_age_aligned.mp4"></video>
<img src="figures/01_trajectories.png"><img src="figures/02_outlet_and_size.png"></html>'''
    (HERE/'OPEN_RESULTS.html').write_text(html)
    print(json.dumps(manifest),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preview',action='store_true');p.add_argument('--animation-smoke',action='store_true')
    a=p.parse_args();main(a.preview or a.animation_smoke,a.animation_smoke)
