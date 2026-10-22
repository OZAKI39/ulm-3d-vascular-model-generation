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
    fig.suptitle(f'新流场有限尺寸微泡轨迹｜{len(rows)} 条｜dt = 0.5 ms',color='white',fontsize=17)
    labels=[Line2D([0],[0],color=COLORS[k],lw=2,label=f'{k}: {len(v)}') for k,v in groups.items()]
    fig.subplots_adjust(bottom=.14)
    fig.legend(handles=labels,loc='lower center',bbox_to_anchor=(.5,.055),ncol=3,labelcolor='white',frameon=False,fontsize=10)
    fig.text(.5,.025,'颜色表示最终出口／终态；原始三维坐标，保留近壁接触及未到达出口的轨迹',ha='center',color='#c0cad8',fontsize=10)
    fig.savefig(out/'01_trajectories.png',dpi=300,facecolor=fig.get_facecolor())
    fig.savefig(out/'01_trajectories.pdf',facecolor=fig.get_facecolor());plt.close(fig)
    plot.close()
    # Outcome and diameter statistics use every trajectory, not display subsamples.
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    keys=list(groups);counts=[len(groups[k]) for k in keys]
    bars=axes[0].bar(keys,counts,color=[COLORS[k] for k in keys]);axes[0].bar_label(bars)
    axes[0].set(ylabel='轨迹数',title='实际出口与终态（分母为全部轨迹）')
    axes[0].tick_params(axis='x',labelrotation=20)
    bins=np.linspace(.75,4,27)
    for key in keys:
        values=[r['diameter_um'] for r in rows if (r['outlet'] or r['status'])==key]
        axes[1].hist(values,bins=bins,histtype='step',lw=1.8,color=COLORS[key],label=key)
    axes[1].set(xlabel='微泡直径（μm）',ylabel='轨迹数',title='固定入口队列的直径与去向');axes[1].legend(fontsize=8)
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
        actor=plot.add_mesh(cloud,color=COLORS[key],point_size=7,render_points_as_spheres=True,lighting=True)
        clouds[key]=cloud;selection[key]=ids
    plot.show(auto_close=False)
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
            draw.text((40,25),f'新流场微泡轨迹 | {len(rows)} 条 | dt = 0.5 ms | 轨迹年龄 {age*1000:.1f} ms',font=title_font,fill='white')
            draw.text((40,1018),'独立轨迹按年龄对齐回放，非真实同时注入；标记大小仅用于显示，实际微泡直径见数据',font=small_font,fill='#cad4e2')
            x=40
            for key in groups:
                draw.text((x,70),f'{key}: {len(groups[key])}',font=small_font,fill=COLORS[key]);x+=320
            writer.append_data(np.asarray(pixels))
    plot.close()
    manifest=dict(count=len(rows),render_backend=gpu_strings,projected_ports=projected,
        video_frames=frames,video_fps=fps,video_seconds=frames/fps,age_end_s=end,
        video_semantics='INDEPENDENT_TRAJECTORY_AGE_ALIGNED_NOT_SIMULTANEOUS_PHYSICAL_INFUSION',
        marker_size='7 pixels display only; physical diameter in cohort',all_tracks_used=True,
        encoder='libx264_CPU',ffmpeg_executable='/usr/bin/ffmpeg',smoke_only=animation_smoke,
        nvenc_available=False,nvenc_failure_log='logs/animation_smoke_retry.log')
    if animation_smoke:
        (out/'animation_smoke_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
        print(json.dumps(manifest),flush=True);return
    (HERE/'data/render_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    html='''<!doctype html><meta charset="utf-8"><title>新流场微泡轨迹</title>
<style>body{max-width:1300px;margin:35px auto;background:#101824;color:#edf2fa;font:18px sans-serif}img,video{width:100%}a{color:#50d9ef}</style>
<h1>新流场微泡轨迹：1500 条，dt = 0.5 ms</h1>
<p>独立有限尺寸微泡积分；动画按轨迹年龄对齐，不能视为真实同时注入。</p>
<p><a href="MICROBUBBLE_RESULTS_ZH.md">中文结果说明</a> · <a href="data/trajectory_catalog.csv">逐轨迹 CSV</a> · <a href="data/outlet_summary.csv">出口统计</a> · <a href="data/trajectories_dt0p5ms.npz">0.5 ms 位置数据</a></p>
<video controls src="animations/microbubble_age_aligned.mp4"></video>
<img src="figures/01_trajectories.png"><img src="figures/02_outlet_and_size.png">'''
    (HERE/'OPEN_RESULTS.html').write_text(html)
    print(json.dumps(manifest),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--preview',action='store_true');p.add_argument('--animation-smoke',action='store_true')
    a=p.parse_args();main(a.preview or a.animation_smoke,a.animation_smoke)
