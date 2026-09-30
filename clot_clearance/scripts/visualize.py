"""Scientific PNG/PDF, offline interactive HTML and animation from saved states."""
import json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import imageio.v2 as imageio
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT=Path(__file__).resolve().parents[1]


def main():
    run=ROOT/'results/straight_pipe_demo';out=ROOT/'visualization';out.mkdir(exist_ok=True)
    figs=out/'figures';figs.mkdir(exist_ok=True);frames=out/'frames';frames.mkdir(exist_ok=True)
    z=np.load(run/'states.npz');c=json.loads((run/'CONFIG.json').read_text());h=json.loads((run/'history.json').read_text())
    x=z['x']*1e3;X=z['X']*1e3;damage=z['damage'];cycles=z['cycles'];vmax=max(.01,float(damage.max()))
    bubble=np.array(c['streaming']['bubble_center_m'])*1e3
    dark='#07101e';white='#e8f0fa';blue='#51b9ff'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    # Report figure: actual scale, common damage range across all snapshots.
    fig=plt.figure(figsize=(13,9),layout='constrained')
    ax=fig.add_subplot(221,projection='3d')
    xx=np.linspace(-3,3,25);theta=np.linspace(0,2*np.pi,48);xx,tt=np.meshgrid(xx,theta);R=c['pipe']['radius_m']*1e3
    ax.plot_surface(xx,R*np.cos(tt),R*np.sin(tt),alpha=.08,color='#879aae',linewidth=0)
    ax.scatter(*x[-1].T,c=damage[-1],cmap='inferno',vmin=0,vmax=vmax,s=12)
    ax.scatter(*bubble,color=blue,s=25);ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',title='Straight pipe / fixed clot wall patch')
    ax.set_box_aspect((2.4,1,1));ax.view_init(elev=23,azim=-62)
    ax=fig.add_subplot(222,projection='3d')
    sc=ax.scatter(*x[-1].T,c=damage[-1],cmap='inferno',vmin=0,vmax=vmax,s=20)
    base=X[z['fixed']];ax.scatter(*base.T,facecolors='none',edgecolors='#3182bd',s=29,linewidth=.6)
    ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',title='Final particle damage (true displacement scale)')
    ax.set_box_aspect((3,2,1));ax.view_init(elev=27,azim=-63);fig.colorbar(sc,ax=ax,shrink=.65,label='Mean incident bond damage')
    ax=fig.add_subplot(223)
    ax.plot(cycles,[r['mean_particle_damage'] for r in h],'o-',label='Mean particle damage')
    ax.plot(cycles,[r['maximum_particle_damage'] for r in h],'s-',label='Maximum particle damage')
    ax.set(xlabel='Represented carrier cycles N',ylabel='Damage',title='Empirical cycle-jump response');ax.grid(alpha=.2);ax.legend()
    ax=fig.add_subplot(224)
    ax.plot(cycles,[r['maximum_displacement_um'] for r in h],'o-',color='#155e9d')
    ax.set(xlabel='Represented carrier cycles N',ylabel='Maximum saved displacement (um)',title='Deformation; final fragments = 1, broken bonds = 0');ax.grid(alpha=.2)
    fig.suptitle('One-way streaming → NOSB-PD clot prototype\nSynthetic verification load; no clinical thrombolysis prediction',fontsize=16)
    fig.savefig(figs/'workflow_results.png',dpi=300);fig.savefig(figs/'workflow_results.pdf');plt.close(fig)
    # Evidence panel for the deliberately imposed cleavage test.
    import pyvista as pv
    cleavage=pv.read(ROOT/'verification/F_artificial_cleavage/vtk/particles_0000.vtp')
    fig=plt.figure(figsize=(8,5),layout='constrained');ax=fig.add_subplot(projection='3d')
    ax.scatter(*(cleavage.points*1e3).T,c=cleavage['fragment_id'],cmap='tab10',s=38,vmin=0,vmax=9)
    ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',title='Artificial kinematic cleavage: 2 connected components\nSoftware verification only; not produced by the pipe load')
    fig.savefig(figs/'artificial_fragmentation.png',dpi=300);fig.savefig(figs/'artificial_fragmentation.pdf');plt.close(fig)
    # Same saved macro states, no synthetic interpolation of physical trajectories.
    images=[]
    for i in range(len(cycles)):
        fig=plt.figure(figsize=(12.8,7.2),facecolor=dark);ax=fig.add_axes([.01,.12,.62,.75],projection='3d',facecolor=dark,computed_zorder=False)
        sc=ax.scatter(*x[i].T,c=damage[i],cmap='inferno',vmin=0,vmax=vmax,s=38,depthshade=False,edgecolors='#8b9aac',linewidths=.3)
        ax.scatter(*bubble,color=blue,s=85,edgecolors=white,linewidths=.6,depthshade=False,zorder=20)
        ax.plot([bubble[0],bubble[0]],[bubble[1],bubble[1]],[bubble[2]-.04,-.59],color=blue,alpha=.5,ls='--')
        ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',xlim=(-.8,.8),ylim=(-.55,.55),zlim=(-1.12,-.39))
        ax.set_box_aspect((2.4,1.6,1));ax.view_init(elev=24,azim=-64)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.label.set_color(white);axis.set_pane_color((.04,.08,.13,1));axis.line.set_color('#526276');axis._axinfo['grid']['color']='#304153'
        ax.tick_params(colors=white)
        fig.text(.055,.94,'STRAIGHT PIPE · CLOT DAMAGE',color=white,fontsize=23,weight='bold')
        fig.text(.76,.80,f'N = {int(cycles[i]):,} cycles',color=blue,fontsize=20,va='top')
        fig.text(.76,.70,f'Mean damage  {damage[i].mean():.4f}\nMax damage    {damage[i].max():.4f}\n\nFragments      {h[i]["number_of_fragments"]}\nBroken bonds   {h[i]["broken_bond_fraction"]:.1%}',color=white,fontsize=14,linespacing=1.65,va='top')
        fig.text(.76,.31,'Fixed base\nTrue displacement scale\nBlue marker: prescribed\nbubble center (enlarged)',color='#a8b8cd',fontsize=12,linespacing=1.5,va='top')
        cb=fig.colorbar(sc,cax=fig.add_axes([.675,.32,.013,.39]));cb.set_label('Particle damage',color=white);cb.ax.tick_params(colors=white)
        fig.text(.05,.035,'Synthetic load / empirical cycle jump · 9 saved macro states · not a resolved ultrasound simulation',color='#a8b8cd',fontsize=12)
        path=frames/f'frame_{i:03d}.png';fig.savefig(path,dpi=100,facecolor=dark);plt.close(fig);images.append(imageio.imread(path))
    import imageio_ffmpeg,os
    os.environ['IMAGEIO_FFMPEG_EXE']=imageio_ffmpeg.get_ffmpeg_exe()
    with imageio.get_writer(out/'clot_damage.mp4',fps=12,codec='libx264',quality=8,macro_block_size=1) as w:
        for im in images:
            for _ in range(12):w.append_data(im[:,:,:3])
    imageio.mimsave(out/'clot_damage.gif',images,duration=700,loop=0)
    # Offline rotatable 3-D view with a synchronized timeline and literal saved frames.
    fig=make_subplots(rows=1,cols=2,specs=[[{'type':'scene'},{'type':'xy'}]],column_widths=[.68,.32],horizontal_spacing=.05)
    fig.add_trace(go.Scatter3d(x=x[0,:,0],y=x[0,:,1],z=x[0,:,2],mode='markers',name='Clot particles',
        marker=dict(size=5,color=damage[0],colorscale='Inferno',cmin=0,cmax=vmax,colorbar=dict(title='Damage',x=.61,len=.6)),
        customdata=np.arange(len(X)),hovertemplate='ID %{customdata}<br>Damage %{marker.color:.5f}<extra></extra>'),row=1,col=1)
    fig.add_trace(go.Scatter3d(x=[bubble[0]],y=[bubble[1]],z=[bubble[2]],mode='markers',name='Prescribed bubble center (marker enlarged)',marker=dict(color=blue,size=7)),row=1,col=1)
    fig.add_trace(go.Scatter(x=cycles,y=[r['mean_particle_damage'] for r in h],mode='lines+markers',name='Mean damage',line=dict(color=blue)),row=1,col=2)
    fig.add_trace(go.Scatter(x=cycles,y=[r['maximum_particle_damage'] for r in h],mode='lines+markers',name='Max damage',line=dict(color='#fcba68')),row=1,col=2)
    fig.add_trace(go.Scatter(x=[cycles[0]],y=[damage[0].mean()],mode='markers',showlegend=False,marker=dict(size=15,color='#ffffff')),row=1,col=2)
    fig.frames=[go.Frame(name=str(i),data=[go.Scatter3d(x=x[i,:,0],y=x[i,:,1],z=x[i,:,2],marker=dict(color=damage[i])),
        go.Scatter(x=[cycles[i]],y=[damage[i].mean()])],traces=[0,4]) for i in range(len(cycles))]
    fig.update_layout(template='plotly_dark',paper_bgcolor=dark,plot_bgcolor=dark,height=660,
        scene=dict(xaxis_title='X (mm)',yaxis_title='Y (mm)',zaxis_title='Z (mm)',aspectmode='data',
                   xaxis=dict(range=[-.85,.85]),yaxis=dict(range=[-.55,.55]),zaxis=dict(range=[-1.15,-.4]),camera=dict(eye=dict(x=1.5,y=-2,z=1.2))),
        xaxis_title='Represented carrier cycles N',yaxis_title='Particle damage',legend=dict(orientation='h',y=1.1),
        margin=dict(l=20,r=30,t=70,b=110),
        updatemenus=[dict(type='buttons',x=.03,y=-.12,buttons=[dict(label='Play',method='animate',args=[None,dict(frame=dict(duration=600,redraw=True),transition=dict(duration=0),fromcurrent=True)]),dict(label='Pause',method='animate',args=[[None],dict(mode='immediate',frame=dict(duration=0,redraw=False))])])],
        sliders=[dict(x=.2,len=.75,y=-.08,currentvalue=dict(prefix='Saved cycle count: '),steps=[dict(label=str(int(cycles[i])),method='animate',args=[[str(i)],dict(mode='immediate',frame=dict(duration=0,redraw=True),transition=dict(duration=0))]) for i in range(len(cycles))])])
    fig.write_html(out/'interactive.html',include_plotlyjs=True,auto_play=False)
    html='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Clot clearance · 直管原型</title>
<style>body{margin:0;background:#07101e;color:#e8f0fa;font:16px/1.65 system-ui,sans-serif}main{max-width:1280px;margin:auto;padding:36px}h1{font-size:36px;margin-bottom:8px}h2{font-size:23px}a{color:#51b9ff}section{margin:30px 0;padding:22px;background:#101e30;border-radius:12px}.kpis{display:flex;gap:35px;flex-wrap:wrap}.kpis b{display:block;font-size:27px;color:#51b9ff}iframe{width:100%;height:690px;border:0}img,video{width:100%;border-radius:8px}.note{color:#b5c6da}code{color:#b4d6ff}</style><main>
<p class="note">BRAVA RESEARCH / ISOLATED STRAIGHT-PIPE EXAMPLE</p><h1>直管中的血栓变形与损伤</h1>
<p>规定局部微流牵引 → NOSB-PD 血栓 → 循环损伤 → 连通分量识别。当前结果用于验证计算流程。</p>
<div class="kpis"><div><b>384</b>血栓粒子</div><div><b>8,000</b>代表循环数</div><div><b>2.681%</b>平均粒子损伤</div><div><b>1</b>最终连通分量</div></div>
<section><strong>结果边界</strong><p>本次主示例发生轻度损伤，未出现断键或碎片，也未模拟血栓清除。蓝色点表示规定的假想微泡中心，标记经过放大；没有求解声场或微泡振荡。500 Hz 的代表性机械载荷与 1 MHz 的循环计数是不同时间尺度，不能视为实际超声疗效预测。</p></section>
<h2>交互结果</h2><p class="note">拖动旋转、滚轮缩放，使用下方时间条查看已保存宏步。位移采用真实比例；颜色范围在所有帧中一致。</p><iframe src="visualization/interactive.html" title="Interactive saved states"></iframe>
<section><h2>动画与报告图</h2><video controls preload="metadata" poster="visualization/frames/frame_008.png" src="visualization/clot_damage.mp4"></video><p><a href="visualization/clot_damage.mp4">MP4</a> · <a href="visualization/clot_damage.gif">GIF</a> · <a href="visualization/figures/workflow_results.pdf">PDF 图</a></p><img src="visualization/figures/workflow_results.png" alt="直管、血栓损伤和位移曲线"></section>
<section><h2>断裂功能验证</h2><p>另一个人工张开测试产生两个连通分量。它验证键断裂和碎片编号功能；张开位移由测试规定，不能解释为主示例的流动导致断裂。</p><img src="visualization/figures/artificial_fragmentation.png" alt="人工张开后的两个连通分量"></section>
<section><h2>可复查文件</h2><p><a href="REPORT_ZH.md">中文报告</a> · <a href="README.md">复现命令</a> · <a href="MODEL_ASSUMPTIONS.md">模型与限制</a> · <a href="references/PAPER_READING.md">论文对应说明</a></p><p><a href="results/straight_pipe_demo/particles.pvd">ParaView 粒子序列</a> · <a href="results/straight_pipe_demo/bonds.pvd">键序列</a> · <a href="results/straight_pipe_demo/traction.pvd">牵引序列</a> · <a href="results/straight_pipe_demo/summary.csv">统计 CSV</a> · <a href="verification/ACCEPTANCE.json">验证记录</a></p></section></main></html>'''
    (ROOT/'OPEN_RESULTS.html').write_text(html)
    print('Saved figures, MP4/GIF and offline HTML:',ROOT/'OPEN_RESULTS.html')


if __name__=='__main__':main()
