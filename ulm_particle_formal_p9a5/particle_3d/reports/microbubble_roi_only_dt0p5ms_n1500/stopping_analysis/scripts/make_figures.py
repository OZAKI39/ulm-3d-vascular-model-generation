"""Existing-data figures only: stop locations, sizes, and velocity evidence."""
from pathlib import Path
import os,csv,json
import numpy as np
HERE=Path(__file__).resolve().parents[1];RUN=HERE.parent
os.environ.setdefault('MPLCONFIGDIR',str(HERE/'logs/mplconfig'))
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt,font_manager
from matplotlib.collections import PolyCollection
FONT=RUN/'assets/NotoSansCJKsc-Regular.otf';font_manager.fontManager.addfont(str(FONT))
plt.rcParams.update({'font.family':font_manager.FontProperties(fname=str(FONT)).get_name(),
    'font.size':11,'axes.unicode_minus':False,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})

def rows(name):return list(csv.DictReader((HERE/'data'/name).open()))
def save(fig,name):
    fig.savefig(HERE/'figures'/f'{name}.png',dpi=300,facecolor='white')
    fig.savefig(HERE/'figures'/f'{name}.pdf',dpi=300,facecolor='white');plt.close(fig)

def main():
    allrows=rows('all_tracks.csv');stops=rows('stopped_tracks.csv');ts=rows('example_timeseries.csv')
    summary=json.loads((HERE/'data/analysis_summary.json').read_text());a=np.load(RUN/'data/gpu_mesh_input.npz')
    tri=a['wall_triangles']*1e6
    xyz=np.array([[float(r[f'final_{axis}_um']) for axis in 'xyz'] for r in stops]);diam=np.array([float(r['diameter_um']) for r in stops])
    fig,axes=plt.subplots(1,2,figsize=(13,7),gridspec_kw={'width_ratios':[1,1]})
    for ax,(i,j) in zip(axes,[(0,2),(1,2)]):
        ax.add_collection(PolyCollection(tri[:,:,[i,j]],facecolor='#dce3e9',edgecolor='none',rasterized=True))
        sc=ax.scatter(xyz[:,i],xyz[:,j],c=diam,cmap='plasma',s=22,linewidths=.25,edgecolors='white',zorder=4)
        for label,p in [('J1',[92,49,111]),('J2',[130.04,82.04,87.18])]:
            ax.plot(p[i],p[j],'+',color='#172b42',ms=10,mew=1.5,zorder=5)
            ax.annotate(label,(p[i],p[j]),xytext=(18,12),textcoords='offset points',arrowprops={'arrowstyle':'-','color':'#172b42'},zorder=6)
        ax.autoscale_view();ax.set_aspect('equal');ax.set(xlabel='XYZ'[i]+' (μm)',ylabel='XYZ'[j]+' (μm)',title=f"{'XYZ'[i]}–{'XYZ'[j]} 投影：真实壁面与停止中心")
        ax.grid(alpha=.12)
    cax=fig.add_axes([.885,.245,.014,.55])
    fig.colorbar(sc,cax=cax,label='停止微泡直径 (μm)')
    fig.suptitle('最新流场：185 / 1500 条轨迹的中途停止位置',fontsize=18,y=.975)
    counts={r['region']:r['count'] for r in summary['terminal_regions']}
    fig.text(.08,.065,f"入口侧：{counts['INLET_to_J1']} 条    J1 半径 5 μm：{counts['J1_r5um']} 条    J1→J2 主干：{counts['J1_to_J2']} 条",fontsize=12)
    fig.text(.08,.026,'圆点表示微泡中心；点大小不代表物理直径。区域按中心位置与已核实管网路径归属，不代表接触面所在边界。',fontsize=10,color='#4b5563')
    fig.subplots_adjust(left=.07,right=.84,bottom=.14,top=.90,wspace=.30)
    save(fig,'01_stop_locations')

    fig,axes=plt.subplots(2,2,figsize=(13,9));fig.subplots_adjust(left=.09,right=.97,bottom=.095,top=.89,hspace=.48,wspace=.32)
    ax=axes[0,0]
    completed=np.array([float(r['diameter_um']) for r in allrows if r['status']=='COMPLETED'])
    bins=np.arange(.8,2.8501,.05)
    ax.hist(completed,bins=bins,color='#447ca8',label=f'已穿出口：{len(completed)}',alpha=.8)
    ax.hist(diam,bins=bins,color='#d96635',label=f'接触支持静止：{len(diam)}',alpha=.9)
    ax.axvspan(completed.max(),diam.min(),color='black',alpha=.2)
    ax.set(xlabel='微泡直径 (μm)',ylabel='条数',title='A  本批结果随粒径明显分离');ax.legend(fontsize=9)
    ax.text(.02,.78,f'通过最大：{completed.max():.4f} μm\n停止最小：{diam.min():.4f} μm',transform=ax.transAxes,fontsize=9)
    ax=axes[0,1];rr=[r for r in ts if int(r['particle_id'])==2]
    t=np.array([float(r['age_ms']) for r in rr]);vp=np.array([float(r['particle_speed_mm_s']) for r in rr]);vf=np.array([float(r['fluid_speed_mm_s']) for r in rr])
    ex=next(r for r in stops if int(r['particle_id'])==2)
    ax.plot(t,vf,color='#447ca8',label='同位置流体速度',lw=1.7)
    # Saved particle velocity is held over the PRECEDING accepted interval.
    ax.step(t,vp,where='pre',color='#d96635',label='微泡平移速度',lw=1.4)
    ax.axvline(float(ex['near_zero_translation_onset_ms']),ls=':',color='#666',lw=1)
    ax.set(xlabel='轨迹年龄 (ms)',ylabel='速度 (mm/s)',title=f"B  微泡 #2，直径 {float(ex['diameter_um']):.3f} μm")
    ax.legend(fontsize=9);ax.text(.44,.49,'流体仍在流动，\n微泡中心停止平移',transform=ax.transAxes,fontsize=11)
    ax=axes[1,0]
    groups=[('流场速度','local_fluid_speed_mm_s','#447ca8'),('近壁阻力后\n接触约束前','unconstrained_hydro_speed_mm_s','#45a185'),('接触约束后','constrained_particle_speed_mm_s','#d96635')]
    for k,(label,key,color) in enumerate(groups):
        v=np.array([float(r[key]) for r in stops]);jit=np.linspace(-.17,.17,len(v))
        ax.scatter(k+jit,v,s=8,alpha=.45,color=color)
        ax.plot([k-.23,k+.23],[np.median(v)]*2,color='black',lw=1.7)
    ax.set_yscale('log');ax.set_xticks(range(3),[r[0] for r in groups]);ax.set(ylabel='平移速度 (mm/s，对数轴)',title='C  全部 185 条停止轨迹的末次状态')
    ax.grid(axis='y',alpha=.15)
    ax=axes[1,1]
    for region,label,color in [('INLET_to_J1','入口→J1','#8055b4'),('J1_r5um','J1 半径 5 μm','#d96635'),('J1_to_J2','J1→J2','#447ca8')]:
        r=[r for r in stops if r['region']==region]
        ax.scatter([float(x['diameter_um']) for x in r],[float(x['distance_to_J1_um']) for x in r],s=18,label=f'{label} ({len(r)})',color=color,alpha=.75)
    ax.axhline(5,color='#999',ls=':',lw=1)
    ax.set(xlabel='微泡直径 (μm)',ylabel='停止中心距 J1 (μm)',title='D  不同粒径停止于不同局部位置');ax.legend(fontsize=9)
    fig.suptitle('停止原因证据：有限尺寸、近壁作用与多面接触约束',fontsize=18,y=.97)
    fig.text(.09,.025,'数值来自原始轨迹和同位置静态复核；未推进任何新轨迹。约 2.02 μm 的分界仅描述本批样本，不是普适通行阈值。',fontsize=10,color='#4b5563')
    save(fig,'02_stop_mechanism')
    print('FIGURES_READY')

if __name__=='__main__':main()
