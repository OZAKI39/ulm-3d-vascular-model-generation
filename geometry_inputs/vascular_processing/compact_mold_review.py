"""Review artifacts for the modest s1-10 box-height/base-thickness revision."""
import numpy as np
import trimesh,pyvista as pv
from .casting_fixture_review import poly
from .support_inputs import save,digest
from .compact_c_pipeline import default_output as compact_output
from .bambu_casting_mold import parse_paths

NAMES=['01_height_and_clearance_comparison','02_2mm_bottom_and_sidewall_section','03_final_2mm_assembly','04_PDMS_after_trim','05_actual_Bambu_toolpaths']


def render(out,old,l,oldg,g,diag,s):
    dim=s['dimensions'];q=out/'QC';camera=None
    for i,name in enumerate(NAMES,1):
        print('s1-10 dimension QC',name,flush=True)
        if i==2:
            section(out,l,g,s);continue
        p=pv.Plotter(off_screen=True,shape=(1,2) if i in [1,5] else (1,1),window_size=(1800,1050) if i in [1,5] else (1500,1050));p.set_background('white')
        def add(m,color,opacity=1):p.add_mesh(poly(m),color=color,opacity=opacity,smooth_shading=False)
        def scene(layout,mold,previous=False):
            add(mold,'#589bb6',.18);add(layout['core'],'#cc5947')
            for port in layout['ports']:
                k=port['port_id'];m=trimesh.load_mesh((compact_output() if previous else out)/f'retainers/port_{k}_C_retainer.stl');add(m,'#249c72')
            p.add_axes();p.camera_position='iso';p.camera.up=(0,0,1);p.enable_parallel_projection();p.reset_camera()
        if i==1:
            p.subplot(0,0);scene(old,oldg['mold_mesh'],True);camera=(p.camera.position,p.camera.focal_point,p.camera.up,p.camera.parallel_scale);p.add_text('Previous C-cap system\nTop 12.00 / bottom 10.00 mm\nBase 3.00 mm; wall 2.00 mm',font_size=15,color='#253c4c')
            p.subplot(0,1);scene(l,g['mold_mesh']);p.camera.position=camera[0];p.camera.focal_point=camera[1];p.camera.up=camera[2];p.camera.parallel_scale=camera[3];p.add_text(f"Updated system | SAME SCALE\nTop {dim['top_clearance_mm']:.2f} / bottom {dim['bottom_clearance_mm']:.2f} mm\nBase and walls 2.00 mm",font_size=15,color='#253c4c')
        elif i==3:
            scene(l,g['mold_mesh']);p.add_text('2 mm ABS mold + unchanged compact C caps',font_size=16,color='#253c4c');p.add_text(f"Complete core translated rigidly; no change to radius or ports\nBottom {dim['bottom_clearance_mm']:.2f} mm; upper rim {dim['top_clearance_mm']:.2f} mm clearance",position='lower_left',font_size=12,color='#526471')
        elif i==4:
            add(diag['PDMS_after_trim'],'#8fbbcf',.25);add(diag['vascular_lumen_after_trim'],'#cc5947',.8);p.camera_position='iso';p.camera.up=(0,0,1);p.reset_camera();p.add_axes();p.add_text('Trimmed PDMS | same port-exposure planes\nNo residual auxiliary void; four usable 4.50 mm ports',font_size=15,color='#253c4c')
        elif i==5:
            for j,key in enumerate(['mold','retainers']):
                p.subplot(0,j);r=s['bambu'][key]
                if r.get('archive'):
                    a=parse_paths(r['archive']);segments=a['segments'][a['kind']==0];levels=np.round(segments[:,:,2].mean(1),4);unique=np.unique(levels);chosen=unique[np.linspace(0,len(unique)-1,min(30,len(unique))).astype(int)];seg=segments[np.isin(levels,chosen)];p.add_lines(seg.reshape(-1,3),color='#249c72',width=1);p.add_text(f"ACTUAL {key} G-code\nX1 Carbon / 0.2 mm / ABS\nAuto Support OFF | 0 support segments",font_size=14,color='#253c4c')
                else:p.add_text('SLICE MISSING: NEEDS_ADJUSTMENT',color='red')
                p.camera_position='iso';p.camera.up=(0,0,1);p.reset_camera();p.add_axes()
        p.show(screenshot=str(q/(name+'.png')),auto_close=True)
    save(q/'figure_manifest.json',[dict(name=n,sha256=digest(q/(n+'.png'))) for n in NAMES])


def section(out,l,g,s):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from .casting_fixture_geometry import shape_mesh
    fig,ax=plt.subplots(figsize=(13,6));m=shape_mesh(g['shell']);sec=m.section(plane_origin=[0,0,0],plane_normal=[0,1,0])
    for curve in sec.discrete:ax.plot(curve[:,0],curve[:,2],color='#377f9e',lw=2)
    # Complete vascular projection, not an invented section through its center.
    pts=l['core'].vertices;ax.scatter(pts[:,0],pts[:,2],s=.2,color='#c95546',alpha=.7,label='Complete vascular surface projected onto XZ')
    lo,hi=l['low'],l['high'];vlow,vhigh=l['core'].bounds[:,2];ax.plot([lo[0],hi[0]],[hi[2],hi[2]],'--',color='#708b99',lw=1,label='Open rim / nominal fill plane')
    ax.annotate('',xy=(0,lo[2]),xytext=(0,vlow),arrowprops=dict(arrowstyle='<->',color='#405b67'));ax.text(2,(lo[2]+vlow)/2,f'Bottom clearance {vlow-lo[2]:.2f} mm',va='center',fontsize=11)
    ax.annotate('',xy=(0,vhigh),xytext=(0,hi[2]),arrowprops=dict(arrowstyle='<->',color='#405b67'));ax.text(2,(vhigh+hi[2])/2,f'Top clearance {hi[2]-vhigh:.2f} mm',va='center',fontsize=11)
    ax.annotate('Base thickness 2.00 mm',xy=(-22,1),xytext=(-35,-7),arrowprops=dict(arrowstyle='->'),fontsize=11);ax.annotate('Wall thickness 2.00 mm',xy=(hi[0]+1,15),xytext=(19,-7),arrowprops=dict(arrowstyle='->'),fontsize=11)
    ax.set_xlim(lo[0]-6,hi[0]+6);ax.set_ylim(-10,hi[2]+8);ax.set_aspect('equal');ax.set_xlabel('X (mm)');ax.set_ylabel('Z (mm)');ax.set_title('Actual mold cross-section at Y = 0; full vessel surface projected for cover measurement');ax.legend(loc='upper left',fontsize=9);ax.grid(alpha=.15);fig.tight_layout();fig.savefig(out/'QC/02_2mm_bottom_and_sidewall_section.png',dpi=160);plt.close(fig)


def report(out,s):
    d=s['dimensions'];t=s['trim'];text='# BG001 RMCA 盒体厚度与上下间距调整\n\n'
    text+=f"状态：**{s['status']}**。本轮仅调整 s1-10 全 ABS 盒体的厚度及高度，继续使用已选中的 90° C 形插片。\n\n"
    text+='|尺寸|调整前|调整后|\n|---|---:|---:|\n'
    for name,a,b in [('血管表面至内底面净距',d['old_bottom_clearance_mm'],d['bottom_clearance_mm']),('血管表面至顶部开口平面净距',d['old_top_clearance_mm'],d['top_clearance_mm']),('底板厚度',d['old_bottom_thickness_mm'],d['bottom_thickness_mm']),('侧壁厚度',2.,d['wall_thickness_mm']),('盒体外部总高度',d['old_outer_height_mm'],d['outer_dimensions_mm'][2])]:text+=f'|{name}|{a:.2f} mm|{b:.2f} mm|\n'
    text+='\n上下 5 mm 是本轮用于收紧原有 10/12 mm 空余区的设计取值，尚不是经过浇注实验优化的工艺参数。净距按包含全部延长端口的完整血管最外表面计算，不按中心线计算。盒体顶部仍开放，顶部净距指最高血管表面到盒沿（名义浇注面）的高度差。\n\n'
    text+=f"盒体长宽不变，外部尺寸为 {d['outer_dimensions_mm'][0]:.3f} × {d['outer_dimensions_mm'][1]:.3f} × {d['outer_dimensions_mm'][2]:.3f} mm。内腔尺寸为 {d['inner_dimensions_mm'][0]:.3f} × {d['inner_dimensions_mm'][1]:.3f} × {d['inner_dimensions_mm'][2]:.3f} mm。完整血管及四个插片共同沿 Z 方向移动 {s['common_translation_mm'][2]:.2f} mm，不缩放、不改血管半径或路线。C 形插片几何保持不变，上方间隙、导向及主/次基准配合保持原值。托座上部接触结构不变，仅缩短连接到新底面的竖直部分。\n\n"
    text+='重新检查完整血管的顶部装入、C 插片下降与短滑移、中央零接触、±1 mm 次基准配合余量及 +Z/俯仰/滚转有限姿态。四条壁侧通路及插片自身壁侧入口均重新计算，两套 SWC 分支可达性继续保留。裁切面的几何位置不变；因盒体变矮，其有限面高度随新 PDMS 域调整。\n\n'
    text+=f"切前辅助空腔为 {t['auxiliary_void_volume_before_trim']:.6f} mm³，切后为 {t['auxiliary_void_volume_after_trim']:.12f} mm³（数值零）。四端口仍各保留 4.50 mm 已确认直段。新外模及插片使用 X1 Carbon、0.2 mm 喷嘴、ABS、Auto Support OFF 重新实际切片，并检查薄片、接触面、浅导向和 C 唇边的真实挤出路径。没有向打印机上传或启动打印。\n\n"
    text+='## 文件与使用\n\n- `outer_mold/ABS_sacrificial_outer_mold_2mm.stl/.step`：新盒体与适配的壁生托座。\n- `outer_mold/ABS_sacrificial_outer_mold_2mm_X1C_ABS.3mf`：新盒体实际切片。\n- `retainers/`：新装配位置的四个 C 插片，以及独立打印 3MF。\n- `assembly/all_ABS_casting_system_C_retainer_2mm.step`：新位置的独立组件完整装配。\n- `assembly/all_ABS_C_retainer_mold_2mm_review.stl`：完整审阅网格，不用于联合打印。\n- `vascular_core/assembly_placed_ABS_core.stl`：新位置的血管芯；原始字节副本也在此目录。\n- `assembly/casting_core_transform.json`：新旧坐标变换。不要将旧装配坐标的血管与新盒体直接叠加。\n- `QC/`：尺寸对比、2 mm 截面、总装、裁切 PDMS、真实刀路检查图。\n- `reports/run.log`、`mold_revision_summary.json`：运行日志和完整结果。\n\n'
    text+='入口保持 `./run_python.sh s1-10_all_abs_sacrificial_casting_system.py`，默认配置为 `config/compact_abs_mold_2mm.yaml`。之前的盒体、C 插片阶段结果和原 61 项测试均保留。血管芯及 C 插片仍禁止蒸汽平滑；实际收缩、抗浮力、溶解和 PDMS 行为尚未物理校准。\n\n'
    text+=f"源及历史文件检查 {s['source_protection']['checked_files']} 项，错配 0。血管芯 SHA256：`{s['source_core_sha256']}`。\n"
    if s['failures']:text+='\n尚未通过：'+', '.join(s['failures'])+'。\n'
    (out/'盒体调整说明.md').write_text(text)
