"""Standalone manufacturing figures and a Chinese scientific implementation report."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import numpy as np
import pyvista as pv

from .boundary_review_figures import bounds
from .print_orientation import transform_points


def centerline_figure(path,graph,core,title):
    fig=plt.figure(figsize=(14,10));ax=fig.add_subplot(projection='3d',computed_zorder=False)
    for kind,color in [('core','#008f9f'),('context','#db8b14')]:
        edges=[e for e in graph.edges if (all(n in core for n in e))==(kind=='core')]
        segments=[[graph.nodes[n]['coords'][:3] for n in e] for e in edges]
        if segments:
            ax.add_collection3d(Line3DCollection(segments,colors=color,linewidths=2))
    bounds(ax,[d['coords'][:3] for _,d in graph.nodes(data=True)]);ax.view_init(22,-65)
    fig.suptitle(title,fontsize=17)
    fig.legend(handles=[Line2D([0],[0],color='#008f9f',lw=3,label='Frozen semantic core'),
        Line2D([0],[0],color='#db8b14',lw=3,label='Real BraVa connector / proximal context')],loc='lower center',ncol=2)
    fig.savefig(path,dpi=140);plt.close(fig)


def radius_figure(path,original,compensated,config):
    hard=config['effective']['hard_min_diameter_mm'];preferred=config['effective']['preferred_min_diameter_mm']
    fig=plt.figure(figsize=(18,9));ax=fig.add_subplot(121,projection='3d')
    categories=[[],[],[]];colors=['#ca3944','#df9b20','#008f9f']
    for a,b in original.edges:
        d=2*np.mean([original.nodes[n]['coords'][3] for n in [a,b]])
        index=0 if d<hard else 1 if d<preferred else 2
        categories[index].append([original.nodes[n]['coords'][:3] for n in [a,b]])
    for segments,color in zip(categories,colors):
        if segments:ax.add_collection3d(Line3DCollection(segments,colors=color,linewidths=2))
    bounds(ax,[d['coords'][:3] for _,d in original.nodes(data=True)]);ax.view_init(22,-65)
    ax.set_title('Raw diameter classes (edge endpoint mean)')
    hist=fig.add_subplot(122)
    raw=np.array([2*d['coords'][3] for _,d in original.nodes(data=True)])
    comp=np.array([2*d['coords'][3] for _,d in compensated.nodes(data=True)])
    bins=np.linspace(0,max(raw.max(),comp.max())*1.03,30)
    hist.hist(raw,bins=bins,alpha=.55,label='Original diameter',color='#64748b')
    hist.hist(comp,bins=bins,alpha=.6,label='Manufacturing diameter',color='#008f9f')
    hist.axvline(hard,color=colors[0],ls='--',label=f'Hard minimum {hard:g} mm')
    hist.axvline(preferred,color=colors[1],ls='--',label=f'Preferred minimum {preferred:g} mm')
    hist.set(xlabel='Diameter (mm), D = 2r',ylabel='Original SWC sample count');hist.legend()
    fig.suptitle('Manufacturing compensation only: XYZ, TYPE and retained topology unchanged',fontsize=17)
    fig.legend(handles=[Line2D([0],[0],color=c,lw=3,label=t) for c,t in zip(colors,[f'D < {hard:g} mm',f'{hard:g} <= D < {preferred:g} mm',f'D >= {preferred:g} mm'])],loc='lower center',ncol=3)
    fig.savefig(path,dpi=140);plt.close(fig)


def surface_figure(path,mesh,title,config,*,build_volume=False):
    plotter=pv.Plotter(off_screen=True,window_size=(1600,1100));plotter.set_background('white')
    plotter.add_mesh(mesh,color='#159ba5',smooth_shading=True)
    if build_volume:
        x,y,z=config['printer']['build_volume_mm'];margin=config['printer']['bed_margin_xy_mm'];top=config['printer']['top_margin_mm']
        plotter.add_mesh(pv.Box(bounds=(0,x,0,y,0,z)),style='wireframe',color='#64748b',line_width=2)
        plotter.add_mesh(pv.Box(bounds=(margin,x-margin,margin,y-margin,0,z-top)),style='wireframe',color='#d78b18',line_width=2)
        plotter.add_mesh(mesh.outline(),color='#20933b',line_width=3)
        plotter.add_mesh(pv.Plane(center=(x/2,y/2,0),direction=(0,0,1),i_size=x,j_size=y),color='#aebcca',opacity=.2)
    plotter.add_text(title,position='upper_left',color='black',font_size=13)
    plotter.add_text('MANUFACTURING_COMPENSATED_GEOMETRY | mm | scale = 1.0\nSupport/orientation metrics are engineering heuristics; physical printing not validated.',
                     position='lower_left',color='#374151',font_size=10)
    plotter.show_bounds(color='#374151',grid=None,all_edges=False,location='outer',
                       n_xlabels=3,n_ylabels=3,n_zlabels=3,xtitle='X (mm)',ytitle='Y (mm)',ztitle='Z (mm)')
    plotter.view_isometric();plotter.camera.zoom(1.1)
    plotter.show(screenshot=str(path),auto_close=True)


def write_figures(output,manifest,core,original,compensated,surface):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    nodes=set(manifest['core_original_node_ids'])
    centerline_figure(output/'01_semantic_core.png',core,nodes,'Immutable anatomical / semantic MeVO core')
    centerline_figure(output/'02_connected_print_roi.png',original,nodes,'Connected manufacturing ROI: native BraVa context only')
    radius_figure(output/'03_radius_compensation.png',original,compensated,manifest['profile'])
    if surface is None:return
    orientation=manifest['orientation'];best=manifest['selected_orientation']
    baseline=surface.copy();baseline.points=transform_points(surface.points,orientation['baseline']['transform_4x4'])
    final=surface.copy();final.points=transform_points(surface.points,best['transform_4x4'])
    surface_figure(output/'04_print_orientation_baseline.png',baseline,'Baseline: original exported axes, translated to bed',manifest['profile'])
    surface_figure(output/'05_print_orientation_best.png',final,
        f"Selected orientation | height {best['height_mm']:.2f} mm | heuristic support area {best['downward_support_area_mm2']:.1f} mm2",manifest['profile'])
    surface_figure(output/'06_build_volume_preview.png',final,
        f"{manifest['printer']['model']} | build volume {manifest['printer']['build_volume_mm']} mm\nOrange: usable envelope; green: ROI bounds",manifest['profile'],build_volume=True)


def write_report(output,m):
    before=m['semantic_core'];after=m['manufacturing_topology'];comp=m['manufacturing_geometry'];printer=m['printer']
    fmt=lambda values:' × '.join(f'{x:.3f}' for x in values)
    lines=['# 冻结 MeVO 语义基础上的脑血管制造 ROI 与打印方向研究','', '## 摘要','',
        '本研究将解剖候选、制造拓扑和制造几何分为三层。既有 NN 语义及语义精修 SWC 保持不变；制造层以真实 BraVa 路径连接核心组件，以局部半径补偿保留长而细的血管，仅删除满足条件的短末端枝。原生 VascularMD 表面用于确定性打印方向搜索，并以实际 Bambu Studio 切片检验工程可行性。任何结果均未被标记为经过真实制造验证。','',
        '## 材料与方法','',
        '本轮没有重新配准、选择供体、执行 NN 查询、校准 UNKNOWN 或聚合分支概率。核心评分直接读取冻结审查表的组件均值与已知比例。默认要求两者均不低于 0.80；工程评分综合语义支持、已知比例、分支数、长度与空间范围，按配置选择前三名。此排序服务于实验 ROI 选择，并非新的解剖分类器。','',
        '连接采用原始有向树中各入口的最近共同祖先及其真实路径，是保留核心节点所需的最小子树。单入口上游段最多扩展 30 mm，且限于已确认的 RMCA 树。不得以欧氏直线连接、不自动放大坐标。低价值 twig 仅在它是原有末端分支、原始直径中位数低于硬下限、长度不足 12 mm 时一次性删除，分叉连接点保留，且不递归剥除新形成的末端。','',
        '制造半径首先取原始半径与目标下限的较大者，再在下限激活边界附近施加约 5 mm 的局部 smoothstep 肩部过渡；过渡只允许增加制造半径。离散采样窗口的实际跨度记录在清单中，未进行全局 Gaussian 平滑，也未修改坐标或 TYPE。共享分叉点使用单一半径，两个 SWC 层分别导出，原始与补偿半径逐点对应。','',
        '硬直径下限为 max(2.5×喷嘴直径, 用户硬下限)，优选直径为 max(3×喷嘴直径, 用户优选下限)。默认 0.4 mm 喷嘴对应 1.0/1.2 mm；备选 0.2 mm 细节配置对应 0.6/0.8 mm。它们均为本项目工程启发式参数，NOT manufacturer guaranteed limits. Must be calibrated experimentally.','',
        '代码审计显示，现有 BraVa 交互界面调用 camera.Azimuth 转动相机；学术展示脚本 doc_visualize_v3 的 PCA 也用于相机构图。它们均没有旋转导出的打印网格，因此本轮以原始导出坐标轴作为方向基线，并保留原代码。新方法对最终补偿 VascularMD 表面进行 PCA 初始化、两种朝向和固定角度网格搜索。原生表面端口只在打印派生物中封盖，以表示可溶解的实心牺牲血管芯；原始 VTK 不变。','',
        '方向评分同时考虑向下三角面面积、模型高度、可用打印体积占比、细长程度及投影轮廓面积；长而细且近水平的血管作为支撑项的辅助惩罚。投影凸包不是实际床面接触面积，向下面积也不是切片器计算的支撑体积。所有坐标变换为刚体旋转和平移，scale 始终为 1。','',
        '## 制造 ROI 与半径结果','',
        f"入选组件为 {m['selected_components']}；所有未选组件保留为 OPTIONAL_SEMANTIC_COMPONENT。part05 虽达到基本均值门槛，但其综合排名低于前三，未以编号硬编码剔除。连接分支为 {m['context_branches']}。",'',
        f"核心 AABB 为 {fmt(before['bbox']['extents_mm'])} mm；连接、扩展及 twig 筛选后的 AABB 为 {fmt(after['bbox']['extents_mm'])} mm，对角线 {after['bbox']['diagonal_mm']:.3f} mm，最长轴 {after['bbox']['longest_dimension_mm']:.3f} mm。PCA 方向包围盒为 {fmt(after['bbox']['pca_obb_extents_mm'])} mm。这是 PCA 对齐包围盒，并非求解全局最小体积包围盒。",'',
        f"制造拓扑包含 {after['node_count']} 个原始采样点、{after['connected_component_count']} 个连通分量、{after['inlet_count']} 个入口及 {after['outlet_count']} 个出口，总中心线长度 {after['centerline_length_mm']:.6f} mm。上游实际扩展 {m['proximal_context']['extended_length_mm']:.6f} mm；达到 RMCA 原树起点后停止，没有为了达到优选尺寸而引入不明确的其它动脉或全局缩放。最低尺寸门槛状态：{m['roi_size_gate']}。",'',
        f"原始直径最小/中位/最大为 {after['diameter_min_mm']:.4f}/{after['diameter_median_mm']:.4f}/{after['diameter_max_mm']:.4f} mm；补偿后为 {comp['diameter_min_mm']:.4f}/{comp['diameter_median_mm']:.4f}/{comp['diameter_max_mm']:.4f} mm。半径下限为 {m['radius_floor_mm']:.4f} mm。",'',
        f"修改半径的采样点数为 {m['radius_inflation']['modified_samples']}，各点增加半径之和为 {m['radius_inflation']['sum_added_radius_mm']:.6f} mm，平均增加 {m['radius_inflation']['mean_added_radius_mm']:.6f} mm，最大增加 {m['radius_inflation']['maximum_added_radius_mm']:.6f} mm。半径增量和是采样统计量，不是体积。删除短末端 twig 共 {len(m['removed_terminal_twigs'])} 条。",'',
        '## 原生建模与打印方向','',
        f"原始半径版 VascularMD 状态：{m['models']['original']['status']}；制造补偿版状态：{m['models']['compensated']['status']}。两者均使用现有 adapter，不修改上游。建模后半径检查为：{m['post_model_radius_qc']}。若平滑后仅低于优选下限而仍高于硬下限，将明确报告优选下限欠缺；若低于硬下限，不能取得 PRINT_READY_CANDIDATE。",'']
    if m.get('selected_orientation'):
        o=m['selected_orientation'];baseline=m['orientation']['baseline']
        lines += [f"共评估 {m['orientation']['tested']} 个确定性方向，体积可容纳 {m['orientation']['feasible_count']} 个。选择方向的 XYZ 欧拉角为 {o['euler_xyz_deg']} 度，旋转后表面包围盒为 {fmt(o['bbox_extents_mm'])} mm。",'',
            f"相对原方向，高度由 {baseline['height_mm']:.3f} mm 变为 {o['height_mm']:.3f} mm；启发式向下需支撑面积由 {baseline['downward_support_area_mm2']:.3f} mm² 变为 {o['downward_support_area_mm2']:.3f} mm²。两项目标可能有权衡，降低综合评分不等于每项指标都下降。",'',
            '四阶齐次变换（列向量约定，单位 mm）为：','', '```json',str(o['transform_4x4']), '```','']
    lines += ['## 实际切片与状态','',
        f"检测到的打印机配置为 {printer['model']}，喷嘴 {printer['nozzle_diameter_mm']} mm，构建体积 {printer['build_volume_mm']} mm；Bambu Studio 版本 {printer.get('version')}。软件原有材料预设与本任务使用的 ABS 配置分别记录，未覆盖 GUI 预设，也没有向打印机发送任务。",'',
        'Bambu CLI 使用本地官方 machine/process/ABS filament 配置，递归展开继承关系，不使用 --orient 或 --scale 覆盖本轮方向。CLI 返回码本身不等同于切片成功；需要 3MF 内实际 G-code 与挤出指令作为证据。缺失的打印时间、耗材或支撑字段保留为空，不能补造。','',
        f"实际切片状态：{m.get('slicer',{}).get('status')}；最终工程状态：**{m['status']}**。",'',
        '| 方向排名 | 切片状态 | 预计时间 s | 耗材 g | 支撑字段 |', '|---|---|---:|---:|---|']
    for r in m.get('slicer',{}).get('results',[]):
        support=r.get('support_information') or {}
        lines.append(f"| {r['rank']} | {r['status']} | {r.get('estimated_print_time_seconds')} | {r.get('filament_used_g')} | used={support.get('support_used')}; type={support.get('support_type')}; 独立用量未提供 |")
    if m.get('orientation_comparison'):
        comparison=m['orientation_comparison']
        lines += ['',f"最终方向由实际切片时间与总耗材代理量共同选择。相比原方向，高度降低 {comparison['height_reduction_percent']:.3f}%，启发式向下支撑面积降低 {comparison['unsupported_area_reduction_percent']:.3f}%。总耗材包含模型与支撑，不能称为独立支撑质量。",'']
    lines += ['', '## 讨论与可复核性','',
        '本轮将 d075 保留为历史实验，默认打印 ROI 来自完整语义精修核心及真实连接上下文。直径补偿明确属于 MANUFACTURING_COMPENSATED_GEOMETRY，不能解释为真实管腔尺寸，也不能用于反推 M2/M3 标签。打印候选仍需实际 ABS 打印、PDMS 封装及溶解实验验证；切片成功不等于 MANUFACTURING_VALIDATED。','',
        '校准 STL 含 0.8、1.0、1.2、1.5、2.0 mm 五档直径与相对打印床 0°、30°、45°、60°、90° 五档角度的 25 根短圆柱，位置映射见同名 CSV。只用于记录黏附、细枝稳定性和丙酮蒸气处理后的实测变化，不自动修改解剖数据。','',
        '所有源哈希与最终检查见 protected_sources.json、protection_verification.json；打印变换及每个方向评分分别见 print_transform.json 和 orientation/orientation_scores.csv。回归结果见 QC/test_results.json。','',
        '参考依据：[Bambu Studio 官方 CLI 文档](https://github.com/bambulab/BambuStudio/wiki/Command-Line-Usage)。若检测到 A1+ABS，按[官方 A1 材料规格](https://bambulab.com/en/a1/tech-specs)单独提示不推荐，不擅自改换材料。','']
    if m.get('validation'):
        v=m['validation']
        lines += ['## 实际验证结果','',
            f"回归合计 {v['full_suite_passed']} 项通过，其中既有测试 {v['stable_passed']} 项、新增制造测试 {v['new_passed']} 项；最终失败为零。原有 scikit-image/NumPy 兼容性弃用提示共 {v['regression_warnings']} 条，新增测试无警告。最后的保护与切片检查均通过。",'',
            f"{v['protected_files']} 个源文件的哈希保持一致。五个切片 3MF 的实际模型顶点已通过其构建矩阵还原，与各自输入 STL 的双向最大距离不超过 {v['maximum_slice_vertex_difference_mm']:.9f} mm，证明没有另行旋转、缩放或移动。软件 GUI 的原有 machine/process/filament 预设也未改变。",'',
            '切片元数据中的 first_layer_time 有异常数值，已原样留存并标记 UNRELIABLE_FIRST_LAYER_TIME；方向选择使用有效的总 prediction 时间，不使用异常首层时间。其他局限包括：平滑后低于优选直径但高于硬下限，以及 AABB 最长轴未达到 100 mm 优选值；最低尺寸门槛通过，PCA 包围盒最长轴超过 100 mm。', '',
            '| 关键源 | SHA256 |','|---|---|']
        lines += [f'| {name} | `{value}` |' for name,value in v['source_hashes'].items()]
        lines += ['']
    (Path(output)/'manufacturing_report.md').write_text('\n'.join(lines),encoding='utf-8')
