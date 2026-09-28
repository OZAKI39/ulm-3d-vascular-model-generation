"""Generate evidence tables and a Chinese delivery report from completed records."""
from pathlib import Path
import csv,json,hashlib,xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports'
def read(p):return json.loads(Path(p).read_text())
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |',*['| '+' | '.join(map(str,row))+' |' for row in rows]])

def main():
    active=read(OUT/'ACTIVE_FLOW.json');case=ROOT/active['case']
    sections=read(OUT/'internal_sections/INTERNAL_SECTION_CHECK.json')
    divergence=read(OUT/'internal_sections/DIVERGENCE_THEOREM_CHECK.json')
    streamlines=read(ROOT/'visualization/data/PREPARATION.json')
    assert sections['source_flow_sha256']==divergence['source_flow_sha256']==streamlines['source_flow_sha256']==active['flow_sha256']
    assert divergence['maximum_closure_error_fraction_inlet']<1e-10
    section_table=table(['边界','向内距离 mm','截面流量 μL/min','与端口差/该端口流量 %','与端口差/总入口流量 %'],
        [[r['boundary'],r['inward_offset_mm'],f"{r['flow_magnitude_uL_min']:.3f}",f"{100*r['relative_deviation_from_port']:.4f}",f"{100*r['deviation_as_fraction_total_inlet']:.4f}"] for r in sections['rows']])
    execution=read(case/'reports/execution.json');cal=ROOT/'cases/calibration_gpu8_production'
    calibration=read(cal/'postprocess/INDEPENDENT_CHECK.json')
    checks=read(ROOT/'microbubble/LOCAL_VERIFICATION.json')
    flow_integrity=read(OUT/'LOCAL_FLOW_VERIFICATION.json')
    assert flow_integrity['PASS'] and flow_integrity['flow_sha256']==active['flow_sha256']
    mb=read(ROOT/'microbubble/data/final_summary.json')
    assert active['status']=='PASS' and execution['status']=='PASS' and checks['PASS'] and mb['numerical_gate']
    geometry=read(ROOT/'inputs/geometry.json');mesh=read(OUT/'mesh_validity.json');quality=read(OUT/'mesh_quality.json');qc=read(OUT/'geometry_qc.json')
    backend=read(ROOT/'gpu_solver_fix/backend_equivalence.json')
    assert backend['passed']
    parity_rows=[dict(step=r['step'],**r['errors']) for r in backend['comparisons']]
    with (OUT/'CPU_GPU_backend_equivalence.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(parity_rows[0]));writer.writeheader();writer.writerows(parity_rows)
    eq=ET.parse(case/'run/solver.xml').find('.//Add_equation')
    prescribed={bc.get('name'):float(bc.findtext('Value')) for bc in eq.findall('Add_BC')}
    m=active['measurements'];w=active['raw_WSS']
    rows=[]
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        inlet=role=='INLET';q=m['Q_in_m3_s'] if inlet else m['outlet_flows_m3_s'][role]
        rows.append(dict(boundary=role,Q_m3_s=q,Q_uL_min=q*6e10,flow_percent=100 if inlet else m['outlet_fractions'][role]*100,
                         prescribed_pressure_Pa='' if inlet else prescribed[role],area_mean_pressure_Pa=m['area_average_pressure_pa'][role],
                         microbubble_count=1500 if inlet else mb['outlets']['O'+str(int(role[-2:]))]))
    with (OUT/'final_boundary_and_microbubble_summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    boundary_table=table(['边界','流量 μL/min','分流 %','设定出口压力 Pa','实测面积平均压力 Pa','微泡数'],
                        [[r['boundary'],f"{r['Q_uL_min']:.6f}",f"{r['flow_percent']:.6f}",('未指定' if r['prescribed_pressure_Pa']=='' else f"{r['prescribed_pressure_Pa']:.6f}"),f"{r['area_mean_pressure_Pa']:.6f}",r['microbubble_count']] for r in rows])
    port_table=table(['端口','实际边界编号','源打印坐标中心 mm'],[[k,v['face_id'],', '.join(f'{x:.6f}' for x in v['center_source_mm'])] for k,v in geometry['ports'].items()])
    case_rows=[]
    for folder in sorted((ROOT/'cases').iterdir()):
        if not folder.is_dir():continue
        marker=folder/'reports/execution.json';interruption=folder/'reports/execution_interruption.json'
        run=read(marker) if marker.exists() else {}
        interrupted=read(interruption) if interruption.exists() else {}
        purpose=read(folder/'purpose.json') if (folder/'purpose.json').exists() else {}
        case_rows.append(dict(case=folder.name,active_physical_flow=folder==case,
                              execution_status=run.get('status',interrupted.get('status','NOT_COMPLETED_OR_NOT_STARTED')),
                              independent_check=(read(folder/'postprocess/INDEPENDENT_CHECK.json')['status'] if (folder/'postprocess/INDEPENDENT_CHECK.json').exists() else 'NOT_AVAILABLE'),
                              backend_diagnostic_only=folder.name in backend['records'],
                              final_step=run.get('final_step',''),mode=purpose.get('mode',''),
                              reason=json.dumps(interrupted or run.get('health_failures',[]),ensure_ascii=False)))
    with (OUT/'case_registry.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(case_rows[0]));writer.writeheader();writer.writerows(case_rows)
    resources=[]
    for folder in [cal,case]:
        e=read(folder/'reports/execution.json');samples=np.array([[float(x.strip()) for x in r['gpu'].split(',')] for r in e['GPU_samples']])
        resources.append(dict(case=folder.name,elapsed_s=e['wall_time_s'],MPI_ranks=e['MPI_ranks'],device_gpu_utilization_mean_percent=float(samples[:,0].mean()),
                              device_gpu_utilization_max_percent=float(samples[:,0].max()),device_gpu_memory_max_MiB=float(samples[:,1].max()),peak_tree_RSS_MiB=e['peak_tree_RSS_MiB_sampled'],
                              linear_solves=len(e['history']['linear_solves']),failed_linear_solves=e['history']['failed_linear_solves']))
    with (OUT/'resource_usage.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(resources[0]));writer.writeheader();writer.writerows(resources)
    convergence_rows=[]
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for folder in [cal,case]:
        e=read(folder/'reports/execution.json');label='Equal-flow calibration' if folder==cal else 'Final pressure outlets'
        history=e['history']['linear_solves'];last={r['step']:r for r in history};v=list(last.values())
        axes[0,0].semilogy([r['step'] for r in v],[min(r['nonlinear_Ri_over_R0'],r['nonlinear_Ri_over_R1']) for r in v],label=label)
        intervals=e['intervals'];axes[0,1].semilogy([r.get('step',j) for j,r in enumerate(intervals)],[r['E_u'] for r in intervals],label=label)
        states=e['states']
        for s in states:
            convergence_rows.append(dict(case=folder.name,step=s['step'],time_s=s['time_s'],mass_error=s['epsilon_mass'],inlet_error=s['epsilon_Q'],velocity_max_m_s=s['velocity_max_m_s'],
                                         **{r+'_percent':x*100 for r,x in s['outlet_fractions'].items()}))
        axes[1,0].semilogy([s['step'] for s in states],[max(s['epsilon_mass'],1e-17) for s in states],label=label)
        if folder==case:
            for role in ['OUTLET_01','OUTLET_02','OUTLET_03']:axes[1,1].plot([s['step'] for s in states],[s['outlet_fractions'][role]*100 for s in states],label=role)
    axes[0,0].axhline(1e-10,color='k',ls='--',lw=.8);axes[0,1].axhline(1e-5,color='k',ls='--',lw=.8);axes[1,0].axhline(1e-6,color='k',ls='--',lw=.8);axes[1,1].axhline(100/3,color='k',ls='--',lw=.8)
    for ax,title,y in zip(axes.flat,['Completed-step nonlinear residual','Saved-interval velocity change','Boundary mass balance','Final pressure-outlet flow split'],['Relative residual','Relative volume L2 change','Absolute imbalance / inlet target','Outlet flow (%)']):
        ax.set(title=title,xlabel='Step' if ax is not axes[0,1] else 'Saved interval index',ylabel=y);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.savefig(OUT/'convergence_and_flow_split.png',dpi=300);fig.savefig(OUT/'convergence_and_flow_split.pdf');plt.close(fig)
    with (OUT/'convergence_summary.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(convergence_rows[0]));writer.writeheader();writer.writerows(convergence_rows)
    report=f'''# BraVa candidate 0 / 15：18 mL/min 流场与微泡交付报告

## 实际完成的内容

使用用户选择的 BALANCED 四端口血管芯，包含人工接管、排除盒体。candidate 0 和 15 是同一通道的两个刚性打印姿态；当前无重力的模型只求解一次，并正确旋转显示位置和矢量。入口固定 18,000 μL/min。正式压力出口的最大分流偏差为 {active['max_split_deviation_percentage_points']:.6f} 个百分点（相对各 33.333333%）。已生成两种姿态的流场动画，以及 1500 条、dt=0.5 ms 的有限尺寸微泡轨迹和英文动画。

同时确认了一项精度限制：人工接管内部截面流量相对各自端口存在约2%–3%的差异，已由独立散度积分核对。因此“边界分流均匀”不等于“局部速度和微泡输运已充分验证”；详细证据及使用范围如下。

## 几何、单位与模型

- 来源：`vascular_printing/开发与使用指南.md` 中 BG001 / RMCA 的已验收四端口制造模型。完整源路径、SHA256 和两姿态矩阵见 `inputs/geometry.json`。
- 源坐标 mm；求解坐标 m。转换：`x_SI = (x_source_mm - origin_source_mm) × 0.001`；源原点为 {geometry['origin_source_mm']} mm。显示坐标按两姿态矩阵转换回 mm。
- 官方 SimVascular/TetGen 网格：{mesh['vertices']} 节点、{mesh['tetra']} 四面体；翻转 {mesh['inverted']}、退化 {mesh['degenerate']}；体网格连通且封闭。最小 SICN={quality['q_min']:.6g}，P5={quality['P5']:.6g}。
- 网格相对源闭合体积差 {qc['volume_relative_error']*100:.6f}%；双向采样壁面距离 P95={qc['wall']['P95_m']*1e6:.3f} μm、最大 {qc['wall']['max_m']*1e6:.3f} μm。名义尺寸 0.2 mm，无专门棱柱边界层。
- svMultiPhysics：不可压缩牛顿流体、刚壁无滑移，P1/P1 稳定化四面体，包含惯性；密度 1056 kg/m³，动力黏度 0.00345312 Pa·s，运动黏度 3.27×10⁻⁶ m²/s。未加入重力、RBC 耦合、血管壁变形或非牛顿本构。
- 这是工程打印通道的血液参数模型，并非已测得的实际灌注液参数，也不是未经制造修正的原生生理血管。

{port_table}

端口来自制造身份映射，不按数字大小推断。源坐标为定位参考，求解实际端口面积和中心见 `reports/geometry_qc.json`。

## 边界与实际分流

先用两个出口的定流量与一个参考压力完成等流量校准，读取各出口面积平均压力，减去共同最小值作为统一参考平移；再将三个出口全部改为压力牵引边界，重新求解并测量真实分流。校准阶段的强制等流量不作为正式自然分流证据。入口为归一化平坦速度轮廓、边缘无滑移，总流量固定；壁面速度为零。

{boundary_table}

“设定压力”是出口自然牵引条件的参数，不保证端口节点压力处处相等，因此与面积平均压力分别列出。压力以 Pa 记录；共同参考平移不会改变出口间压差，三个出口的差异是实现此工程分流目标的设计量，不宣称是人体实测下游压力。

## 求解与独立检查

- CFD 隐式稳态推进 dt=0.01 s；微泡 dt=0.0005 s，二者不能混同。正式解保存于第 {active['final_step']} 步（推进时间 {active['physical_time_s']:.3f} s）。
- 线性 FGMRES 和非线性残差要求均保留 1e-10；连续五个保存间隔满足速度变化≤1e-5、流量变化≤1e-6；入口与整体质量误差≤1e-6。每步最多24次非线性迭代；增加预算不等于放宽标准。
- 正式入口误差 {m['epsilon_Q']:.6g}，质量误差 {m['epsilon_mass']:.6g}；壁面最大速度 {m['wall_velocity_max_m_s']:.6g} m/s。
- 最后两个保存场的相对变化：速度 {active['last_interval_changes']['velocity_relative_L2']:.6g}，压力 {active['last_interval_changes']['pressure_relative_L2']:.6g}，壁面 WSS 面积 L2 {active['last_interval_changes']['wss_area_relative_L2']:.6g}。
- 节点坐标/顺序与网格完全一致，四面体对应同一行；求解器局部顶点置换不用于重新定义网格。后处理使用原始体网格、真实节点速度和相邻四面体梯度。
- 正式流场 SHA256：`{active['flow_sha256']}`。实际 XML、PETSc 选项、网格哈希、原始日志和独立检查均保存在 `{active['case']}/`。

![收敛与分流](convergence_and_flow_split.png)

### 内部截面检查：已确认的局部守恒限制

边界流量平衡接近机器精度，**并不意味着离散速度场在内部逐点无散度**。在四条已核对的人工直接管内，分别取距端口2、5、10 mm的平面，选择对应连通管段，对真实P1速度按三角形精确积分：

{section_table}

最大内部截面偏差为总入口流量的 {100*sections['maximum_deviation_fraction_total_inlet']:.6f}%。另将每个端口与2 mm截面之间的体单元裁成薄段，对P1单元的常值 `div(u)` 独立体积积分；所得值与 `Q_port - Q_section` 一致，最大闭合差/总入口为 {divergence['maximum_closure_error_fraction_inlet']:.6g}。因此这些差异属于当前解析出的P1速度场的局部守恒缺陷，不能归咎于切面积分或色标。稳定化P1/P1方法弱施加不可压缩约束；本项没有证明缺陷的唯一成因，也没有运行加密对照来量化真实离散误差。

原场未重标定、未做速度修补。当前结果可展示该离散模型下的流动和微泡轨迹，边界分流接近均匀；但不能将近壁停留、内部路径比例或局部剪切数值作为已经充分验证的物理预测。截面VTP、数值CSV和散度独立核验见 `internal_sections/`；复现脚本为 `scripts/check_internal_sections.py` 与 `scripts/check_internal_divergence.py`。

![内部截面流量及独立散度核验](internal_sections/internal_flux_and_divergence.png)

## WSS 与可视化

从实际 FEM 节点速度在每个 P1 四面体内计算完整梯度 G，然后用相邻体单元的梯度和单位外法向 n：

`t_mu = mu (G + G^T) n`，`tau_w = t_mu - dot(t_mu,n) n`，`WSS_raw_Pa = norm(tau_w)`。

仅对固体壁面计算，入口出口截面不计入壁面。原始 WSS 是壁面三角形值；另存面积加权节点显示值，用于复用既有连续色彩图。两种数据分别提供动画，不把节点显示平滑当作原始数值更准确。

原始壁面 WSS：面积平均 {w['area_mean_Pa']:.6f} Pa；P5/P50/P95 = {w['P5_Pa']:.6f}/{w['P50_Pa']:.6f}/{w['P95_Pa']:.6f} Pa；范围 [{w['min_Pa']:.6f}, {w['max_Pa']:.6f}] Pa。单位法向投影后的最大绝对法向残量 {w['maximum_abs_normal_traction_Pa']:.6g} Pa。

精确零 WSS 面片 {w['exact_zero_faces']} 个，占壁面积 {100*w['exact_zero_area_fraction']:.6f}%；其中 {w['exact_zero_faces_with_all_wall_tetra_nodes']} 个的相邻四面体四节点全部属于无滑移壁面。对这些 P1 单元，四个节点速度均为零，梯度和原始 WSS 因离散空间而为零，不能自动解释为真实流动滞止。它们保留在原始场中，没有填值、删除或把显示插值作为修复。

复用原来的黑色背景、Turbo 配色、放大固定 Z 轴旋转、旁置箭头和渐变标签；色标数值根据新数据设置，坐标改为 mm。每姿态各5段 1920×1080、24 fps、18 s 动画：流线、局部速度矢量、压力、显示 WSS、原始面片 WSS。流线只是稳态流线，不是微泡轨迹；其数量不代表分流比。全旋转使用同一色标和固定缩放，并保存4K PNG/PDF。

流线使用192个入口通量加权候选点，真实终态计数为 `{streamlines['candidate_outcomes']}`；其中未达出口的55条达到30 s跟踪时限，记录全部保留，没有人为接到出口。动画使用最先自然完成的96条路径，因而不能代表全体入口样本的停留或完成比例。完成路径中抽查9条，将最大步长25 μm减半至12.5 μm，最大路径差 {max(r['maximum_path_difference_m'] for r in streamlines['refinement_comparisons'])*1e6:.6f} μm且出口身份不变；这只是抽样路径的积分检查，不证明所有近壁路径已收敛，更不能证明停留的物理原因。完整记录见 `visualization/data/`。

## CPU/GPU 实际使用与修复证据

服务器 RTX 4090，CPU 配额约7.68核；采用8个 MPI 进程负责装配与局部 LU，CUDA 稀疏矩阵/向量执行 Krylov 运算（PETSc 3.25.5）。不将 CPU 部分称为 GPU 计算。按5秒采样的资源统计见 `resource_usage.csv`；GPU指标是全设备采样，可能包含本任务同期几何预览，不能当作独立 CFD 内核的专属占用。稀疏混合求解不一定持续达到满 GPU 利用率。

首次多进程 CUDA 测试暴露了设备解向量到主机 ghost-local 视图未同步的问题：线性求解可收敛，但导出的修正量错误。只在独立求解器副本中增加 VecGetArrayRead/VecRestoreArrayRead 同步，未修改原求解器、FEM 公式、几何或边界。随后使用相同真实网格、参数和两步时间推进，对照原 CPU8 与修正 GPU8 的完整速度、压力及 WSS，误差见 `gpu_solver_fix/backend_equivalence.json` 和比较表；容许误差1e-8。此验证证明这次后端对照一致，不代替网格或稳态验证，也不推断历史其他配置全部存在同样问题。

所有失败、停止及仅用于后端核对的短算例保留在 `cases/`，不进入正式数据接口。逐目录状态见 `case_registry.csv`；两步后端诊断不满足稳态门槛，其后端比较通过与完整 CFD 通过是两件事。构建清单、单文件补丁、源文件前后哈希和独立可执行文件哈希均保留。

## 有限尺寸微泡

- 原微泡动力学源文件保持字节不变；核对 {checks['protected_source_files_verified']} 个保护文件。只适配新几何、流场身份和入口采样，名义 dt=0.5 ms，必要时内部细分。
- 取入口通量加权、原 SonoVue 条件粒径分布的前1500个有效出生事件，不按出口结果筛选。保留原单向冻结流场、有限尺寸、近壁水动力及事件安全机制，无微泡–微泡或 RBC 相互作用。
- 实际终态：`{mb['statuses']}`；三个出口计数见上表。全部事件完成标记、轨迹文件、审计和源身份哈希通过检查。安全计数：`{mb['safety_totals']}`。
- CPU8 并行积分。实际累计 CPU 时间 {mb['CPU_seconds']:.1f} s，等效占用 {mb['CPU_core_equivalents']:.3f} 核。编译几何核与原 Python 核先进行实际三角形查询、完整单条轨迹及审计的一致性验证；CUDA FP64 单独核对全体梯度及轨迹统计，NVIDIA EGL 渲染，CPU 编码 MP4。
- 数值样本到达率1000/s只是生成有限轨迹集合，不是人体剂量估计。视频按轨迹年龄对齐；不是1500个微泡同时注射的物理录像。英文文字，沿用既有风格，主要元素放大1.30倍。
- 两姿态共用同一组物理轨迹，仅作对应刚性显示变换；原始 SI 轨迹不覆盖。退出计数存在有限样本、有限尺寸和积分离散效应，不必严格等于流量比例。

## 可用范围与限制

已确认当前输入身份、求解残差、稳态变化、边界质量守恒、无滑移、CPU/GPU 后端一致性、原始 WSS 链路、流场与轨迹来源以及媒体完整性。本轮只有一档人脑打印通道网格，未完成空间/时间步收敛研究；局部 WSS 尤其是接管连接、分叉和尖角附近不能宣称网格无关。打印实物尺寸、灌注液、端口阻力和实验分流尚未实测校准。模型没有重力，所以姿态复用适用于本模型；不能推出真实实验改变姿态一定没有影响。

本次额外核验已确认内部截面约百分之几的局部速度通量缺陷；边界残差通过不能覆盖这一限制。微泡的数值安全检查验证的是此冻结离散场中的积分和几何约束，没有证明流场局部精度或停滞的生物物理真实性。

微泡粒径0.75–4 μm，而本次名义网格尺寸200 μm、壁面相对源几何的P95距离约 {qc['wall']['P95_m']*1e6:.3f} μm；这些是不同尺度。防穿壁证书针对当前三角形壁面成立，并不证明微米尺度实际打印壁面的接触位置准确。尤其不能将多面片接触约束导致的停留直接解释为真实血管或打印实验中的微泡黏附。

## 打开与复现

- 本地总入口：`../OPEN_RESULTS.html`。
- 流场：`../visualization/candidate_0/OPEN_RESULTS.html`、`../visualization/candidate_15/OPEN_RESULTS.html`。
- 微泡：`../microbubble/candidate_0/OPEN_RESULTS.html`、`../microbubble/candidate_15/OPEN_RESULTS.html`；详细说明 `../microbubble/MICROBUBBLE_RESULTS_ZH.md`。
- 原始求解、冻结流场和边界表：`../{active['case']}/`。
- 两姿态各自的 SI 流场副本：`../{active['case']}/frozen_flow/candidate_0/flow.vtu` 与 `candidate_15/flow.vtu`；保持标量不变、正确旋转速度与壁面矢量，附变换往返误差检查。
- 服务器：`vast4090:/workspace/brava_flow_roi_18mlmin_20260928/`。使用独立 supervisor；不控制其他项目服务。
- 复现命令与依赖见 `../REPRODUCE.md`。本报告由实际机器记录自动生成，不能用作新算例跳过检查的凭证。
'''
    (OUT/'BRAVA_FLOW_MICROBUBBLE_REPORT_ZH.md').write_text(report)
    html='<!doctype html><meta charset="utf-8"><title>BraVa flow and microbubbles</title><style>body{max-width:1200px;margin:40px auto;background:#101824;color:#edf2fa;font:20px sans-serif}a{color:#50d9ef}video{width:100%}</style><h1>BraVa BG001 RMCA | 18 mL/min</h1><p>Same four-port core, two rigid print poses. 1500 microbubble tracks, dt = 0.5 ms.</p><a href="reports/BRAVA_FLOW_MICROBUBBLE_REPORT_ZH.md">中文结果报告</a>'
    for pose in ['0','15']:
        html+=f'<h2>Candidate {pose}</h2><p><a href="visualization/candidate_{pose}/OPEN_RESULTS.html">Flow movies and figures</a> · <a href="microbubble/candidate_{pose}/OPEN_RESULTS.html">Microbubble movie and statistics</a></p><video controls preload="metadata" src="microbubble/candidate_{pose}/animations/microbubble_age_aligned.mp4"></video>'
    (ROOT/'OPEN_RESULTS.html').write_text(html)
    readme=ROOT/'README.md'
    text=readme.read_text().replace('新微泡批次，待正式流场及动画完成后执行。','已完成的1500条微泡批次、两姿态动画与完整性检查。')
    marker='本轮交付已完成：'
    if marker not in text:
        text=text.replace('本目录是独立新算例，不修改原打印件、小鼠流场或微泡核心代码。',
            '本目录是独立新算例，不修改原打印件、小鼠流场或微泡核心代码。\n\n本轮交付已完成：[总预览](OPEN_RESULTS.html) · [中文结果报告](reports/BRAVA_FLOW_MICROBUBBLE_REPORT_ZH.md) · [边界及微泡汇总](reports/final_boundary_and_microbubble_summary.csv)。唯一正式场由 `reports/ACTIVE_FLOW.json` 指定。')
    readme.write_text(text)
    print('REPORT_COMPLETE',OUT/'BRAVA_FLOW_MICROBUBBLE_REPORT_ZH.md')

if __name__=='__main__':main()
