"""Snapshot-only navigation and factual state summary; no simulation calls."""
from pathlib import Path
import hashlib,json
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260929';R=D/'brava_flow_roi_18mlmin'
active=json.loads((R/'reports/ACTIVE_FLOW.json').read_text());stop=json.loads((R/'microbubble/data/USER_STOP.json').read_text())
m=active['measurements'];case=R/active['case'];p={'OUTLET_01':556.9706555488448,'OUTLET_02':0.,'OUTLET_03':17.13742715356202}
import xml.etree.ElementTree as ET
eq=ET.parse(case/'run/solver.xml').find('.//Add_equation')
p={bc.get('name'):float(bc.findtext('Value')) for bc in eq.findall('Add_BC') if bc.get('name').startswith('OUTLET')}
table='\n'.join(f"| {role} | {q*6e10:.6f} | {100*m['outlet_fractions'][role]:.6f} | {p[role]:.6f} |" for role,q in m['outlet_flows_m3_s'].items())
text=f'''# BraVa 最新结果与停止状态

本次 GitHub 同步只读取、保存已有结果，没有运行 CFD 或微泡积分。

## 已完成的流场

几何是用户选择的 BG001 RMCA BALANCED 四端口血管芯，含人工接管、无盒体。candidate 0 与15是同一几何的刚性打印姿态；本模型不含重力，复用同一物理解并旋转位置与矢量。入口为18,000 μL/min（18 mL/min），CFD隐式推进dt=0.01 s。

采用svMultiPhysics P1/P1稳定化牛顿流体模型，ρ=1056 kg/m³、μ=0.00345312 Pa·s，刚性无滑移壁面。网格135,852节点、692,313四面体，未设专门边界层。先进行等流量校准，再以三个压力出口重新求解；表内是正式压力出口求解的实际分流。

| 边界 | 流量 μL/min | 分流 % | 设定牵引压力 Pa |
|---|---:|---:|---:|
{table}

正式第{active['final_step']}步达到已设定稳态标准，独立检查通过。流场SHA256：`{active['flow_sha256']}`。

- 实际算例：[balanced_pressure_final](cases/balanced_pressure_final/)。
- [candidate 0 流场动画](visualization/candidate_0/OPEN_RESULTS.html) / [candidate 15 流场动画](visualization/candidate_15/OPEN_RESULTS.html)。每姿态各5段，1920×1080、24 fps、432帧，含流线、局部矢量、压力、显示WSS和原始面片WSS；同时保留4K PNG/PDF。
- 原始WSS均值{active['raw_WSS']['area_mean_Pa']:.6f} Pa，P5/P50/P95={active['raw_WSS']['P5_Pa']:.6f}/{active['raw_WSS']['P50_Pa']:.6f}/{active['raw_WSS']['P95_Pa']:.6f} Pa。
- [CPU/GPU等价核验](gpu_solver_fix/backend_equivalence.json)：独立副本修复设备解到主机视图同步后，同网格两步场及WSS相对差约1e-13量级。GPU负责CUDA稀疏Krylov，CPU8负责装配与局部LU；这是后端一致性检查，不是网格独立性验证。

## 已确认的限制

内部人工接管截面与各自端口流量存在约2%–3%的差异；独立体积散度积分与截面通量差一致。边界整体守恒不能代替局部精度证明。详见 [截面CSV](reports/internal_sections/section_flux.csv)、[独立散度CSV](reports/internal_sections/divergence_theorem.csv) 和 [对比图](reports/internal_sections/internal_flux_and_divergence.png)。本轮只有一档人脑网格，不能宣称网格无关；近壁停留和局部WSS不应视为已充分验证的物理预测。

## 用户已停止微泡批次

计划1500条、dt=0.5 ms，保留原SonoVue条件粒径分布。入口1500个样本、CUDA梯度复核，以及原/编译几何核的单条12 s轨迹和审计一致性检查已经完成。该单条为长停留时间截断，不能当作1500条批次结果。

正式批次完成条数 **{len(stop['complete_track_folders'])}**，保留 **{len(stop['incomplete_track_folders'])}** 条未完成轨迹的记录。任务及独立缓存验证均按用户要求停止：`CANCELLED_BY_USER`；自动重启关闭，主流程检查`USER_STOP.json`后退出。没有本批微泡动画，没有完整批次汇总，没有最终完整交付报告。`write_results_report.py`是待满足交付条件的脚本，不代表报告已生成。

新增`query_cache.py`是未完成验证的草稿；停止时批量进程使用的是缓存修改前已经载入的入口代码。不能把此草稿写成已通过数值验证或已用于正式结果。

停止证据：[USER_STOP.json](microbubble/data/USER_STOP.json)、[说明](microbubble/data/CANCELLED_EXECUTION_NOTES.json)、[原始目录](microbubble/)。所有中间记录原样保存，不补造完成标记。

服务器：`vast4090:/workspace/brava_flow_roi_18mlmin_20260928/`。复现说明见[REPRODUCE.md](REPRODUCE.md)；其中计算命令仅供以后明确授权恢复时参考，本次同步不得执行。GPU输入缓存因体积排除，恢复现有文件的方法见仓库同步[排除说明](../sync_metadata/current_20260929/EXCLUSIONS.md)。
'''
(R/'CURRENT_RESULTS_ZH.md').write_text(text)
old=(D/'README.md').read_text();(M/'README_before_sync.md').write_text(old)
(D/'README.md').write_text('''# 血管几何、3D流场、微泡与RBC工作流

2026-09-29增量同步至 `sync/roi-flow-mb1500-wss-stop-analysis-20260928`。原有小鼠流程完整保留；最新新增BraVa人脑打印通道的18 mL/min流场及其动画。**BraVa微泡批次已按用户要求停止，未完成1500条，未生成微泡动画。**

| 内容 | 入口 |
|---|---|
| BraVa最新流场、数值限制与微泡停止状态 | [当前结果说明](brava_flow_roi_18mlmin/CURRENT_RESULTS_ZH.md) |
| BraVa真实网格、压力边界及求解日志 | [算例](brava_flow_roi_18mlmin/cases/balanced_pressure_final/) |
| candidate 0 / 15 的流场动画与4K图 | [candidate 0](brava_flow_roi_18mlmin/visualization/candidate_0/OPEN_RESULTS.html) · [candidate 15](brava_flow_roi_18mlmin/visualization/candidate_15/OPEN_RESULTS.html) |
| BraVa几何、ROI与四端口芯的生成来源 | [开发与使用指南](vascular_printing/开发与使用指南.md) |
| 小鼠血管A、ROI、1D/0D及边界设计 | [ulm_3D_vascular](ulm_3D_vascular/) |
| 已完成的小鼠ROI-only 3D流场 | [算例](ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1/) |
| 已完成的小鼠1500条微泡、四视角及停止原因分析 | [动画](ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/OPEN_RESULTS.html) · [停止分析](ulm_particle_formal_p9a5/particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/stopping_analysis/MICROBUBBLE_STOP_ANALYSIS_ZH.md) |
| 网格、WSS审计、会议疑问图与旋转绘图代码 | [formal_3D_flow_solver](formal_3D_flow_solver/README.md) |
| RBC代码、轨迹及独立模型展示 | [保留结果索引](ulm_particle_formal_p9a5/CURRENT_RESULTS.md) |
| 当前服务器路径 | [路径表](CURRENT_SERVER_PATHS.md) |
| 此次同步范围、排除项、哈希核验 | [同步报告](sync_metadata/current_20260929/SYNC_REPORT_ZH.md) |

HTML需下载后在浏览器中打开。人脑与小鼠输入、压力、粒径样本和完成状态分别记录，不能混用；RBC结果不表示已完成人脑新场的RBC耦合计算。

本轮只收集和核验，不运行CFD或轨迹。未上传运行环境、可重建二进制和不重要的大型缓存；路径、大小、哈希与恢复方式见[排除说明](sync_metadata/current_20260929/EXCLUSIONS.md)。上游缺失的LFS载荷仅记录指针，不能当作实际数据。部分入口仍绑定原WSL/服务器路径和科学身份契约，快照不是可在任意机器直接启动的安装包。

BraVa流场的边界分流接近各三分之一；已确认内部截面存在约百分之几的局部速度通量缺陷。单网格、边界守恒与残差通过均不证明局部WSS或微泡输运精度。旧细档真实血管CFD、出口压力敏感性及本次BraVa微泡均保持各自的用户取消状态。

上次快照详细入口保存在[2026-09-28说明](sync_metadata/current_20260929/README_before_sync.md)。本次以`ca0ae424d01717f6a231bbea91405cbea8f095a8`为父提交，正常追加到指定分支，不改写历史。
''')
transforms=[]
for name in ['ACTIVE_VASCULAR_WORKFLOW.md','CURRENT_SERVER_PATHS.md']:
    pth=D/name;before=pth.read_text();(M/('original_'+name)).write_text(before)
    if name=='ACTIVE_VASCULAR_WORKFLOW.md':
        prefix='# 2026-09-29 新增 BraVa 独立流程\n\n最新人脑打印通道结果：[BraVa状态与结果](brava_flow_roi_18mlmin/CURRENT_RESULTS_ZH.md)。18 mL/min CFD及两姿态流场动画已完成；微泡批次已按用户要求停止，不能与下列已完成的小鼠1500条混同。\n\n以下保留小鼠工作流：\n\n'
        before=before.replace('temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/','ulm_flow_mean_2p0_mmps/')
    else:
        prefix='# 2026-09-29 新增 BraVa 服务器路径\n\n| 内容 | 路径 |\n|---|---|\n| BraVa当前独立算例 | `vast4090:/workspace/brava_flow_roi_18mlmin_20260928/` |\n| 正式压力出口CFD | 该目录下 `cases/balanced_pressure_final/` |\n| candidate 0 / 15流场动画 | 该目录下 `visualization/candidate_0/`、`candidate_15/` |\n| 已停止的微泡批次 | 该目录下 `microbubble/`；`data/USER_STOP.json`为取消标记 |\n| 隔离的GPU同步修复求解器 | 该目录下 `gpu_solver_fix/svmultiphysics` |\n\n微泡服务为STOPPED、autostart=false、autorestart=false；本次同步不恢复计算。下列小鼠/RBC路径为独立保留流程：\n\n'
    pth.write_text(prefix+before)
    transforms.append(dict(path=name,original_copy=str((M/('original_'+name)).relative_to(D)),snapshot_sha256=hashlib.sha256(pth.read_bytes()).hexdigest(),reason='Navigation/status only; scientific source/data unchanged'))
(M/'documentation_transforms.json').write_text(json.dumps(transforms,ensure_ascii=False,indent=2)+'\n')
print('Snapshot documentation written')
