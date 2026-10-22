"""Record measured outcomes and a minimal append-only research progress update."""
import difflib
import re
import subprocess
import numpy as np
from .common import *

def main():
    s=read(BASE/'analysis/SUMMARY.json'); verify(); groups=s['groups']; runs=s['runs']; c=s['campaign']
    resources=read(next(BASE.glob('*/resources.json')))
    probe=[next(r for r in runs.values() if r['N']==n and r['action']=='PAIR' and abs(r['requested_s']-1e-4)<1e-12) for n in [32,64,128]]
    lines=['# F10Z 补充：64、128 泡有限规模仿真', '',
        '在既有 F10Z 完整结果上新增泡数试验，原声场、APECSS、逐泡声力、润滑、Sequential QP 和壁面查询源码全部保持原样。新增的是可变泡数的试验输入、轨迹保存和统计；没有用团级代理代替逐泡积分。', '',
        '保持 `EXPLORATORY`、`HOMOGENEOUS_WATER_CONTROL_PROXY`、`PHYSICS_UNCALIBRATED`、`ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED`、`training_ready=false`。', '',
        '## 预登记设置', '',
        '- 泡径：半径 1.78 µm；频率：2.5 MHz；背景流严格为零。',
        '- 复用原四路均匀水圆片声场缓存。X−、X+ 各 1036.201658 Pa-equivalent，相位分别 0°、210°；Y 两枚关闭。该幅值不是电压，也不是泡所在处声压。',
        '- N=32/64/128；随机种子 18002/18064/18128。初始最小中心距 5 µm，均为完全分散初态。',
        '- 初始采样立方体半宽按 N 的立方根增大，名义数密度固定为 6.85871×10¹⁴ m⁻³。每种泡数仅一个初态，布局及声场覆盖范围随之改变，不能据此单独归因于泡数或声称统计普适。',
        '- 原 dt=1 µs，对照 dt=0.5 µs；声场缓存 3/1.5 µm。先各做 0.1 ms 测速，再按预登记成本门槛选择 2/1/0.5 ms。总远程预算 1200 s，每条 240 s，没有提高数值容差。', '',
        '## 实际结果', '',
        '|泡数|实际窗口 ms|Rg 变化 µm|Rg 变化 %|连通分量|最大局部组泡数|COM 位移 µm|',
        '|---:|---:|---:|---:|---:|---:|---:|']
    for n in [32,64,128]:
        g=groups[str(n)]
        if 'result' not in g:
            lines.append(f'|{n}|NOT_RUN|—|—|—|—|—|');continue
        r=g['result'];co=r['components_initial_final'];la=r['largest_initial_final'];d=np.linalg.norm(r['displacement_m'])*1e6
        lines.append(f'|{n}|{r["completed_s"]*1e3:g}|{r["Rg_change_m"]*1e6:.6f}|{r["Rg_change_percent"]:.4f}|{co[0]} → {co[1]}|{la[0]} → {la[1]}|{d:.6f}|')
    lines+=['','32 泡 2 ms 一行复用已完成 F10Z 数据；本轮独立复算其前 0.1 ms，与历史逐泡位置严格核对。64、128 泡均为新的完整逐泡积分。所有位置与时间来自实际接受状态，未插值、填充缺失帧或缩放群内间距。','',
        '连通分量依据中心距 ≤2.1R 的几何邻近图；最大局部组不是物理融合团。Rg 减少表示去掉质心平移后的三维收缩，不能仅凭它认定全部微泡已形成稳定聚集团。XY 图中重叠的投影可以来自不同 Z 位置，不代表三维球体重叠。不同 N 的观察窗口不同，不能直接以终点收缩百分比排列其聚集能力。','',
        '## 收敛、实际压力和失败边界','']
    for n in [64,128]:
        g=groups[str(n)]
        if 'result' not in g:continue
        r=g['result'];off=g['off'];comps=s['comparisons'][str(n)]
        err=max((v['errors']['fixed_id_m'] for v in comps if 'errors' in v),default=float('nan'))
        rg=g['Rg_difference_bound_m'];pressure=r['pressure_range_Pa']
        lines.append(f'- {n} 泡：OFF 最大固定 ID 位移 {off["max_fixed_ID_motion_m"]:.3g} m；dt/cache 连续指标准入 {g["numeric_admitted"]}，事件拓扑对照 {g["topology_agrees"]}。最大逐泡位置对照差 {err*1e9:.6g} nm；最大 Rg 对照差 {(rg or 0)*1e9:.6g} nm。实际泡位置声压范围 {pressure[0]:.6g}–{pressure[1]:.6g} Pa。')
    lines += ['', '本次 64 泡轨迹出现非零接触反力，最大 3 个活跃接触约束、2 次 Sequential QP；128 泡在所跑 1 ms 内尚未激活硬接触约束，最小中心距约 3.56130 µm。因此 128 泡结果证明的是这一分散到邻近过程的可运行性，尚未验证 128 泡密集接触的求解性能。128 泡 2 ms 依据运行前登记的成本门槛未运行，并非数值求解失败。', '',
        '每步都保留原硬几何、力平衡、互补、接触修正和零墙修正门槛。非零声压必须在响应表 0.1–40000 Pa 内，无静默裁剪。未改变 1 nm 硬几何门槛，也未改变 17.8 nm 总重叠预算、178 nm 接触修正预算。详细逐步统计、比对误差与未运行原因见 analysis/SUMMARY.json 和 scaling_v1/DECISIONS.json。', '',
        '## 真实计算成本', '',
        f'本轮 {len(runs)} 条积分记录，PASS={sum(r["status"]=="PASS" for r in runs.values())}，FAIL={sum(r["status"]!="PASS" for r in runs.values())}；远端完整进程墙钟 {c["wall_s"]:.3f} s（{c["wall_s"]/60:.2f} min），逐条计算及保存墙钟合计 {sum(r["wall_s"] for r in runs.values()):.3f} s。此处不包含代码开发和本地绘图时间。', '',
        '|泡数|0.1 ms 测速每步 ms|相对 32 泡|', '|---:|---:|---:|']
    for r in probe: lines.append(f'|{r["N"]}|{r["mean_step_s"]*1e3:.3f}|{r["mean_step_s"]/probe[0]["mean_step_s"]:.2f}×|')
    lines+=['',f'监测 GPU 显存峰值 {max((r.get("gpu_MiB",0) for r in resources),default=0):.0f} MiB，GPU 利用率峰值 {max((r.get("gpu_percent",0) for r in resources),default=0):.0f}%；实际为 GPU 声力/壁面 + CPU 稠密阻力/接触求解。服务器上其他 CPU 项目保持运行，本次数据不是严格独占 CPU 的性能基准。','',
        '全粒子对数由 496 增至 2016、8128，约束矩阵列数也随 3N 增大。更多泡并不能简单按显存容量估算速度；接触出现后比初始分散短测更贵。未在本轮进行千泡或秒级外推。','',
        '## 产物与复现','',
        '- figures/Fig_particle_count_morphology.png / .pdf：真实半径的 XY 初始/终点投影，同一坐标比例。',
        '- figures/Fig_particle_count_response.png / .pdf：Rg、质心、分量比例和同窗口短测耗时。',
        '- scaling_v1/runs/：原始每个接受步逐泡坐标、真实时间、固定 ID、逐步数值指标、受力抽样、输入和耗时。',
        '- PREREGISTRATION.json、INITIAL_STATE_AUDIT.json、RUN_SCHEDULES.json：运行前输入、种子、控制、成本停止规则及哈希。',
        '- evidence/：Git 执行身份、传输验证、视觉检查及科研文档 before/after/diff。',
        '- 运行：`PYTHONPATH=src:. .venv/bin/python -m scripts.straight_tube_exp1.f10z_scaling.remote --label <unique_label>`；仅用于新的独立执行，不能覆盖已有标签。',
        '- 后处理：同样环境执行 `-m scripts.straight_tube_exp1.f10z_scaling.analyze`。', '',
        f'积分 Git 快照：`{c["identity"]["execution_commit"]}`；分支：`dev/f10z-bubble-count-scaling-20261018`。历史结果未覆盖，本轮未推送 GitHub。','',
        '## 下一步','',
        '优先对同一泡数增加独立初态复算，再研究保持同一物理范围的密度变化。若继续扩大到 256 泡，应先处理稠密接触矩阵的计算成本并保留原数值验收；不通过减少相互作用、放宽接触容差或增大时间步长来制造加速结论。先前未校准的近场声力、润滑、接触和真实 PDMS 声场限制仍然存在。','']
    (BASE/'REPORT.md').write_text('\n'.join(lines))
    progress=Path('/home/lzy/projects/microbubble_aggregation.md')
    backup=progress.with_name(progress.name+'.bak_F10Z_NSCALING_20261018')
    before=progress.read_text(); marker='【F10Z 泡数扩展，2026-10-18】'
    paragraph=marker+'本轮基于原逐泡引擎补做 64、128 泡，保持半径 1.78 µm、2.5 MHz、零背景流及原数值容差；初态随 N 的立方根扩大，维持近似相同数密度。'
    for n in [64,128]:
        g=groups[str(n)]
        if 'result' not in g: paragraph+=f'{n} 泡正式轨迹 NOT_RUN。';continue
        r=g['result'];paragraph+=f'{n} 泡实际 {r["completed_s"]*1e3:g} ms，Rg 变化 {r["Rg_change_percent"]:.4f}%，几何分量 {r["components_initial_final"][0]}→{r["components_initial_final"][1]}，最大局部组 {r["largest_initial_final"][1]} 泡；dt/cache 连续量准入 {g["numeric_admitted"]}，OFF 完全静止。'
    paragraph+=f'本轮远端总墙钟 {c["wall_s"]:.3f} s；较大泡数按预登记成本选择短窗口，没有冒充完整 2 ms 或秒级结果。每种 N 仅一套布局，不能声称纯泡数因果或物理融合。新结果和两组 PNG/PDF 位于原 F10Z 工程的 results/straight_tube_exp1/f10z_scaling，保留 EXPLORATORY、PHYSICS_UNCALIBRATED、training_ready=false。'
    evidence=BASE/'evidence';evidence.mkdir(exist_ok=True)
    if marker not in before:
        assert not backup.exists(), 'Do not overwrite progress backup'
        backup.write_bytes(progress.read_bytes())
        anchor='### 尾声：故事里还没讲完的部分'
        assert before.count(anchor)==1
        after=before.replace(anchor,paragraph+'\n\n'+anchor)
        assert re.findall(r'^#+ .+$',before,re.M)==re.findall(r'^#+ .+$',after,re.M)
        assert after.replace(paragraph+'\n\n','',1)==before
        progress.write_text(after)
        (evidence/'progress.before.md').write_text(before);(evidence/'progress.after.md').write_text(after)
        (evidence/'progress.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='before',tofile='after')))
        write(evidence/'PROGRESS_AUDIT.json',dict(backup=str(backup),before_sha256=sha(backup),after_sha256=sha(progress),
            headings_unchanged=True,old_content_preserved=True))
    else: assert paragraph in before, 'Do not silently rewrite previous scientific progress'
    write(evidence/'DELIVERY_IDENTITY.json',identity())
    print('Report and progress document saved')

if __name__=='__main__':main()
