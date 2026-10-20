"""Evidence-driven report and static local index, no physical run."""
from common import *
import shutil,html
def main():
 result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());q=json.loads((O/'SCIENTIFIC_OUTCOME.json').read_text())
 geo=json.loads((O/'IRREGULAR_REFERENCE_GEOMETRY.json').read_text());elastic=json.loads((O/'ELASTIC_CHECK.json').read_text());audit=json.loads((O/'TERMINAL_RESTART_AUDIT.json').read_text())
 old=q['statistics']['Y regular'];new=q['statistics']['Y V4']
 comparisons='\n'.join(f'| {label} | {old[key]:.6g} | {new[key]:.6g} |' for label,key in [('游离材料点数','free_nodes'),('游离分量数','free_components'),('最大游离分量节点数','largest_nodes'),('时程最大游离分量节点数','peak_largest_nodes'),('20–200 点目标块末态数','target_20_200'),('同一窗口历史最多目标块','peak_target_count'),('目标块质量 / kg','target_mass_kg'),('1–2 点细屑占初始质量','debris_1_2_mass_fraction')])
 lines=f"""# 原 Y 形血管血块上的 V4 不规则 LPS 应用

任务 `y_clot_irregular_current_method_v4`；RUN `{RUN}`。**{result['status']}**；实际接受到 **{result['accepted_time_s']*1e3:.8f} ms / 3.6 ms**。本轮只运行一条新物理轨迹，先到 2 ms，通过后从同一原生检查点继续；历史 Y 规则点阵只读复用。旋转视频仅改变相机，未增加物理工况。

## 迁移内容与守恒基准

继续原 Y 方法 v2 的工程构造三瓣壁附血块，原 1927 材料点、1561 锚点、原材料 ID、强弱相位与附着面保留。它不是病人几何重建。原 10 个强区的平移裁剪在此几何内留下 8 个区号；882 强相点、1045 弱相点。相位标签不因扰动重新采样，强相积分比例保持 45.7706279%。

本案例守恒 **0.240875 mm³ PD 积分体积 / 0.240875 mg 质量**。原连续几何约 0.271577 mm³ 是另一种口径；长方体 V4 的 1.015625 mm³ 不适用于本案例。

生成参考点云时，1515 个边界/界面材料点不施加初始扰动，只扰动 412 个内部点。正式动力学仍让全部 1927 材料点运动，仅固定 1561 个锚点。固定种子 20261019，最大扰动 12.5 µm；最小间距约 28.4943 µm。共享映射控制体积经六个辅助四面体积分，全部体积为正；没有事后缩放 Vi 来凑总量。原始参考键 157637 → 153986（本体 78618、界面 22041、锚层支架 53327）。新连接按 150.5 µm horizon 和原连续几何可见性构建，之后断键不重连。

边界为保持原形状而固定，故外表面仍显示规则排布；不能把图中的内部扰动解释成真实纤维网络。V4 同时改变内部位置、参考邻接和求积，不是仅改变点的画法。

## 材料、载荷及求解

复用已冻结非均匀体积原生变分 LPS：K=100000/3 Pa，G=940 Pa，ρ=1000 kg/m³，逐节点质量 ρVi、邻点 Vj 和重新计算的 mi。强弱相阈值和原 Y 界面/支架阈值保持。完整力来自同一离散势，κ 和 d 不可逆；碎块由存活原始材料图确定，不用空间聚类。

指定 +X 当前表面切向滑移代理，驻止表面标称峰值 1000 Pa；0–2 ms smoothstep 上升，2–2.4 ms 平台，2.4–2.8 ms 撤载，此后保留被动相对拖曳至 3.6 ms。μ=1 mPa·s，h_eff=25 µm，β=40 Pa·s/m。该指定外层速度上限 25 m/s **不是已求解声流速度**。载荷按当前位置、法向、当前面积与速度求积，附着背面在未脱离时排除；脱离后还需壁间隙与真实管腔检查。

L1 为存活支撑三角面；L2 为真实游离分量的映射控制体面，保留跨面全部原始桥键与遮挡检查。它仍是 **PARTIAL_CURRENT_SURFACE_APPROXIMATE_COUPLING**。壁接触复用原 Y 几何、20 µm 有效半径、50 N/m 刚度与 0.5 阻尼比，阻尼系数用实际 ρVi；壁力与流体力分开记录。

本轮没有新解 CFD、声场、微泡群、双向 FSI 或长期输运。材料与异质性均为 **EXPLORATORY_UNCALIBRATED**，载荷为 **MECHANICAL_SHEAR_SURROGATE**。

## 数值核验及计算范围

原 Y 同节点小仿射诊断中，体积响应相对变化 {elastic['relative_modulus_changes']['hydro']:.5g}，xy/xz/yz 剪切响应相对变化 {100*elastic['relative_modulus_changes']['xy']:.5f}% / {100*elastic['relative_modulus_changes']['xz']:.5f}% / {100*elastic['relative_modulus_changes']['yz']:.5f}%。这是离散一致性诊断，不是宏观模量或各向同性标定。

原生力/势与独立完整变分模型、能量梯度、面力合力/力矩/虚功、壁接触势梯度、非均匀质量重启续步检查通过。正式 BASE=0.25 µs，最小事件子步 0.03125 µs，预测 Δd≤0.1；不设 4 GiB / 1800 s 停机预算。

接受 {result['accepted_steps']} 步，预测拒绝 {result['rejected_predictions']} 次；最大归一化能量残差 **{100*result['max_relative_R']:.7f}%**，门槛 1%。壁中心最小净距 {q['minimum_wall_clearance_m']*1e6:.6g} µm。末态重启审计 **{audit['status']}**，审计检查点是否等于正式末态：{audit['matches_formal_accepted_end']}。全部真实失败诊断（若有）保留，最后接受状态与失败试步分开。

能量账本采用 Wdrive+Wdrag+Wwall=ΔU+ΔK+D+ΔUc−Jc+R。墙面势能是诊断列，已记壁外功时不重复加到右端；R 不并入耗散。单 rank、单线程，无 GPU；累计墙钟 {result['wall_s']/3600:.4f} h，进程峰值 RSS {result['peak_RSS_bytes']/2**30:.4f} GiB。

## 与历史 Y 规则点阵的实际结果

| 指标 | 历史 Y 规则点阵 | 当前 Y V4 |
|---|---:|---:|
{comparisons}

两者所用端态真实时间分别为 {old['time_s']*1e3:.8f} / {new['time_s']*1e3:.8f} ms。历史时程选择不晚于新轨迹接受末态的保存状态。碎块质量按各点 Vi 加总；20–200 是原始点数区间，不等同固定物理直径。游离图分量、节点标记和连续组织外壳也不等同。

逐成员集合的出生、分裂和末态右删失保存在 `FRAGMENT_SURVIVAL.json`。这是一条受约束点云样本，不能据末态块数宣称更真实、更稳定或已经验证超声碎栓机制。

## 可视化与数据入口

- [静态 HTML 入口](../../results/{TASK}/{RUN}/OPEN_RESULTS.html)
- `figures/FIG_01_Y_VESSEL_LOCATION.png`：原 Y 血管内实际位置。
- `figures/FIG_02_INITIAL_CLOT.png` 与 `FIG_02_FINAL_CLOT.png`：实际初态与末态；`FIG_02_LARGE_FRAGMENT_CLOT.png` 展示保存态中的最大游离分量阶段。
- `figures/FIG_03_ACTUAL_SURFACE_LOAD.png`：实际保存状态重放的当前表面载荷、真实力方向。
- `figures/FIG_04_LOAD_AND_FRAGMENT_HISTORY.png`、`FIG_05_ENERGY_LEDGER.png`、`FIG_06_FRAGMENT_SIZE_DISTRIBUTION.png`：数值曲线与历史对照。
- `videos/VIDEO_01_Y_CLOT_V4_FRAGMENTATION.mp4`：接受态碎裂演化，1080p，Z 轴竖直，红色点云，位移×1，无中间状态插值。
- `rotation/VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_4K.mp4` 及 `_1080P.mp4`：原三视图内容合成的一圈相机环绕，18 s / 30 fps，世界 Z 竖直，保留一份图例及完整坐标网格；物理时间始终 t=0。

载荷场完整 FP64 数据为 `ACTUAL_LOAD_FIELD.npz`。显示箭头为实际面力的确定性子集，统一长度比例；完整求积面不抽样。所有节点球形标记体积等于 Vi，仅是积分点表达，没有生成或补洞的连续碎块外壳。图表、视频来源哈希和完整解码结果见 `DELIVERY_VALIDATION.json` 与 `ARTIFACT_MANIFEST.json`。

本轮代码与小型证据在本地独立分支归档；大原生状态/检查点完整保留服务器。没有自动推送 GitHub，也没有重跑历史工况。
"""
 name='REVIEW_Y_CLOT_IRREGULAR_CURRENT_METHOD_V4.md'
 (D/name).write_text(lines);(O/name).write_text(lines);shutil.copy2(D/'REPRODUCE.md',O/'REPRODUCE.md')
 figs='\n'.join(f'<figure><img src="figures/{p.name}" loading="lazy"><figcaption>{html.escape(p.stem)}</figcaption></figure>' for p in sorted((F/'figures').glob('*.png')))
 page=f'''<!doctype html><html lang="zh"><meta charset="utf-8"><title>Y 血块 V4 实算与旋转</title><style>body{{max-width:1300px;margin:30px auto;font:18px/1.7 sans-serif;color:#263642}}img,video{{width:100%}}figure{{margin:24px 0}}a{{color:#8f1738}}</style><h1>原 Y 血块上的 V4 不规则 LPS</h1><p>{result['status']}；实际接受 {result['accepted_time_s']*1e3:.8f} ms。1927 材料点 +1561 锚点；412 内部点扰动，原形状边界保持。工程构造血块、未标定材料、指定机械剪切代理。</p><p><a href="../../../analysis/{TASK}/{name}">中文开发与科研报告</a> · <a href="../../../verification/{TASK}/{RUN}/SCIENTIFIC_OUTCOME.json">结果数据</a></p><h2>实际碎裂时程</h2><video controls src="videos/VIDEO_01_Y_CLOT_V4_FRAGMENTATION.mp4"></video><p>实际接受态，位移×1；Z 竖直。球形标记表示 Vi，不是连续外壳。</p><h2>初态相机绕 Z 轴旋转</h2><video controls loop src="rotation/VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_1080P.mp4"></video><p><a href="rotation/VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_4K.mp4">4K 主视频</a> · <a href="rotation/VIDEO_CLOT_DISCRETIZATION_Z_ROTATION_1080P.mp4">1080p PPT 版</a>。18 s / 30 fps，物理时间 t=0，只有相机旋转。</p>{figs}</html>'''
 (F/'OPEN_RESULTS.html').write_text(page)
 final=dict(task=TASK,run=RUN,status=result['status'],accepted_time_s=result['accepted_time_s'],formal_trajectories=1,statistics=q['statistics'],report=f'analysis/{TASK}/{name}',results=f'results/{TASK}/{RUN}/OPEN_RESULTS.html',delivery_validation='PENDING',scientific_labels=[HET_LABEL,LABEL,'PARTIAL_CURRENT_SURFACE_APPROXIMATE_COUPLING'])
 dump(O/'FINAL_STATUS.json',final);dump(D/'FINAL_STATUS.json',final)
 (D/'RESUME_STATUS.md').write_text(f'''# Y clot V4 — terminal physical result

RUN {RUN}: {result['status']}, accepted {result['accepted_time_s']*1e3:.8f} ms. The one authorized physical trajectory has finished. Do not relaunch the formal service or repeat the pilot.

Report: REVIEW_Y_CLOT_IRREGULAR_CURRENT_METHOD_V4.md. Results: ../../results/{TASK}/{RUN}/OPEN_RESULTS.html. Delivery validation runs next; FINAL_STATUS.json records its result. Server native archives and immutable physical source hashes remain in verification/{TASK}/{RUN}.

The separate rotation is t=0 camera motion only. Any later renderer correction is postprocessing and must preserve physical states.
''')
 print('REPORT_AND_INDEX_COMPLETE',flush=True)
if __name__=='__main__':main()
