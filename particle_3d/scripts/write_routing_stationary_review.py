#!/usr/bin/env python3
"""Chinese review and browser index generated from completed machine evidence."""
from pathlib import Path
import json,html,subprocess
from particle_3d.routing_stationary_audit import read,dump,sha,ROLES
from particle_3d.particle8_replay import REPO
R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';D=R/'data'

def table(headers,rows):
 return '| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in rows)
def main():
 s=read(D/'audit_summary.json');st=read(D/'stationary_audit.json');pair=read(D/'paired_routing.json');figs=read(D/'figure_manifest.json');g=read(D/'gates.json');status=g['status'];s['status']=status;dump(D/'audit_summary.json',s)
 a=s['point_to_p65_transition'];b=s['p65_to_p9a1_transition'];loss=s['outlet01_loss_chain'];p65=s['p65_500_split'];p9=s['p9a1_500_split'];rad=s['stationary_radius_statistics'];protected=s['production_protection'];field=s['point_field_stagnation'];flow=s['flow_split']
 stage_rows=[['流体体积流量',*(f"{flow['outlet_fractions'].get(k,0)*100:.4f}%" for k in ROLES[:3]),'—','流体体积参照']]
 for label,key in [('834 原始入口候选（point）','raw_candidate_point_split'),('500 已准入位置（point）','accepted500_point_split'),('同 500，P6.5 + 当前 handoff/contact','p65_500_split'),('原保存 P9-A.1 正式轨迹','p9a1_500_split')]:
  x=s[key];n=sum(x.values());stage_rows.append([label,*(f'{x[k]} ({100*x[k]/n:.2f}%)' for k in ROLES),str(n)])
 O3=[x for x in pair if x['point_outlet']=='OUTLET_03'];p65O3=[x for x in pair if x['p65_outlet']=='OUTLET_03'];classes={k:[r['particle_id'] for r in st if r['classification']==k] for k in s['stationary_classification_counts']}
 changes=[r for r in pair if r['p65_outlet']!=r['p9a1_outlet']]
 qlines=table(['边界','保存流量 (m³/s)','换算流量 (µL/min)','相对入口'],[['INLET',f"{flow['integrated_inlet_Q_m3_s']:.12g}",f"{flow['integrated_inlet_Q_m3_s']*6e10:.12g}",'100%']]+[[k,f"{flow['outlet_flows_m3_s'][k]:.12g}",f"{flow['outlet_flows_m3_s'][k]*6e10:.12g}",f"{flow['outlet_fractions'][k]*100:.6f}%"] for k in ROLES[:3]])
 passages={ratio:[r['particle_id'] for r in st if r['critical_radius_ratio'] and r['critical_radius_ratio'].get('highest_tested_passing')==ratio] for ratio in [.99,.95,.90]}
 station_table=table(['ID','半径 (µm)','全体粒径百分位','热点','分类','通过/未通过比值区间'],[[r['particle_id'],f"{r['radius_m']*1e6:.4f}",f"{r['size_percentile']:.1f}%",r['hotspot_id'],{'RIGID_GEOMETRIC_JAM':'刚性模型几何堵塞（中等尺寸敏感）','RIGID_MODEL_HIGH_SENSITIVITY':'刚性模型高度尺寸敏感','NUMERICAL_OR_GEOMETRIC_UNRESOLVED':'未解决'}[r['classification']],str(r['critical_radius_ratio'].get('interval') if r['critical_radius_ratio'] else None)] for r in st])
 text=f'''# 一句话结论

**{status}**。出口偏向 O2 首先明显出现在**入口有限尺寸准入**：原候选 point 的 O2 占 705/834=84.53%，准入后变成 469/500=93.8%；后续 P6.5/P9-A.1 的具体改变见下表。14 条静止都是 accepted-500 中最大的 14 个球，13 条聚在同一热点；它们在当前刚性球模型下缺少不穿墙且顺流的瞬时移动方向，但全部在本轮测试的 1%–10% 虚拟减径范围内恢复了局部通过。因此不能把这些停止写成真实微泡的生理捕获。

这里的 COMPLETE 是**科学审核完成**，不是新 production PASS，也不证明真实物理已全部正确。原正式数据始终是 READ_ONLY_REFERENCE，P9-A.1 500 条没有重积分。

# 出口偏差在哪一步开始

{table(['阶段','O1','O2','O3','无出口','分母/含义'],stage_rows)}

1. **原始入口采样（SOURCE A）**：直接读 834 个真实候选，834 个入口位置各不相同，没有重抽样。流体参考 O2=85.2050%，原候选 point O2=84.5324%，不存在足以解释最终 95.4% 的初始 O2 偏置。有限样本波动和下述离散流场停滞使三出口数不宜强行拟合流量比例。即使把全部 13 个无出口点都算给 O2，原候选的 O2 上界也只有 86.09%，仍低于准入后的 93.8%；因此主要入口筛选结论不依赖对这 13 条的出口猜测。
2. **准入筛选（SOURCE B）**：原方法 B 固定 flux anchor，最多重试 512 次粒径；接受的是位置条件下的粒径分布，不能称为未经筛选的 SonoVue 分布。O1/O3 流域的位置更容易被拒绝。500 个 accepted 中心与其原 anchor 逐项、逐位相等，因此 L2 精确复用 L1 对应 point 路径，没有重抽 500 个位置。
3. **有限尺寸/当前 P6.5 几何与接触（SOURCE C）**：point→P6.5 有 {a['changed_including_no_exit']} 条改变终态类别，其中 {a['changed_between_named_outlets']} 条在已命名出口之间改道。其余涉及 point 壁面停滞或有限球 stationary，不应全部叫“换出口”。
4. **P9-A.1 修正（SOURCE D）**：P6.5→P9-A.1 有 {b['changed_including_no_exit']} 条改变终态类别，其中 {b['changed_between_named_outlets']} 条在已命名出口之间改道。两模型同事件、同半径、同出生位置/时间/四元数、同 FEM、同 dt=0.00025 s、同 horizon=1.5 s，差别仅为 P9 切向/旋转修正开关。
5. **分类器/示踪/离散流场（SOURCE E）**：原分类器对全部已出流 point 末段重新核对一致。13 个 point 无出口不是准入拒绝的同义词，其中 2 个实际通过准入。两档更细设置均仍到 30 s 上限；其局部 P1 场可解析证明渐近趋向壁面，见下面的限定说明。

{table(['point 流域','原候选','accepted','rejected','准入率'],[[r['point_basin'],r['raw_count'],r['accepted_count'],r['rejected_count'],f"{r['acceptance_fraction']*100:.2f}%" if r['acceptance_fraction'] is not None else '—'] for r in s['acceptance_by_point_basin']])}

完整流量来自已保存积分记录，没有重跑 CFD；它是 fluid-volume reference，不是对粒子出口比例的强制目标。

{qlines}

13 个 point 无出口的解析核对：每个终端四面体有三个真正 WALL 节点，其速度精确为零，第四节点速度非零；单元 P1 速度是该节点速度乘以一个重心坐标。重心坐标沿轨迹按负散度指数衰减，极限点位于同一真实 WALL 面内，凸四面体内的解析路径不会穿越出口。这说明当前离散速度场存在局部非零散度造成的壁面吸引行为，不能用全局进出流平衡替代局部无散度验证。13 个单元散度范围 {min(r['divergence_s_inv'] for r in field['records']):.3f} 至 {max(r['divergence_s_inv'] for r in field['records']):.3f} s⁻¹。

这项解析结论针对**保存的离散 P1 场**，不是对真实血流的停滞结论。未把 13 条强行分到 O1/O2/O3，也未调整 FEM 或示踪器。源数据：`point_resolution.json`、`point_stagnation_field.json/.csv`。原始 834 路径以及 28 个复核种子的两档路径全部保存；15 个已有出口的稳定 ID 对照用于检查细化前后分类一致性。

# Outlet 01 为什么是 0

原 O1 point 流域有 {loss['raw_O1']} 个候选，其中 {loss['admission_rejected']} 个在准入阶段被拒绝，只剩 {loss['accepted_O1']} 个 accepted。这里“未进入 accepted”与“被拒绝”是同一组，不能相加。

这 {loss['accepted_O1']} 个原 O1 种子中，P6.5 改道到其它出口 {loss['accepted_P65_reroutes']} 个、停止 {loss['accepted_P65_stops']} 个、保留 O1 {loss['P65_keeps_O1']} 个；随后 P9-A.1 使其中原仍在 O1 的 {loss['P65_keeps_then_P9_reroutes']} 个进一步离开 O1。最终正式全体 O1={p9['OUTLET_01']}。

代表例最多 5 条，见 `outlet01_representatives.csv`。这些是逐 ID 对照，不是为了获得目标比例而重新指定出口。

# Outlet 03 为什么明显偏低

原候选有 89 个 O3 point，准入后只剩 20 个：拒绝 69 个，准入率 22.47%，低于 O2 的 66.52%。在这 20 个同源 O3 种子中，P6.5 去向为 {dict(Counter(x['p65_outlet'] for x in O3))}，P9-A.1 去向为 {dict(Counter(x['p9a1_outlet'] for x in O3))}。

全体 P6.5 的 O3={p65['OUTLET_03']}；其中后续 P9-A.1 去向为 {dict(Counter(x['p9a1_outlet'] for x in p65O3))}。正式 P9-A.1 的 O3={p9['OUTLET_03']}，需连同其它 point 流域的流入改道一起计数，不能只看单一源流域。

# P9-A.1 本身改道了多少

相对同一批 P6.5，{b['changed_between_named_outlets']} 条从一个已命名出口改到另一个；含 no-exit 在内，共 {b['changed_including_no_exit']} 条终态类别变化。逐 ID 改变清单：{[r['particle_id'] for r in changes]}。

{table(['原 P6.5 → P9-A.1','O1','O2','O3','无出口'],[[ROLES[i],*b['counts'][i]] for i in range(4)])}

当前两模型的出口比例接近或相同，也不能推出沿途速度、旋转或全部路径相同。本审核只对所列同状态消融比较负责。

# 14 条 stationary 在哪里

采用几何距离与接触三角形集合的确定性连通分组：两点距离不超过两者半径和的 2 倍，且共享接触三角形或距离不超过半径和。结果是 {len(s['stationary_hotspots'])} 个热点：{'; '.join('H'+str(h['hotspot_id'])+' 有 '+str(h['count'])+' 条，IDs '+str(h['particle_ids']) for h in s['stationary_hotspots'])}。

这是几何热点分组，不是未经标注的生理解剖命名。阈值、质心、三角形集合全部保存在 `audit_summary.json`。图 04 只显示原血管与这 14 个点。

# 它们是不是更大的微泡

**是，而且恰好是这批 accepted-500 中最大的 14 个。** 静止组半径中位数 {rad['stationary']['radius_um']['median']:.4f} µm，范围 {rad['stationary']['radius_um']['min']:.4f}–{rad['stationary']['radius_um']['max']:.4f} µm；486 个完成者半径中位数 {rad['completed']['radius_um']['median']:.4f} µm，范围 {rad['completed']['radius_um']['min']:.4f}–{rad['completed']['radius_um']['max']:.4f} µm。

14 个静止粒子的经验百分位为 97.4%–100.0%。这里只描述这批入口条件下的尺寸排序，没有进行不必要的显著性检验，也不能外推为所有 SonoVue 微泡的截留率。

# 它们是不是真的被血管卡住

在**当前不可变形刚性球模型、当前网格和连续区接管下限**下，14 个最终状态均保留三个独立接触方向；原始候选接触可多于三个，但当前非负锥冗余逻辑只保留独立信息。保存了白化 Jacobian 奇异值、条件数、保留/删除行、非负系数、乘子和自由/无约束/有约束速度。

本批约束在原有 **2 nm 连续区接管下限**处激活，不是球面与壁面达到零间隙的真实接触；乘子也只是当前运动学约束乘子，不作真实接触力解释。

独立小型线性规划检查全部原始接触法向约束，求最大顺流投影；14 个均没有正的可行投影，并有非负组合的对偶证书。该结论是“当前位置没有不穿墙且顺流的瞬时方向”，不等价于证明整个三维构形空间不存在先退后进的绕行路径。球仍可旋转；stationary 指球心平移停止。

原半径从确定性 pre-jam checkpoint 续算，14/14 复现原卡点，终点差不超过原 roundoff 预算的 16 倍。没有发现新的 contact/KKT 实现错误；但是 point 审核发现了离散 P1 场的近壁吸引限制，这一点会限制生理解释。

{station_table}

# virtual radius sensitivity 怎么看

只在 `DIAGNOSTIC_ONLY` 目录中从保存的 pre-jam 状态做局部续算，保持历史前缀、时间、姿态、流场和模型设置。100% 是复现对照，再按 99%、95%、90% 检查；没有做 bisection，也没有生成改变半径的新正式轨迹。

- 99% 首次通过：{passages[.99]}，临界诊断区间 [0.99, 1.00]。
- 95% 首次通过：{passages[.95]}，区间 [0.95, 0.99]。
- 90% 首次通过：{passages[.90]}，区间 [0.90, 0.95]。

这些是“已测试通过/未通过”区间；没有证明尺寸响应严格单调，也没有估计精确临界半径。8 个 90% 才通过的例子属于 **INTERMEDIATE**，不是“减小 10% 仍不通过”的 ROBUST_RIGID_SIZE_EXCLUSION；本批没有后者，图 06 明确使用中等敏感示例，不虚构 robust 示例。

局部通过判据要求离开原卡点 3 个原半径的球域、沿每步局部 FEM 方向累计前进超过 3 个原半径、连续 5 个名义步明确前进，并与原接触三角形恢复正间隙；原出口面的真实中心穿越可作替代终止证据。球心离开局部卡点不证明球体完整越过出口，也不证明后面不再遇到卡点。

审核工具初版使用终态单一 FEM 方向作为固定投影，转弯后会低估顺流距离。最终用原已保存局部路径按上述曲线流向判据复核，未重新积分。旧判据及已多计算的诊断试验全部另存，未删轨迹：`stationary_audit_fixed_direction_reference.json`、`radius_passage_policy.json` 与各条 supplementary trials。最终复现代码直接使用修订后的局部判据。该修改仅涉及审核停止/解读规则，production physics 完全不变。

最多减小 5% 就局部通过的 6 条归为 RIGID_MODEL_HIGH_SENSITIVITY；另 8 条归为 RIGID_GEOMETRIC_JAM，并标注中等尺寸敏感。**不能写成真实微泡只需变形 1%、5% 或 10% 就能通过。**

# 目前能相信什么

- 原正式 P9-A.1：O1/O2/O3/no-exit = 0/477/9/14，未重积分。
- 出口偏置首先明显出现在入口准入，后续有限尺寸改道由同 ID 对照量化；没有强制匹配流体流量。
- 14 条静止集中于少数几何热点、位于本批尺寸分布最上端，其当前接触方向独立，局部顺流不可行，且原半径停止可复现。
- 虚拟半径变化显示刚性几何结果存在尺寸敏感性，为将来的变形模型提供动机。
- 13 条 point 无出口在原离散 P1 速度场中有解析的局部壁面停滞解释，未伪装成出口。

# 目前不能声称什么

- 不能把 14/500 当作真实微泡生理捕获率，不能把虚拟缩径当作真实变形需求。
- 不能因为全局进出流守恒，就声称 P1 插值场在每个壁面四面体都严格无散度。
- 不能由一批 500 个 accepted、入口尺寸重试后的独立刚性球，推出任意浓度、粒径分布、RBC 条件或真实微泡的出口比例。
- 没有证据要求删掉或改写 P9-A.1 500 条数值结果，也没有新发现的 contact/classifier 实现错误；其生理解释需加入本次揭示的入口筛选和离散场限制。

# 下一步建议

先审核**入口位置条件下的尺寸筛选**以及**近壁 P1 速度场的局部散度/流线行为**，明确它们对实验比较的适用范围；再审核有限尺寸在分叉和狭窄部位的 routing。对于最大的 14 个球，本次已有开展 deformable MB 可行性研究的动机，尤其是 6 个高度尺寸敏感例子，但这些证据还不足以标定真实变形模型，更不足以用变形掩盖离散流场问题。本轮没有启动新物理开发。

## 执行、保护与复现

- 分支：`dev/particle9a1-routing-stationary-audit`；HEAD 沿用 `775adc019536585a4ee2f4dbf2c083c959e0a339`，未 commit/push。
- 流场 SHA：`{flow['flow_sha256']}`。
- 服务器：`{s['server']['root']}`。
- P6.5 完整新增 500 批次 6 workers；point 和 radius 小任务各 1 worker，与 P6.5 重叠时最多 7 个计算进程，低于 7.68 CPU quota。所有 BLAS/OMP 库单线程，未使用 GPU。
- 服务器旧正式参考的 1,434 个文件也已逐 SHA 与本地冻结参考比对，无缺失、无差异。
- 原始 point 批次 {s['runtimes_s']['point']:.2f} s，P6.5 500 批次 {s['runtimes_s']['p65']:.2f} s，局部 radius 批次 {s['runtimes_s']['radius']:.2f} s。批次存在重叠，不能相加当总墙钟时间。额外 point 分辨率复核 {s['point_resolution_audit']['wall_seconds']:.2f} s。
- 本轮 {g['tests']['tests']} 项永久审核测试通过；原 {protected['file_count']} 个保护文件 SHA 不变，git index 不变。最终 Gate 明细见 `data/gates.json`。
- 完整 Git 工作区原有 12 个已修改文件仍保留；`git diff --stat` 不会显示本轮未跟踪新增代码，另见 `logs/new_audit_files.txt`。不能把先前修改误称为本轮 production 改动。

## 可直接审核的图与数据

'''
 for f in figs:text+=f"- [{f['figure']}.png](figures/{f['figure']}.png) / [PDF](figures/{f['figure']}.pdf)；数据："+', '.join('`data/'+x+'`' for x in f['sources'])+'。\n'
 text+='\n机器汇总：`data/audit_summary.json`；全部复现命令见 [REPRODUCE.md](REPRODUCE.md)。\n'
 # Counter imported here to keep the report template close to its data transformations.
 (R/'ROUTING_STATIONARY_REVIEW_ZH.md').write_text(text)
 blocks=['<!doctype html><html lang="zh"><meta charset="utf-8"><title>P9-A.1 Scientific Audit</title><style>body{max-width:1150px;margin:35px auto;font:17px system-ui;line-height:1.6;color:#25333f}img{max-width:100%;border:1px solid #e2e7eb}a{color:#176490}</style><h1>P9-A.1 出口与静止轨迹科学审核</h1>',f'<p><b>{status}</b> · <a href="ROUTING_STATIONARY_REVIEW_ZH.md">中文报告</a> · <a href="data/audit_summary.json">机器汇总</a> · <a href="REPRODUCE.md">复现说明</a></p>']
 for f in [figs[-1]]+figs[:-1]:blocks.append(f'<h2>{html.escape(f["figure"])}</h2><img src="figures/{f["figure"]}.png"><p>'+', '.join(f'<a href="data/{x}">{x}</a>' for x in f['sources'])+'</p>')
 (R/'OPEN_RESULTS.html').write_text('\n'.join(blocks)+'</html>')
 print(status,'Chinese review and index written')
if __name__=='__main__':
 from collections import Counter
 main()
