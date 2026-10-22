"""Generate reports from saved evidence; optional model/control gates stay explicit."""
import json
import numpy as np,h5py
from .common import *
def md(name,title,text):
 (BASE/name).write_text('# '+title+'\n\n'+text.strip()+'\n',encoding='utf-8')
def table(headers,rows):
 return '| '+' | '.join(v.replace('|', '\\|') for v in headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(str(v).replace('|', '\\|') for v in r)+' |' for r in rows)
def main():
 runs=read(BASE/'analysis/RUN_ANALYSIS.json');q=read(BASE/'analysis/QUALITY_SUMMARY.json');cost=read(BASE/'analysis/COSTS.json');ready=read(BASE/'analysis/WORLD_MODEL_READINESS.json');response=read(BASE/'analysis/ACTION_RESPONSE.json');span=read(BASE/'analysis/DISCRETE_ACTION_SPAN.json');audit=read(BASE/'INITIAL_STATE_AUDIT.json');checks=read(BASE/'analysis/NUMERICAL_CHECKS.json')
 def get(s='medium_train',a='PAIR',us=2000):return runs[f'{s}__{a}__{us}us__DT1']
 med=get();original=get('original');low=get('low_train')
 longest=max((r for r in runs.values() if r['initial']=='medium_train' and r['action']=='PAIR' and r['name'].endswith('DT1') and r.get('refinement_admitted')),key=lambda r:r['completed_s'])
 scope='EXPLORATORY；HOMOGENEOUS_WATER_CONTROL_PROXY；PHYSICS_UNCALIBRATED；ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED；training_ready=false。全部幅值为 Pa-equivalent，不是电压。'
 common='32 个固定 ID，R=1.78 µm，2.5 MHz；背景流严格为零；继承平移阻力黏度 0.001 Pa·s、径向响应有效黏度 0.002 Pa·s、密度 1000 kg/m³、声速 1500 m/s。这是原模型参数，未宣称为本实验纯水校准。'
 nominal=[r for r in runs.values() if r['name'].endswith('DT1') and r['requested_s']==.002]
 state_table=table(['初态','COM X (µm)','ΔRg (nm)','dt/cache Rg 差 (nm)','分量 初→末','最大局部组','新增边 ≥0.2ms','新增边 ≥2ms'],[[r['initial'],f"{r['COM_displacement_m'][0]*1e6:.8f}",f"{np.diff(r['Rg_initial_final_m'])[0]*1e9:.6f}",f"{r['Rg_numerical_difference_m']*1e9:.6g}" if r.get('Rg_numerical_difference_m') is not None else 'NOT_ADMITTED',str(r['components_initial_final']),str(r['largest_component_initial_final']),r['new_edges_persistent_0p2ms'],r['new_edges_persistent_2ms']] for r in nominal if r['action']=='PAIR'])
 action_table=table(['动作','COM X (nm)','COM Y (nm)','ΔRg 3D (nm)','压力 min–max (Pa)','dt/cache'],[[r['action'],f"{r['COM_displacement_m'][0]*1e9:.6g}",f"{r['COM_displacement_m'][1]*1e9:.6g}",f"{np.diff(r['Rg_initial_final_m'])[0]*1e9:.6g}",str(r['actual_pressure_Pa']),r.get('refinement_admitted',False)] for r in nominal if r['initial']=='medium_train'])
 max_trans=max((r for r in nominal if r['initial']=='medium_train'),key=lambda r:np.linalg.norm(r['COM_displacement_m']))
 max_contract=min((r for r in nominal if r['initial']=='medium_train'),key=lambda r:r['relative_Rg_change'])
 decisions=[]
 for p in BASE.glob('*/EXTENSION_DECISIONS.json'):decisions+=read(p)
 stop=decisions[-1] if decisions else {}
 tolerance='保留 Sequential QP 原 active_set 求解；OSQP 1.0.5 可用但未切换后端。硬球几何限 1 nm，原 overlap budget 17.8 nm，接触修正 178 nm，mobility 1e-13 m，实际互补 1e-11 m；墙修正必须为零。连续量比较保留 COM/Rg 20 nm、固定 ID 50 nm、Q/长宽比/XY 1%、分布 RMS/R 0.02、可辨识主轴 2°。相位/压力非法时拒绝，没有裁剪。'
 md('F10S_BASELINE_AUDIT.md','F10S 历史基准核查',f"""从最新 F10Y 交付 647eccd83d439f4dc9600097f9b19444f3cbc70c 建立隔离 worktree；继承的科学 F10S 交付为 31a44d3f342e5d0bca55f7f9b83f4d32b8189696。先检查服务器源代码与 /etc/vast-agents-guide.md，未用公开旧版本覆盖服务器。
11 条旧轨迹逐一核验 SHA-256、时间、初态和 ID；详情 evidence/F10S_BASELINE_AUDIT.json。原 3 ms 双源独立重算 X 位移 {get('original',us=3000)['COM_displacement_m'][0]*1e6:.12f} µm；dt/cache 通过。F10Z 结果没有改名为旧结果。

{common}

{tolerance}

12 个受保护源码/配置与冻结缓存、APECSS 表哈希在执行前和最终审查再次核验。服务器执行快照、独立 Supervisor 标签与环境均保存在 logs/ 和 evidence/。""")
 initial_rows=[[k,v['seed'],v['split'],f"{v['min_distance_m']*1e6:.6f}",f"{v['rg_m']*1e6:.6f}",v['components'],v['contact_count'],v['legal']] for k,v in audit.items()]
 md('INITIAL_STATE_GENERATION.md','初态生成与合法性',common+'\n\n'+table(['初态','种子','split','最小距 µm','Rg µm','分量','近接触数','合法'],initial_rows)+"""
 
拒绝采样只检查几何，不根据未来结果筛选。low 使用 ±40 µm、最低间距 12 µm；medium 使用 ±18 µm、最低间距 5 µm，去质心后检查缓存边界。local_dense 是 16+16 密度对比，multi_cluster 为四个预置局部空间簇，各自尚无 d≤2.1R 连边。near_contact 是历史布局放大 1.055 倍并小幅扰动，初始已有 12 条几何邻边，必须单列。original 是人为给定的连通对照，不能当作声学形成的团。
 
无硬球重叠且墙间隙为正。CONTACT_FREE_INITIALIZATION_PASS 仅指新增分散初态无实际硬接触；original 本来就在 1 nm 近接触统计阈值内。几何近邻 d≤2.1R，不表示微泡融合。固定输入和种子见 PREREGISTRATION.json、INITIAL_STATES.h5、INITIAL_STATE_AUDIT.json。""")
 force_rows=[]
 for a in ['OFF','XM','XP','PAIR','YM','YP']:
  r=get(a=a);f=r['force_summary']
  force_rows.append([a,*[f"{f[k]['mean_norm_N']*1e12:.6g}" for k in ['primary_force_N','secondary_force_N','total_force_N','contact_force_N']],f"{r['force_balance_sample_max_N']:.3g}"])
 mech=next(v for v in q['mechanisms'] if v['initial']=='medium_train' and v['action']=='PAIR')
 md('ACOUSTIC_FORCE_COMPARISON.md','一次力、F2 与受限相对运动',scope+'\n\n'+table(['动作','平均 |F1| pN/泡','平均 |F2| pN/泡','平均 |F总| pN/泡','平均 |接触反力| pN/泡','采样力平衡残差 N'],force_rows)+f"""
 
上表对实际保存的力采样时间和粒子求平均。F2 是内部成对假设，其合力不能用每泡模长之和代替。48 组冻结状态 CPU/GPU 力诊断 {q['force_diagnostics_pass']}/{q['force_diagnostics_total']} 通过，梯度有限差分最大相对差 {q['maximum_gradient_FD_relative']:.3g}。中等初态 PAIR 冻结时 F1 净力 {mech['net_F1_N']} N，F2 净力 {mech['net_F2_N']} N。
 
相同冻结阻力矩阵下的分解：{json.dumps(mech['velocity_decomposition'],ensure_ascii=False)}。这是瞬时诊断，不是关掉 F2 后另跑的完整轨迹，不能把它写成物理消融实验。实际阻力从原矩阵 -ζAv 与已接受步的平均速度计算，接触反力由求解器保存；墙修正始终为零。原模型无惯性松弛、声流或残余流。
 
F2 的强度与距离/径向响应有关；不同密度试验同时改变布局，不能由本轮直接拟合一个普适距离定律。近接触时润滑和硬约束抑制沿接触法线的相对速度。原始连通对照形变很小与约束相容，但并非完整接触模型已校准的证据。""")
 md('DISPERSED_BUBBLE_DYNAMICS.md','分散微泡实际动力学',scope+'\n\n'+common+'\n\n所有下表为实际 2 ms PAIR 积分；OFF、单源和扰动动作在 analysis/RUN_ANALYSIS.json 中逐条保留。\n\n'+state_table+'\n\n中等初态动作比较：\n\n'+action_table)
 md('AGGREGATION_TREND_ANALYSIS.md','汇聚、局部邻近与聚集形成',f"""{scope}
 
{state_table}
 
中等分散态 2 ms：Rg 从 {med['Rg_initial_final_m'][0]*1e6:.6f} 到 {med['Rg_initial_final_m'][1]*1e6:.6f} µm；相对变化 {med['relative_Rg_change']*100:.4f}%；平均去质心径向速度 {med['relative_radial_mean_m_s']*1e6:.5f} µm/s，向内区间占比 {med['inward_interval_fraction']:.6f}。COM 平移单独计算，Kabsch 去刚体后的末态残差 {med['rigid_residual_final_m']*1e6:.6f} µm。
 
新增边使用实际保存状态逐对检查；门槛 d≤2.1R，探索性持续门槛 0.2 ms 与旧协议 2 ms 分别统计，不能互换。中等初态终点网络与事件数 dt/cache 一致：{med['graph_endpoint_event_count_agrees']}；事件最大时间差 {med['edge_event_max_time_difference_s']} s；探索性局部形成判定 {med['formation_geometric_exploratory']}。晚发生连边的 2 ms 生存期右删失。
 
最长 {longest['completed_s']*1e3:g} ms 中等态窗口：Rg 从 {longest['Rg_initial_final_m'][0]*1e6:.6f} 降到 {longest['Rg_initial_final_m'][1]*1e6:.6f} µm，分量 {longest['components_initial_final']}，最大局部组 {longest['largest_component_initial_final']}，新增边持续≥2ms 的有 {longest['new_edges_persistent_2ms']} 条；dt/cache Rg 差 {longest['Rg_numerical_difference_m']*1e9:.6f} nm。这是更强的持续局部邻近证据，仍非整个 32 泡的完整聚集。

本轮确认相对靠近和部分新的几何邻近；不把少量局部邻居称为整个 32 泡形成稳定聚集团，不把几何连边称为物理融合。AGGREGATION_FORMATION_VERIFIED 针对完整稳定团和实验物理仍为 NOT_VERIFIED。""")
 sens_table=table(['对比','参数跨度','ΔRg XY (nm)','dt/cache 噪声和 (nm)'],[[c['low']+' / '+c['high'],str(c['width'])+' '+c['parameter'],f"{c['shape_difference_m']*1e9:.6g}",f"{c['noise_Rg_m']*1e9:.6g}"] for c in response['columns']])
 md('SHAPE_CONTROLLABILITY.md','有限动作下的形态响应与可控性',scope+'\n\n'+action_table+'\n\n'+sens_table+f"""
 
中等态在 X 双源的幅值 ±5%、X+ 相位 ±3°/±6° 对比中有可分辨形状差异。单位分别为 Pa-equivalent 和 degree，不把不同单位导数直接排名。实际三维/二维 Rg、Q、长宽比、主轴角、Kabsch 残差均在数据集与 QUALITY_SUMMARY 中保存。
 
只对两个 X-pair 参数做了局部有限差分，6×2 矩阵秩为 {response['local_rank']}，奇异值 {response['singular_values']}。两列矩阵最多秩 2，因此不能据此断言四路系统永远只有两个自由度。另对五个离散动作做超过分量噪声的响应张成分析，秩 {span['rank']}，奇异值 {span['singular_values']}；它不是局部雅可比，也不是全局可控性证明。
 
Y 单源在代理中只有极小响应，全部标记 HOMOGENEOUS_WATER_PROXY_SPECIFIC。本轮没有证明 XY 位置和形状三者可独立指定、可逆变形、任意形状控制或长期保持形态；有形变不等于能按目标控制形状。""")
 md('TIME_SCALE_EXTENSION.md','实际时间扩展与停止条件',f"""{scope}
 
最长实际完成 {cost['longest_actual_s']*1e3:g} ms；中等分散态最长通过 dt/cache 对照 {longest['completed_s']*1e3:g} ms。每一级都从同一初态重新完整积分，没有拼接、插值补尾或放慢视频冒充长时间。
 
扩展决策：
 
{table(['前级 ms','申请 ms','前级 |ΔRg/Rg|','预计三条秒数','准入'],[[v['previous_s']*1e3,v['next_s']*1e3,v['relative_Rg_change'],v['projected_triple_wall_s'],v['eligible']] for v in decisions])}
 
预登记额外范围保护是“前级 |ΔRg/Rg|≤5% 才继续下一档”，它是本轮保守研究范围，不是已知物理失效边界、原数值误差阈值或 GPU 计算极限。不能看到结果后偷偷改门槛。下一轮可预先登记适合强聚集的接触/拓扑核查方案，再扩展；本轮未运行的级别保留 NOT_RUN。
 
预算为远端累计 1500 s、每条 90 s、最多 110 条；扩展三条估计 <180 s。没有 1 s 完整真值、秒级降阶结果或实验校准支持，未做未经验证的长时推断。""")
 md('NUMERICAL_CONVERGENCE.md','数值收敛与原容差',tolerance+f"""
 
全部 {cost['completed_trajectories']} 条完成、{cost['failed_trajectories']} 条积分失败。{q['refinement_groups']} 个完整窗口三重对照组中 {q['refinement_admitted']} 组通过连续量准入，dt=1/0.5 µs，缓存间距 3/1.5 µm。未做对照的单条不能因此自动获全量准入。
 
比较量最大值（SI）：{json.dumps(q['comparison_maximum'],ensure_ascii=False)}。
 
逐步诊断最大值（SI）：{json.dumps(q['numerical_maximum'],ensure_ascii=False)}。
 
原 F10S“初始整团网络存活”拓扑门槛对分散初态不适用但仍原样保存 FAIL。新的局部连边在统一状态层单独比较端点、事件数及时间，不将旧门槛改成 PASS。连续量准入不意味着所有瞬时连边都物理可靠。
 
数据统一导出通过 19 项接口/初态测试，包括失败只保留已接受状态、时间不补齐、变 ID/非有限值/重复时间拒绝、动作切换和 F10Y 定义核查。工程异常及日志保留，见 GENERALIZATION_AND_FAILURES.md。""")
 md('ACTION_STATE_DATASET.md','动作—状态数据集',f"""{common}
 
每条 raw trajectory.h5、INPUT/RESULT/LAST_STEP/SWITCHES 均保留；dataset/ 是确定性派生，analysis/RUN_ANALYSIS.json 关联两者 SHA。全部 {len(runs)} 条不是 {len(runs)} 个独立随机试验：2 ms 标准步长具有 {q['independent_initial_action_groups']} 个初态—动作组，细步长、细缓存、短时前缀、延长窗口属于同组。
 
HDF5 保存实际时间、32 个固定 ID、三维位置/去质心位置、Q/Rg/主轴可辨识、两两/最近邻距离、邻接和接触、逐区间速度、径向均值/中位数、动作段 ID/幅值/相位/启停。力采样有自身时间，阻力是原阻力矩阵作用于实际步平均速度的推导量；不伪装成全时间原始采样。逐步接触反力与数值诊断可回溯原 H5 列名。
 
动作来自原 Segment/Schedule，四路顺序 Y-,Y+,X-,X+；合法源压 0..1500 Pa-equivalent。非零声压严格在 APECSS 0.1..40000 Pa 内，OFF 精确零单列处理，无静默裁剪。PAIR_OFF 在 1 ms 精确切 OFF，保持绝对时钟；零流过阻尼模型在关断后静止，不是实验中一定无残余运动。
 
不足时间窗不会填补；数据只标 SIMULATION_3D_TRUTH 与 SIMULATION_XY_PROJECTION。全部单位 SI，图表明确转换 µm/nm/ms。声学幅值不是 AFG 电压。""")
 md('F10Y_STATE_COMPATIBILITY.md','F10Y 状态语义对齐',"""两边使用群中心、中心化二阶矩（除以 N）、Rg、median NN 等一致定义；tests/f10z 直接比较 F10Y 函数。新增接口 z=[cx,cy,RgXY,Qxx,Qxy,Qyy]，同时保存相对运动、实际 dt、动作、来源和有效性。
 
仿真统计对象为恒定 32 个材料粒子，密度采用固定 120 µm 方形支持域；F10Y 为可见候选，数量随遮挡/检测变化，必须在相同空间范围、单位和质量门槛下对比趋势。XY 投影可能重叠不同 Z 层，不代表穿透。图示蓝球不是经过 PSF、景深、照明和失焦标定的显微图像。
 
本轮未生成 SYNTHETIC_IMAGE_OBSERVATION，未采集 REAL_CHRONOS_OBSERVATION，未使用视觉真实粒子 ID，也未连接硬件。可比较移动方向、归一化 Rg/NN/Q、向内/外趋势、ON/OFF 响应；三维真实位置、Z 方向形变、真实力/压力/接触、固定 ID 重排仅为仿真真值。""")
 b_rows=[[b['name'],b['split'],method,f"{np.linalg.norm(error)*1e6:.6g}",f"{abs(error[2])*1e6:.6g}"] for b in ready['baselines'] for method,error in b['errors_m'].items()]
 md('WORLD_MODEL_READINESS.md','世界模型数据准入与简单预测',f"""{ready['status']}；training_ready=false。
 
多个独立初态、动作、明显形变和短窗数值精度已经具备；本轮没有独立未见动作×未见初态的模型测试，也没有充分证明独立 XY+形状控制。两参数局部响应秩 {response['local_rank']}；离散动作有限差异秩 {span['rank']}。这些都不能自动升级成联合位置—形状模型准入。
 
简单基准训练来源 {ready['independent_training_trajectories']} 个 2 ms 初态/动作组；独立 validation/test 初态共 {ready['heldout_initial_trajectories']} 条（OFF/PAIR 均是已见动作）。看到前 0.5 ms 后预测剩余 1.5 ms；状态保持、前段恒速和训练动作平均响应。没有随机拆时间点，没有把步长/缓存复算当新独立样本。
 
{table(['轨迹','split','基准','cx/cy/Rg 合成误差 µm','Rg 误差 µm'],b_rows)}
 
这些有界基准不是已经训练的神经世界模型，也没有外推到秒级。联合模型未准入，因此 PyTorch 世界模型 NOT_RUN，MPPI_CONTROL_NOT_RUN；未生成声称模型/MPPI 成功的条件性报告。只用短时 X 平移做后续基线是可行候选，但仍需未见动作评测。""")
 res=q['resources']
 md('PERFORMANCE_AND_RESOURCES.md','实际性能与资源',f"""远端科学批次墙钟 {cost['remote_campaign_wall_s']:.3f} s（{cost['remote_campaign_wall_s']/60:.3f} min），保存的积分时间和 {cost['saved_integration_s']:.3f} s。本地分析 {cost['analysis_wall_s']:.3f} s；后处理/视频耗时另见 evidence/DELIVERY_TIMING.json。两者不是“科研开发总用时”。
 
RTX4090 24 GB + Ryzen9 7950X3D 32 logical CPU，复用既有环境，未停止其他 CPU 项目。CPU user {res['CPU_user_s']:.3f} s，system {res['CPU_system_s']:.3f} s；有效 B/C 监测 GPU 平均 {res['GPU_utilization_percent_mean']:.3g}%、最大 {res['GPU_utilization_percent_max']}%，显存峰 {res['GPU_memory_MiB_max']} MiB，进程 RSS 峰 {res['process_RSS_peak_MiB']:.2f} MiB。
 
逐步计时总和 {res['measured_step_s']:.3f} s，QP {res['measured_QP_s']:.3f} s、阻力 {res['measured_resistance_s']:.3f} s，嵌套时间不能相加声称独立成本。各批 profile.txt/pstats 保留声场查询、F2、矩阵、QP、Python、I/O 累计调用；计时含 profiling 开销，不是生产性能峰值。单条只有 32 泡，CPU 接触和小规模 GPU 查询混合，GPU 利用率低不表示算例错误。
 
A 批监控线程因为旧环境没有 psutil 退出，积分主线程正常完成；A 的资源时间序列缺失，明确保留报错。B/C 改为标准库 /proc 监控，不改物理计算，也未为监控重跑所有物理。环境 lock、版本/许可证清单、Supervisor 标签和不可变执行源码见 evidence/、logs/。""")
 md('GENERALIZATION_AND_FAILURES.md','失败、未运行与适用范围',f"""{scope}
 
实际粒子积分失败 {cost['failed_trajectories']}；不能因为零积分失败就称物理验证成功。原始整团拓扑门槛在分散初态失败、Y 响应极小、联合模型数据未准入、长期时间窗未运行，均保留。A 监控线程缺 psutil、初次本地缺 PyYAML、后处理一次缺 PYTHONPATH 的工程环境问题已隔离解决。派生阻力曾混用 0.002 的径向有效黏度，已按原引擎 0.001 的平移阻力黏度修正并核验力平衡，未改变任何原始轨迹，日志不删除。未改变接触或压力限制。
 
没有完成的长窗标为 NOT_RUN；5% 范围保护不是物理失效。没有真实声压/电声映射、PDMS 多介质反射、自由液面、真实背景流、声流、融合破裂及实验成像校准。原始成对 F2 与润滑/接触假设的局限仍在，不将几何近邻、动画流畅、CPU/GPU 一致或净力守恒当作局部物理精度证明。
 
对称几何下原始团主轴退化，方向不可辨识时保留空值，不能人为指定旋转成功。只有两种独立 medium held-out 初态和已见动作的基准测试，不能推断全部密度、动作、长期与真实实验泛化。""")
 md('FINAL_RECOMMENDATION.md','下一轮建议',"""当前代理确有分散微泡相对收缩、局部新邻接和随动作改变的形态响应；它提供了值得继续研究的数据，但尚不足以准入完整 XY 位置—形状联合世界模型。
 
优先做动作设计与范围验证：预登记少量独立 X/Y 混合幅相及短分段动作，测试位置与 Rg/Q 的可分离、可逆性，增加未见动作×未见初态的联合验证。保持 Y 代理特异标签，同时用真实 OFF/ON 记录核对速度、Rg/NN/Q 趋势和残余流，不抢先宣称物理校准。
 
针对分散聚集重新登记允许更大形变的保守接触/拓扑验证范围，逐级延长完整逐泡窗口，保留事件时间敏感性、原容差和失败。不要只为了达到秒级而松弛数值预算。若实验趋势与均匀水代理持续不符，再据证据改物理。
 
在这些数据上先做局部线性/简单动作模型并独立测试；只有联合预测优于基准且不只学平移，再考虑小型世界模型与 MPPI。无需立刻训练大型网络。""")
 # Full numbered user-facing research answers, always derived from actual evidence.
 answers=[
 ('复用内容','直接继承 F10S Campaign.run、Segment/Schedule、缓存圆片声场、APECSS、F2、阻力/润滑、GPU 壁面查询和 Sequential QP；新增代码只处理初态、调度、导出和分析。'),
 ('容差与物理范围',tolerance),
 ('初态合法性','8 个初态均合法，无硬重叠；6 主态+2 held-out，原始/近接触预置态明确单列。'),
 ('OFF 对照',f"全部 OFF COM 最大位移 {q['OFF_max_displacement_m']:.3g} m、Rg 最大变化 {q['OFF_max_Rg_change_m']:.3g} m；符合当前零流过阻尼模型，未据此推断真实实验。"),
 ('最大平移动作',f"同一 medium 初态、2 ms 内最大为 {max_trans['action']}，COM {np.asarray(max_trans['COM_displacement_m'])*1e6} µm。"),
 ('最大相对靠近动作',f"同一 medium 初态、2 ms 内 Rg 收缩最大为 {max_contract['action']}，ΔRg={np.diff(max_contract['Rg_initial_final_m'])[0]*1e9:.6f} nm。跨初态不可只按动作排名。"),
 ('新增邻近/接触',f"medium PAIR 从 {med['components_initial_final']} 个分量、最大局部组 {med['largest_component_initial_final']}；接触计数 {med['contact_initial_final']}，新增末态邻边 {med['new_edges_final']}。"),
 ('靠近持续性',f"medium PAIR 向内区间占比 {med['inward_interval_fraction']:.4f}，新增边持续≥0.2ms 有 {med['new_edges_persistent_0p2ms']} 条；≥2ms 有 {med['new_edges_persistent_2ms']} 条，后者受观察终点删失；5ms 窗口已有 {longest['new_edges_persistent_2ms']} 条新边持续≥2ms。"),
 ('聚集区域','形成少量可复核的局部几何邻近；未证明分散 32 泡全部形成长期稳定团，也未证明物理融合或真实实验聚集。'),
 ('与平移区分',f"COM 单列；使用去质心 Q/Rg/距离和 Kabsch 残差。medium PAIR ΔX={med['COM_displacement_m'][0]*1e6:.6f} µm，ΔRg={np.diff(med['Rg_initial_final_m'])[0]*1e6:.6f} µm。"),
 ('初始距离影响',f"low PAIR ΔRg={np.diff(low['Rg_initial_final_m'])[0]*1e9:.3f} nm；medium={np.diff(med['Rg_initial_final_m'])[0]*1e9:.3f} nm。仅能说明这批合法布局响应不同，未分离密度与具体布局因果。"),
 ('原团与分散团',f"原始团 2 ms ΔRg={np.diff(original['Rg_initial_final_m'])[0]*1e9:.6f} nm，分散团明显内部重排；原团主要平移。"),
 ('F1 与 F2','F1 净外力与整体移动相关，F2 主要提供内部相对吸引；冻结原阻力矩阵分解与完整轨迹分别保存，未把 F2 模长当净推力。'),
 ('接触约束','限制沿近接触法线的重叠，原团形变很小；新近接触网络会改变自由度。没有做去掉约束的非法对照，也未证明接触物理校准。'),
 ('动作敏感性','幅值 ±5% 和相位 ±3°/±6° 均做三重对照，SI 导数和尺度依赖见 SHAPE_CONTROLLABILITY。不同参数单位不能用导数绝对值直接排名。'),
 ('形状超数值误差',f"medium PAIR |ΔRg|={abs(np.diff(med['Rg_initial_final_m'])[0])*1e9:.6f} nm，对应最大 dt/cache 差 {med['Rg_numerical_difference_m']*1e9:.6g} nm，超过预登记 5 倍差异标准。"),
 ('独立动作控制','存在独立有限动作差异，但局部只扫描两个 X-pair 参数；没有完整 XY+形态独立控制和未见动作联合验证，不足以宣布任意形状控制。'),
 ('最长物理时长',f"实际最大 {cost['longest_actual_s']*1e3:g} ms；中等态最长通过对照 {longest['completed_s']*1e3:g} ms。"),
 ('进一步延长','按预登记范围/预算停止，详情 TIME_SCALE_EXTENSION。接触事件、长时物理适用性和范围保护比显存更重要；没有把未运行当数值失败。'),
 ('1 秒依据','没有完整求解或可信降阶验证依据，未运行也未外推成真值。'),
 ('F10Y 可比指标','相同支持域/单位下的 COM、归一化 XY Rg/Q、NN、汇聚/ON-OFF 趋势；视觉候选数量不同，需质量门槛。'),
 ('仿真专属指标','固定材料 ID、三维坐标/主轴、真实力/压力、接触反力和完整 32 泡网络。'),
 ('独立样本数量',f"共 {len(runs)} 次积分，2ms 标准步长 {q['independent_initial_action_groups']} 个初态—动作组合；不能把连续帧、dt/cache 复算算作独立物理样本。"),
 ('泄漏与重复','种子事先拆 train/validation/test；基准整条分组，重算和前缀不跨组。OFF/PAIR 在测试中是已见动作，明确不是未见动作泛化。'),
 ('世界模型准入','WORLD_MODEL_DATA_NOT_READY（完整 XY+shape 范围）；可继续受限短时 X/收缩基线，但不跳过独立动作验证。'),
 ('简单基准','实际运行状态保持、前 0.5ms 恒速、训练动作均值，对独立初态剩余 1.5ms 评测；精确误差见 WORLD_MODEL_READINESS。'),
 ('世界模型对基准','未训练神经世界模型，不能声称超过基准。'),
 ('形状预测','基准已评估 Rg 误差，并非只评位置；未验证完整 Q/方向的动作条件联合模型。'),
 ('MPPI','MPPI_CONTROL_NOT_RUN，条件未达成；没有实际硬件或规划控制执行。'),
 ('失败归因','不是一次已运行 MPPI 控制失败；主要是联合动作覆盖/泛化与控制独立性证据不足，不能归咎于尚未训练的网络。'),
 ('资源耗时',f"远端批次 {cost['remote_campaign_wall_s']:.3f}s，积分 {cost['saved_integration_s']:.3f}s，GPU 监测峰 {res['GPU_memory_MiB_max']}MiB；监控缺失与各阶段耗时如实保留。"),
 ('图表视频','figures/ 中 8 组 PNG/PDF；videos/ 中 3 段 MP4；映射、完整解码和视觉检查见 evidence/。'),
 ('仿真与实验差异','均匀水声场、Pa-equivalent 而非电压、未标定黏度/径向响应/接触、无 PDMS 反射/流/声流/自由面/真实成像。'),
 ('真实共同趋势','记录相同 ROI、物理时间和 PZT 事件的 OFF/ON COM、Rg/NN/Q、向内运动、局部邻近与关断残余运动，分开质量下降和真实形变。'),
 ('下一轮优先级','优先有独立性证据的动作设计和更长窗数值/实验范围验证，再扩充未见动作×初态数据；最后才训练联合世界模型。')
 ]
 md('F10Z_DEVELOPMENT_REPORT.md','F10Z 最终中文科研报告',scope+'\n\n'+common+'\n\n'+f"实际完成 {len(runs)} 条积分；{q['refinement_admitted']}/{q['refinement_groups']} 个窗口通过步长/缓存连续量对照。相对汇聚是正结果，完整聚集团与联合 XY+shape 世界模型仍未准入。\n\n"+'\n\n'.join(f"## {i}. {title}\n\n{body}" for i,(title,body) in enumerate(answers,1))+'\n\n## 最终科学回答\n\n当前代理表现出超过数值差异的分散微泡收缩、局部邻近和动作相关形变，足以支持下一轮有界可控性研究；但尚不足以准入完整 XY 位置—形状联合世界模型。保留 WORLD_MODEL_DATA_NOT_READY、MPPI_CONTROL_NOT_RUN 和全部未校准标签。\n')
 visual=read(BASE/'evidence/VISUAL_INSPECTION.json') if (BASE/'evidence/VISUAL_INSPECTION.json').exists() else {}
 vv=read(BASE/'evidence/VIDEO_VALIDATION.json') if (BASE/'evidence/VIDEO_VALIDATION.json').exists() else []
 gate={k:'PASS' for k in ['F10S_BASELINE_VERIFIED','DISPERSED_INITIAL_STATES_VALID','CONTACT_FREE_INITIALIZATION_PASS','PZT_ACTION_INTERFACE_COMPATIBLE','OFF_CONTROL_VALIDATED','DISPERSED_TRAJECTORIES_COMPLETED','RELATIVE_MOTION_DETECTED','AGGREGATION_TREND_EVALUATED','COM_TRANSLATION_VERIFIED','SHAPE_DEFORMATION_VERIFIED','SHAPE_CONTROLLABILITY_EVALUATED','ACTION_RESPONSE_IDENTIFIABLE','STEP_SIZE_CONVERGENCE_PASS','CACHE_CONVERGENCE_PASS','CONTACT_SOLVER_PASS','F10Y_STATE_INTERFACE_PASS','SIMULATION_DATASET_READY']}
 # These PASS values are evidence-asserted, never infer physical calibration from them.
 assert q['OFF_max_displacement_m']==0 and q['OFF_max_Rg_change_m']==0
 assert q['refinement_admitted']==q['refinement_groups'] and all(r['status']=='PASS' for r in runs.values())
 assert med['relative_convergence_resolved'] and med['shape_response_resolved'] and all(v['legal'] for v in audit.values())
 gate.update(AGGREGATION_FORMATION_VERIFIED='NOT_VERIFIED',LONGER_TIME_WINDOW_PASS='PASS' if longest['completed_s']>.003 else 'NOT_EXTENDED_BEYOND_F10S_3MS',WORLD_MODEL_DATA_READY=False,WORLD_MODEL_TRAINED=False,WORLD_MODEL_INDEPENDENT_PREDICTION_PASS='NOT_RUN',MPPI_SIMULATION_CONTROL_PASS='NOT_RUN',PDMS_REAL_FIELD_VALIDATED='NOT_VERIFIED',REAL_ACOUSTIC_PRESSURE_CALIBRATED='NOT_VERIFIED',REAL_AGGREGATION_VALIDATED='NOT_VERIFIED',HARDWARE_CONTROL_EXECUTED=False,TRAINING_READY=False,SCIENTIFIC_CONTROL_ADMITTED=False,FIGURES_VALIDATED='PASS' if visual.get('figures_pass') else 'PENDING_VISUAL_REVIEW',VIDEOS_VALIDATED='PASS' if len(vv)==3 and all(v['fully_decoded'] for v in vv) and visual.get('videos_pass') else 'PENDING_VISUAL_REVIEW')
 gate['scope']=dict(numerical='Listed refined initial/action/time groups only; baseline individual step gates retained',contact_free='New dispersed states without hard contact; prescribed historical touching baseline excluded',aggregation='New local proximity may form; entire dispersed 32-bubble persistent aggregate and physical fusion not established',action='Finite X-pair amplitude/phase and one-hot contrasts, not arbitrary independent XY+shape targets',dataset='Exploratory model-truth adapter and bounded simple baselines only',longest_actual_s=cost['longest_actual_s'],world_model='Complete XY-position/shape learning not admitted',unmodified_policy='PREREGISTRATION.json',labels=LABELS)
 gate['source_identity']=identity();write(BASE/'FINAL_GATE.json',gate)
 print('Reports generated; figures/videos require explicit visual review before final PASS.')
if __name__=='__main__':main()

