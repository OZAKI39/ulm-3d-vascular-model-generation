from common import *
def read(n):return json.loads((O/(n+'.json')).read_text())
s=read('FINAL_STATUS');e=read('DYNAMIC_ENERGY_LEDGER');h=read('HISTORICAL_TRAJECTORY_COMPARISON');r=read('RUN_RESOURCES');d=read('DOMAIN_CHANGE_AUDIT');c=read('CONTACT_AND_OVERLAP_AUDIT');p=read('INITIAL_STATE_PARITY');a=read('FULL_700US_DYNAMIC_AUDIT');co=read('INDEPENDENT_COMPONENT_RECONSTRUCTION');lib=read('SOURCE_AND_LIBRARY_IDENTITY');doc=read('WSL_DOCUMENT_UPDATE') if (O/'WSL_DOCUMENT_UPDATE.json').exists() else None
assert s['status']=='PASS'
answers=[
('上轮688µs丢点原因是否确认？','确认。旧域fff、各轴±250µm；节点73在687.75µs的y=249.953987866µm、vy=0.400969997m/s、合力为0，下一原生VV漂移预测y=250.054230365µm，越出旧域；旧原生日志确报340/341。新运行真实688µs状态保留该点并与预测相容。它是数值域删除节点，未发现由此可归因的软化自由能错误。'),
('新域范围？','X/Y/Z均为[−400,+400]µm，boundary f f f。由18个旧完整态的最大绝对坐标249.954µm、全旧接受步最大速度0.792842m/s、余下12.25µs估计9.712µm运动，再加100µm空白余量，向上取整到50µm格点。速度外推是规划估计而非严格未来速度上界；每次候选仍检查实际VV预测包络。'),
('原始节点、键、m和材料保持？','是。341节点、12047参考键、245自由点/96规定运动点；质量ρV、完整固定m、r0、V和taper及G=940Pa、K=33333.333333Pa、ρ=1000kg/m³、dx=50µm、δ=150.5µm、α=0、s0=0.15、sc=0.18均保持。参考文件SHA与上轮相同，原始键ID逐步核对，末态m精确相等。'),
('初始力和储能改变？','没有。新旧域的原始坐标、力、完整储能及接触能等初始物理量逐值一致；无缩放、平移、墙、周期相互作用或额外接触。初始能量/力比较见INITIAL_STATE_PARITY.json。'),
('异常退出能保存最后接受态？','限定测试通过。用真实683µs旧原生检查点，恢复同名fix后扩大域，在候选推进前注入一次绑定异常。最后合法x/v/F/κ/d/U/D/time及原生拓扑、special排除全部读回一致，事件落盘，候选没有覆盖接受态。异常处理仅用已缓存合法状态，不再gather或强写restart。主运行每接受步刷新标量、事件、FP64节点流；每40接受步轮换两槽原生restart，沿用损伤检查点，越界危险前另存。硬进程终止仅保证既有检查点及已落盘数据，不宣称任何故障下均零丢失。'),
('是否从t=0正式重跑？','是。唯一正式物理工况从原始无损伤t=0出发。P0与P1均无物理推进；未拼接旧前缀或续算12.25µs。'),
('到达700µs？',f'是。实际累计时间{s["final_time_s"]:.17g}s，浮点容差内为700µs，升载、保持、卸载和末保持全部完成，无停止原因。'),
('341原始点全部保留？',f'是。2805个合法记录（含初态）逐一检查原生ID恰为1–341；最小数值边界余量{d["minimum_actual_clearance_m"]*1e6:.6f}µm，末态节点73的y=254.865870326µm。'),
('接受/拒绝次数？','2804接受，4次Δd拒绝；2796个0.25µs接受步、8个0.125µs接受步，覆盖2800个原基准区间。4次完整回滚全部通过。冻结最小0.03125µs保持，实际不需要更小步。'),
('最大Δd？',f'{s["max_accepted_dd"]:.17g}≤0.10，门槛未放宽。'),
('最大全程能量余额？',f'{s["max_relative_energy_balance"]*100:.9f}%，低于2%；末态R={e["final"]["R_J"]*1e12:.12f}pJ，末态归一化余额{e["final"]["relative_R"]*100:.9f}%。原归一化和完整能量账本保持。'),
('动态接触保持？',f'是。1.0dx及刚度不变；232.5µs首次触发，{c["active_accepted_steps"]}个接受步有接触、峰值{c["max_pairs"]}对、末态{c["final_pairs"]}对。逐步70374条真实接触对记录与计数完全相符，独立力最大差{c["max_force_error_N"]:.6g}N、接触能差{c["max_U_error_J"]:.6g}J。'),
('末态失效键多少？','4937条完全失效，7110条存活；相对旧687.75µs的4918失效，多19条。失效键原生张力为零，没有重复提交。'),
('末态真实连通分量？','4个，192、147、1、1，合计341。两个大分量最大点间距离分别547.167303µm和433.336158µm，质量分别2.4000×10⁻⁸kg、1.8375×10⁻⁸kg；两个孤立节点73、283各1.25×10⁻¹⁰kg。真实坐标、速度、质心、包络和尺寸逐分量保存，点分量自身最大点间距为0。这里只称离散图分量，未认证实际血栓碎块。'),
('原生与独立拓扑一致？','是。正常结束时直接获取全部原始atom ID与7110个原生存活端点对，和原始12047键中native d<1集合完全一致；原生失败标志、原始键ID及独立图分量也一致。'),
('新旧重叠时程相容？',f'是。0–687.75µs共2756个相同物理时刻的标量比较，18个旧完整FP64状态直接比较，无插值、旧域重跑或轨迹修正。每步失效数、接触数、dt相同；全部可用关键态失效ID集合一致。坐标最大分量差{h["keyframe_max_absolute_errors"]["x"]:.6g}m、速度{h["keyframe_max_absolute_errors"]["velocity"]:.6g}m/s、力{h["keyframe_max_absolute_errors"]["force"]:.6g}N、d差{h["keyframe_max_absolute_errors"]["d"]:.6g}，属于浮点量级。首软化85µs/键3393、首完全失效92.5µs/7键一致。旧逐步x/v和全量个体失效事件未保存，不能声称对它们逐项比较。'),
('CPU/内存/墙钟？',f'主积分及在线核查墙钟{r["main_wall_s"]:.6f}s、CPU {r["main_CPU_s"]:.6f}s，峰值RSS {r["peak_RSS_bytes"]/1024**2:.3f}MiB。外部父进程记录含导入与末尾压缩的墙钟{r["external_process"]["wall_s"]:.6f}s、CPU {r["external_process"]["user_CPU_s"]+r["external_process"]["system_CPU_s"]:.6f}s。1 MPI×1线程，无GPU、无重新编译。原生候选积分{r["native_candidate_wall_s"]:.6f}s；其余包括独立科学检查和持久化，不能将主耗时全部称原生积分成本。'),
('显著接触重叠仍存在？',f'是。最大{c["max_overlap_dx"]:.12f}dx≈{c["max_overlap_m"]*1e6:.6f}µm；末态{c["final_overlap_dx"]:.9f}dx。保持原数值接触核算并不认证软罚力的物理不可穿透性。'),
('可进入小型非均匀裂纹研究？','可作为下一轮受控研究的准入基础；本轮不执行新的试样，不自动部署V3或扩大计算。'),
('为何不能报告实际Gc？','完整动力学准入不提供已认证新生裂面面积、全局时间/空间收敛或与加载及接触无关的裂纹耗散；D包含全试样分布损伤，不能直接除以任意截面。宏观强度和实验血栓材料也未标定。历史G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²是原瞬断LPS理想标称值，不是新软化模型的实际Gc。'),
('WSL文档实际更新读回？',f'是。实际路径/home/lzy/projects/clot_rupture.md，before/after/diff及读回记录完整；修改前SHA {doc["before_sha256"]}，修改后SHA {doc["after_sha256"]}；移除本轮插入段可逐字恢复全部原文。' if doc else '待本轮实际WSL写入及读回后生成最终报告。'),
('下一轮唯一建议？','LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1。先限定小型非均匀裂纹自然演化与面积/耗能可辨识性，不自动增加材料组合、完整V3、CFD、微泡或AI任务。')]
# Derive singleton IDs from actual topology, never from recollection.
singletons=[q['atom_ids'][0] for q in co['components'] if q['node_count']==1];answers[13]=(answers[13][0],answers[13][1].replace('73、283','、'.join(map(str,singletons))))
body=f'''# LPS_NATIVE_SOFTENING_DOMAIN_RETENTION_AND_FAILURE_SAFE_OUTPUT_V1

Run `{RUN}` · 2026-10-14 · 唯一正式工况 · 当前服务器实际隔离库与源文件为准。

**结论：保持原材料、本构、加载协议与事件分辨算法不变，扩大固定数值域并改进异常持久化后，LAMMPS已安全完成341节点、700µs完整渐进软化动力学。**

| 准入项 | 结果 |
|---|---|
| DOMAIN_RETENTION_PASS | PASS，341原始ID全部保留 |
| FULL_700US_PASS | PASS，原始t=0至700µs，2804接受/4拒绝 |
| ENERGY_PASS | PASS，最大相对余额0.372426823%≤2% |
| TOPOLOGY_PASS | PASS，原生端点集合与独立重建一致 |
| 小型非均匀裂纹 | 下一轮可研究，本轮NOT_RUN |
| 材料Gc、宏观强度、全局时间收敛、接触物理精度、完整V3 | 均未认证 |

## 22项直接回答

'''+''.join(f'**{i}. {q}** {answer}\n\n' for i,(q,answer) in enumerate(answers,1))
body+='## 700µs完整能量账本\n\n保持 `W = ΔU + ΔK + D + ΔUc − Jc + R`，初始U/K/Uc均为0，D仍为同一新坐标的完整U_old−U_new。带符号材料拓扑−D不重复计入。\n\n| 项 | 末态 |\n|---|---:|\n'
for k in ['W_J','U_J','K_J','D_J','Uc_J','Jc_J']:body+=f'| {k[:-2]} | {e["final"][k]*1e9:.12f} nJ |\n'
body+=f'| R | {e["final"]["R_J"]*1e12:.12f} pJ |\n\n损伤力切换功累计{e["final"]["switch_work_J"]*1e12:.12f}pJ，仍为有限时间离散余额的主要来源；余额没有藏入D。独立分解最大余量{a["max_errors"]["decomposition_J"]:.6g}J。\n\n'
body+='## 原生改动和异常证据\n\n只修改Python输入的数值域及持久化驱动，没有C++改动或重新编译。复用库SHA：`'+lib['library_sha256']+'`。冻结驱动身份见DOMAIN_PROTOCOL_FREEZE.json；生产库、env.sh、独立模型及旧manifest均保护。\n\n新读data后、fix建立前执行不带remap的change_box。读取旧restart时先用同名fix恢复per-atom记录并run 0，然后扩大旧box并run 0；新域restart无需再次改box。这遵守锁定change_box.cpp中对restart per-atom fix数据的限制。实际输入和旧/新差分已保存。\n\n异常日志区分接受态、当前候选和最近原生检查点；每步节点流保留x/v/总力/材料力，损伤变化和实际接触对另存。最终滚动检查点时间为700µs。P1只做一次定向绑定异常测试，未故意重跑旧盒子到lost atom，也未开展故障矩阵。服务器缺少/usr/bin/time，改用外部Python父进程记录真实子进程资源；这一步未推进物理模拟。首次图3绘制因初态失效键集合为空而中断，仅修正空集合绘图逻辑，没有重算动力学。\n\n'
body+='## 文件、可视化与复现\n\n五组PNG/SVG及[图表入口](OPEN_RESULTS.html)均来自真实FP64数据。逐步节点流布局见full/STREAM_LAYOUT.json；完整关键态、原始事件、实际子步、checkpoint和全历史比较均在本目录。完整输出字节数以封存RUN_RESOURCES.json及PAYLOAD_MANIFEST.json为准。复现入口见REPRODUCE.txt，封存目录禁止重跑覆盖。\n\n对照资料：[LAMMPS change_box](https://docs.lammps.org/change_box.html)、[固定非周期边界](https://docs.lammps.org/boundary.html)、[restart](https://docs.lammps.org/read_restart.html)。实际兼容性依据本目录保存的锁定源码，不升级上游。\n\n本轮开发、验证、图表和文档交付后停止，等待科研审查。\n'
(O/('REVIEW_'+TASK.upper()+'.md')).write_text(body)
fragment=f'''
**本轮计算域保留与异常保存审查（{RUN}，2026-10-14）**：在保持341个实际节点、12047原始参考键、245自由点/96规定运动点、完整m/r0/V/taper、G/K/ρ、dx/δ、s₀=0.15、s_c=0.18、α=0、半余弦60µm加载与原事件子步算法不变的前提下，**首次从原始无损伤t=0安全完成700µs完整动力学**。复用上轮隔离库，SHA仍为`{lib['library_sha256']}`，没有C++修改、重新编译或生产环境替换。

旧fff数值域各轴±250µm，节点73在687.75µs的y=249.953987866µm，下一688µs漂移预测250.054230365µm；这是固定数值盒子删除节点，旧日志实报340/341。本轮根据旧18个完整态的坐标范围、全时程最大速度、余下12.25µs运动估计及100µm空白余量，冻结**各轴±400µm、boundary f f f**。不使用remap，不加墙或周期性；数值盒子不是物理血管壁。初态节点、原始键、质量、m及力/完整能量逐值一致。全程真实坐标包络最终Y最大254.865870326µm，新域最小余量145.134129674µm，341点全部保留。

| 本轮核查 | 实际结果与科学边界 |
|---|---|
| 四项准入 | **DOMAIN_RETENTION_PASS / FULL_700US_PASS / ENERGY_PASS / TOPOLOGY_PASS**全部通过；单独从t=0新运行，没有旧前缀拼接。 |
| 原生事件积分 | 2804接受、4次Δd拒绝并完整回滚；2796个0.25µs步、8个0.125µs步覆盖2800个基准区间。最小允许0.03125µs及Δd≤0.10门槛不变，实际最大Δd=0.09786763438279739。无κ/d回退、重复损伤提交或非有限值。 |
| 异常保存 | 真实683µs旧检查点上的一次绑定异常测试通过；完整合法态及原生restart读回一致，候选与接受态分开。每步刷新标量/事件/FP64节点流，两滚动槽每40接受步更新并保留原损伤检查点；异常处理不再gather或强写restart。硬终止只保留已有落盘数据，不宣称全故障零丢失。 |
| 新旧轨迹 | 对旧0–687.75µs的2756个标量时刻和18个完整态直接比较，dt/失效数/接触数及关键态失效ID一致；坐标最大差1.21973×10⁻¹⁹m、速度3.41914×10⁻¹⁵m/s、力3.76930×10⁻²⁰N，属于浮点量级。旧逐步x/v缺失，不虚构比较。 |
| 末态拓扑 | **4937失效、7110存活**，正常原生gather全部341点和存活端点，与独立原图重建完全一致；实际分量**192、147、1、1**，两个大分量最大点间距547.167303/433.336158µm，质量2.4000×10⁻⁸/1.8375×10⁻⁸kg。它们是离散图分量，尚非认证的真实血栓碎块。 |
| 接触 | 1.0dx及刚度不变；1818个接受步实际接触，峰值87对，末态67对。最大重叠仍0.246649196dx≈12.33246µm，接触物理精度未解决。 |
| 成本 | 1 MPI×1线程，主计算墙钟54.473s、CPU35.475s、峰值RSS659.688MiB；外部父进程含导入/压缩墙钟55.215s。无GPU、无额外物理工况。 |

700µs末态完整账本仍为 $W=\\Delta U+\\Delta K+D+\\Delta U_c-J_c+R$：**W=4.575844829448nJ、U=0.525683244103nJ、K=0.674093509273nJ、D=3.272865781378nJ、Uc=0.126299310063nJ、Jc=0、R=−23.097015368980pJ**。全程最大归一化余额**0.372426823%≤2%**；末态0.347156805%。D定义与累计尺度未改变，损伤力切换功约23.094596565pJ仍是主要数值余额来源；没有把余额补入材料耗能。

结论是当前341点试样通过完整动力学准入，可作为下一轮小型非均匀裂纹研究的基础；**全局时间收敛、宏观抗拉强度、材料Gc、真实新裂面与接触物理精度均未认证**。D不是已识别裂面单位面积耗能。历史 **G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²**仍为原瞬断LPS理想标称值，全部旧记录及173.156623/230.776817nJ保持。五组PNG/SVG见[本轮图表索引](temp_storage/clot_github_sync_20261011/clot_clearance/results/{TASK}/{RUN}/OPEN_RESULTS.html)，[22项中文报告](temp_storage/clot_github_sync_20261011/clot_clearance/verification/{TASK}/{RUN}/REVIEW_{TASK.upper()}.md)及原始FP64/事件可追溯。下一轮唯一建议**LPS_NATIVE_SOFTENING_SMALL_NONAFFINE_CRACK_VALIDATION_V1**；本轮不自动执行P4、V3、CFD、微泡或AI，交付后停止等待科研审查。

'''
(O/'WSL_DOCUMENT_INSERT.md').write_text(fragment)
repro=f'''This archive is a completed single-case experiment. Never rerun into {O}.
Server: ssh -p 15283 root@95.3.33.46
Environment: source /workspace/y_clot/env.sh
Isolated library: {LIB}
Library SHA256: {lib['library_sha256']}
Production library and environment were NOT changed. C++ recompiles: 0.

Actual execution sequence (already executed once):
export LD_LIBRARY_PATH={LIB.parent}:$LD_LIBRARY_PATH
export PYTHONPATH={PREV}/source/python:$PYTHONPATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python {A}/prepare.py
python {A}/launch.py preflight
python {A}/freeze.py
python {A}/launch.py full
python {A}/analyze.py
python {A}/figures.py
python {A}/report.py

Any future independent reproduction requires a NEW isolated driver copy and NEW RUN identifier in common.py. prepare.py and launch.py refuse existing outputs. Preserve all model/reference/gate identities from DOMAIN_PROTOCOL_FREEZE.json. Do not restart or concatenate the archived prefix and call it a fresh full path.
Postprocessing is read-only in relation to physical states. No simulation is needed to regenerate plots from existing outputs. report.py is rerun after the actual WSL document readback JSON has been synchronized.
Native output units and binary stream layout: full/STREAM_LAYOUT.json. Actual force is N; complete energy is J; material force-density normalization is inherited unchanged.
Latest normal native restart and its matching FP64 state: full/LATEST_NATIVE_CHECKPOINT.json (700 us). Original native build source/python wrapper is the prior immutable isolated archive, whose manifest was verified before this task. No new duplicate full source tree or library bundle was created.
'''
(O/'REPRODUCE.txt').write_text(repro)
print('REPORT',bool(doc),'SINGLETON_IDS',singletons)
