from common import *
import html,shutil
def main():
 f=json.loads((O/'PLOT_AND_REPORT_FACTS.json').read_text());st=f['status'];lo=f['local'];en=f['ledger']['final'];res=f['resources'];co=f['contact'];top=f['topology'];ident=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text())
 doc=json.loads((O/'WSL_DOCUMENT_UPDATE.json').read_text()) if (O/'WSL_DOCUMENT_UPDATE.json').exists() else {}
 st.update(damage_force_energy_numerical_consistency='PASS_ON_ACCEPTED_0_TO_687p75us_ONLY',failed_incomplete_trial_time_s=.000688,full_700us='FAILED_ATOM_LOSS',next_task='LPS_NATIVE_SOFTENING_DOMAIN_RETENTION_AND_FAILURE_SAFE_OUTPUT_V1',durable_output_limitations='No valid final restart, no full all-bond sparse event stream, exact online error maxima and full peak RSS unavailable after exception. Last accepted full FP64 and every accepted scalar row preserved.')
 dump(O/'FINAL_STATUS.json',st);f['status']=st;dump(O/'PLOT_AND_REPORT_FACTS.json',f)
 design='''# 事件分辨原生积分：设计、冻结与实际范围

冻结文件：EVENT_INTEGRATOR_FREEZE.json。只改变研究库中 fix_move::reset_dt；新增 event_state_audit.cpp 只读输出完整 special 邻接表。原生 bond_bpm_peri_variational、动态协调 fix、pair、fix_nve、Verlet 和自由能源码不变。编译器和锁定上游均沿用，不替换生产库。

每个基准区间长度 0.25 µs。真实原生候选步先完成第一半踢、位置漂移、解析夹持运动、POST_INTEGRATE κ/d 与拓扑提交、新力计算及 POST_FORCE 单位转换、第二半踢。Python读取实际候选Δd；超过0.10时整步拒绝。关闭试算实例，以接受态的完整原生restart重建实例和同名fix，再从恢复的x/v/F重新积分两段h/2。后半段使用前半段真正接受的新力和新历史。绝不剪裁d或将一次提交人为拆成多次。

采用只读的VV漂移与冻结损伤公式预测值≥0.075，或当前已处于子步时，才保存完整检查点。预测只控制I/O，不决定正式接受；实际Δd始终是唯一损伤门槛。预测坐标与候选坐标差须≤2e−18m；意外未保存检查点的越界只能停止。本次预测最大坐标误差为0。只有35个完整运行检查点，未在2755个接受步逐步写restart。

拒绝后核对所有点的x/v/F、材料力、κ/d、U/D/lastdrop/事件数/提交数/token、K/Uc、步号、物理时间、实际存活键端点及完整1–2/1–3/1–4排除表。4次完整运行回滚和4次局部回滚全部逐值一致。force-evaluation/cache诊断计数在重建实例时重新开始，它们不属于材料历史；全局trial编号和拒绝计数由驱动单独持久记录。Python账本只在候选接受之后累积，拒绝候选不进入W/D/Jc/R。

二分下限预注册为0.03125 µs；实际仅用到0.125 µs，没有测试/宣称更小步长。每个基准区间通过队列完整覆盖后才恢复0.25 µs；覆盖误差容差1e−15 s，仅处理浮点累计舍入。实际4个相邻基准区间二分，640.25–641.25 µs有8个真实半步；随后恢复基准步长。

锁定版本 Input::timestep 先update_time再改变dt并调用各fix的reset_dt；Update保存atime/atimestep，restart序列化ATIME/ATIMESTEP/TIMESTEP/NTIMESTEP。新fix_move允许六个完全规定的equal变量使用可变dt，其他样式仍拒绝。边界变量读取真实time；没有step×0.25µs换算。原生fix_nve的reset_dt同步两次半踢系数。

物理账本仍为 W=ΔU+ΔK+D+ΔUc−Jc+R。W逐接受步用真实夹持Δx的梯形反力功，加夹持动能变化；动量/VV等式使用实际h；D用同一新坐标的U_old−U_new。负D仅显示，不重复累计。力切换误差独立列出。归一化累计尺度与原轮完全相同，未借调分母掩盖误差。局部窗口仅为对比继承旧0–640.25µs标量账本；正式新路径从t=0重新积分，没有轨迹拼接。

局部真实restart来自最小2561步确定性重放；与旧640.25µs完整FP64逐值一致。局部候选也精确复现旧640.5µs现场。新完整路径在687.75µs之前通过逐步门槛；下一候选发生固定盒子丢点，停止而未扩大域、改材料、放宽门槛或继续计算。

本算法仅在观测到的4次损伤拒绝上认证完整回滚；不承诺任意原生异常可恢复。原生丢点异常暴露驱动finally仍调用完整gather/restart的缺陷，导致部分内存统计未落盘。未用NPZ伪造最终restart。此限制及丢点属于下一轮优先问题。

参考说明以本地锁定源码为准。官方文档用于核对一般接口：[NVE/Velocity-Verlet](https://docs.lammps.org/fix_nve.html)、[时间步命令](https://docs.lammps.org/timestep.html)、[规定运动](https://docs.lammps.org/fix_move.html)、[restart重建及同名fix](https://docs.lammps.org/read_restart.html)、[Fix接口](https://docs.lammps.org/Developer_write_fix.html)。标准[dt/reset](https://docs.lammps.org/fix_dt_reset.html)并不提供本项目的损伤增量事务准入，本轮没有采用它。
'''
 (O/'EVENT_INTEGRATOR_DESIGN.md').write_text(design)
 report=f'''# LPS_NATIVE_SOFTENING_EVENT_RESOLVED_INTEGRATION_V1

运行标识：{RUN}；审查日期：2026-10-14（Europe/Zurich）。

**结论：局部事件分辨与完整回滚通过；完整700µs协议未通过。** 原0.11363469损伤跳跃已用真实0.125µs子步分辨，材料、自由能及0.10门槛均未改变。唯一新完整路径最后接受到687.75µs；下一688µs候选在固定计算盒子发生丢点。不能写成可信完整700µs，也不准入P4。本轮按停止条件结束，没有扩大域续算。

## 24项回答

1. **实现方式。** Trial→Check→Accept/Reject→Subdivide。候选由真正LAMMPS Velocity-Verlet产生；只读预测器选择检查点，实际原生Δd决定接受。0.25µs基准，二分至预注册下限0.03125µs，实际最小0.125µs。

2. **完整回滚与重积分。** 是。原生restart恢复x/v/F/κ/d/失效图/special排除/能量/事件计数/时间；重新创建同名fix。局部及完整运行各4次回滚，全部物理数值逐值一致、端点和排除表精确一致。每段真正重新漂移与半踢，不是分拆损伤提交。

3. **累计时间边界。** 是。仅隔离fix_move的六equal变量规定运动允许reset_dt，变量仍读累计time。最大已保存关键态夹持坐标差{co['max_loading_x_error_m']:.4e}m、速度差{co['max_loading_v_error_m_s']:.4e}m/s。每个接受步在线检查；全时程逐步精确最大夹持误差因异常收尾未落盘，不能把关键态最大值冒充全部步最大值。

4. **bond1813演化。** 节点24/245。640.25µs κ=0.15000702139906827，d=0.0002808428166676996。旧候选640.5µs Δd=0.11363469187762688精确复现后被拒绝。两个真实半步Δd为0.05736735957307094与0.056267200299590134；640.5µs新κ=0.15290300115822433、d=0.11391540268932877。之后三个相邻基准步也因同一键拒绝，均二分；该键在643µs实际达到完全失效。

5. **拒绝污染。** 未观察到。局部/正式共8次事务回滚的物理状态、损伤事件数、提交token、实际键和排除表均精确恢复。候选临时提交仅属于随后丢弃的分支；正式标量账本不累计被拒绝分支。诊断求力次数重启是已声明行为，不能混作物理事件数。

6. **最小实际dt。** 0.125µs；0.0625与0.03125µs未使用、未认证。正式时程有8个0.125µs接受步，其余2747个为0.25µs。

7. **0.113635跃迁解决情况。** 解决了该局部损伤时间分辨问题。两段真实积分后的合计损伤仍约0.11363456，允许跨两个物理步累计；每段自身均满足0.10。坐标与速度也发生可测改变，不是仅将d数字分拆。

8. **局部力与能量改善。** 同一640.25→640.5µs窗口，力切换功从1.1231398166fJ降至0.5597610589fJ；余额增量从+0.4120480272fJ变为−0.1581792339fJ，绝对值下降约61.6%。VV缺陷绝对值约降为1/4。相对旧候选，实际末态全体节点最大坐标分量差7.72048pm、速度差1.55935×10⁻⁵m/s、总力差1.34497×10⁻¹¹N；两路径该时刻都4834条完全失效键。原生/独立局部材料力最大差1.05879×10⁻²¹N。该局部证据不证明全局收敛。

9. **是否完整新运行。** 确实从原始无损伤t=0新启动唯一正式路径，无旧前缀拼接；但没有完成700µs。2755个接受子步覆盖687.75µs；下一688µs候选中止。

10. **卸载与最终保持。** 0–300升载、300–350保持、350–650卸载均完成；650–687.75µs末保持完成37.75µs，尚缺12.25µs。650µs实际状态已存FP64并绘图，700µs未绘制。

11. **子步数。** 正式接受2755、损伤拒绝4、另有1次未能完成的丢点候选；共2760次原生积分尝试。局部接受15、拒绝4。开发阶段另一次t0首步token检查错误的启动试算保留在development_attempt_02，未伪称正式接受或额外完整路径。

12. **最大接受Δd。** {st['max_accepted_dd']:.14f}，小于0.10。局部与完整路径的最大值一致。原始0.25µs大候选均保留，不混入接受曲线。

13. **能量相对余额。** 接受区间最大{100*st['max_relative_energy_balance']:.9f}%，通过原2%门槛；末态R={en['R_J']*1e12:.9f}pJ，末态相对值{en['relative_R']*100:.9f}%。这是0–687.75µs的结果，不是完整700µs准入。

14. **损伤力切换误差。** 原局部事件约下降50.2%；完整累计切换功为{en['switch_work_J']*1e12:.9f}pJ。没有旧700µs固定步长对照，不能声称全局误差下降；旧路径仅到640.5µs。仍未证明全局时间收敛或断裂积分二阶精度。

15. **拓扑。** 首软化仍85µs/键3393；首次完全失效92.5µs/7键。687.75µs共有4918失效、7129存活；由保存的完整d及原始图重建连通分量为**192、147、1、1**，共4个。650µs为340+1。最终原生gather_bonds因为后续丢点未能执行，不能声称完成最终直接端点对照；每个接受步已核对总键数与原始ID。最终完整d/κ与图均保存。

16. **动态接触。** 逐接受步在线核对实际排除、接触力与Uc；每次拒绝后排除表逐值一致。1769个接受步有接触，最多87对，末态85对。已保存关键态独立接触力最大差{co['contact_force_error_N']:.4e}N、能量差{co['contact_energy_error_J']:.4e}J；局部每步的精确最大误差另存。接触源、刚度与1.0dx均不变。

17. **重叠限制。** 最大{co['max_penetration_dx']:.9f}dx，即{co['max_penetration_m']*1e6:.6f}µm，物理接触精度仍未认证。此问题与本次固定盒子删除自由点是不同问题，未同时标定接触。

18. **成本。** 正式区间等价2751个基准步，新增4个接受子步及4次损伤拒绝，正常完成的原生试算增加8次（约0.291%）；另1次异常未完成。检查点35次，写入约{res['transaction_full_wall_s']['checkpoint_wall_s']:.3f}s，回滚约{res['transaction_full_wall_s']['restore_wall_s']:.3f}s，预测约{res['transaction_full_wall_s']['predict_wall_s']:.3f}s。正式接受区间墙钟由锁文件/逐步CSV时间戳估算{res['full']['main_wall_s']:.3f}s；相对上轮每步耗时外推同区间约增加{(res['measured_full_to_estimated_fixed_ratio']-1)*100:.1f}%，不是另跑的固定步长对照。局部{res['local']['main_wall_s']:.3f}s，峰值{res['local']['peak_RSS_bytes']/1024**2:.2f}MiB；正式峰值RAM因异常收尾未保存，明确为UNAVAILABLE，不用局部RAM代替。交付总字节数见manifest与LOCAL_SYNC_RECEIPT。

19. **实现错误与新瓶颈。** 本构C++未发现新错误。锁定fix_move原本不允许变dt，已作限定修正。开发时发现Python不导出special表、初始token=-1的断言错误，均保留证据并修正。正式运行的新停止原因是盒子仍为±250µm、边界fff：节点73最后y=249.953987866µm、vy=0.400969997m/s、F=0；下一VV漂移预测y=250.054230365µm。原生日志同时报告340而非341点。由这些证据判断为固定域丢点，而非MPI重复ID。异常后完整gather与write_restart也失败，暴露收尾容错缺陷；逐步CSV、18个完整关键状态与事务数据已只读恢复，未伪造restart或重算填补缺失统计。

20. **材料与能量不变。** G=940Pa、K=100000/3Pa、ρ=1000kg/m³、dx=50µm、δ=150.5µm、s0=.15、sc=.18、α=0；341/12047原图、完整m、r0、V、taper和夹持组不变。bond/协调fix/pair/nve/Verlet源码SHA逐一一致。

21. **P4条件。** 未具备本轮要求的完整700µs数值条件，P4非均匀小裂纹仍NOT_RUN/NOT_ADMITTED。局部损伤分辨通过不能抵消全程原生身份失败。

22. **为什么不能报Gc。** 虽然出现多个图连通分量，尚无经认证的新生裂面面积、准静态传播功、接触精度与全时程收敛；D=3.257717822nJ不能直接除以任意名义面积宣称材料Gc。G0_THEORY=0.07170070833333J/m²、G0_DISCRETE=0.06969623734234J/m²仍属于旧瞬断理想标称量。

23. **WSL文档。** 已修改并实际读回/home/lzy/projects/clot_rupture.md，原全文及①–⑨保留。状态{doc.get('status','PENDING')}；修改后与读回SHA为`{doc.get('after_sha256','PENDING')}`。before/after/diff及修改前SHA见WSL_DOCUMENT_UPDATE.json与wsl_document/。文档明确写明700µs未通过。

24. **下一轮。** 优先LPS_NATIVE_SOFTENING_DOMAIN_RETENTION_AND_FAILURE_SAFE_OUTPUT_V1：明确保持材料和物理载荷，审查计算域如何容纳真实分离点，并使异常退出先保存可用接受态与资源/事件日志，再评估同协议准入。只有完整协议通过后再审查小型非均匀裂纹。本轮不自动执行。

## 末态能量与复现边界

|量|最后接受687.75µs|
|---|---:|
|W|{en['W_J']*1e9:.12f} nJ|
|U|{en['U_J']*1e9:.12f} nJ|
|K|{en['K_J']*1e9:.12f} nJ|
|D|{en['D_J']*1e9:.12f} nJ|
|Uc|{en['Uc_J']*1e9:.12f} nJ|
|Jc|0 J|
|R|{en['R_J']*1e12:.12f} pJ|

本轮库SHA：`{ident['new_sha256']}`；生产库保持原SHA，未修改env.sh。全部本轮材料、几何、门槛及源文件身份在冻结文件和SOURCE_AND_LIBRARY_IDENTITY.json。五组PNG/SVG见OPEN_RESULTS.html。正式全时程的全键稀疏事件流及精确在线误差最大值在异常发生时尚在内存、未成功写出，报告明确该数据缺口；已有逐步事件计数及关键时刻完整ID不受影响。

**对最终问题的直接回答：机制确实解决了原局部损伤跳跃，但没有使本模型完成可信的700µs加载—卸载协议。新的实际瓶颈是688µs固定计算域丢点及异常收尾，不能用局部成功替代全程通过。**
'''
 (O/'REVIEW_LPS_NATIVE_SOFTENING_EVENT_RESOLVED_INTEGRATION_V1.md').write_text(report)
 index='''<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>事件分辨原生积分结果</title><style>body{max-width:1400px;margin:32px auto;padding:0 24px;font:18px/1.6 system-ui;color:#203340;background:#fbfcfd}h1{font-size:30px}.status{padding:18px;border-left:6px solid #c34049;background:#fff0ed}figure{margin:40px 0;background:white;padding:16px;border:1px solid #ccd5dc}img{width:100%;height:auto}a{color:#08728c}figcaption{margin-top:12px}</style><h1>原生软化损伤事件分辨积分</h1><p>2026-10-14 · 341节点 / 12047原始键 · 1 MPI</p><div class="status"><b>局部事件 PASS；完整700µs FAILED_ATOM_LOSS。</b><br>2755个接受子步、4次损伤拒绝，最小0.125µs；最大Δd=0.09786763，接受区间最大能量余额0.37243%。最后接受687.75µs，下一688µs候选在固定盒子丢失节点73。材料与门槛均未改变。</div>'''
 names=[('Fig1_event_substeps','图1：接受/拒绝、局部二分与恢复基准步长'),('Fig2_bond1813_history','图2：真实bond1813与节点24/245演化'),('Fig3_actual_dynamic_protocol','图3：实际动力学至687.75µs；红色尾区未完成'),('Fig4_energy_ledger','图4：完整自由能账本与局部力切换误差'),('Fig5_actual_material_points','图5：固定相机下的真实341材料点；没有700µs伪造末态')]
 for name,label in names:index+=f'<figure><a href="{name}.png"><img src="{name}.png" alt="{label}"></a><figcaption>{label} · <a href="{name}.svg">SVG矢量图</a></figcaption></figure>'
 rel=f'../../../verification/{TASK}/{RUN}/'
 index+=f'<p><a href="{rel}REVIEW_LPS_NATIVE_SOFTENING_EVENT_RESOLVED_INTEGRATION_V1.md">24项中文报告</a> · <a href="{rel}FINAL_STATUS.json">状态JSON</a> · <a href="{rel}NATIVE_ATOM_LOSS_DIAGNOSIS.json">实际丢点诊断</a> · <a href="{rel}EVENT_INTEGRATOR_DESIGN.md">算法设计</a></p><p>未认证：完整700µs、全局时间收敛、接触物理精度、材料Gc及P4准入。正式峰值RAM及部分内存事件流未能在异常后写出，详见报告。</p></html>'
 (F/'OPEN_RESULTS.html').write_text(index)
 (O/'OPEN_RESULTS.html').write_text(f'<!doctype html><meta charset="utf-8"><title>结果索引</title><p><a href="../../../results/{TASK}/{RUN}/OPEN_RESULTS.html">打开五组科研结果图与中文报告</a></p>')
 reproduce=f'''READ-ONLY FINAL EVIDENCE: {O}
This run is sealed and failed at attempted 688us. Do not rerun in place.
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LD_LIBRARY_PATH={O}/build:$LD_LIBRARY_PATH
export PYTHONPATH={O}/source/python:$PYTHONPATH
Actual workflow: prepare_build.py (source clone + restricted fix_move patch), incremental read-only event_state_audit.cpp compilation, freeze.py, replay.py (2561 minimal steps), run_protocol.py local, local_admission.py, run_protocol.py full (one new accepted path, stops at fixed-box atom loss), recover_stop_evidence.py, analyze.py, figures.py, report.py.
The two pre-formal driver issues are archived in development_attempt_01/02 and DEVELOPMENT_CORRECTIONS.json. Changing only the initial-token assertion did not alter the local positive-step arithmetic. Final executed source hashes are in EVENT_INTEGRATOR_FREEZE.json.
An independent reproduction must use a NEW RUN id and NEW output directory in a copied common.py, same geometry/material and actual locked source/compiler; never edit a sealed output. No reproduction or continuation was performed after the atom-loss stop.
Latest replay input is genuine native accepted_640p25.restart, not visual NPZ. read_restart must recreate identical fix IDs. Final accepted687.75 FP64 is not a native restart and must not be used as if it were one.
Local state/trajectory and full scalar ledger are FP64. Full all-bond sparse event stream and exact global online maxima were not flushed after fatal native exception. recover_stop_evidence.py only processes already saved outputs; its next-drift prediction is explicitly not an accepted/native saved candidate.
New library SHA256: {ident['new_sha256']}
Old library SHA256: {ident['old_sha256']}
Build flags: CMAKE_BUILD_TYPE=Release BUILD_MPI=ON BUILD_OMP=ON BUILD_SHARED_LIBS=ON PKG_BPM=ON PKG_MOLECULE=ON MPI_CXX_COMPILER=/workspace/y_clot/deps/openmpi/bin/mpicxx; cmake --build build -j12. Same upstream 8de817dd79bfe4525d5d39246a212d833e6dee07. No production deployment.
'''
 (O/'REPRODUCE.txt').write_text(reproduce)
 print('REPORT_WRITTEN')
if __name__=='__main__':main()
