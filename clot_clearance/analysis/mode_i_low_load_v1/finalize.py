from common import *
import time,resource,csv
T=time.perf_counter()
def read(name):return json.loads((O/name).read_text())
a=read('STATIC_EQUILIBRIUM_AUDIT.json');n=read('NATIVE_FORCE_GRADIENT_CROSSCHECK.json');s=read('INDEPENDENT_EQUILIBRIUM_CANDIDATE.json');m=read('MAIN_SOLVER_EXIT.json');rec=read('RECOVERY_SOLVER_EXIT.json');imp=read('IMPLEMENTATION_CHECKS.json');r=read('NOTCH_VS_LOADING_END_STRETCH.json');di=read('UNACCEPTED_DISPLACEMENT_DIAGNOSTIC.json');pf=read('PREFLIGHT.json');doc=read('WSL_DOCUMENT_UPDATE.json');tangent=read('READ_ONLY_TANGENT_AUDIT.json');old=read('before_recovery_INDEPENDENT_EQUILIBRIUM_CANDIDATE.json')
before=read('HISTORICAL_PROTECTION_BEFORE.json');changed=[];missing=[]
for p,v in before.items():
 if not Path(p).is_file():missing.append(p)
 elif sha(p)!=v['sha256']:changed.append(p)
protection=dict(status='PASS' if not changed and not missing else 'FAIL',files_checked=len(before),changed=changed,missing=missing,production_library_sha256=sha(PROD),diagnostic_library_sha256=sha(B/'liblammps.so'),env_sha256=sha('/workspace/y_clot/env.sh'))
dump(O/'HISTORICAL_PROTECTION_AFTER.json',protection)
resources=dict(server_cpu_cgroup=pf['cgroup'],server_affinity=pf['affinity'],python=pf['python'],scipy_version=pf['scipy_version'],threads=1,GPU_used=False,native_MPI_ranks=1,main_solve=dict(wall_s=m['wall_s'],accepted_iterations=988,accepted_state_files_including_initial=989,energy_force_evaluations=m['evaluations'],CPU_user_s=old['CPU_user_s'],CPU_system_s=old['CPU_system_s'],peak_RSS_MiB=old['peak_RSS_KiB']/1024,exit='MAIN_WALL_BUDGET'),recovery=dict(wall_s=rec['wall_s'],invocations=1,accepted_corrections=0,HVP_calls=rec['HVP_calls'],CG_iterations=600,CPU_user_s=s['CPU_user_s'],CPU_system_s=s['CPU_system_s'],peak_RSS_MiB=s['peak_RSS_KiB']/1024,exit='CG_STOP_600'),native_check=dict(wall_s=n['wall_s'],peak_RSS_MiB=n['peak_RSS_KiB']/1024,coordinate_evaluations=3,reference_run0=1),recorded_energy_force_evaluations=imp['raw_energy_force_evaluations'],main_plus_recovery_wall_s=m['wall_s']+rec['wall_s'],solver_and_native_wall_s=m['wall_s']+rec['wall_s']+n['wall_s'],physical_cases=1,formal_dynamics=0,all_other_scans_or_physics=0,figure_resources=read('FIGURE_RESOURCES.json'),preparation_resources=read('PREPARATION_RESOURCES.json'),postprocess_resources=read('POSTPROCESS_RESOURCES.json'),wall_scope='Stage wall times exclude interactive development, SSH, hashing and WSL transfer. Main cutoff checked at evaluation boundaries (240.162 s actual); recovery remained within pre-reserved120 s. No expanded budget.',cost_scope='Full accepted coordinates and residuals, contact every evaluation, stretch diagnostics and compression cost more than old CSV-only 24.9s precheck. No measured timing split; no claim that entire increase is solver arithmetic.')
paths=[p for base in [O,F,R/'analysis/mode_i_low_load_v1'] for p in base.rglob('*') if p.is_file() and '__pycache__' not in p.parts];resources.update(output_snapshot_files=len(paths),output_snapshot_bytes=sum(p.stat().st_size for p in paths));dump(O/'RUN_RESOURCES.json',resources)
st={}
def status(name,val,scope,evidence):st[name]=dict(status=val,scope=scope,evidence=evidence)
for name in ['SPECIMEN_IDENTITY','PRISTINE_REFERENCE_M','PRECRACK_IDENTITY','PERIODIC_REFERENCE_TOPOLOGY']:status(name,'PASS','Exact original 16000 nodes, 909392 bonds, 16416 cuts, fixed full m and periodic images','SPECIMEN_IDENTITY_AUDIT.json')
status('LOW_LOAD_PHYSICAL_INPUT','PASS','Actual 2 Pa native force added once after one V conversion','LOW_LOAD_FORCE_APPLICATION_AUDIT.json')
status('LOAD_FORCE_MOMENT_CONSISTENCY','REFERENCE_LOAD_PASS_CURRENT_EQUILIBRIUM_FAIL','Reference total load/moment zero. Current residual torque fails; full-vector periodic scope and Y component reported.','STATIC_EQUILIBRIUM_AUDIT.json')
status('SOLVER_STOP_CAUSE','DIAGNOSED','Old20Pa ftol; new2Pa main wall budget, reserved recovery CG600 cap','SOLVER_FAILURE_ROOT_CAUSE_AUDIT.md')
status('RIGID_MODE_HANDLING','PASS_SCOPE_JUSTIFIED','Only 3 translations removed. Weak Y rotation is not removed. Final unprojected forces used.','RIGID_MODE_AND_GAUGE_AUDIT.json')
status('SCALING_AND_GRADIENT_VALIDITY','PASS_NUMERICAL_IDENTITY_ONLY','Consistent q/Pi scaling; native gradient and read-only HVP FD pass, not convergence','IMPLEMENTATION_CHECKS.json')
for name,key in [('ABSOLUTE_FORCE_CONVERGENCE','absolute_force'),('RELATIVE_FORCE_CONVERGENCE','relative_force'),('ENERGY_PLATEAU','energy_plateau'),('DISPLACEMENT_PLATEAU','displacement_plateau'),('SMALL_STRAIN_ADMISSION','small_strain')]:status(name,'PASS' if a['checks'][key] else 'FAIL','Pre-frozen gates; plateau refers to last3 accepted states. Last10 retained separately.','STATIC_EQUILIBRIUM_AUDIT.json')
status('NATIVE_FORCE_ENERGY_CONSISTENCY','PASS_FORCE_GRADIENT','Native force matches full independent gradient; native bond PE is incomplete and not used as U','NATIVE_FORCE_GRADIENT_CROSSCHECK.json')
status('NOTCH_LOADING_END_ELASTIC_RESPONSE','NOT_CERTIFIED_UNACCEPTED_STATE','Candidate ratio only; does not establish credible elastic ranking or first failure','NOTCH_VS_LOADING_END_STRETCH.json')
status('LOW_LOAD_SECANT_COMPLIANCE','UNAVAILABLE','No accepted loaded endpoint','LOW_LOAD_COMPLIANCE_AND_ENERGY.json')
status('LOADING_TIMESCALE_SCOPE','WAVE_SCALES_ONLY_QUASISTATIC_NOT_CERTIFIED','H/W crossing 1.031/5.157ms; no accepted static-shape slow mode, no24ms dynamics','LOADING_TIMESCALE_AUDIT.md')
status('NEW_BOND_BREAKS_ZERO','PASS_STATIC_FIXED_TOPOLOGY','No new cuts; fixed diagnostic break no, actual s00 remains0.15','NATIVE_FORCE_GRADIENT_CROSSCHECK.json')
status('STATIC_EQUILIBRIUM_ACCEPTED','FAIL','Absolute/relative force and torque not admitted; optimizer success not substituted','STATIC_EQUILIBRIUM_AUDIT.json')
status('GC_EFFECTIVE_NOT_MEASURED','NOT_MEASURED','No natural fracture, new area or loading-work path; null rather than0','LOW_LOAD_COMPLIANCE_AND_ENERGY.json')
for name in ['PRODUCTION_LIBRARY_UNCHANGED','HISTORICAL_RESULTS_PROTECTED']:status(name,protection['status'],f'{len(before)} protected files checked, no historical replacement','HISTORICAL_PROTECTION_AFTER.json')
status('WSL_DOCUMENT_UPDATE',doc['status'],'Actual WSL section⑤ negative evidence update, all old text retained, SHA readback','WSL_DOCUMENT_UPDATE.json')
next_task='对同一保存2 Pa状态预注册一次弱方向/预条件静态收敛专项，保持物理载荷和真实自由度，先通过原生残差验收；不再追加更大载荷或夹具试错。'
final=dict(task=C['task'],run_id=C['run'],outcome='STATIC_STATE_NOT_ACCEPTED_NUMERICAL_CONVERGENCE_LIMIT',credible_small_strain_static_baseline=False,credible_notch_over_loading_end_ranking=False,scientific_answer='2Pa已处于冻结的小伸长范围，原生载荷及梯度正确，但在主路径和一次有界恢复后仍未满足静态力与力矩门槛；不能认证弹性基线、柔度或Gc。',statuses=st,physical_load_cases=1,main_static_solves=1,reserved_recovery_invocations=1,accepted_recovery_steps=0,native_static_instances=1,formal_fracture_cases=0,new_natural_breaks=0,GC_EFFECTIVE=None,G_SEPARATION_PROTOCOL=None,G0_THEORY=.07170070833333332,G0_DISCRETE=.06969623734233767,production_changed=False,next_unique_priority=next_task,stopped=True,completed_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
dump(O/'FINAL_STATUS.json',final)
front=r['regions']['notch_front'];end=r['regions']['loading_ends'];U=s['last']['U_LPS_J'];W=s['last']['Fdotu_J']
qa=[
('是否使用原始16000点？','是。逐ID坐标及参考NPZ SHA ef145b6d6db6c17584dce6cbea1fa1b9d95847d71a69265a73df6334d1695c03一致；没有重建近似试样。'),
('原参考键、预切和m是否一致？',f"909392完整参考键，初始预切16416、存活892976；131196条参考跨Y缝键、38992存活跨z=0键。先完整初始化再预切，原生m逐位保持，相对独立值最大差{n['reference_m_max_relative']:.5g}。"),
('2 Pa节点力与面积正确吗？','上下各800真实原ID，每点面积dx²=2.5e-9m²，力±5e-9N；每面2e-6m²、合力±4e-6N。原生已施加，不只是写在独立势能里。'),
('合力、力矩自平衡吗？',f"参考外载合力/力矩为零。最终总残余力归一化{a['global_force_relative']:.6g}通过；当前残余力矩={a['current_residual_torque_N_m']} N·m，归一化{a['global_torque_relative']:.6g}高于1e-10。固定周期胞X/Z取矩有边界力矩限制；相容Y分量约−4.712e-18N·m也超过冻结2e-18N·m尺度，不能忽略。"),
('20 Pa提前停止直接原因？','旧L-BFGS-B满足ftol目标下降条件，而非绝对力收敛。旧1132函数评估日志不是接受步序列，线搜索及条件数不能补造。见专门只读诊断。'),
('平移零模或其他病态方向？',f"3个平移HVP精确为零。固定Y胞仅Y转动保持内部长度，但dead-load势不保持，未投影。保存2Pa未接受态的随机方向与Y转动切向Rayleigh商约0.1361855与1.11644e-7N/m，差约1.22e6倍；支持弱方向尺度差异，不是精确谱条件数或已定位全部病态原因。"),
('优化数值尺度？','q=u/(1e-6m)，Pi/(1e-11J)，梯度=-P(F)×1e5；ftol=0、gtol=1e-11，对应单分量1e-16N。缩放不自动改变Hessian条件数。'),
('是否规约，有无改变约束？','仅均匀质心平移规约，无固定点、支撑弹簧或转动投影。平移探针保留伸长、储能与开度；最终使用未投影原生节点力，规约反力为零。'),
('求解方法？','一次L-BFGS-B，从真实零载A0开始；Pi=完整固定m独立LPS U+现有保守接触U−ΣF·u。接触每次评估，未用不完整原生bond PE。'),
('需要恢复吗？','是。预留的一次解析HVP Newton-CG同状态恢复，共601 HVP调用；600次CG未返回收敛，因此没有线搜索接受修正，立即停止。'),
('新增几次静态求解？','1个物理2Pa工况，1次主路径、1次同状态恢复调用；1个原生静态实例含完整参考run0和3次无积分器坐标/力评价。正式动力断裂0。'),
('优化器为何停止？','本轮主路径在240s预算检查停止，没有SciPy成功返回；恢复因CG_STOP_600停止。原驱动错误地排除所有异常后的恢复，已只修正调度启用预先保留的120s额度；没有增加总预算或第二恢复。旧停止记录与所有坐标保留。'),
('最大自由力？',f"原生{a['native_max_force_N']:.12g}N=1.187075pN，是1e-13N门槛的约11.87倍，FAIL。"),
('相对残余通过吗？',f"{a['native_relative_force']:.12g}，高于1e-5，FAIL。"),
('能量和位移平台？',f"末3个不同接受态：相对能量跨度{a['energy_span_last3_relative']:.12g}≤1e-8，位置跨度{a['position_span_last3_m']:.12g}m≤1e-10m。窗口3在运行前已冻结，上轮20Pa无显式平台门槛；旧wide用10。额外报告末10能量跨度{a['energy_span_last10_relative']:.12g}未通过同一阈值；末10位置{a['position_span_last10_m']:.12g}m通过。不能用短窗口平台掩盖力学失败。"),
('最大伸长低于0.03吗？',f"是，{a['max_positive_stretch']:.12g}。本轮失败来自数值收敛，非已证实载荷超小应变范围。"),
('原生力与梯度一致？',f"相对差{n['native_vs_independent_force_relative']:.12g}≤1e-9；最大分量差{n['max_abs_internal_force_difference_N']:.12g}N。HVP与1nm方向差分相对差{imp['HVP_FD_relative']:.12g}≤1e-6。通过表示算子一致，不表示残差收敛。"),
('是否保留原始键与裂纹？','是，参考全键保留，固定拓扑诊断只屏蔽原16416预切键；活键892976准确，未从删键邻域重算m。'),
('新增断键？','0。实际所有活键类型1，s00=0.15不变；break no用于静态诊断，最大s同时低于0.03。无自然损伤事件。'),
('裂尖最大伸长？',f"仅未接受候选：{front['max_positive_stretch']:.12g}，s/s00={front['max_s_over_s00']:.9g}；原键{front['reference_bond_id']}，atom {front['atom_ids']}，参考中点{front['reference_midpoint_m']}m。"),
('加载端最大伸长？',f"仅未接受候选：{end['max_positive_stretch']:.12g}，s/s00={end['max_s_over_s00']:.9g}；原键{end['reference_bond_id']}，atom {end['atom_ids']}，参考中点{end['reference_midpoint_m']}m。"),
('可信低载态裂尖仍领先吗？',f"未建立可信低载平衡，不能回答为已认证领先。保存候选比值{r['notch_to_loading_end_max_stretch_ratio']:.9g}、约高3.42%，仅空间诊断；不预测真实动态首断。"),
('最大位移与平均开度？',f"仅未接受候选：最大{di['max_displacement_m']*1e6:.9f}µm，atom {di['max_displacement_atom_id']}；平均端面开度{di['mean_opening_m']*1e6:.9f}µm；平均位移模长{di['mean_displacement_magnitude_m']*1e6:.9f}µm。全部16000点u已保存，不认证为弹性基线。"),
('单点割线柔度？','UNAVAILABLE/null。由于受载态未验收，本轮不报告可信割线柔度，也不画已测连续载荷—开度曲线。'),
('独立储能增量？',f"零载U0=0；候选U及增量{U:.12g}J={U*1e12:.9f}pJ，仅未接受静态坐标的可回收能。"),
('终点外力位移收缩？',f"候选ΣF·u={W:.12g}J={W*1e12:.9f}pJ；Pi={s['last']['Pi_J']:.12g}J。不是路径外功。"),
('二者合理对应吗？','不作已验收线性弹性一致性结论。状态未平衡，不能把近似半F·u关系当作能量/力验收或材料实验验证。'),
('剪切传播尺度？','cs=sqrt940/1000=0.969535971m/s；H穿越1.031421ms，W穿越5.157106ms。'),
('24ms已准静态吗？','没有证据。24ms约4.654个W穿越，裂臂慢弯曲及裂纹传播尺度未认证；没有对未接受位移做已认证模态或静态柔度推断，也没有运行24ms动力学。'),
('为何不能报告Gc？','无自然扩展、新生面积和物理加载路径外功，且静态准入失败。G0_THEORY=0.07170070833333、G0_DISCRETE=0.06969623734234J/m²只保留为历史标称参考。实际Gc、分离协议能和带符号拓扑释放均null。'),
('需要改原生LPS或材料吗？','本轮没有证据要求修改。生产库、G/K/s00/alpha、horizon、点距、接触及显式theta语义均不变。失败不能直接归咎于材料太脆。'),
('资源和墙钟？',f"现场CPU cgroup={pf['cgroup']}；单进程/单线程，GPU未用。主路径{m['wall_s']:.6f}s、峰值{old['peak_RSS_KiB']/1024:.3f}MiB；恢复{rec['wall_s']:.6f}s；原生{n['wall_s']:.6f}s、峰值{n['peak_RSS_KiB']/1024:.3f}MiB。主路径记录1036次储能/力评估、988接受迭代；恢复1次额外状态评估和601HVP。末次显示评估计数从接受点续接少计2个旧试探，原日志保留，按真实1037行修正资源统计。全部989接受状态坐标/残余保留，最终规模以manifest为准。"),
('哪些科研图？','5组PNG/SVG：实际点阵/存活键/2Pa载荷，主求解与一次失败恢复，候选×1变形，分区伸长及未认证柔度，分项科学准入。所有数据来自实际状态，标STATIC_STATE_NOT_ACCEPTED；无医学外壳、虚构断裂或动画。'),
('保护历史及生产吗？',f"{protection['status']}，逐文件复核{len(before)}个保护资产；变化{len(changed)}、缺失{len(missing)}。未reset/clean/stash/commit/push/替换生产库或env.sh。"),
('WSL文档实际更新？',f"{doc['status']}。真实{doc['actual_path']}仅追加⑤负面证据和结果链接，①–⑨及原有历史文本全部保留。"),
('备份diff SHA在哪？',f"本run的wsl_document/保存before、after、逐行diff和审计。before {doc['before_sha256']}；after/readback {doc['after_sha256']}。WSL_DOCUMENT_UPDATE.json为入口。"),
('下一轮唯一优先项？',next_task+' 本轮已停止。')]
report=f'''# LPS_MODE_I_LOW_LOAD_ADMISSION_AND_SOLVER_RECOVERY_V1\n\nRun `{C['run']}`。\n\n**尚未建立可信的小变形静态弹性基线。** 2 Pa的实际载荷、原始试样、原生力/独立梯度及小伸长范围通过；主路径240秒和一次600步CG恢复后，绝对/相对残余力及冻结力矩门槛仍未通过。没有放宽门槛、重试另一优化器或启动断裂。\n\n| 核心量 | 实际值 | 裁定 |\n|---|---:|---|\n| 最大伸长 | 0.01537935146 | ≤0.03，通过范围检查 |\n| 原生最大残余力 | 1.187074951e-12 N | >1e-13，失败 |\n| 原生相对残余 | 2.390067922e-5 | >1e-5，失败 |\n| 原生/独立梯度差 | 7.239375529e-12 | ≤1e-9，通过算子核对 |\n| 正式动力断裂 | 0 | 未运行 |\n| 可信柔度 / 实际Gc | UNAVAILABLE | 不填0 |\n\n裂尖与加载端的保存候选伸长比约1.03421，但该状态未验收，**不能据此确认可信低载态中裂尖仍领先，也不能预测真实首断位置**。\n\n'''
for i,(question,answer) in enumerate(qa,1):report+=f'## {i}. {question}\n\n{answer}\n\n'
report+='## 逐项科学状态\n\n| 项目 | 状态 | 范围 |\n|---|---|---|\n'+''.join(f'| {k} | {v["status"]} | {v["scope"]} |\n' for k,v in st.items())
report+='\n## 复现与限制\n\n结果入口[OPEN_RESULTS.html](OPEN_RESULTS.html)，完整轨迹和残差见STATIC_SOLVER_TRACE.csv、accepted_optimizer_states/；没有物理时间轴。REPRODUCE.txt区分只读复核与需要下一轮授权的新求解。\n\n[SciPy官方L-BFGS-B停止条件](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html)支持旧ftol/gtol解释，实际算法和环境身份以本run源码与PREFLIGHT为准。未重新进行大规模论文检索，也未把公开文档当作本地原生实现的替代依据。\n'
(O/'REVIEW_LPS_MODE_I_LOW_LOAD_ADMISSION_AND_SOLVER_RECOVERY_V1.md').write_text(report)
(O/'REPRODUCE.txt').write_text('''This task is CLOSED. Do not rerun solve.py or recover_saved.py in this run directory.
Read-only review:
  source /workspace/y_clot/env.sh
  cd /workspace/y_clot/project
  python analysis/mode_i_low_load_v1/figures.py   # regeneration only, overwrites this run's figures
  Review FINAL_STATUS.json, STATIC_EQUILIBRIUM_AUDIT.json and the full report.
Historical execution order (not an instruction to launch more work):
  prepare.py -> frozen geometry/load/criteria/plan
  supervised solve.py, one main path from A0, stopped at main wall budget
  supervised recover_saved.py, one reserved recovery, CG_STOP_600, zero accepted steps
  LD_LIBRARY_PATH=<final-double-fix build>:$LD_LIBRARY_PATH native_check.py
  postprocess.py -> checks.py -> figures.py
  WSL update_document.py with before/diff/SHA -> finalize.py -> manifest and transfer checks
Environment threads: OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1.
Native library SHA: 00b246ebe4224e648c27b74afb64fc22b9cea4ba6ffb23f7c40d96c37ca22ae9.
The immutable old 20Pa state is read only. Every accepted new state has x and full residual in SI.
No integrator, natural fracture, external-work path, new crack area or Gc exists this turn.
Any further solve must use a new run ID and separately approved scientific protocol.
''')
dump(O/'SOURCE_INPUT_SHA256.json',{str(p):sha(p) for p in [*Path(__file__).parent.glob('*.py'),Path(__file__).with_name('config.json'),O/'LOW_LOAD_DESIGN_FREEZE.json',O/'STATIC_SOLVER_RECOVERY_PLAN.json',O/'reference_topology.npz',O/'load_identity.npz',O/'native_static/body.data',O/'native_static/in.lammps',O/'native_static/actual_load_assignments.in']})
print('FINALIZED',len(st),'independent statuses; history',protection['status'],flush=True)
