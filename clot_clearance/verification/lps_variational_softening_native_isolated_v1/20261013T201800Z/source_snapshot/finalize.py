from common import *
import os,sys,time,shutil,subprocess,datetime,re
start=time.perf_counter()
def read(n):return json.loads((O/n).read_text())
identity=read('SOURCE_AND_LIBRARY_IDENTITY.json');binary=read('ISOLATED_BINARY_IDENTITY.json')
assert sha(LIB)==binary['sha256']
protected=read('HISTORICAL_PROTECTION_BEFORE.json');bad=[]
for rel,v in protected.items():
 p=R/rel
 if not p.is_file() or p.stat().st_size!=v['bytes'] or sha(p)!=v['sha256']:bad.append(rel)
base_bad=[]
for p,v in identity['files'].items():
 if not Path(p).is_file() or sha(p)!=v['sha256']:base_bad.append(p)
assert not bad and not base_bad,(bad,base_bad)
base_manifest={}
for p in (OLD/'source_owner_fix/src').rglob('*'):
 if p.is_file():
  rel=p.relative_to(OLD/'source_owner_fix/src');v=dict(sha256=sha(p),bytes=p.stat().st_size);assert sha(O/'source/src'/rel)==v['sha256'];base_manifest[str(rel)]=v
assert len(base_manifest)==binary['original_src_files_verified']
dump(O/'BASE_SOURCE_TREE_MANIFEST.json',dict(base=str(OLD/'source_owner_fix/src'),files=base_manifest))
for n,h in binary['new_sources'].items():assert sha(A/n)==h and sha(O/'source/src/BPM'/n)==h
binding=identity['lammps_python_binding'];assert sha(binding['path'])==binding['sha256']
after=dict(status='PASS',files_checked=len(protected),bytes_checked=sum(v['bytes'] for v in protected.values()),historical_mismatches=bad,source_library_environment_files_checked=len(identity['files']),source_library_environment_mismatches=base_bad,old_source_files_unchanged=len(base_manifest),production_library_sha256=sha(R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so'),diagnostic_library_sha256=sha(OLD/'build_owner_fix/liblammps.so'),env_sha256=sha('/workspace/y_clot/env.sh'),binding_unchanged=True)
dump(O/'HISTORICAL_PROTECTION_AFTER.json',after)
# Freeze exact new code, retain baseline full source separately on the server.
(O/'source_snapshot').mkdir(exist_ok=True)
for p in A.iterdir():
 if p.is_file() and p.suffix in ('.py','.h','.cpp'):shutil.copy2(p,O/'source_snapshot'/p.name)
(O/'isolated_library').mkdir(exist_ok=True);shutil.copy2(LIB,O/'isolated_library/liblammps.so');assert sha(O/'isolated_library/liblammps.so')==binary['sha256']
ldd=subprocess.run(['ldd',str(LIB)],capture_output=True,text=True,check=True);(O/'ISOLATED_LIBRARY_LDD.txt').write_text(ldd.stdout);assert 'not found' not in ldd.stdout
visual=read('FIGURE_VISUAL_REVIEW.json');assert visual['status']=='PASS'
for v in visual['figures']:
 assert sha(F/(v['name']+'.png'))==v['png_sha256'] and sha(F/(v['name']+'.svg'))==v['svg_sha256']
npz_audit=[]
for p in sorted(O.rglob('*.npz')):
 if p.is_relative_to(O/'source') or p.is_relative_to(O/'build'):continue
 with np.load(p,allow_pickle=False) as z:
  for k in z.files:
   a=z[k]
   if a.dtype.kind=='f':assert a.dtype.itemsize==8 and np.isfinite(a).all(),(p,k)
  npz_audit.append(dict(file=str(p.relative_to(O)),arrays=len(z.files),sha256=sha(p),status='PASS_FINITE_FP64_FLOAT_ARRAYS'))
dump(O/'RAW_FP64_DATA_VERIFICATION.json',dict(status='PASS',files=npz_audit))
doc=read('wsl_document/DOCUMENT_UPDATE_AUDIT.json');assert doc['status']=='PASS'
assert sha(O/'wsl_document/clot_rupture.before.md')==doc['before_sha256'] and sha(O/'wsl_document/clot_rupture.after.md')==doc['after_sha256']
status=read('FINAL_STATUS.json');status.update(PRODUCTION_LIBRARY_UNCHANGED='PASS_SHA256_UNCHANGED',HISTORICAL_RESULTS_PROTECTED=f'PASS_{len(protected)}_FILE_SHA256',WSL_DOCUMENT_UPDATE='PASS_ACTUAL_WSL_READBACK_WITH_BACKUP_DIFF_SHA',FIGURE_VISUAL_REVIEW='PASS_5_PNG_REVIEWED_SVG_VALIDATED')
dump(O/'FINAL_STATUS.json',status)
checks=read('CHECK_RESOURCES.json');paths=read('PATH_RESOURCES.json');pa=read('NATIVE_LOADING_UNLOADING_RELOADING.json')['records'];mpi=read('MPI_DAMAGE_HISTORY_AUDIT.json');restart=read('RESTART_DAMAGE_HISTORY_AUDIT.json');contact=read('POSTBREAK_CONTACT_AUDIT.json');sort=read('ATOM_SORT_DAMAGE_IDENTITY.json');b1=read('implementation_history/pre_safety_cleanup/BUILD_RESOURCES.json');b2=read('BUILD_RESOURCES.json')
resources=dict(scope='Measured individual task stages; not elapsed conversational/development/transfer time. Native RSS is Python+LAMMPS process, not standalone kernel allocation.',initial_build=b1,final_incremental_build=b2,final_mechanical_checks=checks,final_controlled_paths=paths,final_two_MPI_wall_s=mpi['wall_s'],final_two_MPI_RSS_KiB_per_rank=mpi['RSS_KiB_per_rank'],restart_wall_s=restart['wall_s'],restart_RSS_KiB=restart['RSS_KiB'],contact_wall_s=contact['wall_s'],contact_RSS_KiB=contact['RSS_KiB'],final_stage_launches=read('FINAL_LAUNCH_RESOURCES.json'),first_identity_launches_including_not_exercised_MPI=read('IDENTITY_LAUNCH_RESOURCES.json'),first_paths=read('implementation_history/pre_safety_cleanup/PATH_RESOURCES.json'),figure_resources=read('FIGURE_RESOURCES.json'),initial_figure_resources=read('implementation_history/figure_layout_before_review/FIGURE_RESOURCES.json'),peak_build_sampled_CPU_cores=max(s['sampled_peak_CPU_cores'] for b in [b1,b2] for s in b['stages']),peak_build_sampled_process_tree_RSS_bytes=max(s['sampled_peak_process_tree_RSS_bytes'] for b in [b1,b2] for s in b['stages']),compile_job_limit=12,BLAS_OMP_threads=1,MPI_ranks=2,GPU_used=False,cgroup=identity['cgroup'],protocol_counts=dict(main_graphs=1,formal_sc_cases=2,states_per_case=354,total_full_paths_including_pre_cleanup_revalidation=4,mechanical_gate_batches=2,FD_energy_evaluations_per_batch=32,component_batches=2,sort_batches=2,MPI_batches=2,MPI_batches_with_actual_migration=1,restart_scenario=1,restart_native_instances=3,contact_scenario=1,contract_rejection_cases=4,CPP_correction_cycles=1),per_case_solver={n:dict(native_evaluations=v['native_evaluations'],endpoint_solves_per_evaluation_per_rank=341,total_endpoint_solves_per_rank=v['endpoint_star_solves_per_rank'],internal_solver_wall_s=v['native_internal_solver_wall_s'],mean_internal_solver_wall_ms=1e3*v['native_internal_solver_wall_s']/v['native_evaluations'],accepted_states=v['states'],wall_s=v['wall_s'],wall_ms_per_accepted_state_including_Gauss_work_and_output=1e3*v['wall_s']/v['states']) for n,v in pa.items()},large_runs=0,full_V3_runs=0,CFD_runs=0,AI_runs=0)
resources['server_storage_at_finalization']={name:int(subprocess.check_output(['du','-s','-B1',str(p)],text=True).split()[0]) for name,p in [('isolated_source',O/'source'),('build',O/'build'),('verification_including_source_build',O),('figures',F),('new_analysis',A)]}
resources['payload_byte_accounting']='Exact deliverable bytes and count in PAYLOAD_MANIFEST.json after seal; source/build full trees retained server-side, not duplicated in WSL.'
dump(O/'RUN_RESOURCES.json',resources)
repro='''LPS_VARIATIONAL_SOFTENING_NATIVE_ISOLATED_V1 — reproduction and immutable evidence
Accepted run: 20261013T201800Z. Never rerun into this accepted directory.
Server project is flat /workspace/y_clot/project. First read /etc/vast-agents-guide.md,
/workspace/y_clot/project/AGENTS.md and /workspace/y_clot/README_MIGRATION_ZH.txt.
Use the existing isolated clot environment; no global installation or env.sh edits.

1. Verify PAYLOAD_MANIFEST.json hashes (paths relative to project root),
ISOLATED_BINARY_IDENTITY.json, SOURCE_AND_LIBRARY_IDENTITY.json and the previous
accepted 127-file manifest. Baseline prior-run asset SHA is the hard admission gate.
Complete baseline source remains at the path in BASE_SOURCE_TREE_MANIFEST.json.
This artifact set intentionally does not copy the multi-package full source/build to WSL.
The final library is also preserved at isolated_library/liblammps.so.

2. Fresh reproduction: copy analysis/native_softening_v1 to a new analysis directory.
In that copy's common.py assign a NEW RUN and A to that NEW analysis directory.
Leave PREV/OLD and all material constants unchanged. prepare.py rejects an existing O.
The fresh-build wrapper tolerates absent SETUP_RECOVERY.json; the originally executed
wrapper is archived in implementation_history/build_recipe_before_portability_fix.
This wrapper-only portability cleanup does not change C++ or compiled-library identity.

3. In a dedicated shell (replace paths with the new run/copy):
source /workspace/y_clot/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python "$task_analysis/prepare.py"
python "$task_analysis/build.py"
export LD_LIBRARY_PATH="$task_output/build:$LD_LIBRARY_PATH"
export LAMMPS_LIBRARY_PATH="$task_output/build"
python "$task_analysis/checks.py"
python "$task_analysis/paths.py"
python "$task_analysis/components.py"
python "$task_analysis/identity.py" sort
mpirun --allow-run-as-root -np 2 python "$task_analysis/identity.py" mpi
python "$task_analysis/restart.py"
python "$task_analysis/contact.py"
python "$task_analysis/contract_checks.py"
python "$task_analysis/postprocess.py"
python "$task_analysis/figures.py"
Run each command only if its predecessor succeeded; use `set -e` in scripts.
The wrapper checks actual /proc/self/maps library SHA; do not fall back to production.
No Velocity-Verlet integration or physical dt is used: accepted-state kinematic protocol.

4. New-run document update, visual review and sealing require new real before/after
checks. Do not reuse this run's WSL/visual PASS record in another run. finalize.py
assembles this run's evidence and is not a universal fresh-run scientific driver.
Binary bitwise identity across another compiler/build path is not promised;
compare declared source/build identity and FP64 mechanical gates.

5. Native interface: new bpm/peri/variational + versioned reference.native;
`lpsvar_action`: op1 diagnostic state import before any accepted commit, op2 explicit
collective accepted-state commit between runs; repeated identical token is a no-op.
`lpsvar_output`: read-only FP64 force/energy/history. See native.py and implementation
audit for units (N/m3 native f, multiply V for N; complete ebond U/V, exported U in J).
All 2MPI API operations must be called collectively; this is not a production API.
'''
(O/'REPRODUCE.txt').write_text(repro)
u=read('UNDAMAGED_NATIVE_PARITY.json')['records'];partial=read('PARTIAL_DAMAGE_FORCE_ENERGY_PARITY.json')['records'];fd=read('FIXED_DAMAGE_GRADIENT_AUDIT.json');binarydiff=read('BINARY_LIMIT_DIFFERENCE.json');p18=pa['soft_sc018'];p24=pa['soft_sc024'];pro=read('OFFLINE_MODEL_PROVENANCE.json')
# The numbered questions below correspond one-to-one to the user's required 40.
qa=[
('是否完全复用了上一轮正式离线软化模型？',f"是。先核验上一轮 127 项封存资产，manifest SHA `{pro['manifest_sha256']}`；正式 model.py SHA `{pro['formal_model_sha256']}`。完整自由能、独立端点状态、双线性 κ→d、固定 m/V/r/taper 逐项继承。公式来自完整冻结源码，不来自摘要。"),
('是否避免恢复旧版两端共享内部状态的错误？','是。每一无向原始键共享 κ/d，但两个端点有各自的邻域内部变量 v。完全断开 A 块移动后，B 块非接触力变化严格为 0；未恢复上一轮初稿的共享 v。'),
('新原生 C++ 修改了哪些文件？','仅新增 `src/BPM/bond_bpm_peri_variational.h` 和 `.cpp`，新样式 `bpm/peri/variational`；原 4,729 个 src 文件 SHA 相同。真实新增源码在 `source_snapshot/`，完整新增补丁 `NATIVE_SOURCE_PATCH.diff`。Python 驱动在 `analysis/native_softening_v1/`，不属于生产路径。'),
('实际编译出的诊断库 SHA 是多少？',f"最终隔离库 SHA256：`{binary['sha256']}`，{binary['size']:,} 字节。服务器 `build/liblammps.so`，交付副本 `isolated_library/liblammps.so`。修正前库另存 implementation_history，不能混作最终版本。"),
('旧生产库是否保持不变？',f"是。结束重新计算 SHA：`{after['production_library_sha256']}`。旧双修复诊断库仍为 `{after['diagnostic_library_sha256']}`；env.sh 与 Python 绑定未改。"),
('未损伤原生力与旧 LPS 是否一致？',f"是。实际运行新库中原样保留的最终双修复 bpm/peri 与新样式，对照逐 atom ID 三维力；旧样式源文件逐字节与基线一致。参考态力 0；小仿射相对 L2 差 {u[1]['legacy_force_relative_L2']:.9e}，最大分量差 {u[1]['legacy_force_max_abs_N']:.9e} N；非仿射相对差 {u[2]['legacy_force_relative_L2']:.9e}。"),
('未损伤完整储能是否一致？',f"是，与独立冻结完整能量比较。小仿射绝对差 {u[1]['energy_abs_J']:.9e} J，非仿射 {u[2]['energy_abs_J']:.9e} J；参考态 U=0。本轮新 style 的完整 native tally 乘 V 后一致。旧 bpm/peri 不完整 PE 没有被用作完整储能参考。"),
('当前局部内部状态采用何种求解方法？','精确对角加秩一消元，等价于每个端点正定二次能量极小化；不是迭代优化、伪逆或经验缩放。6 个端点另以稠密 FP64 系统独立核对，最大 v 差 1.0587911841e−21 m。'),
('是否利用了低秩结构提高计算效率？','是。每个端点以标量 p 消去秩一耦合，每次求力仅对有向边做线性遍历。341 个 endpoint solves / 次 / rank；每条正式路径 2,131 次能量力求值、726,671 个端点求解。每次内部求解约 0.43 ms；这不包括全部 Python、通信和输出开销。'),
('部分损伤的原生力与离线参考最大差是多少？',f"同一实际坐标和历史：sc=.18 固定状态最大分量差 {partial[0]['force_max_abs_N']:.9e} N、完整 U 差 {partial[0]['energy_abs_J']:.9e} J；sc=.24 分别 {partial[1]['force_max_abs_N']:.9e} N、{partial[1]['energy_abs_J']:.9e} J。两个状态 κ/d 逐元素一致；另核对非均匀损伤。相对 L2 差分别 {partial[0]['force_relative_L2']:.9e}、{partial[1]['force_relative_L2']:.9e}。"),
('固定损伤下力是否等于自由能负梯度？',f"在已测试自由度内通过。8 个自由度、2 个步长、正负中心扰动共 32 次 native 能量计算；最大绝对差 {fd['max_error']['abs_N']:.9e} N，相对实际全局分量尺度差 {fd['max_error']['relative_global']:.9e}，真实分母 {fd['force_scale_N']:.9e} N。未设置人为分母下限，扰动期间提交次数 0。"),
('新模型是否保留固定参考 m 和原始积分权重？','是。341 个参考 m、V=dx³、原始 r 与 taper 均来自冻结点阵，损伤不会重新归一化 m。排序、MPI、restart 均核查固定 m；全程未换点阵或 horizon。'),
('是否正确处理 d=0、0<d<1 和 d=1？','是。d=0 恢复完整弹性，部分损伤使用冻结串联柔度，d=1 释放端点约束且该键非接触力严格为零；不加 epsilon 残余刚度。注意新模型二值极限不等于旧直接删键法，另见下方单列差异。'),
('κ 与 d 是否具有稳定原始键身份？','是。从 body.data 恢复真实 12,047 个原始 bond ID，以 bond ID + 两端 immutable atom tags + type + 零镜像标识；不使用当前 bond loop 下标。即使物理 bond slot 被删除，完整历史仍在版本化原始图中。重复端点对、周期/三斜几何明确不支持，重复身份已实测拒绝。'),
('是否有任何卸载愈合？','没有。已测试加载、卸载、压缩、再加载全路径 κ/d 不下降；完全失效 restart 与后续闭合也没有愈合。只对本轮已执行状态作此结论。'),
('是否出现负损伤耗散？','没有。记录原始有符号 U_old−U_new，不取绝对值、不截负值。两条正式路径最小损伤增量为 0 J，所有已测增量非负。'),
('原生 sc=0.18 和 0.24 的耗能分别是多少？',f"sc=.18：{p18['total_native_D_J']:.16e} J = {p18['total_native_D_J']*1e9:.12f} nJ；sc=.24：{p24['total_native_D_J']:.16e} J = {p24['total_native_D_J']*1e9:.12f} nJ。为同一运动学协议的离散累计耗能，不是单位面积材料 Gc。"),
('是否复现已接受离线耗能和广义反力？',f"是。累计 D 绝对差分别 {p18['native_offline_D_abs_J']:.9e}、{p24['native_offline_D_abs_J']:.9e} J。两个 sc 的 Q 峰值均 38.3625 mN，初始广义刚度 5115 N/m。旧瞬断 143.859377 nJ 仅复用已接受离线对照，未虚构原生部分损伤旧曲线。"),
('实际力学与功账本是否闭合？',f"是。沿冻结的分段仿射路径，以旧 d 两点 Gauss 外功积分，再在接受坐标提交损伤；W−ΔU−D 最大相对残差分别 {p18['work_balance_max_relative']:.9e}、{p24['work_balance_max_relative']:.9e}。无动能/阻尼项，没有把 run0 次数称为物理时间。连续解析积分与离散更新的历史差异没有被当作 C++ 错误。"),
('完全失效后断块之间是否仍存在虚假弹性耦合？','本轮未发现：2,470 个跨面键完全失效，A 块独立位移时对应非接触跨键力=0，B 块力变化=0。整体路径最终 12,047 个键全失效，原生 bonds=0、最终节点力和 U 均严格为 0。'),
('部分损伤键是否与接触重复施力？',f"未出现。原生 1–2 special exclusion 在 0<d<1 时保留。部分损伤闭合的附加 pair 力最大 {contact['partial_extra_pair_force_max_N']:.9e} N，为 FP64 相减舍入，不能写成数学上的严格零。"),
('完全失效后的原生接触是否实际经过验证？',f"是。完全失效提交移除原生键与对应 1–2 屏蔽，正常邻居重建后闭合出现 49 对真实接触；保留原 1.0dx 罚力。逐物理对仅计一次的参考与 native 最大分量差 {contact['contact_force_error_max_N']:.9e} N，相对差 {contact['contact_force_relative']:.9e}；闭合 κ/d 不变。大开口只为触发诊断断面，不代表有效材料应变区间，未接入壁面/V3。"),
('是否真实触发原子排序？','是。标准 atom_modify sort + run0 触发 341 个 owned 数组槽非恒等改变，保存排序前后实际 ID 序列。未使用手动直调 Atom::sort。'),
('排序前后损伤历史是否一致？','是。原始 ID、κ、d、m、节点力、完整 U 的排序差均为 0。排序后新接受状态产生预期 217 个失效事件且只提交一次。'),
('是否真实触发 MPI ownership 迁移？','是。最终 2MPI 分区实际迁移 49 个原子，ghost 数从 [145,196] 变为 [194,147]。第一次分区位于同一空隙，没有迁移，明确保留 MPI_MIGRATION_NOT_EXERCISED，不计通过。最终分区 −40 µm 按原坐标层位置确定。'),
('跨 rank 损伤更新是否正确？','本轮 2MPI 用例通过。2,470 条跨 rank 原始键，更新前 49 条部分损伤，更新后 217 条跨 rank 失效；所有 rank 历史/账本一致，一个 accepted commit 对应 217 个唯一事件。实际单 rank 与两 rank 的力及 U 差为 0，迁移前后 κ/d/m/F/U 也均相同。每 rank 复制小图并只对 owned atom 施力，不能外推生产伸缩性。'),
('restart 后是否保持完整 κ、d 及参考 m？','是。frame91 的 12,047 个部分损伤键写实际 LAMMPS restart，重新构造实例读回后 κ/d/m/IDs/F/U/D 均一致。继续后续 262 个接受状态，与不中断原生路径的力、U、d、D 差均为 0。完全失效 restart 也实测无愈合。'),
('是否重复累计历史耗能或失效事件？','没有。两个主路径最终各 354 次 accepted commit、12,047 次首次失效；每条路径 4 个重复 token 检查无新增 D/事件。restart 重放已保存 token 无新增，排序/MPI/接触重复提交也不增加。相同 token 配不同坐标返回 −5 并保持历史；旧 style 接口返回 −99。'),
('是否发现新的原生 Bug？','代码审查发现第一版一个未被科学路径调用的辅助 lpsvar_sort 入口不安全：直接 Atom::sort 绕过 exchange/borders 生命周期且缺少 sortfreq 条件。已用唯一一次最小 C++ 修正删除，数学核未改，最终库重做必要门槛与两路径。未发现力学核差错。首次 MPI 未覆盖迁移是诊断分区问题，修正并保留记录。拷贝符号链接和缺少 /usr/bin/time 的准备错误另存 SETUP_RECOVERY，不伪装成力学失败。'),
('本轮是否修改材料 G/K/s0/sc？','没有。G=940 Pa、K=100000/3 Pa、ρ=1000 kg/m³、dx=50 µm、horizon=150.5 µm、alpha=0、s0=.15、sc=.18/.24 固定。未扫描、拟合或扩大 sc。'),
('当前是否能保证宏观峰值强度不变？','不能。只证明同一仿射受控协议下广义 Q 峰值相同；38.3625 mN 不是宏观抗拉强度，5115 N/m 不是材料模量。'),
('是否能够报告实际材料 Gc？','不能。本轮 D 单位是 J，尚无自然裂纹传播、有效裂面及几何/尺度收敛认证。material_Gc_J_m2=null。'),
('能否说明 6–8 J/m² 目标已经实现？','不能。没有实现或认证实验断裂阻力，也未通过调参追求该值。历史 G0_THEORY=0.07170070833333 J/m²、G0_DISCRETE=0.06969623734234 J/m² 保留其历史语境，不能替代本轮材料 Gc。'),
('实际新增多少正式小型测试？','只有一个真实 341 节点主图、两个预注册 sc 的正式 354 状态路径。机械门槛组包括 3 未损伤状态、历史 5 键冻结状态、2 部分损伤状态、1 非均匀固定损伤状态、6 个端点稠密核对和 32 次差分能量求值；随后各一个断块、排序、2MPI、restart、接触场景，另 4 项接口拒绝检查。因唯一安全修正，机械门槛/两条路径/断块/排序各执行两批，完整主路径实际共 4 条；MPI 执行 2 次，其中仅 1 次真正迁移通过。restart 1 场景含 3 个原生实例。未把这些求值计成不同材料试验。'),
('CPU、内存、墙钟和文件成本是多少？',f"首次配置+编译 {b1['wall_s']:.6f} s（编译段 {b1['stages'][1]['wall_s']:.6f} s），最终增量 {b2['wall_s']:.6f} s。编译限 12 jobs，采样峰值 {resources['peak_build_sampled_CPU_cores']:.3f} 核、进程树 RSS {resources['peak_build_sampled_process_tree_RSS_bytes']/2**20:.2f} MiB。最终单 rank 门槛 {checks['wall_s']:.6f} s、两路径 {paths['wall_s']:.6f} s（CPU {paths['CPU_s']:.6f} s、RSS {paths['RSS_KiB']/1024:.2f} MiB）；2MPI {mpi['wall_s']:.6f} s，每 rank 峰值 RSS {mpi['RSS_KiB_per_rank']} KiB；restart {restart['wall_s']:.6f} s，接触 {contact['wall_s']:.6f} s。GPU 未使用。逐段完整成本、首次复验、服务器源码/build 磁盘占用在 RUN_RESOURCES.json；精确交付字节数在 PAYLOAD_MANIFEST.json。编译 CPU/RSS 为每 0.5 s 采样，不声称捕获短瞬态；MPI 各 rank 峰值不能相加称同时峰值；开发、传输耗时未计入上述求解墙钟。"),
('生成哪些实际科研可视化？','5 组各 PNG/SVG：原生/离线逐 ID 力与 U；6 状态原生点阵损伤；Q/U/W/D/残差；原始跨面键、部分损伤、真实分离、实际闭合接触；排序/MPI/restart 历史、力/能量误差及成本。OPEN_RESULTS.html 静态浏览。所有新曲线源自保存的实际 C++ 状态，旧瞬断离线对照显式标注。5 张 PNG 已实际打开审图，调整标签重叠后复核。'),
('生产原生库及历史数据是否受到保护？',f"是。结束核验 {len(protected):,} 个历史文件、{after['bytes_checked']:,} 字节，SHA/字节数零变化；另 {len(identity['files'])} 个源/库/环境资产与 4,729 个原源码文件一致。无 reset/clean/stash/commit/push、无 rsync --delete、无生产部署，最终 Git 分支/HEAD 仍单列记录。"),
('WSL clot_rupture.md 是否实际更新？','是，实际 `/home/lzy/projects/clot_rupture.md` 已读取、备份、并发 SHA 核验、插入科研结果并读回一致。保留所有旧文字、章节、参数及离线历史结论。'),
('备份、diff 和 SHA 在哪里？',f"本次结果目录 `wsl_document/`：clot_rupture.before.md、clot_rupture.after.md、clot_rupture.diff、DOCUMENT_UPDATE_AUDIT.json。before SHA `{doc['before_sha256']}`；after SHA `{doc['after_sha256']}`。包含实际 WSL 路径和读回检查，不只是待应用补丁。"),
('是否具备进入下一轮非均匀裂纹验证的条件？','具备进入隔离的小型非均匀裂纹研究的条件，非生产部署准入。下一轮需另行研究非仿射裂纹演化、宏观强度、真实分离能与尺度/时间层收敛。本轮到此停止，未自动启动下一轮、V3、CFD、微泡或 AI。')]
assert len(qa)==40
report=f'''# LPS_VARIATIONAL_SOFTENING_NATIVE_ISOLATED_V1 科研验收报告

Run ID：`{RUN}`。服务器项目：`{R}`。现场分支 `{identity['git_branch']['stdout'].strip()}`，HEAD `{identity['git_head']['stdout'].strip()}`。

**核心结论：已经在真正的 LAMMPS 原生 C++ 中复现上一轮冻结模型的节点力、完整储能、损伤历史及不可逆耗能。该结论限于 341 节点 / 12,047 原始键的规定运动学研究与本轮有限身份/接触诊断，不能推广为 V3、宏观强度或实验材料 Gc 认证。**

## 主要数值证据

| 同一规定加载协议 | sc=0.18 | sc=0.24 |
|---|---:|---:|
| 实际 native 接受状态数 | 354 | 354 |
| 最终不可逆 D / nJ | {p18['total_native_D_J']*1e9:.12f} | {p24['total_native_D_J']*1e9:.12f} |
| 与冻结离线 D 的绝对差 / J | {p18['native_offline_D_abs_J']:.9e} | {p24['native_offline_D_abs_J']:.9e} |
| 全路径节点力最大分量差 / N | {p18['max_force_abs_N']:.9e} | {p24['max_force_abs_N']:.9e} |
| 全路径完整 U 最大差 / J | {p18['max_energy_abs_J']:.9e} | {p24['max_energy_abs_J']:.9e} |
| κ/d 回退 | 无 | 无 |
| 最小有符号耗散增量 / J | 0 | 0 |
| 功账本最大相对残差 | {p18['work_balance_max_relative']:.9e} | {p24['work_balance_max_relative']:.9e} |
| 最終有效原生键 / 首次失效事件 | 0 / 12047 | 0 / 12047 |

关键独立证据：固定损伤梯度相对差 {fd['max_error']['relative_global']:.9e}；完全失效跨块响应严格 0；真实排序 341 槽；真实 2MPI 迁移 49 atom；restart 后 262 帧继续一致；闭合产生 49 对原生接触。

## 冻结模型与实施边界

完整推导及变量对应见 [C++ 审查](VARIATIONAL_CPP_IMPLEMENTATION_AUDIT.md)、[模型冻结](NATIVE_SOFTENING_DESIGN_FREEZE.json)、[逐文件身份](SOURCE_AND_LIBRARY_IDENTITY.json)。新 style 是对最终双修复源码的两文件独立扩展。使用低秩精确局部消元、独立端点内部状态、共享不可逆 κ/d、固定参考 m 和积分权重。损伤在显式接受 token、固定真实坐标处一次性提交；trial/run0/差分不改变历史。完整 U 与力出自相同变分能量。

所有新力由真实 LAMMPS atom force 按已固定体积换算 N，完整能量导出 J；旧 BPM force-density 和 bond tally 约定明确区分。全节点规定坐标协议不是自由动力学或准静态平衡，无物理 dt、动能和阻尼项。没有使用优化器求值次数定义物理时间。

历史五键断裂：新模型对自身冻结参考正确，但相对旧删键 LPS 的节点力向量 L2 差为 **{binarydiff['old_vs_new_force_difference_N']:.12e} N**。旧保存原生力与独立旧二值参考差 {binarydiff['old_binary_definition_preserved_by_offline_baseline']:.9e} N。这是原始冻结新自由能对失效端点内部约束松弛的本构差异，不是 sc 调整或移植误差；未强行宣称二者相等。

状态按稳定原始键身份存储；每 rank 复制整个小图、汇集 owned 坐标并仅向 owned atom 施力。因此通过 2MPI 的实际迁移，不代表大规模 MPI 伸缩性。未实际触发原子 nmax 扩容，该项保留 NOT_EXERCISED；多体积、多类型、周期/重复端点、任意复杂 topology 及动态积分未认证。restart 使用独立 LPSVAR_RESTART_V1，错误版本实测拒绝。

## 40 项明确回答

'''+''.join(f'### {i}. {q}\n\n{a}\n\n' for i,(q,a) in enumerate(qa,1))+'''## 分项验收状态

不以单一 PASS 替代边界。下表与 FINAL_STATUS.json 一致：

| 验收项 | 结论 |
|---|---|
'''+''.join(f'| `{k}` | `{v}` |\n' for k,v in status.items())+'''
## 复现与审计入口

- [五组图](OPEN_RESULTS.html)、[FP64 数值比较](NATIVE_SOFTENING_COMPARISON.csv)、[资源记录](RUN_RESOURCES.json)。
- [复现步骤](REPRODUCE.txt)、[最终库身份](ISOLATED_BINARY_IDENTITY.json)、[新增源码补丁](NATIVE_SOURCE_PATCH.diff)。
- [唯一 C++ 安全修正](MINIMAL_REPAIR_RECORD.json)、[补丁](MINIMAL_REPAIR.diff)、`implementation_history/pre_safety_cleanup/` 保留初版库、源码、真实计算与首次未触发迁移的结果。
- [历史保护结束检查](HISTORICAL_PROTECTION_AFTER.json)、[文档更新](wsl_document/DOCUMENT_UPDATE_AUDIT.json)、[图像审阅](FIGURE_VISUAL_REVIEW.json)。
- `PAYLOAD_MANIFEST.json` 封存可交付文件；`DELIVERY_VERIFICATION.json` 在 WSL 校验实际下载结果后生成。服务器保留完整隔离 source/build，本地交付实际新增源码、最终库副本和必要证据，不同步整个冗余构建树。

本轮不自动进入下一轮研究。宏观峰值强度、材料 Gc、6–8 J/m² 目标、实验断裂验证、完整 V3 和 AI 材料卡均未认证。
'''
assert len(re.findall(r'^### \d+\.',report,re.M))==40
(O/'REVIEW_LPS_VARIATIONAL_SOFTENING_NATIVE_ISOLATED_V1.md').write_text(report)
endgit={n:subprocess.check_output(c,cwd=R,text=True).strip() for n,c in [('branch',['git','branch','--show-current']),('head',['git','rev-parse','HEAD']),('status',['git','status','--short'])]}
assert endgit['branch']==identity['git_branch']['stdout'].strip() and endgit['head']==identity['git_head']['stdout'].strip();dump(O/'ENDING_GIT_STATE.json',endgit)
dump(O/'FINALIZATION_AUDIT.json',dict(status='PASS',UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),wall_s=time.perf_counter()-start,report_questions=40,baseline_protected=True,source_snapshot_files=len(list((O/'source_snapshot').iterdir())),figure_reviewed=True,document_evidence_verified=True,build_wrapper_portability_fix='Only handle missing historical SETUP_RECOVERY on fresh replay; no C++/binary changes. Archived original executed recipe.',production_deployed=False,next_phase_started=False))
# Hash all delivered evidence and code, excluding large server-only source/build and seal receipts.
files={}
for root in [A,O,F]:
 for p in sorted(root.rglob('*')):
  rel=p.relative_to(R)
  if not p.is_file() or '__pycache__' in p.parts or 'upstream_review' in p.parts:continue
  if p.is_relative_to(O/'source') or p.is_relative_to(O/'build'):continue
  if p.name in ['PAYLOAD_MANIFEST.json','PAYLOAD_MANIFEST.sha256','DELIVERY_VERIFICATION.json','SERVER_DELIVERY_RECEIPT.json','finalize_launcher.log']:continue
  files[str(rel)]=dict(sha256=sha(p),bytes=p.stat().st_size)
manifest=dict(run_id=RUN,files=files,file_count=len(files),bytes=sum(v['bytes'] for v in files.values()),excluded='Raw full source/build trees retained server; caches and self/transfer receipts excluded',source_baseline_manifest_sha256=sha(O/'BASE_SOURCE_TREE_MANIFEST.json'),new_library_sha256=sha(LIB))
dump(O/'PAYLOAD_MANIFEST.json',manifest);(O/'PAYLOAD_MANIFEST.sha256').write_text(sha(O/'PAYLOAD_MANIFEST.json')+'  PAYLOAD_MANIFEST.json\n')
print(json.dumps(dict(status='FINALIZED',files=len(files),bytes=manifest['bytes'],protected=len(protected),manifest_sha256=sha(O/'PAYLOAD_MANIFEST.json'),library_sha256=sha(LIB)),indent=2),flush=True)
