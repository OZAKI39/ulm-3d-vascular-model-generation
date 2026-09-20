"""Generate human handoff documents only from validated frozen manifests."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def read(p):
    return json.loads((ROOT/p).read_text())


def put(name, text):
    (ROOT/name).write_text(text.strip()+'\n')


def main():
    flow = read('frozen_reference/flow/flow_field_manifest.json')
    mesh = read('frozen_reference/mesh_manifest.json')
    boundary = read('frozen_reference/boundary_manifest.json')['boundaries']
    base = read('frozen_reference/baseline_summary.json')
    runtime = read('frozen_reference/run/runtime_manifest.json')
    meta = read('sync_metadata/main_base.json')
    physics, measured = base['physics'], base['measurement']
    stage = read('reports/sv1_3q/pytest_summary.json')
    historical_failures = '\n'.join('- `'+s+'`' for s in stage['historical_failures'])
    records = [('背景流', flow), ('正式体网格', mesh)] + [(k, v) for k,v in boundary.items()]
    table = '| 用途 | 相对此目录的路径 | SHA256 |\n|---|---|---|\n' + '\n'.join(
        f'| {label} | [{d["path"]}]({d["path"]}) | `{d["sha256"]}` |' for label,d in records)
    faces = '| 语义 | 实际 face ID | 面积（m²） |\n|---|---:|---:|\n'+'\n'.join(
        f'| {k} | {v["sv_face_id"]} | {v["area_m2"]:.16g} |' for k,v in boundary.items())
    freeze = '''Mesh convergence: **NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT**。

Time-step sensitivity: **NOT PERFORMED BY USER DECISION / NOT PLANNED BEFORE PARTICLE DEVELOPMENT**。

CPU/GPU field-level equivalence: **DEFERRED BY USER DECISION**。

这些是 USER-ACCEPTED PROJECT DECISIONS，不是通过的验证。FEM performance tuning: frozen after Stage Q。Stage Q is the frozen particle-development background-flow baseline；不能解释为网格收敛、时间步收敛、CPU/GPU 科学等价或新的 FEM production approval。'''
    put('PARTICLE_HANDOFF.md', f'''
# Frozen FEM → Particle-0 Handoff

FEM DEVELOPMENT: **FROZEN**。Frozen FEM performance baseline: **Stage SV1.3Q**，状态 **{base['stage_status']}**。

Next work: **Particle-0 → Particle-8**。从现在开始 FEM 是只读背景流；本分支完成冻结交接，尚未开始任何粒子代码开发。

## 项目已接受的验证边界

{freeze}

## 单向耦合

第一版仅允许 FEM → RBC / MB。RBC / MB → FEM 不允许，粒子运动不会重新求流场。不要重跑 FEM、调整网格、时间步、PETSc 或 ILU 来启动 Particle-0。

## Frozen Field Contract

以下路径全部相对于本文件所在目录，都是分支中的实际文件，无需原 WSL 或 GPU 服务器。

{table}

正式体网格：{mesh['nodes']} 节点，{mesh['tetra']} 四面体，每单元 {mesh['nodes_per_cell']} 个节点（VTK type {mesh['vtk_cell_type']}）。外表面见 [{read('frozen_reference/boundary_manifest.json')['exterior']['path']}]({read('frozen_reference/boundary_manifest.json')['exterior']['path']})；其哈希及全部边界见 [boundary manifest](frozen_reference/boundary_manifest.json)。

坐标使用原始 Cartesian x/y/z（右手系），单位 m；未平移、旋转或重新缩放，未赋予解剖轴标签。边界框 `[xmin,xmax,ymin,ymax,zmin,zmax]` = `{mesh['bounding_box_m']}` m。流场和正式网格节点坐标及顺序完全一致；逐行比较确认 tetra cell rows 也一一对应，差异只在单元内部的局部顶点排列（从 volume 到 flow 的排列为 `{mesh['flow_local_vertex_permutation_from_volume']}`）。`tetra_id` 约定为**正式 volume mesh 的零起始 cell index**；其 `GlobalElementID` 是从 1 开始。本冻结文件对的 cell index 对应关系已验证，未来替换文件时必须重新核对，不能只凭相同单元数推断。

冻结时刻为 step **{flow['step']}**，`TimeValue = {flow['TimeValue']:.17g}` s；step {flow['first_formal_steady_step']} 首次通过连续 5 区间稳态，完成当前时间步后在 step {flow['safe_stop_step']} 安全停止。背景流按稳态固定使用；该时间是 FEM 状态的生成时刻，不是粒子必须采用的初始时刻。

`Velocity`：POINT association，Float64，3 分量，m/s。`Pressure`：POINT association，Float64，1 分量，Pa。VTU 无 cell arrays。导出的是四节点四面体上的节点表示；这说明导出数据格式，不额外断言求解器所有内部空间的离散阶次。完整信息见 [flow manifest](frozen_reference/flow/flow_field_manifest.json)。

velocity gradient、vorticity、strain rate：**NOT STORED IN FROZEN VTU**。Particle-0 必须从 tetra velocity field 推导梯度，并明确单元内插值、跨面梯度和边界处的取值规则。本次未实现这些功能。

## FlowSampler contract（未来接口，尚未实现）

输入 `position[3]`：原坐标系中的 SI 位置，单位 m。

| 返回值 | 来源与约定 |
|---|---|
| `velocity[3]` | 由节点 `Velocity` 在所定位四面体内插值，m/s |
| `pressure` | 由节点 `Pressure` 插值，Pa |
| `velocity_gradient[3,3]` | Particle-0 推导，约定 G[i,j] = ∂u_i/∂x_j，s⁻¹ |
| `vorticity[3]` | Particle-0 推导 curl(u) = (G[2,1]−G[1,2], G[0,2]−G[2,0], G[1,0]−G[0,1])，s⁻¹ |
| `strain_rate[3,3]` | Particle-0 推导 (G+Gᵀ)/2，s⁻¹ |
| `inside_lumen` | Particle-0 的点定位结果；边界容差必须显式定义 |
| `tetra_id` | 正式体网格零起始 cell index；域外建议为 -1 |

接口建议：域外明确返回 `inside_lumen=false`，不得默默外推并伪装为域内结果；其余域外返回值及异常规则由 Particle-0 固定并测试。不得把节点数组读取成功等同于已完成 FEM sampling 验证。

## 单位与物理条件

geometry/length = m；time = s；velocity = m/s；pressure = Pa；density = kg/m³；dynamic viscosity = Pa s；Q = m³/s。

实际配置：ρ = {physics['density_kg_m3']:.16g} kg/m³，μ = {physics['dynamic_viscosity_pa_s']:.16g} Pa s，Q target = {physics['inlet_volume_flow_m3_s']:.17g} m³/s，dt = {base['dt_s']:.17g} s。此 Q 是参考数值条件，不是实验泵流量。运行 XML 的入口有符号值为 {base['BCs']['INLET']['Value']}，保留了既有入口离散归一化校正，不能替换为手抄的 Q target。

建议 3D Particle solver 内部核心计算全部保持 SI：length m、time s、velocity m/s、force N、torque N m、pressure Pa、viscosity Pa s。若显示 µm，只在 I/O / visualization 层显式用 `1 m = 10^6 µm` 转换，禁止隐式混用。

## Boundary semantics 与法向

{faces}

所有 accepted VTP 的三角形绕序已逐面核对：`n_out = cross(x1-x0,x2-x0)/|cross|` 指向流体域外，验证依据是所属四面体中心与面中心。入口 inward normal = `-n_out`。`∫u·n_out dA` 在入口为负、出口为正，正入口流量 Qin 是入口有符号积分的负值。壁面使用各三角形局部法向；面积加权平均法向不能代替局部壁面法向。

VTP 未存储 normals 数组，Particle-3 / Particle-7 应按上述约定计算；不要猜测方向。边界和外表面的三角形集合构成完整、不重叠的外边界。编号来自 [face_map.json](configs/face_map.json)，与实际 XML 和 Stage Q 输入哈希交叉核对。

## Solver provenance

PETSc **{runtime['PETSc_version']}**（commit `{runtime['PETSc_commit']}`），CUDA **{runtime['CUDA_version']}**，Open MPI **{runtime['MPI_version']}**（Fortran bindings={runtime['MPI_bindings']}），1 RTX 4090、MPI ranks={runtime['MPI_ranks']}、OMP={runtime['OMP_NUM_THREADS']}。svMultiPhysics upstream commit `{runtime['upstream']['upstream_commit']}`。

正式配置为 GMRES/right、restart100、rtol1e-10、atol1e-24、max_it2000、ASM overlap2、ILU(2)。RA：首次 fresh ILU 健康求解迭代数作为 I_ref；严格超过 1.5× 标记下次重建，age≥5 强制重建。实际完整运行中增长触发与失败恢复均未触发，不能据此证明 adaptive trigger 本身比固定 R5 更快。

实际 XML 中有旧 `petsc-jacobi` 初始标签，运行时 PETSC_OPTIONS 将其覆盖为 ASM/ILU；以 [实际 options](frozen_reference/run/PETSC_OPTIONS.txt) 和已核对的 [runtime manifest](frozen_reference/run/runtime_manifest.json) 为准。配置没有被“清理”成另一份新配置。环境证据见 [GPU manifest](frozen_reference/run/gpu_environment_manifest.json)。

源码在 [vendor/svMultiPhysics_stage_q](vendor/svMultiPhysics_stage_q)，生命周期、safe-stop、Stage P/Q ILU 修改链见 [PATCH_PROVENANCE.md](patches/PATCH_PROVENANCE.md)。未提交 solver binary、PETSc/CUDA/MPI 安装树或 build cache。最终 checkpoint 已提交，仅用于 FEM 可复现证据，不是 Particle-0 背景流输入。

## 新聊天第一件事

读取本文件和 [FROZEN_FEM_BASELINE.md](FROZEN_FEM_BASELINE.md)，运行下述只读完整性检查，再为 Particle-0 的采样接口、点定位和梯度验证制定实现计划；粒子实现留到新聊天。

在此 FEM 子目录创建自己的 Python≥3.11 环境并 `python -m pip install -e .`。已验证的读取依赖版本见 [python_reader_versions.json](sync_metadata/python_reader_versions.json)。

```bash
python -B scripts/fem_freeze_sync/validate_frozen.py
python -B scripts/fem_freeze_sync/verify_manifest.py
python -B -m pytest -q tests/test_fem_sync_manifest.py tests/test_frozen_*.py tests/test_particle_handoff_contract.py tests/test_no_particle_implementation_yet.py tests/test_main_history_preserved.py
```

这些命令不运行 CFD，不连接远端；新校验脚本从自己的路径定位仓库，运行目录无需是原 WSL 路径。历史开发脚本的移植限制、Stage Q 原测试命令及历史失败见 [PORTABILITY_NOTES.md](PORTABILITY_NOTES.md)。

## 2D reference（只读链接）

参考分支：[sync/ulm-microbubble-traj-gen-2d-20260917](https://github.com/OZAKI39/ulm-3d-vascular-model-generation/tree/sync/ulm-microbubble-traj-gen-2d-20260917)。WSL provenance：`/home/lzy/projects/ulm_microbubble_traj_gen_2D`。

可参考 `particle_field_sampling.py`、`particle_mobility.py`、`particle_mobility_transport.py`、`particle_collisions.py`、`particle_constrained_step.py`、`particle_inlet_flux.py`、`particle_counter_rng.py`。`red_blood_cell_transport.py` 是旧 reduced-order RBC model，Particle-0→8 不继续把它作为正式 RBC 动力学。此次没有复制、修改或 merge 2D 分支。

后续路线见 [PARTICLE_RESEARCH_ROADMAP.md](PARTICLE_RESEARCH_ROADMAP.md)。
''')
    metrics = {'rho (kg/m³)':physics['density_kg_m3'], 'mu (Pa s)':physics['dynamic_viscosity_pa_s'],
        'Q target (m³/s)':physics['inlet_volume_flow_m3_s'], 'dt (s)':base['dt_s'],
        '首次正式稳态 step':base['first_formal_steady_step'], 'safe-stop step':base['safe_stop_step'],
        '最终 TimeValue (s)':flow['TimeValue'], 'wall time (s)':base['wall_time_s'],
        '最大速度 (m/s)':flow['max_velocity_m_s'], '最低压力 (Pa)':flow['pressure_range_pa'][0], '最高压力 (Pa)':flow['pressure_range_pa'][1],
        'Qin (m³/s)':measured['Q_in_m3_s'], **{k+' (m³/s)':v for k,v in measured['outlet_flows_m3_s'].items()},
        'mass error |Qout−Qin|/Qtarget':measured['epsilon_mass'], 'inlet error':measured['epsilon_Q'],
        '最大壁面速度 (m/s)':measured['wall_velocity_max_m_s'], 'KSP attempts':base['statistics']['KSP_solves'],
        'KSP iterations':base['statistics']['total_iterations'], 'ILU numerical rebuilds':base['ILU_rebuild_count'],
        'ILU reuse':base['ILU_reuse_count'], 'linear unrecovered failures':base['linear_unrecovered_failures'], 'nonlinear failures':base['nonlinear_failures']}
    metrics_table='| 项目 | 实际值 |\n|---|---:|\n'+'\n'.join(f'| {k} | {v} |' for k,v in metrics.items())
    put('FROZEN_FEM_BASELINE.md', f'''
# Frozen SimVascular FEM Baseline — Stage SV1.3Q

状态：**{base['stage_status']}**。FEM DEVELOPMENT: FROZEN。下表由本分支 frozen artifact fresh-read 和 accepted Stage Q JSON 自动生成，生成器是 `scripts/fem_freeze_sync/write_handoff.py`。

{metrics_table}

正式 SimVascular 网格有 {mesh['nodes']} 节点、{mesh['tetra']} 四面体和 {mesh['surface_facets']} 外表面三角形，单位 m。不是旧 DOLFINx 网格。网格 SHA256：`{mesh['sha256']}`。

BC：INLET 为 steady Dirichlet flat profile、impose_flux、zero_out_perimeter；原 XML 入口值 `{base['BCs']['INLET']['Value']}`；三个出口为零 Neumann，WALL 为零 Dirichlet。ρ、μ、Q、dt 和 generalized-alpha spectral radius {base['spectral_radius']} 保持实际运行值。

稳态门槛：每 {base['production_policy']['save_interval_steps']} 步保存，连续 {base['production_policy']['steady_last_intervals']} 区间 E_u≤{base['production_policy']['velocity_change_limit']}、E_Q≤{base['production_policy']['flow_change_limit']}，mass gate≤{base['production_policy']['mass_limit']}。step {base['first_formal_steady_step']} 触发停止请求；step {base['safe_stop_step']} 是允许完成的 in-flight timestep，final VTU 与 native checkpoint 时间相同。这里的 mass error=0 是当前浮点积分输出，不是物理误差恰好为零的断言。

背景流：[{flow['path']}]({flow['path']})，{flow['size_bytes']} bytes，SHA256 `{flow['sha256']}`。`Velocity`/`Pressure` 是 Float64 节点数据，均 finite，wall no-slip PASS；速度梯度、涡量和应变率没有存储。

{freeze}

[接口入口](PARTICLE_HANDOFF.md) · [Stage Q 原报告](reports/sv1_3q/REPORT.md) · [fresh-read 结构化摘要](frozen_reference/baseline_summary.json) · [原始验收](reports/sv1_3q/REAL_VASCULAR_GPU_ILU_REUSE_WINNER_acceptance.json)
''')
    put('PARTICLE_RESEARCH_ROADMAP.md', '''
# Particle-0 → Particle-8 研究路线

FEM 已冻结。以下仅是研究路线，本次没有实现任何一阶段；统一采用第一版单向耦合 FEM → particles。

- **Particle-0 — Frozen FEM sampling**：定义冻结场读取、四面体点定位、插值和速度梯度接口，验证单位、边界和域外行为。
- **Particle-1 — single MB**：研究单个微泡在冻结背景流中的运动与基本力学，建立可检查的单粒子案例。
- **Particle-2 — single rigid RBC / orientation**：建立单个刚性 RBC 的平动、转动与取向模型，不延用旧 reduced-order RBC 模型作为正式动力学。
- **Particle-3 — wall interaction**：研究最近壁面、法向与粒子壁面相互作用，明确接触与穿透处理。
- **Particle-4 — particle contact**：研究粒子之间的接触与约束，并验证多体接触一致性。
- **Particle-5 — near-field hydrodynamic resistance**：研究近场流体阻力及其与接触模型的配合，避免重复计入作用。
- **Particle-6 — LAMMPS bridge**：研究与 LAMMPS 的数据及积分接口，明确职责、单位和同步方式；本阶段不安装它。
- **Particle-7 — inlet/outlet + continuous injection**：实现阶段将研究入口通量、持续注入、出口离域及可复现随机性。
- **Particle-8 — full RBC + MB suspension**：组合经验证的模块研究 RBC 与 MB 悬浮体系，报告适用范围和数值限制。
''')
    put('PORTABILITY_NOTES.md', f'''
# 可移植性与历史证据

冻结读取入口、manifest 校验和新测试只使用当前 checkout 相对路径。它们不需要 `/home/lzy`、原 FEM Git 历史、远端 SSH 或 GPU。Python 读取包见 `pyproject.toml` 与 [版本记录](sync_metadata/python_reader_versions.json)。

历史 `scripts/`、`configs/`、`reports/` 和 `tests/` 保留来源字节，不据此承诺原开发机自动化命令可在任意机器直接重跑。尤其 `remote.py`/`invoke.py` 依赖旧 SSH transport 和旧本地 Git blob，build/runner 脚本引用旧 `/workspace/` 安装前缀，部分旧报告依赖 `/home/lzy/` 及未打包历史输出。完整扫描见 [portability_inventory.json](sync_metadata/portability_inventory.json)，每项列路径和行号。不要为运行本阶段校验去执行这些历史命令。

原 Stage Q 20 项测试保留不改，可在当前分支运行：

```bash
python -B scripts/fem_freeze_sync/prepare_test_paths.py
python -B -m pytest -q tests/test_sv13q_*.py
```

准备脚本只在被 Git 忽略的 `external/sv13q/` 生成可删的源码副本，恢复原测试期待的路径及上游 `.gitattributes` 原字节，不执行源码或下载依赖。已提交 vendor 的 `.gitattributes` 仅关闭上游 LFS filter；原文件位置列于 `sync_metadata/upstream_attribute_map.json`（根规则保存在 `sync_metadata/upstream_gitattributes.txt`）。上游 950 个历史示例/安装包本来就是未下载的 LFS pointer；其原始指针文本、OID、逻辑大小和 SHA 全部归档在 `sync_metadata/upstream_lfs_pointer_inventory.json`，不把 raw pointer blob 交给 GitHub 触发未知 LFS 对象检查，也不下载这些示例。准备脚本可在 ignored 源码副本里恢复原指针字节；正式冻结 flow/mesh/boundary/checkpoint 均为实际完整文件，不是 LFS pointer。此分支不启用 LFS。

Stage Q 当时的完整历史 suite：{stage['full_passed']} passed、{stage['full_failures']} failed、{stage['full_skipped']} skipped；原 Stage Q 本阶段测试 {stage['stage_passed']} passed、{stage['stage_failures']} failed。下面是仍保留的早期失败证据，不被改写为 PASS，也未转换成 pytest xfail：

{historical_failures}

本次同步只运行新增完整性测试及可闭合的 Stage Q 测试；不声称精简分支的整个历史 pytest 都能运行。原始 suite 的很多测试需要被排除的依赖树、旧输出或原开发 Git 历史，失败原因不能混同于新的科学失败。[旧完整测试记录](reports/sv1_3q/pytest_summary.json) 原样保留。

`source_manifest.csv` 是本次同步文件清单，不等同于旧阶段 `delivery_manifest.json`：后者完整记录了当时开发机交付范围，包含本次明确不提交的 build/dependencies/intermediates。遗漏清单 [omitted_artifacts.json](sync_metadata/omitted_artifacts.json) 提供原 WSL 路径、size、SHA256；无必要 Particle-0 文件依赖这些遗漏。

Frozen solver XML 保留原相对路径，在 `frozen_reference/run/solver.xml` 下可完整解析到相邻 `SV_MESH/`。`frozen_reference/run/STOP_SIM` 与 `sv13q_stop_ack.json` 是完成时的停止证据；此目录用于只读交接，不是新启动目录，不要拿留下的 STOP_SIM=0 作为初始运行输入。Native checkpoint 保留最终 step71，随仓库提交，不需 pointer。历史运行时 binary/library 的 SHA 保留为 provenance，不上传二进制安装树。
''')
    put('patches/PATCH_PROVENANCE.md', f'''
# Frozen Stage Q 源码来源

svMultiPhysics upstream：{runtime['upstream']['upstream_remote']}，commit `{runtime['upstream']['upstream_commit']}`。

`vendor/svMultiPhysics_stage_q/` 保存 WSL accepted Stage Q 的可编译源码及许可证。原 source snapshot 包含的 950 个未下载上游示例/安装包指针归档为 JSON provenance，不作为 raw LFS pointer blob 提交。全部 {runtime['upstream']['verified_files']} 个原始文件按 `reports/sv1_3q/source_patch.json` 的 after hashes 核对；仅发布包装中的两份 `.gitattributes` 用无 LFS 规则替换，原字节归档在 `sync_metadata/upstream_gitattributes.txt`。这不改变 solver source。

继承顺序：

1. Stage N：[PETSc lifecycle adapter](sv1_3n/svmp_petsc325_compat.patch) 与 [来源说明](sv1_3n/PATCH_PROVENANCE.md)，在 MPI_Finalize 前完成 PETSc cleanup。
2. Stage O：[final-output-on-stop](sv1_3o/final_output_on_stop.patch)，native stop 时保留 final VTU/checkpoint。
3. Stage P：[within-timestep ILU reuse](sv1_3p/reuse_within_timestep.patch)。
4. Stage Q：[跨步复用、adaptive recovery 与最终 stop timestamp](sv1_3q/ilu_rebuild_policy.patch)。单独的 `stop_ack_timestamp.patch` 是历史中间补丁，不要重复应用在最终完整补丁上。

从原 upstream 到最终 Q 的累计补丁：[upstream_to_frozen_stage_q.patch](upstream_to_frozen_stage_q.patch)。它包含新增 `sv13q_reuse.h`。恢复源码时选“vendor source + prepare_test_paths 恢复原快照占位指针”或“upstream + 累计补丁”，不要同时应用累计补丁和分阶段补丁。累计补丁已经在临时 pristine upstream 树中应用并逐文件比对，见 `sync_metadata/patch_reconstruction.json`；没有编译或运行 solver。

Stage M 的 PETSc 3.19 CUDA ghost backport 及 repair 历史也保留在 `sv1_3m/`，它属于历史路径，当前 PETSc 3.25.5 不需要把旧 PETSc patch 重新套上去。

实际 build provenance：`reports/sv1_3q/remote/svmp_reuse_build.json`。最终 executable SHA256 `{runtime['solver_sha256']}`，PETSc library SHA256 `{runtime['PETSc_library_sha256']}`。未提交 executable、CMake cache 或 toolkit。旧 wrapper/build 脚本保留供审计，移植限制见 `../PORTABILITY_NOTES.md`。
''')
    put('reports/HISTORY.md', '''
# FEM 历史与证据范围

当前背景流唯一入口是 Stage SV1.3Q 的 frozen_reference，不能用旧流场替代。

- [SV1](sv1/REPORT.md) / [SV1.1](sv1_1/REPORT.md)：真实 SimVascular meshing 与线性求解稳定化。旧短算例质量守恒失败仍是失败证据。
- [SV1.2](sv1_2/REPORT.md)：完整 transient-to-steady 验证。
- [SV1.3](sv1_3/REPORT.md)：CPU early-stop；保留历史 production 身份，不在本阶段运行 CPU。
- [G](sv1_3g/REPORT.md)、[H](sv1_3h/REPORT.md)、[J](sv1_3j/REPORT.md)：GPU runtime、CUDA/PETSc 旧版兼容与矩阵测试；部分路线失败，被后续升级取代。
- [L](sv1_3l/REPORT.md)：修复缺失的 MPI Fortran predefined datatypes，继续暴露 CUDA ghost vector 兼容问题。
- [M](sv1_3m/REPORT.md)：旧 PETSc CUDA ghost backport，继续暴露 teardown/lifecycle 问题。
- [N](sv1_3n/REPORT.md)：现代 PETSc 3.25.5 + CUDA13.2，完成 lifecycle repair 和真实 GPU flow bring-up；历史 G/H/J/L/M 失败不因该结果被改写成成功。
- [O](sv1_3o/REPORT.md)：GPU 简单调优及 native final-output-on-stop。
- [P](sv1_3p/REPORT.md)：ASM2/ILU2、同时间步复用，保留 fallback。
- [Q](sv1_3q/REPORT.md)：跨时间步 reuse 比较及最终 RA 完整稳态，作为冻结的粒子开发背景流。

关键原报告、现有测试与小型证据保留；巨型机器 inventory、build 和重复中间输出不迁入 Git。遗漏路径与 SHA 在 `../sync_metadata/omitted_artifacts.json`。原 Stage Q 完整 suite 的 8 项历史失败见 `../PORTABILITY_NOTES.md`；本次未重跑 CFD 或将旧失败标成 PASS。
''')
    put('README.md', '''
# Frozen 3D FEM / SimVascular baseline

FEM DEVELOPMENT: **FROZEN after Stage SV1.3Q**。本目录保存粒子开发需要的源码、accepted mesh/flow/boundaries、配置、证据与校验；尚未实现 Particle-0。

从 [PARTICLE_HANDOFF.md](PARTICLE_HANDOFF.md) 开始，再读 [FROZEN_FEM_BASELINE.md](FROZEN_FEM_BASELINE.md) 和 [PARTICLE_RESEARCH_ROADMAP.md](PARTICLE_RESEARCH_ROADMAP.md)。同步详情见 [SYNC_REPORT_FEM_PARTICLE_HANDOFF.md](SYNC_REPORT_FEM_PARTICLE_HANDOFF.md)，历史及可移植性见 [PORTABILITY_NOTES.md](PORTABILITY_NOTES.md)。

新环境：Python≥3.11，`python -m pip install -e .`。只读验证：

```bash
python -B scripts/fem_freeze_sync/validate_frozen.py
python -B scripts/fem_freeze_sync/verify_manifest.py
```

科学数值源于实际 Stage Q artifact，不源于聊天手抄。Mesh convergence / time-step sensitivity 未做且不计划在粒子开发前补做；CPU/GPU field equivalence deferred，均为用户接受的项目决定，不是 PASS。禁止把本分支当作开始进一步 FEM 调优的授权。
''')
    put('SYNC_REPORT_FEM_PARTICLE_HANDOFF.md', f'''
# Stage FEM-FREEZE-SYNC 同步报告

## 从哪里同步？

唯一 source of truth：`{meta['source_FEM_path']}`。所有 accepted 资料来自 WSL；没有从 GPU 服务器重新下载，没有运行 CFD。代码/测试/配置逐文件复制并校验，原 WSL Git 历史和工作树保持原状。

## 同步到哪里？

GitHub：`{meta['repository']}`，分支 `{meta['branch']}`。WSL 独立克隆生成 commit，实际 remote 为 `{meta['remote_url']}`。

## main 有没有被修改？

**NO**。分支由 fetch 后最新 `origin/main` 的 `{meta['base_commit']}` 创建；main commit date `{meta['base_commit_date']}`。仅向新 sync 分支 push，无 force。差异限定为本 FEM 子目录及 root README 的短入口链接，原 main 代码不改写；最终审计见 `sync_metadata/main_history_audit.json`。

## 冻结 FEM 是哪一阶段？

Stage SV1.3Q，`{base['stage_status']}`。RA 完整稳态 {base['wall_time_s']:.6f} s，首次正式稳态 step{base['first_formal_steady_step']}、安全停止 step{base['safe_stop_step']}。这只是 frozen particle-development background-flow baseline，未扩大 scientific approval。

## 最终背景流是哪一个文件？

[{flow['path']}]({flow['path']})，{flow['size_bytes']} bytes，SHA256 `{flow['sha256']}`。

## Particle-0 需要哪些文件？

{table}

同时读取 `frozen_reference/mesh_manifest.json`、`frozen_reference/boundary_manifest.json`、`frozen_reference/flow/flow_field_manifest.json` 和 `frozen_reference/baseline_summary.json`。单位、字段、法向、cell ordering 约定均在交接文档。Checkpoint 是额外可复现证据，Particle-0 不必读取它。

## 哪些内容没有同步？

PETSc/MPI/CUDA source/install/build trees、solver executable、CMake cache、编译中间文件、Python 环境、缓存、core、下载包、远端 home、SSH key/token/credential，以及巨大机器 inventory/重复中间运行输出。上游 svMultiPhysics 源码及其许可证保留在 vendor；950 个未下载历史示例/安装包的 LFS pointer 原文保存在 `upstream_lfs_pointer_inventory.json`，没有启用或下载 LFS。

完整候选扫描、复制映射及遗漏清单见 `sync_metadata/copied_from_wsl.json`、`sync_metadata/source_manifest.csv`、`sync_metadata/omitted_artifacts.json`。Stage Q full solver log、steady history、adaptive trace、stop request/ack、原始验收与环境摘要都保留。旧环境 inventory 用 hash pointer 替代，不将 65–82 MiB 的机器清单放入普通 Git。

## 有没有大文件被 pointer 替代？

必要 scientific artifact **没有**：最终 flow、正式 mesh、全部边界和最终 checkpoint 都是完整 Git 文件。历史机器清单、依赖 archive、重复中间输出的 WSL path/size/SHA 记录在 `omitted_artifacts.json`；取回时从其记录的 source 路径复制并重新核对 SHA，禁止把未取回文件声称为已验证。上游已有的 950 个未下载示例/安装包只有指针：其原文/OID/逻辑大小另列于 `sync_metadata/upstream_lfs_pointer_inventory.json`，这些从来不是正式冻结 scientific artifact。没有 zip 科学二进制绕过大小限制。文件大小检查见 `sync_metadata/file_size_audit.json`。

## 有哪些 validation 明确没有做？

{freeze}

本次执行 frozen artifact fresh-read、同 run config/time/hash 一致性、边界法向与流量复核、来源 manifest、无 secret、大文件及 main 历史审计；运行同步测试和原 Stage Q artifact tests，均不求解 CFD。测试记录见 `sync_metadata/test_summary.json`；原完整历史失败仍见 Stage Q 原 pytest 报告。

## 下一个聊天从哪里开始？

[PARTICLE_HANDOFF.md](PARTICLE_HANDOFF.md) → Particle-0。此次未写 FrozenFEMField、point locator、interpolation、RBC/MB dynamics、LAMMPS 或 resistance solver。

## 发布验收

本阶段最终状态由 `sync_metadata/stage_status.json` 记录；远端 push 与 GitHub 内容重新读取的证据在 `sync_metadata/remote_verification.json`。远端验证按已发布 commit 锚定，最终 commit 的 local/remote HEAD 一致性另由交付终端记录验证，避免把文件自己的 commit SHA 递归写入自身。
''')
    print('Generated handoff, baseline, roadmap, portability, provenance, history and sync report from validated artifacts.')


if __name__ == '__main__':
    main()
