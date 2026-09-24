from pathlib import Path
import json,shutil,collections

A=Path(__file__).resolve().parents[1];C=json.loads((A/'context.json').read_text());D=Path(C['destination'])
M=D/'sync_metadata'/C['snapshot'];SV='formal_3D_flow_solver/FEM_SimVascular';R='particle_3d/reports/network_derived_flow_mb_validation_v1'
edits={}
def put(rel,body,why):
    p=D/rel
    if p.exists():
        saved=M/'original_docs'/rel;saved.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,saved)
    p.write_text(body.strip()+'\n');edits[rel]=why

put('README.md',f'''# Network A/H0 → 3D FEM → Microbubble

这是截至 2026-09-25（柏林时间）的工作流同步分支，包含代码瘦身后的本地实现、必要输入、冻结流场、轨迹、可视化、正式报告，以及服务器源码和执行记录。

当前结果为 **NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY**：同一组 30 个初始状态，OLD/NEW 各 29 条完成、1 条 stationary、0 条数值失败。500 条新流场生产计算尚未启动。

| 内容 | 仓库入口 |
| --- | --- |
| 血管 A、ROI、端口、Network H0 | [vascular_network](vascular_network/README.md) |
| H0 3D FEM 算例、求解、校验 | [FEM_SimVascular]({SV}/README.md) |
| 当前 OLD/NEW 配对轨迹与科学代码 | [particle_3d](particle_3d/README.md) |
| 配对验证正式报告 | [NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md]({R}/NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md) |
| 压力、WSS、流线、局部矢量原可视化接口 | [rotate_visualization]({SV}/rotate_visualization/OPEN_RESULTS.html) |
| 服务器保留源码、配置、日志、执行证据 | [server_evidence](server_evidence/README.md) |
| SonoVue 固定分布 | [sonovue_size_distribution_v0](sonovue_size_distribution_v0/README.md) |
| 历史 FEM / RBC / 2D 来源 | `formal_3D_flow_solver/FEM/`、`particle_3d/reports/`、`source_dependencies/` |
| 同步范围、清理记录、逐文件来源与排除原因 | [同步报告](sync_metadata/{C['snapshot']}/SYNC_REPORT_ZH.md) |

## 下载与检查

```bash
git clone --single-branch --branch {C['branch']} \\
  https://github.com/{C['repo']}.git
cd ulm-3d-vascular-model-generation
python3 sync_metadata/{C['snapshot']}/verify_snapshot.py
```

完整性检查只需 Python 标准库；本分支使用普通 Git 对象。使用具备 numpy/scipy/VTK/pyvista/pytest/Pillow/matplotlib 等依赖的 Python 3.11+ 环境后运行：

```bash
python sync_metadata/{C['snapshot']}/verify_snapshot.py --scientific-inputs
python sync_metadata/{C['snapshot']}/run_checks.py --output /tmp/ulm-workflow-checks
```

发布前在临时迁移目录中完成 **115 项回归，全部通过**（Particle 24、Network 49、流场 42）。测试使用仓库内文件；只在临时副本中迁移 5 份 JSON 的历史绝对路径，源码、数组、公式、参数及测试断言保持原字节。原科学清单中的 92 个解释器缓存条目在临时验证清单中排除，全部 114 个实际科学源文件仍校验原 SHA。原报告和原清单保留。

## 流场角色

| 路径 | 角色 |
| --- | --- |
| `{SV}/frozen_reference/` | Particle 默认使用的历史 2 mm/s 场 |
| `{SV}/flow_cases/mean-2p0-mmps/frozen_flow/` | OLD，SHA `129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d` |
| `{SV}/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/frozen_flow/` | NEW/H0，SHA `064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4` |
| `{SV}/upstream_stage_q_reference/` | Flow 原工作树的早期 Stage Q 冻结交付，用于历史流场回归 |
| `{SV}/legacy_inputs/` | Particle 迁移前输入的历史角色 |

流场求解使用固定的 svMultiPhysics/PETSc GPU 工具链；微泡轨迹由 CPU 进程按泡并行。必要构建脚本、源补丁、版本与服务器记录已保留；已安装环境与二进制可按记录重建。

原始全脑数据、虚拟环境/构建树、多余渲染帧、中间 CFD 状态、重复部署包及非必要大文件未纳入。依赖当前回归的约 29.5 MiB 流线候选池作为重要文件保留。历史全流程脚本仍有原机器路径，重新计算前需按[迁移说明](sync_metadata/{C['snapshot']}/PORTABILITY_ZH.md)配置；本次没有重算 CFD 或轨迹。
''','New branch-relative workflow entry and verified input roles.')

put('particle_3d/README.md',f'''# 当前 Particle 科学代码

当前任务是 Network H0 新旧流场的 30 泡配对验证，状态 `NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY`。OLD/NEW 各 29 条完成、1 条 stationary、0 条数值失败。

- 运行入口：[runner.py](reports/network_derived_flow_mb_validation_v1/scripts/runner.py)。
- 正式报告：[NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md](reports/network_derived_flow_mb_validation_v1/NETWORK_DERIVED_FLOW_MB_VALIDATION_ZH.md)。
- 科学实现：`src/particle_3d/`；原科学快照与固定输入：`reports/network_derived_flow_mb_validation_v1/server_bundle/`。
- 核心调用：`Particle9AStepper → particle82a_integration.integrate_admitted → particle6_stepper / particle65_motion → 阻力与接触求解`。轨迹在 CPU 上按泡并行。

`FrozenFEMField` 进行四面体场查询，壁面/阻力/接触共享模块及历史数学回归均保留。点示踪对照使用 `NativePointTracer`。旧阶段名称仍包含实际依赖，未按编号删除科学核心。

从仓库根目录运行 [run_checks.py](../sync_metadata/{C['snapshot']}/run_checks.py)，可完成包含当前 Particle 24 项的完整 115 项回归；它在临时目录迁移历史路径，不修改本仓库数据。

新计算必须指定新的输出名。CPU 配对复现可使用：

```bash
env PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VTK_SMP_MAX_THREADS=1 \\
python particle_3d/reports/network_derived_flow_mb_validation_v1/scripts/runner.py \\
  --root particle_3d/reports/network_derived_flow_mb_validation_v1/server_bundle \\
  --label NEW --workers 6 --ids all --output replay_new_001
```

该命令会实际重新积分，入口会验证科学快照哈希并拒绝覆盖已有输出。本次同步仅运行回归与完整性检查。依赖版本见 `requirements-validation-lock.txt` 和各项目 `pyproject.toml`。

各阶段历史数据、RBC 示例与独立升力审核保留其原科学适用范围。跨仓库入口见[主 README](../README.md)。
''','Portable Particle entry; original science implementation retained byte-for-byte.')
put('particle_3d/CURRENT_WORKFLOW_ZH.md','''# 当前 Particle 工作流

当前为 Network H0 的 30 泡 OLD/NEW 配对验证，状态 `NETWORK_FLOW_MB_VALIDATION_PASS_WITH_STATIONARY`，各 29 条完成、1 条 stationary、0 条数值失败。500 条新流场生产尚未启动。
当前入口、调用链和复现说明统一见 [README.md](README.md)。
''','Replace obsolete machine-specific navigation with the current entry.')
put(SV+'/README.md',f'''# 3D FEM / SimVascular

当前是 Network A/H0 出口压力驱动的冻结 3D 血管流场。

- H0 算例生成：[fem_h0_case.py](../../vascular_network/network_1d0d/fem_h0_case.py)。
- 远程求解：[solve_a_h0_fem_remote.py](../../vascular_network/scripts/solve_a_h0_fem_remote.py)。
- 物理校验：[validate_a_h0_fem.py](../../vascular_network/scripts/validate_a_h0_fem.py)。
- 共享实现：`src/sv_validation/`；当前日志解析：`scripts/sv13q/flow_parser.py`。
- OLD 与 H0 算例：[flow_cases](flow_cases/)；当前派生量入口：`scripts/flow_2mmps/`。
- 原旋转可视化接口：[OPEN_RESULTS.html](rotate_visualization/OPEN_RESULTS.html)。
- 固定求解器源码：[vendor/svMultiPhysics_stage_q](vendor/svMultiPhysics_stage_q/)，其三个修改文件与本次服务器实际源文件逐字节相同；运行记录见 [server_evidence](../../server_evidence/README.md)。

`frozen_reference` 供当前 Particle 使用；`upstream_stage_q_reference` 保留 Flow 工作树的旧冻结输入。两者角色不同。独立回归脚本在临时目录选择对应输入，避免改写任何已有冻结目录。
Stage N 的 PETSc/CUDA 与 Stage L 的 MPI/Fortran 构建脚本、启动 wrapper 和 Stage Q reuse patch 仍有实际用途，已保留。

完整性和当前 42 项流场测试包含在[主 README](../../README.md)的统一检查命令中。Taylor–Hood 验证已因资源成本停止，代码/日志作为历史证据保留，不代表已有完整 P2/P1 稳态场。
''','Clarify consolidated frozen-input roles and use branch-relative links.')
put('vascular_network/README.md',f'''# vascular A → ROI → Network H0

本目录保留血管预处理、Network A/H0 的源码、配置、测试和必要源数据。

| 阶段 | 入口 |
| --- | --- |
| 小鼠 ROI | `s1-1_swc_roi_generate_mouse.py`、`utils/schmid_pkl/`、`utils/sampling/` |
| Ultraliser 表面 | `s2_swc_stl_model_generate.py`、`utils/cfd_lumen/` |
| 端口与 CFD 表面 | `s3_cfd_1D_data_preprocess.py`、`s4_cfd_surface_prepare.py` |
| Network A 基础图与精确切口 | `scripts/run_a_network_1d0d_boundary_v1.py`、`network_1d0d/` |
| H0 假设与求解 | `scripts/run_a_network_h0_v2.py`、`network_1d0d/idealized_h0.py` |
| 3D H0 算例与验收 | `network_1d0d/fem_h0_case.py`、`scripts/solve_a_h0_fem_remote.py`、`scripts/validate_a_h0_fem.py` |

H0 仍依赖 v1 图/端口，因此保留 v1 的完整必要输入哈希闭包。只选择本案例的源 SWC、ROI、几何和配置，没有复制完整原始脑数据集。
当前报告见 [A_NETWORK_H0_REVIEW_ZH.md](reports/a_network_1d0d_boundary_v2_idealized/A_NETWORK_H0_REVIEW_ZH.md) 和 [ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md](reports/a_network_1d0d_boundary_v2_idealized/ROI_3D_A_H0_PRESSURE_VALIDATION_ZH.md)。

原 JSON 中的绝对路径按原样保留以保存来源；[统一回归工具](../sync_metadata/{C['snapshot']}/run_checks.py)会在临时目录迁移这些路径并运行全部 49 项 Network 回归。

表面重建仍使用 Ultraliser，`feed_radius_um = source_radius_um * 0.91`，H5 第四列为直径 `2 * feed_radius_um`，CFD 坐标从微米按 `1e-6` 缩放为米。配置与说明在 `configs/`、`docs/`；第三方构建环境按记录重建。
其余 LBM、human/MeVO/TopBrain 的已有源码保留各自范围；并行进行的 TopBrain/BRAVA 大型输出不属于本次血流/Particle 数据快照。
''','Consolidated network subtree with original provenance retained.')

put('server_evidence/README.md','''# 服务器源码与执行证据

各子目录保留原 `/workspace/<目录名>/` 下当前仍存在的源码、配置和选定日志；Git 会复用相同内容的对象。这些历史证据不作为默认开发代码入口。

- `particle_network_flow_mb_validation_v1_20260924T210047Z/`：当前 30 泡配对验证、源快照和执行记录。
- `flow_mean_2p0_mmps_A_H0_20260924T181140Z/`：当前 H0 求解器入口、输入与实际运行证据。
- `flow_mean_2p0_mmps_A_H0_TaylorHood_20260924T194503Z/`：已停止的 P2/P1 验证代码与资源/失败记录。
- `formal_3D_flow_solver_FEM_SimVascular_lzy/`：保留的 Stage L/N/Q 工具链脚本、配置、日志和版本来源。
- `pinned_solver_source/`：实际运行的 svMultiPhysics 修改文件、上游 commit 和完整 tracked diff。

第三方安装目录、CUDA/PETSc/MPI 二进制、解释器缓存、重复部署包及大中间结果未上传。历史脚本中的服务器路径保持原样；重建和重新运行时应按当前机器配置迁移。
已在代码瘦身中删除的 371 个服务器脚本没有恢复到本分支。
''','Explain server snapshots and pinned build provenance.')

# Keep numeric bytes independent of Git LFS or checkout line-ending settings.
for p in D.rglob('.gitattributes'):
    if '.git' in p.parts:continue
    if 'filter=lfs' in p.read_text(errors='replace'):
        rel=str(p.relative_to(D));original=p.read_text();put(rel,original+'\n# Ordinary Git snapshot: retain actual bytes without LFS drivers.\n* -text -filter !diff !merge\n','Disable inherited LFS filters for this ordinary-Git snapshot.')
put('.gitattributes','# Preserve original scientific bytes; this snapshot requires no new Git LFS objects.\n* -text -filter\n','Preserve exact scientific bytes in ordinary Git objects.')
put('.ignore','''# Keep current development code visible; opt into archived evidence explicitly.
server_evidence/
sync_metadata/
**/vendor/
**/outputs/
**/results/
**/logs/
**/input_data/
**/inputs/
**/frozen_reference/
**/upstream_stage_q_reference/
**/legacy_inputs/
**/flow_cases/
**/reports/*
!particle_3d/reports/network_derived_flow_mb_validation_v1/
**/server_bundle/
**/__pycache__/
''','Exclude bulky evidence from default ripgrep context without affecting Git or Python.')

for p in (A/'portable_checks_final').iterdir():shutil.copy2(p,M/p.name)
for name in ['package_local.py','collect_server.py','write_docs.py']:
    q=M/'export_tools'/name;q.parent.mkdir(exist_ok=True);shutil.copy2(A/'scripts'/name,q)
(M/'packaging_edits.json').write_text(json.dumps(edits,ensure_ascii=False,indent=2)+'\n')
print('DOCUMENTATION_AND_VALIDATION_RECORDS_WRITTEN',len(edits))
