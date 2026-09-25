from pathlib import Path
import json,shutil,hashlib,gzip
A=Path(__file__).resolve().parents[1];C=json.loads((A/'context.json').read_text());D=Path(C['destination']);M=D/'sync_metadata'/C['snapshot'];rel=str(M.relative_to(D));P='particle_3d/reports/particle9a4_population_inlet'
changes={}
def edit(name,text):
 p=D/name;q=M/'original_docs'/name;q.parent.mkdir(parents=True,exist_ok=True)
 if p.exists():shutil.copy2(p,q)
 p.write_text(text.strip()+'\n');changes[name]='Publication navigation/packaging only; original saved under sync metadata.'
edit('README.md',f'''# Network H0 → 3D FEM → Microbubble / RBC

当前同步分支：`{C['branch']}`，基于 `abb3ab05fb8dcddcf120765702b682f23a63a71d`。包含当前本地工作文件、P9-A.4新源码和完整审计结果、历史RBC模块、服务器代码/日志/配置和来源记录。原工作树与服务器科学结果保持原位；主分支未合并。

最新 **P9A4_POPULATION_INLET_READY_FOR_LARGE_SAMPLE_REVIEW**：100,000 source proposals → 19,221 accepted births；NEW Network-H0上首30泡为28 completed / 2 supported stationary / 0 solver fail。恒定绝对浓度是明确的模型假设。没有运行正式500。

| 内容 | 仓库入口 |
| --- | --- |
| 新P9-A.4正式报告 | [入口群体审核]({P}/P9A4_POPULATION_INLET_REVIEW_ZH.md) |
| P9-A.4代码 | [continuous_infusion.py](particle_3d/src/particle_3d/continuous_infusion.py)、[population_inlet_p9a4.py](particle_3d/src/particle_3d/population_inlet_p9a4.py) |
| 科学代码逻辑审核 | [CURRENT_CODE_LOGIC_AUDIT_ZH.md]({P}/CURRENT_CODE_LOGIC_AUDIT_ZH.md) |
| 100k全部候选事件/拒绝、births、checkpoint | [data/inlet100k]({P}/data/inlet100k/)；大文本gzip无损保存，先按下文还原 |
| 30泡轨迹及point数据 | [outputs/smoke30]({P}/outputs/smoke30/) |
| 8组英文PNG/PDF图件 | [figures]({P}/figures/)；[HTML总览]({P}/OPEN_RESULTS.html)下载后浏览 |
| vascular A、ROI及Network H0 | [vascular_network](vascular_network/README.md) |
| NEW H0 FEM算例/网格/配置/冻结流场 | [mean-2p0-mmps-A-H0-pressure-v1](formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/) |
| 微泡/RBC完整科学源码、测试、合同 | [particle_3d](particle_3d/README.md) |
| 历史RBC运动、变形代理和共流展示 | [particle2](particle_3d/reports/particle2/)、[particle3](particle_3d/reports/particle3/)、[rbc_mb_flow_rotation](particle_3d/reports/rbc_mb_flow_rotation/) |
| 固定SonoVue分布 | [sonovue_size_distribution_v0](sonovue_size_distribution_v0/) |
| 压力/WSS/流线/速度矢量可视化 | [rotate_visualization](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/) |
| 服务器代码和执行记录 | [server_evidence](server_evidence/README.md) |
| 当前目录映射 | [WORKFLOW_PATHS_ZH.md](WORKFLOW_PATHS_ZH.md) |
| 本次同步范围、排除规则及验证 | [SYNC_REPORT_ZH.md]({rel}/SYNC_REPORT_ZH.md) |

## 克隆、校验和无损还原

```bash
git clone --single-branch --branch {C['branch']} \\
  https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git
cd ulm-3d-vascular-model-generation
python3 {rel}/verify_snapshot.py --scientific-inputs
python3 {rel}/restore_large_artifacts.py
```

普通Git对象，无需Git LFS。还原脚本验证gzip及原始字节SHA，恢复3份科学大文本到原路径；已有同SHA文件可重复执行，遇到不同内容拒绝覆盖。原科学源码、合同及原始报告不因压缩而修改。还原后可直接读取原checkpoint和旧报告所引用的JSON/JSONL。原阶段delivery manifest中绝对机器路径是历史来源记录，不是克隆目录要求。

使用已记录的科学Python依赖运行：

```bash
python {rel}/run_checks.py --output /tmp/ulm-p9a4-checks-new
```

输出目录须为全新目录。验证范围为portable 115项 + P9-A.4 35项 + legacy入口/open-cap 15项，共165项。历史HemoCell/LAMMPS及独立TopBrain开发只保存快照，不因此获得当前主线的验收状态。依赖版本见 [requirements-verified.txt](sync_metadata/network_h0_particle_20260925/requirements-verified.txt)。历史绝对路径迁移规则见[既有迁移说明](sync_metadata/network_h0_particle_20260925/PORTABILITY_ZH.md)。

## 输入版本

NEW H0 SHA：`064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。当前新入口显式验证该流场；Particle默认frozen_reference仍是历史OLD equal-pressure 2 mm/s场，不能混用。更早P8.2A/PPT500和RBC理想化Poiseuille演示保留各自来源；RBC演示不冒充NEW H0 RBC生产轨迹。Taylor-Hood仍为已停止验证阶段。

同步排除环境、构建/安装树、完整上游脑数据、渲染帧、大型重复部署包和非必要中间状态。重要100k科学数据无损保留。原报告中的“未push”是开发阶段历史记录；本分支是其后的发布快照。
''')
edit('particle_3d/README.md',f'''# Particle：最新 P9-A.4 与历史微泡/RBC

当前入口为 [P9-A.4]({P.removeprefix('particle_3d/')}/P9A4_POPULATION_INLET_REVIEW_ZH.md)：固定NEW H0，Poisson source、条件SonoVue D≤4 µm、全入口通量采样、单次真实WALL+handoff筛选，所有拒绝保留。P9-A.1动力学与原B/C源码不变；首30泡28 completed、2 supported stationary、0 solver fail。

源码 `src/particle_3d/continuous_infusion.py` 与 `population_inlet_p9a4.py`；合同 `contracts/P9A4_CONTINUOUS_INFUSION_V1.json`；测试 `tests/particle9a4_population_inlet/`；所有报告/日志/图件/原始30泡结果在 `reports/particle9a4_population_inlet/`。

100k ledger、accepted births及Method B大文本以gzip保存。先从仓库根执行 `python3 {rel}/restore_large_artifacts.py`，恢复逐字节相同的原路径，再使用原审计/重放入口。[发布复现说明](../{rel}/SYNC_REPORT_ZH.md)。

RBC与MB共用该包。`rbc.py`、`rbc_distribution.py`、`rbc_orientation.py`、`rbc_integrator.py`、`rbc_capillary_surrogate.py`与`coflow_rotation.py`均保留；RBC历史结果位于`reports/particle2/`、`particle3/`、`rbc_mb_flow_rotation/`。最后一个是理想化Poiseuille演示。

`reports/network_derived_flow_mb_validation_v1/`保留前一阶段OLD/NEW配对证据。旧P9-A.1/A.2/A.3/A.3B及P8.2A/PPT仍在原相对路径；旧合同与结果不改。默认frozen_reference是历史OLD场，P9-A.4须经专用NEW loader。没有自动正式500。
''')
edit('particle_3d/CURRENT_WORKFLOW_ZH.md',f'''# 当前 Particle 工作流

当前阶段为 P9-A.4。[主报告](reports/particle9a4_population_inlet/P9A4_POPULATION_INLET_REVIEW_ZH.md)、[机器摘要](reports/particle9a4_population_inlet/data/final_summary.json)、[源码](src/particle_3d/continuous_infusion.py)。状态 READY_FOR_P9A4_LARGE_SAMPLE_REVIEW，正式500未运行。

先按仓库根README还原3份gzip科学文本，再运行旧报告里的入口审计/恢复操作。Python/运行环境不打包；[当前同步验证入口](../{rel}/run_checks.py)。

NEW H0与未改P9-A.1动力学；source/accepted群体分别记录。原源码和所有阶段原始报告作为证据保留，本文件只提供当前导航。
''')
edit('server_evidence/README.md',f'''# 服务器源码与执行证据

来源 `root@50.115.148.16:4159`。本目录保留前一版服务器快照，并纳入当前额外代码、配置、日志和关键科学数据。原始服务器目录与源文件没有搬移。

最新P9-A.4服务器目录为 `/workspace/particle9a4_population_inlet_20260924T233708Z`；其源码、154个已回传输出与本次root `particle_3d/`内容相同时，使用manifest引用而不复制。新增部署manifest位于对应服务器证据目录。原始时间戳和绝对路径属于来源记录。

H0求解、Network配对验证、P8/P9历史运行、shear-lift，以及历史HemoCell RBC、LAMMPS/Palabos资料均已检查。历史模拟的代码/日志保留不等于对其物理模型作新的审核。HemoCell/LAMMPS的原生大状态、完整第三方源码/工具链不全量打包。

本轮逐文件映射、去重引用与排除原因：[server_delta_inventory.json.gz](../{rel}/server_delta_inventory.json.gz)；[汇总](../{rel}/server_delta_summary.json)。`REPRESENTED_BY_EXISTING_ARTIFACT`指相同SHA内容已在仓库其他路径；representation明确普通字节或gzip原始字节。所有此类引用在整理后校验。`EXCLUDED_FROM_DELTA`只表示本次不新增，不删除父提交已经保存的历史证据。

固定svMultiPhysics补丁及实际源文件仍在 `pinned_solver_source/`；完整当前源与许可证见根目录`vendor/`。CUDA/PETSc/MPI安装、二进制和Python环境留在服务器，按构建/版本记录复现。
''')
edit('WORKFLOW_PATHS_ZH.md',f'''# 当前代码、数据与服务器路径

这是适合克隆后浏览的入口。原机器上335条关键路径和本次只读目录扫描见 [原始路径索引]({rel}/path_inventory/WORKFLOW_PATHS_20260925_ZH.md)；该记录中的绝对路径用于定位原机器，不是GitHub相对链接。

| 本地/服务器来源 | 仓库内对应 |
| --- | --- |
| `/home/lzy/projects/ulm_3D_vascular/` | [vascular_network/](vascular_network/) |
| `/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/` | [formal_3D_flow_solver/FEM_SimVascular/](formal_3D_flow_solver/FEM_SimVascular/) |
| `/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/` | [particle_3d/](particle_3d/)；受当前新增P9-A.4扩展 |
| `/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/` | [particle_3d/](particle_3d/)；新模块、合同、测试与报告完整保存 |
| `/home/lzy/projects/sonovue_size_distribution_v0/` | [sonovue_size_distribution_v0/](sonovue_size_distribution_v0/) |
| `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/` | [server_evidence/flow_mean_2p0_mmps_A_H0_20260924T181140Z/](server_evidence/flow_mean_2p0_mmps_A_H0_20260924T181140Z/) |
| `/workspace/particle9a4_population_inlet_20260924T233708Z/` | 主体与当前 [P9-A.4报告/输出]({P}/)相同；其余按server manifest定位 |
| `/root/particle8_2_runs/` | [server_evidence/particle8_2_runs/](server_evidence/particle8_2_runs/)及去重映射 |
| `/workspace/hemocell_restore/` | [server_evidence/hemocell_restore/](server_evidence/hemocell_restore/)及去重映射 |
| `/workspace/microbubble_lammps/` | [server_evidence/microbubble_lammps/](server_evidence/microbubble_lammps/)及去重映射 |

代码在src/scripts，模型合同在contracts，科学数据在data/outputs/run，日志在logs及各阶段reports。目录级导航不替代逐文件manifest；被排除的大文件仍留在原机器，理由及路径见同步元数据。
''')
raws=json.loads((M/'packed_artifacts.json').read_text())
edit('.gitignore',(D/'.gitignore').read_text()+'\n# Restored exact scientific text is generated from tracked gzip artifacts.\n'+'\n'.join('/'+r['original_path'] for r in raws)+'\n')
edit('.ignore',(D/'.ignore').read_text()+'\n!particle_3d/reports/particle9a4_population_inlet/\nparticle_3d/reports/particle9a4_population_inlet/data/\nparticle_3d/reports/particle9a4_population_inlet/figures/\nparticle_3d/reports/particle9a4_population_inlet/logs/\nparticle_3d/reports/particle9a4_population_inlet/outputs/\n')
(M/'packaging_edits.json').write_text(json.dumps(changes,indent=2)+'\n')
shutil.copy2(A/'audit/packed_restore_validation.json',M/'packed_restore_validation.json')
shutil.copy2(A/'audit/server_reference_validation.json',M/'server_reference_validation.json')
shutil.copy2(A/'audit/source_git_before.json',M/'source_git_before.json')
for p in (A/'scripts').glob('*.py'):
 if p.name!='server_delta_payload.py':
  q=M/'export_tools'/p.name;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
print('NAVIGATION_COMPLETE')
