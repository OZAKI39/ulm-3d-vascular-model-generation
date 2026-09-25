# Network H0 → 3D FEM → Microbubble / RBC

当前同步分支：`sync/p9a4-flow-particle-rbc-20260925`，基于 `abb3ab05fb8dcddcf120765702b682f23a63a71d`。包含当前本地工作文件、P9-A.4新源码和完整审计结果、历史RBC模块、服务器代码/日志/配置和来源记录。原工作树与服务器科学结果保持原位；主分支未合并。

最新 **P9A4_POPULATION_INLET_READY_FOR_LARGE_SAMPLE_REVIEW**：100,000 source proposals → 19,221 accepted births；NEW Network-H0上首30泡为28 completed / 2 supported stationary / 0 solver fail。恒定绝对浓度是明确的模型假设。没有运行正式500。

| 内容 | 仓库入口 |
| --- | --- |
| 新P9-A.4正式报告 | [入口群体审核](particle_3d/reports/particle9a4_population_inlet/P9A4_POPULATION_INLET_REVIEW_ZH.md) |
| P9-A.4代码 | [continuous_infusion.py](particle_3d/src/particle_3d/continuous_infusion.py)、[population_inlet_p9a4.py](particle_3d/src/particle_3d/population_inlet_p9a4.py) |
| 科学代码逻辑审核 | [CURRENT_CODE_LOGIC_AUDIT_ZH.md](particle_3d/reports/particle9a4_population_inlet/CURRENT_CODE_LOGIC_AUDIT_ZH.md) |
| 100k全部候选事件/拒绝、births、checkpoint | [data/inlet100k](particle_3d/reports/particle9a4_population_inlet/data/inlet100k/)；大文本gzip无损保存，先按下文还原 |
| 30泡轨迹及point数据 | [outputs/smoke30](particle_3d/reports/particle9a4_population_inlet/outputs/smoke30/) |
| 8组英文PNG/PDF图件 | [figures](particle_3d/reports/particle9a4_population_inlet/figures/)；[HTML总览](particle_3d/reports/particle9a4_population_inlet/OPEN_RESULTS.html)下载后浏览 |
| vascular A、ROI及Network H0 | [vascular_network](vascular_network/README.md) |
| NEW H0 FEM算例/网格/配置/冻结流场 | [mean-2p0-mmps-A-H0-pressure-v1](formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/) |
| 微泡/RBC完整科学源码、测试、合同 | [particle_3d](particle_3d/README.md) |
| 历史RBC运动、变形代理和共流展示 | [particle2](particle_3d/reports/particle2/)、[particle3](particle_3d/reports/particle3/)、[rbc_mb_flow_rotation](particle_3d/reports/rbc_mb_flow_rotation/) |
| 固定SonoVue分布 | [sonovue_size_distribution_v0](sonovue_size_distribution_v0/) |
| 压力/WSS/流线/速度矢量可视化 | [rotate_visualization](formal_3D_flow_solver/FEM_SimVascular/rotate_visualization/) |
| 服务器代码和执行记录 | [server_evidence](server_evidence/README.md) |
| 当前目录映射 | [WORKFLOW_PATHS_ZH.md](WORKFLOW_PATHS_ZH.md) |
| 本次同步范围、排除规则及验证 | [SYNC_REPORT_ZH.md](sync_metadata/p9a4_flow_particle_rbc_20260925/SYNC_REPORT_ZH.md) |

## 克隆、校验和无损还原

```bash
git clone --single-branch --branch sync/p9a4-flow-particle-rbc-20260925 \
  https://github.com/OZAKI39/ulm-3d-vascular-model-generation.git
cd ulm-3d-vascular-model-generation
python3 sync_metadata/p9a4_flow_particle_rbc_20260925/verify_snapshot.py --scientific-inputs
python3 sync_metadata/p9a4_flow_particle_rbc_20260925/restore_large_artifacts.py
```

普通Git对象，无需Git LFS。还原脚本验证gzip及原始字节SHA，恢复3份科学大文本到原路径；已有同SHA文件可重复执行，遇到不同内容拒绝覆盖。原科学源码、合同及原始报告不因压缩而修改。还原后可直接读取原checkpoint和旧报告所引用的JSON/JSONL。原阶段delivery manifest中绝对机器路径是历史来源记录，不是克隆目录要求。

使用已记录的科学Python依赖运行：

```bash
python sync_metadata/p9a4_flow_particle_rbc_20260925/run_checks.py --output /tmp/ulm-p9a4-checks-new
```

输出目录须为全新目录。验证范围为portable 115项 + P9-A.4 35项 + legacy入口/open-cap 15项，共165项。历史HemoCell/LAMMPS及独立TopBrain开发只保存快照，不因此获得当前主线的验收状态。依赖版本见 [requirements-verified.txt](sync_metadata/network_h0_particle_20260925/requirements-verified.txt)。历史绝对路径迁移规则见[既有迁移说明](sync_metadata/network_h0_particle_20260925/PORTABILITY_ZH.md)。

## 输入版本

NEW H0 SHA：`064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4`。当前新入口显式验证该流场；Particle默认frozen_reference仍是历史OLD equal-pressure 2 mm/s场，不能混用。更早P8.2A/PPT500和RBC理想化Poiseuille演示保留各自来源；RBC演示不冒充NEW H0 RBC生产轨迹。Taylor-Hood仍为已停止验证阶段。

同步排除环境、构建/安装树、完整上游脑数据、渲染帧、大型重复部署包和非必要中间状态。重要100k科学数据无损保留。原报告中的“未push”是开发阶段历史记录；本分支是其后的发布快照。
