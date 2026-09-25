# 当前代码、数据与服务器路径

这是适合克隆后浏览的入口。原机器上335条关键路径和本次只读目录扫描见 [原始路径索引](sync_metadata/p9a4_flow_particle_rbc_20260925/path_inventory/WORKFLOW_PATHS_20260925_ZH.md)；该记录中的绝对路径用于定位原机器，不是GitHub相对链接。

| 本地/服务器来源 | 仓库内对应 |
| --- | --- |
| `/home/lzy/projects/ulm_3D_vascular/` | [vascular_network/](vascular_network/) |
| `/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/` | [formal_3D_flow_solver/FEM_SimVascular/](formal_3D_flow_solver/FEM_SimVascular/) |
| `/home/lzy/projects/ulm_particle_3d_particle0/particle_3d/` | [particle_3d/](particle_3d/)；受当前新增P9-A.4扩展 |
| `/home/lzy/projects/ulm_particle_population_inlet_p9a4/particle_3d/` | [particle_3d/](particle_3d/)；新模块、合同、测试与报告完整保存 |
| `/home/lzy/projects/sonovue_size_distribution_v0/` | [sonovue_size_distribution_v0/](sonovue_size_distribution_v0/) |
| `/workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/` | [server_evidence/flow_mean_2p0_mmps_A_H0_20260924T181140Z/](server_evidence/flow_mean_2p0_mmps_A_H0_20260924T181140Z/) |
| `/workspace/particle9a4_population_inlet_20260924T233708Z/` | 主体与当前 [P9-A.4报告/输出](particle_3d/reports/particle9a4_population_inlet/)相同；其余按server manifest定位 |
| `/root/particle8_2_runs/` | [server_evidence/particle8_2_runs/](server_evidence/particle8_2_runs/)及去重映射 |
| `/workspace/hemocell_restore/` | [server_evidence/hemocell_restore/](server_evidence/hemocell_restore/)及去重映射 |
| `/workspace/microbubble_lammps/` | [server_evidence/microbubble_lammps/](server_evidence/microbubble_lammps/)及去重映射 |

代码在src/scripts，模型合同在contracts，科学数据在data/outputs/run，日志在logs及各阶段reports。目录级导航不替代逐文件manifest；被排除的大文件仍留在原机器，理由及路径见同步元数据。
