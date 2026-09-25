# Particle：最新 P9-A.4 与历史微泡/RBC

当前入口为 [P9-A.4](reports/particle9a4_population_inlet/P9A4_POPULATION_INLET_REVIEW_ZH.md)：固定NEW H0，Poisson source、条件SonoVue D≤4 µm、全入口通量采样、单次真实WALL+handoff筛选，所有拒绝保留。P9-A.1动力学与原B/C源码不变；首30泡28 completed、2 supported stationary、0 solver fail。

源码 `src/particle_3d/continuous_infusion.py` 与 `population_inlet_p9a4.py`；合同 `contracts/P9A4_CONTINUOUS_INFUSION_V1.json`；测试 `tests/particle9a4_population_inlet/`；所有报告/日志/图件/原始30泡结果在 `reports/particle9a4_population_inlet/`。

100k ledger、accepted births及Method B大文本以gzip保存。先从仓库根执行 `python3 sync_metadata/p9a4_flow_particle_rbc_20260925/restore_large_artifacts.py`，恢复逐字节相同的原路径，再使用原审计/重放入口。[发布复现说明](../sync_metadata/p9a4_flow_particle_rbc_20260925/SYNC_REPORT_ZH.md)。

RBC与MB共用该包。`rbc.py`、`rbc_distribution.py`、`rbc_orientation.py`、`rbc_integrator.py`、`rbc_capillary_surrogate.py`与`coflow_rotation.py`均保留；RBC历史结果位于`reports/particle2/`、`particle3/`、`rbc_mb_flow_rotation/`。最后一个是理想化Poiseuille演示。

`reports/network_derived_flow_mb_validation_v1/`保留前一阶段OLD/NEW配对证据。旧P9-A.1/A.2/A.3/A.3B及P8.2A/PPT仍在原相对路径；旧合同与结果不改。默认frozen_reference是历史OLD场，P9-A.4须经专用NEW loader。没有自动正式500。
