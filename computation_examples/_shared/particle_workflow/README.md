# 当前微泡与 RBC 工作区

最新微泡结果为 **ROI-only-balanced-pressure-v1 新流场、1500 条、名义 dt = 0.5 ms**，已完成计算、CUDA 数据复核及动画交付。入口为 [OPEN_RESULTS.html](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/OPEN_RESULTS.html)；结果及保留基线见 [CURRENT_RESULTS.md](CURRENT_RESULTS.md)。

- 本批计算、后处理及复现入口：[particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/](particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/REPRODUCE.md)。
- 本批固定入口队列：`particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/data/cohort.json`。沿用原有粒径分布及单次准入规则，在新 FEM 场上生成，按来源顺序取首 1500 个有效样本。
- 有限尺寸积分及 RBC 源码：`particle_3d/src/particle_3d/`；原有 114 个受保护源码文件保持不变。
- RBC 代码、轨迹及展示：[RBC 结果索引](CURRENT_RESULTS.md#rbc-代码轨迹与展示)，属于各自保留的模型，本轮未计算新场 RBC 耦合轨迹。
- HemoCell RBC 副本：`server_evidence/hemocell_restore/`。
- 共享 Python：`/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python`。

核验最新批次：

```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python particle_3d/reports/microbubble_roi_only_dt0p5ms_n1500/scripts/verify_collected.py
```

原有 Network-H0、dt=1 ms、CORE500 保留在 `particle_3d/reports/particle9a5_formal_trajectories/`，其入口来源、源码身份及几何/材料依赖继续保留。`scripts/check_current.py` 和 `scripts/verify_current_data.py` 仍针对该 H0 基线；本批使用上面的独立核验入口。原有加载器先校验 H0 依赖，本批 `campaign.py:setup` 再加载并严格校验新流场 SHA；`CURRENT_WORKFLOW.json` 记录当前入口及保留基线。

`frozen_reference` 是网格、材料和 RBC 的复现依赖。部分源码保留早期阶段名，仍被积分器导入或属于 114 个源码身份契约，不能仅凭名称判断为废弃代码。

[完整当前工作流](/home/lzy/projects/ACTIVE_VASCULAR_WORKFLOW.md) · [服务器路径](/home/lzy/projects/CURRENT_SERVER_PATHS.md)
