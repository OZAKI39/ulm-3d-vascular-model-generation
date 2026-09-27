# 当前微泡与 RBC 工作区

当前微泡入口：`particle_3d/reports/particle9a5_formal_trajectories/`，Network-H0 新边界条件流场，名义 **dt = 1.0 ms**。当前计算、数据和结果索引见 [CURRENT_RESULTS.md](CURRENT_RESULTS.md)。

- 计算代码：`particle_3d/src/particle_3d/`。
- 入口出生样本与生成代码：`particle_3d/reports/particle9a4_population_inlet/`；正式队列继续使用已冻结的 accepted births，不重新抽样。
- 运行、收集、分析与渲染：`particle_3d/reports/particle9a5_formal_trajectories/scripts/`。
- RBC：`particle_3d/src/particle_3d/rbc*.py`、`particle_3d/scripts/*rbc*.py` 和 `particle_3d/reports/rbc_mb_*/`；几何与轨迹验证数据在 `reports/particle2/`、`reports/particle3/`。
- HemoCell RBC 的本地代码和结果副本：`server_evidence/hemocell_restore/`。
- 共享 Python：`/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python`。该旧名称目录现仅承载共享环境，没有旧粒子工程代码。

运行当前检查：

```bash
/home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python scripts/check_current.py --output /tmp/ulm-current-checks
```

加 `--full-simulation` 可运行两个需要重新计算完整轨迹的并行/积分器一致性测试。默认检查覆盖当前输入加载、短轨迹积分、队列、入口采样及 RBC 展示数据。已有正式轨迹另由 `scripts/verify_current_data.py` 校验完成标记与每个结果文件的 SHA256。

2026-09-27 已按用户要求删除历史工作树、开发报告、旧轨迹和服务器配置副本。保留的 `network_derived_flow_mb_validation_v1` 是当前 NEW 流场适配器；旧 `frozen_reference` 是网格/黏度/RBC 复现依赖，**不能用作当前微泡的流场入口**。当前微泡加载器仍强制校验 NEW SHA。

源码包内一些模块保留早期阶段名称，但仍被当前积分器/RBC 导入，或属于当前采样器强制校验的 114 个源码文件；删除或改名会破坏现有队列身份校验。历史 `CURRENT_PARTICLE_INPUT.json` 属于旧入口契约，当前入口以本文件与 `CURRENT_WORKFLOW.json` 为准。

[完整当前工作流](/home/lzy/projects/ACTIVE_VASCULAR_WORKFLOW.md) · [服务器路径](/home/lzy/projects/CURRENT_SERVER_PATHS.md)
