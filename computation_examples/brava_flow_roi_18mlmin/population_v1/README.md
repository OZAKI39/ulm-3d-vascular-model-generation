# BraVa 共享物理时间微泡群体模拟 V1

本目录是独立的新任务。既有 CFD、历史单微泡源码及已停止的旧批次均为只读参考。
当前用户配置：C_MB = 4×10⁹ m⁻³；名义 dt = 0.5 ms；每个正式出口累计至少 20 次首次球心穿越后停止灌注，再观察在途微泡。源事件数、接受数、活跃数和完成数分别记录。

正式计算已完成：1565 个源事件全部准入，1.240815660 s 达到每出口至少 20 个并停源，观察至 13.240815660 s。实际出口 O1/O2/O3 为 458/527/437；另有 141 条观察截止、2 条数值失败、0 条支持静止。数量守恒通过，但完整物理轨迹成功标志为 false；不能将失败和删失当作成功输运。

交付入口：[总预览](OPEN_RESULTS.html)、[中文验证报告](POPULATION_VALIDATION_REPORT_ZH.md)、[复现说明](REPRODUCE.md)、[新增文件清单](FILES_CHANGED.md)。正式数据只使用 `runs/realization_0000_contactfix/`，三个中止开发批次保留作为失败证据。

实现位于 `ulm_particle_formal_p9a5/particle_3d/src/particle_3d/population_v1/`；测试位于同工程的 `tests/population_v1/`。本目录保存契约、验证、运行记录、结果与展示。

## 阶段记录

- Phase 1：`reports/ARCHITECTURE_ZH.md`、`reports/PHASE1_SOURCE_REVIEW.json`；审查实际模块、原物理契约及 BraVa 输入身份。
- Phase 2：`test_phase2.py`，`scripts/validate_real_reference.py`；共享时钟、SoA、出口移除。真实入口及近壁状态在相同物理时间分段下与原模块一致。短窗口验证，不代表整条真实轨迹验证完成。
- Phase 3：`test_phase3.py`；精确球间隙、连续扫掠、组件求解及 40 泡稀疏接触链。`validation/phase3_collision.json` 保留所有时间步数据。斜向压缩伸展场的粗时间步具有明显轨迹误差；PASS 指约束/稳定 ID/细化方向验证，不是 dt 精度认证。
- Phase 4：Poisson 来源、精确正通量入口、原粒径分布、一次准入、精确出生切步、删除和守恒。
- Phase 5–7：真实 FEM tetra walking、空间散列、有限三角形壁面、GPU FP64 阻力及接触；困难接触组件明确回退原 CPU 权威解。
- Phase 8：独立 realization 批处理，来源和时钟独立，禁止跨 realization 碰撞。
- Phase 9：26 项最终测试、真实场独立检查、CPU/GPU 分阶段性能、正式群体运行、独立账本/几何复核和英文动画。实际结果在 `validation/server/`，逐阶段说明在 `validation/phase1/` 至 `phase9/`。

## 已执行检查命令

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
export PYTHONPATH="$PWD/source/particle_3d/src"
export POPULATION_VALIDATION_DIR="$PWD/validation"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
env/bin/python -B -m pytest source/particle_3d/tests/population_v1 -q
```

服务器隔离环境为 `/workspace/brava_flow_roi_18mlmin_20260928/population_v1/env`，源码及验证都在同一新目录下；不会启动原来已经停止的批次。真实流场 SHA256 为 `73228d094c253cec75dbb0030476b687171d0d0fef7816a845392dad220285dc`。

## 数值与物理边界

本模型是不可变形球体、冻结单向 CFD、过阻尼阻力与真实有限三角形接触模型。球间 V1 仅含领先阶法向近场阻力及同时接触；没有声学、壳破裂、聚并、RBC 或双向流固耦合。继承的 2 nm 交接尺度是连续模型边界，不是分子接触或真实壳厚。模型支持静止不代表生理栓塞。只有出口实际首次球心穿越才计入用户要求的每出口 20 个。
