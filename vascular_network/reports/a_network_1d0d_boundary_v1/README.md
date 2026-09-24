# A-Network Steady 1D/0D Boundary Model v1

本轮结论与结果：[中文主报告](A_NETWORK_1D0D_REVIEW_ZH.md)。状态 A_NETWORK_SOURCE_AMBIGUOUS；已完成几何来源/拓扑/精确端口/阻力审计与通用库解析测试，未生成 A 水力解或 3D BC。

在项目根目录使用已具备 numpy/scipy/pyvista/matplotlib/PyYAML/pytest 的 Python。此次复用只读的 `/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python`，没有修改 Particle 代码或环境。

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
python scripts/run_a_network_1d0d_boundary_v1.py
python -m pytest -q tests/network_1d0d --junitxml=reports/a_network_1d0d_boundary_v1/logs/pytest.xml
python scripts/render_a_network_1d0d_boundary_v1.py
python scripts/report_a_network_1d0d_boundary_v1.py
```

审计入口针对当前被冻结来源，默认输出仅写本报告目录；可通过 --output 改为新目录。它不运行旧 s1/s2，不写已有生产 case，也不根据旧 assumed-inlet 资料自动越过 source gate。

所有 null 表示未计算，空 hydraulic CSV 只定义 schema。Figure 01 有 PNG/PDF；Figure 02–06 因前置条件失败不生成。data/uncomputed_outputs.json 逐项说明。日后源点确认后仍需添加经过审阅的 case boundary 声明与完整求解驱动，不能仅改 summary 状态。
