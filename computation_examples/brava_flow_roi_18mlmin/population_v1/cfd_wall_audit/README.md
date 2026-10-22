# 只读 CFD 壁面审计

入口：WALL_CFD_AUDIT_REPORT_ZH.md / OPEN_RESULTS.html。所有 CSV 使用 SI，列名包含单位；layer 0/1/2=W0/W1/W2，3=其它体单元；tetra_id 和 wall_triangle_id 均从0开始。全部原始坐标/节点值在 wall_adjacent_cells.csv 与 recomputed_arrays.npz 中保留。

服务器复现（不调用 CFD 或微泡积分器）：

```bash
cd /workspace/brava_flow_roi_18mlmin_20260928/population_v1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PYTHONPATH="$PWD/source/particle_3d/src"
env/bin/python cfd_wall_audit/scripts/audit_field.py --environment data/environment.json --output cfd_wall_audit --commit 51257fb4f
env/bin/python cfd_wall_audit/scripts/validate_recomputation.py --population-root "$PWD"
env/bin/python cfd_wall_audit/scripts/supplement.py --population-root "$PWD"
env/bin/python cfd_wall_audit/scripts/particle_overlap.py --population-root "$PWD"
env/bin/python cfd_wall_audit/scripts/plot_audit.py --population-root "$PWD"
env/bin/python cfd_wall_audit/scripts/build_report.py --population-root "$PWD"
```

再次运行只覆盖本目录派生结果，建议复制到新的审计目录后执行。`inputs/visible_speedup_candidates.csv` 是上一轮加速审查的只读参考副本；不从运动图像臆造新事件。环境依赖 NumPy/SciPy/Pandas/PyVista/VTK/Matplotlib，GPU 距离核对用已安装 PyTorch CUDA 和项目 GPUWall。

包括体积及面积等比重的分母见 AUDIT_CONFIG.json；粒子逐段驻留近似归给起始单元。不要把状态次数当驻留时间或把 WSS 代码一致性当物理精度验证。

簇的壁面覆盖面积仅累计簇内W0直接拥有的正式WALL面；只有W1/W2单元的簇可能对应0直接壁面面积，不等于附近没有壁面。所有单元的最近壁面三角形ID仍保存在recomputed_arrays.npz。
