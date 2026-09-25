# 紧凑 ABS/PDMS 制造候选

默认模型：BALANCED。旧大模型：`../manufacturing_roi/`，角色 FULL_CONTEXT_REFERENCE，原文件不变。

状态：PRINT_READY_CANDIDATE。尺寸目标、原始半径及拟合后直径的实际限制见 `compact_report.md`。

NOT manufacturer guaranteed limits. Must be calibrated experimentally.

- `candidates/{MINI,BALANCED,RICH}/`：原始半径/补偿 SWC、映射、删枝原因、原生 VTK、封口打印 STL。
- `BG001_RMCA_BALANCED_print_candidate.3mf`：实际切片成功时生成，可直接打开 Bambu Studio。
- `orientation/`：仅 BALANCED 的 Top 3 方向，Top1 执行一次切片。
- `QC/`：三候选同视角同比例、旧模型对比、方向对比和测试记录。

从项目根目录运行，使用新的输出目录：

```bash
/home/lzy/projects/ulm_particle_3d_particle0/.venv/bin/python tools/compact_bg001_rmca.py --output outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_run02
```

配置：`config/compact_manufacturing_profile.yaml`。支持 `--stage prepare/model/finish/verify`；默认 all。三个候选共用同一语义组件，默认不运行大模型入口。没有自动全局缩放或打印提交。
