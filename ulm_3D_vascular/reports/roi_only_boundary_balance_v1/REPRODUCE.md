# 复现与证据说明

此目录为独立 ROI-only 数值设计，不替换 full-A 历史工作。

## 执行顺序（完成状态以 progress.json 为准）

1. `optimize_roi_boundary_balance.py --output /home/lzy/projects/temp_storage/github_sync_20260927/ulm_3D_vascular/reports/roi_only_boundary_balance_v1`：先 mixed-BC H0 回归，后 basis/直接/迭代检查。
2. 原有三组加新 ROI 测试：155 passed；日志 `pytest.log`。
3. 同一纯 0D 入口增加 `--freeze`，写入不可覆盖的设计。只重复轻量 0D，无 CFD。
4. 执行下面 prepare 命令；预检后创建唯一 case。
5. rsync 上传 case 及原 runner、新排他 launcher。private supervisor 首次启动 `roi_balance_gpu`。
6. 只运行这一个科学 CFD；无 CFD-based retuning、无 Particle/RBC。
7. CFD 完成后回传全部日志、原生场、checkpoint；执行独立验收及报告。

## 纯 0D 可安全复核

在**新的审计目录**执行，不覆盖已经冻结的本目录：

```bash
cd /home/lzy/projects/temp_storage/github_sync_20260927
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B \
  ulm_3D_vascular/scripts/optimize_roi_boundary_balance.py \
  --output /tmp/roi_only_independent_review
```

## 原始一次性 prepare 命令（已执行，不要重建覆盖）

```bash
/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B ulm_3D_vascular/scripts/prepare_roi_balance_fem.py \
  --frozen-design /home/lzy/projects/temp_storage/github_sync_20260927/ulm_3D_vascular/reports/roi_only_boundary_balance_v1/roi_boundary_design_frozen.yaml \
  --source-case /home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1 \
  --target-case /home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1 \
  --report /home/lzy/projects/temp_storage/github_sync_20260927/ulm_3D_vascular/reports/roi_only_boundary_balance_v1
```

## 完成后的只读数值复核

验收器会写入新 frozen_flow，禁止覆盖；审核者若要重跑，应复制新 case 到独立审核目录、仅在副本内移走已有 frozen_flow，再传副本路径和一个新的 report 输出目录。原 case 保留不变。

```bash
/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B ulm_3D_vascular/scripts/validate_roi_balance_fem.py \
  --case <independent-case-copy> --source-case /home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1 \
  --flow-root /home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular --report <independent-review-output>
```

报告生成仅读完成后的证据（无求解）：

```bash
/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/.venv/bin/python -B ulm_3D_vascular/scripts/report_roi_balance_forward.py \
  --case /home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-ROI-only-balanced-pressure-v1 --source-case /home/lzy/projects/temp_storage/github_sync_20260927/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1 --report /home/lzy/projects/temp_storage/github_sync_20260927/ulm_3D_vascular/reports/roi_only_boundary_balance_v1
```

## 服务器路径和一次性启动证据

- SSH alias `vast4090`。
- 目录 `/workspace/flow_roi_only_balance_20260928`。
- private supervisor 配置 `server/supervisord.conf`；autostart=false，autorestart=false。
- case `mean-2p0-mmps-A-ROI-only-balanced-pressure-v1`。
- 实际唯一启动命令保存在 case `reports/one_scientific_case_dispatch.json`。
- 不再次运行启动器；其 O_EXCL claim 会拒绝已有 dispatch。
- 监视状态可用 `ssh vast4090 supervisorctl -c /workspace/flow_roi_only_balance_20260928/supervisord.conf status`。

## 附件定义

- `evidence/`：关键几何、验证的端口、H0 提取状态、旧不可达证书及实际核心代码副本；来源/哈希见 sources.json。
- `roi_geometry_cache.npz`：仅 ROI 109 节点、108 条边；SI 单位，不含下游分量。
- `roi_boundary_design_frozen.yaml`：纯 0D 冻结设计，校验和保护。
- `roi_fem_handoff.json`：真实延伸段 40 个截面面积、半径、阻力、两次 gauge 来源。
- `configuration_diff.json` / `final_preflight.json`：只改三个 outlet pressure 值的输入审计。
- 原始大场及运行日志放在唯一 case 目录；报告目录无需重复其完整数据。
