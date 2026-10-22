# 复现与附件说明

审计日期：2026-09-27。所有命令从 `wss_audit/` 的上一级 `FEM_SimVascular/` 执行；不会改写正式结果。所有结果的实际判定见主报告，不能将脚本退出 0 等同于血管 WSS 已收敛。

## 本次实际执行

```bash
cd /home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B wss_audit/scripts/audit.py > wss_audit/logs/audit.log 2>&1
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B wss_audit/scripts/analytical_validation.py > wss_audit/logs/analytical_validation.log 2>&1
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B wss_audit/scripts/extra_checks.py > wss_audit/logs/extra_checks.log 2>&1
PYTHONDONTWRITEBYTECODE=1 /home/lzy/projects/temp_storage/ulm_particle_3d_particle0/.venv/bin/python -B wss_audit/scripts/make_figures.py > wss_audit/logs/make_figures.log 2>&1
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -B wss_audit/scripts/verify_bundle.py > wss_audit/logs/verify_bundle.log 2>&1
```

数值环境：Python 3.13.11、NumPy 2.5.3、SciPy 1.18.1、PyVista 0.49.0、VTK 9.7.0；正式绘图还需要 Pillow、Matplotlib、imageio-ffmpeg 0.6.0。首个 Python 环境缺少 imageio-ffmpeg，因此绘图实际用了现有 Particle 环境，没有安装或修改任何依赖。

核心 WSS 文件恢复命令（本次已执行；原文件未恢复到正式目录）：

```bash
git -C /home/lzy/projects/ulm_flow_mean_2p0_mmps show \
  fab4cf0eb3c8f8ad1ce40dfe995181c147cb9a9d:formal_3D_flow_solver/FEM_SimVascular/scripts/flow_2mmps/compute_field_diagnostics.py \
  > wss_audit/evidence/source/compute_field_diagnostics.py
```

### 另一台机器上的独立复核

复制整个 `wss_audit/` 即可运行 `python -B scripts/verify_bundle.py` 和 `python -B scripts/analytical_validation.py`。前者直接用原项目函数重算随附 H0 NPZ，并逐数值比较随附原 WSS VTP；后者调用相同函数执行解析场验证。需安装上述 NumPy/SciPy/PyVista/Matplotlib。它们不依赖项目其他目录、不连接服务器。`audit.py` 的完整时序分析及 `make_figures.py` 的原样式复现则需要原项目和现存时序 VTU。

## 尚未执行的 3D 求解验证：明确与已运行结果区分

服务器 `vast4090` 上已只读核实当前算例、求解器路径和 SHA256，见 `logs/server_paths_check.log`、`logs/server_hashes.log`。本轮没有启动新的 3D 求解器任务。一次原血管求解历史记录约 1675 s，三档网格重算还需重新生成匹配的网格和边界；不能用相同网格的两个 BC 算例冒充网格收敛。

### 固定网格、只改一个出口的后续方法

将本目录的 `prepare_bc_followups.py` 和原 runner 的审计副本上传至新的服务器审计目录后，按下面顺序执行。以下是**后续命令，本次未执行**：

```bash
python -B prepare_bc_followups.py \
  --source-case /workspace/flow_mean_2p0_mmps_A_H0_20260924T181140Z/mean-2p0-mmps-A-H0-pressure-v1 \
  --output /workspace/wss_audit_followup/bc_cases
python -B solve_a_h0_fem_remote.py \
  --case /workspace/wss_audit_followup/bc_cases/O2_plus1pct \
  --reference-root /workspace/flow_mean_2p0_mmps_20260922
python -B solve_a_h0_fem_remote.py \
  --case /workspace/wss_audit_followup/bc_cases/O2_minus1pct \
  --reference-root /workspace/flow_mean_2p0_mmps_20260922
```

首先核对 reference-root 下原 runner 所需的 `src/sv_validation` 和 `scripts/sv13q/flow_parser.py`；服务器旧辅助目录的完整性未在本轮逐项验证，若缺失，需把本地 F 目录的对应模块复制到新审计目录，并令 `--reference-root` 指向该副本。求解器和 build 元数据路径由 runner 使用现存服务器路径；勿改正式 case。扰动为 ±29.3201527108 Pa（当前 O2−O3 压差的 1%），保持 O1/O3、Q、μ、ρ、网格、dt、PETSc 配置一致；比较区域面积平均、P5/P95、低值区面积及流量，而非用图像颜色定优劣。

### 血管三档网格

本轮仅找到一档有效求解网格。后续在独立网格目录使用相同几何及端口标记，设置 `global_edge_size` 为约 0.392、0.294、0.196 μm；如采用局部近壁加密，应单独记录实际壁面法向间距。现有 `generate_tetra_mesh.py` 的 0.8 fallback 仅允许主网格失败后使用，**不能直接拿它充当收敛流程**。复制该脚本到审计目录，参数化输入、输出及尺寸后再生成三档，并复制 `audit_mesh.py` 生成匹配的四面体、边界、节点 ID；不允许沿用旧网格上的速度值作为新解。各档保持同样几何基准、BC、μ、ρ和时间收敛标准，时间步保持共同稳定的小值或单独报告时间步调整并补时间步敏感性。达到稳态后，用同一 WSS 代码及物理区域定义重算；采用几何区域积分、面积分位数及低/高区域位置比较，禁止逐节点硬配不同网格。

### 求解器圆管验证

本轮的圆管 NPZ 是**解析速度在节点上的取样**，不是求解器结果。进一步验证时：R=4 μm、L=40 μm、μ=0.00345312 Pa·s、ρ=1056 kg/m³、平均入口速度 0.002 m/s、Q=1.00530964915e−13 m³/s、刚性无滑移壁、充分发展抛物线入口；选远离两端的中部壁面。边界应与解析充分发展流一致，并分别报告求解器速度/压力误差和 WSS 恢复误差。不能把当前 6.90624 Pa 理论值直接用于原分叉血管。

## 数据字典与范围

- `region_summary.csv`：两个 BC 算例的区域统计；分位数按面片面积加权。`raw_below_2Pa` / `raw_above_30Pa` 的空间掩码固定由当前 H0 决定，旧算例行在同一掩码上比较。其他区域可重叠，不可累加面积。
- `wall_facets.csv`：45,221 个实际壁面三角形，SI 计算后用列名显式标注 μm、μm²、Pa、s⁻¹；ID 为 NPZ 中从 0 开始的原始编号。display_centroid 是显示节点标量的三角形重心插值，不是原始面片 WSS。
- `largest_adjacent_jumps.csv`：共享一条网格边的邻居，按 WSS 差值排序前 30 对；并非主观挑图。
- `outlet2_path_profile.csv`：从 J1 沿 SWC 支路至 O2，1 μm 分箱、面积加权。距路径小于 2.5 倍插值半径；分叉首箱可能覆盖融合区，不能当作精确横截面周向平均。
- `mesh` 质量采用项目已有 minSICN；近壁高度是边界三角形到所属四面体对顶点的垂直距离 3V/A，非棱柱第一层厚度。
- 低值 <1、<2 Pa 和高值 >30 Pa 是审计操作性阈值，不是小鼠生理正常/异常界限。
- `time_convergence.csv`、`residual_solves.csv`、`boundary_values.csv`、`poiseuille_validation.csv`、`extension_scale_check.csv` 各自标明算例、单位及方法限制。
- `evidence/source/`：原核心函数、执行脚本和求解器源码快照；未修改原件。`input_hashes.json` 与 `preservation_check.json` 支持被审计输入未改变的检查。

首次编写审计脚本时修正过两个仅在审计目录中的错误：H0 VTU 文件名应带 `_A_H0`；解析圆管边界索引改为 `x[b][:,:,:2]`。正式源码未改动，最终日志对应修正后的实际成功运行。这些开发错误没有作为生产 WSS 缺陷计入。
