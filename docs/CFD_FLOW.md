# Seeder + Musubi CFD FLOW

本阶段只消费已经通过检查并由人工接受的最终 CFD 表面，不重新运行 SWC、ROI、Ultraliser、全局一维流动、VMTK 延长、重网格、封口或压力修正。正式方法名为 `PROTEUS_COMPATIBLE_SEEDER_MUSUBI_STEADY_LBM_BASELINE`，表示其数据结构与后续 PROTEUS 导入兼容，而不宣称复现 PROTEUS 的全部论文设置。

流程固定为：最终带实体标签的 CFD 表面 → 五个边界分片 → Seeder 均匀笛卡尔/树网格 → Musubi 稳态 LBM → 完整速度与压力 VTU → PROTEUS 导入元数据。

表面中的 `CellEntityIds` 与边界清单共同确定 `wall`、`inlet`、`outlet_01`、`outlet_02` 和 `outlet_03`。五个 STL 仅做 `µm × 10⁻⁶ = m` 的单位换算，不移动、不平滑、也不重建几何。Seeder 的种子点由入口面积中心沿内法向作有限的确定性搜索，并由闭合管腔的内外判定确认。

当前基线只采用 `dx = 0.20 µm`，不进行自动分辨率扫描。Musubi 使用单相 Newtonian 血液、D3Q19、BGK、刚性无滑移壁面和稳态条件。入口由官方 `mfr_eq` 精确保持质量/体积流量；上游要求的抛物线分布会被记录，但本固定版本的有效入口分布如实报告为官方质量流量边界。三个出口使用官方 `pressure_eq`，其压力来自表面延长修正后的 `P_solver_boundary_pa`。共同物理压力参考保证负表压不会变成负的 LBM 绝对压力，同时严格保持三个出口间的压差。

运行前会估算流体格点数和内存，超过当前可用内存 60% 时直接停止，不会自动改粗网格。Seeder 和 Musubi 各最多执行一次；失败后不修改 `dx`、`dt`、松弛参数、边界类型或求解器重跑。只有稳态收敛、质量守恒误差不超过 1%、有限值、实际格子 Mach 数小于 0.05、VTU 字段及 PROTEUS 元数据全部通过时，状态才是 `CFD_FLOW_MUSUBI_BASELINE_PASS_PENDING_GRID_CONVERGENCE`。

## 运行

在项目根目录使用 pmp Python：

```powershell
D:\anaconda3\envs\pmp\python.exe cfd_flow.py
```

也可以显式给出唯一配置文件：

```powershell
D:\anaconda3\envs\pmp\python.exe cfd_flow.py configs/cfd_flow.yaml
```

正式 APES 可执行程序运行在配置指定的 WSL2 Ubuntu 中。需要固定提交的官方 Seeder、Musubi、`mus_harvesting`、MPI 启动器和 Fortran MPI 编译环境。缺失时状态为 `CFD_FLOW_ENVIRONMENT_BLOCKED`，不会替换为其他 CFD 求解器或生成伪造流场。

## 输出

每次运行写入 `outputs/cfd_flow/musubi_anchor003274_<timestamp>/`，仅包含 `input/`、`geometry/`、`seeder/`、`musubi/`、`flow/`、`qc/`、`figures/` 和 `proteus/`。成功结果的主要文件是：

- `flow/flow_field.vtu`：米制坐标、Cartesian hexahedral cells，包含三分量 cellData `velocity_phy`、`pressure_phy` 和 `pressure_gauge_pa`。
- `proteus/proteus_flow_metadata.json`：PROTEUS 导入字段、单位和入口等效直径。
- `figures/`：速度、表压和入口流线三张复核图。
- `qc/`：输入表面、边界分片、Seeder 网格、LBM 缩放、边界条件、收敛、流场、质量守恒、出口压力和 PROTEUS 兼容性证据。

该结果是第一组稳态 Newtonian 基线，不是已经完成网格收敛研究的最终解。成功后的下一步是另行开展 Musubi 网格间距收敛研究。

## D3Q19 显式体积黏度

固定版本 Musubi 的 `mus_load_fluid()` 对 `kind = 'fluid'`、`layout = 'd3q19'`
不会在缺少 `bulk_viscosity` 时采用默认值，而是终止配置加载。因此本基线按照同一
官方源码所附 D3Q19 示例，显式采用 `bulk_viscosity = (2/3) * nu_phy`。当前运动黏度
为 `3.27e-6 m2/s`，由此计算得到体积黏度 `2.18e-6 m2/s`。该值是求解器所需的官方
基线参数，并非小鼠实验测得的生理体积黏度。配置和运行记录分别使用
`MUSUBI_D3Q19_REQUIRED_EXPLICIT_PARAMETER` 与
`OFFICIAL_MUSUBI_BASELINE_TWO_THIRDS_KINEMATIC_VISCOSITY` 说明其来源和选择策略。

Musubi-only recovery 直接引用已经通过 SHA、文件清单、层级和五类边界检查的冻结
Seeder mesh，不调用 Seeder，也不复制体网格。只有 Musubi 稳态收敛后才允许执行
一次 harvester；导出失败时继续使用冻结的 Musubi solution，不重新运行求解器。
