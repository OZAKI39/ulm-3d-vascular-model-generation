# Step 3：真实血管纯流体边界短程验证

此目录是唯一新增源码目录。HemoCell / Palabos core、官方案例及已有 tracked files 均只读；所有构建、运行与报告写入独立 Step3 run directory。没有 RBC、粒子、周期驱动或体力驱动。

物理参数来自当前旧工程配置、常量、已接受的 runtime contract 和推广输出；`prepare_contracts.py` 只读取这些文件，不运行旧求解器。生理入口/出口身份仍为 ASSUMED。

`prepare_numerics.py` 保存逐端口 API 证据、重新推导的 lattice units、预检查、四个延长段测量面的真实 STL 截面及积分权重。它读取 Step2 冻结的 cap triangle IDs，使用原来的 STL reader；不重新识别端口。

`vascularPureFluid.cpp` 复用 Step2 闭合体素化的原生构造步骤以建立 TriangleHash，逐项检查原生体素、变换、cap 三角形坐标与标签。原生 Guo 方法需要的外侧 ghost 是数值支持，不是新的物理管腔。物理质量只统计原 Step2 closed lumen。入口用原生 plug profile，将真实三角面上的积分归一化到指定 Q；它是 **STEP3_NUMERICAL_INLET_PROFILE_ASSUMPTION**，不是实验测得的分布。三个出口分别由原生 DensityNeumann profile 给定固定密度（恒定表压、法向动量外推）。没有自定义 population reconstruction。

CPU/MPI 只使用系统工具。示例命令中的 RUN 必须是已准备且未使用的新输出目录；不要覆盖已有证据。先审查 `contracts/physical_bc_contract.json`、`bc_api_feasibility.md`、`lattice_unit_contract.json`、`execution_contract.json` 和 `diagnostics/lbm_preflight.json`，所有 gate 通过后才运行：

```bash
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
/usr/bin/cmake -S /home/lzy/projects/hemocell_starter/test_code/step3_vascular_pure_fluid -B "$RUN/build" -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DCMAKE_BUILD_TYPE=Release
/usr/bin/cmake --build "$RUN/build" -j 1
cd "$RUN"
OMP_NUM_THREADS=1 /usr/bin/mpirun -n 1 "$RUN/build/vascular_pure_fluid" "$RUN"
/usr/bin/python3 -B /home/lzy/projects/hemocell_starter/test_code/step3_vascular_pure_fluid/export_fields.py "$RUN"
```

0、1、10、100、1000 步连续推进，每一步检查流体和 ghost 的有限值、正密度、Mach；0/10/100/1000 保存场。1000 步不是稳态证明。初始化 rhoLU=1、u=0，前 10 步是明确的数值启动 ramp，随后恢复完整物理目标。

流量在 cap 内侧 4dx 的真实截面上积分实时格点速度。MPI 的 `MultiBlockLattice::get` 返回复用代理，不可缓存其地址；本实现缓存所属 atomic bulk 的格点，并在第 10 步逐格点与即时读取的 MPI API 结果交叉验证。首次实现的该诊断缺陷已单独保留于本轮输出 `attempts/attempt1_invalid_flux`；不能使用其中的流量列验收。

`export_fields.py` 导出原体素的稀疏 VTU 子集，不做重采样。先 Threshold `FluidMask=1` 查看物理场。Step2 的 497 个端口标签位于边界支持格点，保持原值并由 `PhysicalFieldValid=0` 标明；这些标记格点上的零值只是占位。单位详见输出 `field_units.json`。

最终报告保留数值验收与人工 ParaView 验收的区别。本轮不增加稳态长算、网格收敛、RBC 或 Git 提交。
