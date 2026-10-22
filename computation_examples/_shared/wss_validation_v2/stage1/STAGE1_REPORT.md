# 阶段一：工程修复与数值回归

本阶段已运行。所有原冻结输入和图像保持原样；输出写入 stage1/recomputed/。

## 修复内容

1. 将原 `compute_field_diagnostics.py` 的七个核心函数提取为 `solver_support/src/flow_solver_support/wss.py`。原文件 SHA256 为 `efa12e218bcaf99ba07db7194d3b8da79251da709c0a54a97ad6f6179be6e438`，与旧 COMPUTE_VALIDATION.json 对应。六个纯数值函数 AST 完全一致；surface 函数仅增加局部 PyVista 导入，移除该导入后 AST 也一致。模块导入不读取算例、不导入日志解析器或绘图库。
2. `wss_case.py` 从显式 case/run/solver.xml 的 Constant Viscosity 读取 μ，并与同一 case/policy.json 的 ν×ρ 或 μ 字段比较；正值和有限性检查，错误材料会被拒绝。当前 μ=0.00345312 Pa.s，ρ=1056 kg/m3，ν=3.27e-6 m2/s。
3. 生产入口 `rotate_visualization/prepare_surface_data.py` 支持 --case/--config/--flow/--arrays/--mesh/--output。不再定位不存在的旧脚本；默认兼容现有 BUILD_CONTEXT，但本次实际运行显式指定 H0 和独立输出。保存 XML、policy、体网格数组、速度 VTU、NPZ、核心模块及入口哈希。
4. 同时核对 NPZ/VTU 字段、逐单元节点集合、边界及节点顺序，保留无滑移、有限值和总流量检查。没有改变 WSS 公式、显示映射、BC 或色标。
5. 阶段二首次实际求解暴露了输出坐标的精度问题：svMultiPhysics 的 VtkData.cpp 240–260 使用默认 vtkPoints，输出为 float32。新增的严格 double 相等检查因此拒绝了坐标；实际数值恰好等于原输入转 float32。最小修复增加 coordinate_identity，只接受精确原坐标或精确 float32 序列化，并继续检查连接关系；梯度始终使用原输入网格坐标，速度/压力不插值。任意偏移的负测试会失败。H0 在修复后再次回归仍为零差。

## 已运行证据

| 检查 | 结果 |
|---|---:|
| 原始面片 WSS 最大绝对差 | 0 Pa |
| 节点显示 WSS 最大绝对差 | 0 Pa |
| 牵引向量、法向、面积、节点/面片/父单元 ID 差 | 全部为 0 |
| 仿射梯度最大绝对误差 | 1.13687e-13 s^-1 |
| 仿射牵引最大绝对误差 | 9.99201e-16 Pa |
| 刚体旋转最大 WSS | 2.19456e-16 Pa |
| Couette 理论/计算 | 4.3164 / 4.3164 Pa |
| 法向反转后的模长差、向量和 | 0 / 0 Pa |
| μ 与 ν×ρ 不一致 | 按预期抛出 ValueError |
| 任意坐标偏移 | 按预期拒绝 |
| H0 入口相对目标误差 | 2.87046e-10 |
| H0 总边界质量相对误差 | 5.08497e-17 |

验收容差为 1e-9 Pa，仅表示同输入回归一致，不能用来说明物理精度。数据详见 regression.json、data/analytical_tensor_checks.json、recomputed/COMPUTE_VALIDATION.json；命令和原始输出见 ../COMMANDS.md 与 ../logs/stage1_*.log。

源码最小差异使用本轮开始前文件作基线保存于 ../evidence/production_fix.patch，避免把既有清理改动混入本轮差异。原包装器及历史核心完整副本位于 ../evidence/before/。三份生产修复文件已加入本地 Git 暂存；未替用户提交或覆盖其他未提交修改。

追加并已复测的来源绑定：生产入口核对冻结 manifest 中的 input_configuration_sha256；新算例还核对 input_hashes.json 内的 XML、policy、mesh_arrays。后续冻结输出同时保存配置/网格/参数哈希。新增负测试将 μ 与 ν 同时乘2，使材料文件内部一致，仍被真实生产入口以“Configuration differs from frozen solve”拒绝；没有生成 WSS 输出。H0 零差回归与所有张量测试再次通过。详见 data/frozen_configuration_negative_test.json，脚本 frozen_configuration_negative_test.py。
