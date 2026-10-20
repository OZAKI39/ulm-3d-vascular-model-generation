# 声流／剪切多碎块研究复现

本轮目录 `verification/clot_acoustic_streaming_shear_multifragment_v1/20261017T132706Z`。不要在已接受目录重跑或覆盖。新计算先复制本轮分析目录，在 common.py 中设置新的 RUN，并重新生成冻结文件。模型和数据的标签是科学证据的一部分。

服务器激活 `/workspace/y_clot/env.sh`，`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`。沿用 `verification/clot_multifragment_breakup_v1/20261017T110849Z/build/liblammps.so`，SHA `465f8144c6da5c9de3aec4137b9c325f53c26737b199c536fb4a20845511ef77`。通过 LD_LIBRARY_PATH 指向该 build；该库及原始 material law 未修改，也没有安装新依赖。MPI必须1 rank。旧建库配方与上游锁定提交见该前序目录 BUILD_AUDIT.json、analysis/clot_multifragment_breakup_v1/build.py；本轮不需要再次编译。

执行顺序（只在新目录）：

1. `python audit_inputs.py`：读取并哈希实际 F1/F10S/APECSS 文件。优化器在独立环境 `/workspace/y_clot/envs/inverse_opt_v1`，不要求 native 环境安装。
2. `python acoustic_streaming_shear.py`：生成解析模态 Stokes 场，独立有限差分检查。
3. `python shear_loading.py`：原始 top_ids 生成实际节点三角形，冻结非负载荷映射及两例协议。`python v3_mapping_audit.py`：真实V3面片只读映射比较。`python audit_discretization.py`保存节点分摊面积、局部方向/幅值偏差及历史均匀剪切可比性。
4. `python multifragment_objective.py`：冻结独立碎块约束。`python check_coupling.py`：3个非零时刻原生run 0耦合检查，不推进动力学。
5. `python manage.py stokes`：独立Supervisor前台作业。A正式结束且没有目标碎裂后，`python manage.py inverse`，总共两组主要原生动力学。不得自动第三例或优化扫描。
6. `python analyze.py`：原始存活连接、损伤网络、寿命/形态、三维桥接、两张独立物理与机械表和完整能量账本。`python current_moment_audit.py`另输出实际关键态力矩，区分history的名义参考坐标力矩。
7. `python validate_outputs.py --case stokes` 与 `--case inverse`：每例一次全接受态损伤图重放、FP64关键态/加载核对，合并SHA校验。不会重新积分。
8. `python figures.py`、`python report.py`：实际数据PNG/SVG及20问中文报告/HTML；人工查看最终图。
9. 同步到WSL，在本地执行 `update_wsl_document.py`，备份/检查写前SHA/保存diff，回传WSL更新证据；`python update_entrypoints.py`保留旧入口并添加本轮当前入口，`python seal.py`封存，最后本地运行 `local_readback.py`。

声流模型界限：u=(sin(ky)v(z,t),0,0)，周期700µm、间隙25µm、μ1mPa·s、ρ1000kg/m³，上边界4mm/s规定速度，有限边缘与流固反馈省略。声压、APECSS半径响应和质心轨迹都不自动转换为流体剪切。B的1000Pa仅是同模式的机械牵引缩放，不代表高速度Stokes流被物理验证。参考表面固定，脱离后仍加载原始节点，不能称为跟随当前边界的声流。

关键输出：每接受步history.csv、原始ID不可逆damage_events.jsonl.gz、逐接受态CONNECTED_COMPONENT_HISTORY.jsonl.gz及原始成员registry；keys中真实FP64 x/v/force/material_force/κ/d/native能量。剪切方向与实际受力节点由 SHEAR_TRACTION_FIELD.npz 和冻结 ShearLoad 完整确定；history中的时刻足以重建每步力。出生/分裂/内部存活率/最终原成员保持独立保存。物理未知量使用UNAVAILABLE/null；无脱离时质量比分母为0，不填虚假0%。

源代码/依赖SHA和补丁随交付保存；不下载无关项目全量HDF5、上游构建树、二进制库、论文全文或重复缓存。实际原生库留服务器，源码构建来源与哈希可追溯。本轮交付新任务文件并保留历史。
