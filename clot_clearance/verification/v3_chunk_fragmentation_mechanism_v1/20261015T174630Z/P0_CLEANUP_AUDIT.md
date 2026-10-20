# P0 小范围清理与保留审计

实际运行链：launch.py → run_protocol.py → Native/native.py → 固定原始 reference、bpm/peri/variational 与 LPS/VARIATIONAL/DYNAMIC → external.py 的真实图/载荷/拖曳。Integrator 执行 trial/check/reject/subdivide；mechanics.py 独立复核完整力、能量、接触与 VV；journal.py 保存接受标量、稀疏损伤/接触和节点帧。analyze/figures/report 是归档后处理，prepare/freeze/repair_registry_build 是构建与冻结依据，archive_rejections/local_readback 是无损数据与交付核验入口。

逐文件 Python AST/导入/函数及 SHA 见 P0_CALL_CHAIN.json。四个隔离 C++ 文件与旧 endpoint-star Python 参考模型均由正式库/冻结清单引用；本轮不删除或重写。多个旧轮次中重复的 Contact、Native、图统计服务各自冻结复现；在没有证明其无历史价值时，不把“没有当前调用”误当作废弃。

本轮组织调整：新增独立目录；用 common.py 统一只读路径、SI 常量、SHA、原始状态读取，并直接导入已验证的完整 endpoint-star 能量实现。后续拓扑、候选、接触和小试样共用单个 graph_ops.py，避免重复实现材料/锚点可达性。未发现可证明安全删除且无历史价值的旧科研源码，因此旧源码删除数为 0；不为清理数量破坏正式验证或原生构建依赖。未作无关重构。
