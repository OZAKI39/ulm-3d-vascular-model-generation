# 最佳可实现分流：交付索引

主报告为 `BEST_FEASIBLE_BALANCE_FORWARD_REPORT_ZH.md`。本目录仅对应本轮有限范围设计与唯一 CFD，不替换旧 equal-split 失败记录，也不替换轨迹生产输入。

- `final_preflight.json`：启动前冻结的配置、设计、哈希、测试及单次 CFD 约束；不能因运行完成而改写。
- `scientific_case_registry.json`、`delivery_status.json`：实际执行/交付状态；原生 dispatch 和 completion 在新 case 的 `reports/`。
- `best_feasible_fem_handoff.json`、`handoff_H0_regression.json`、`configuration_diff.json`：延伸段传递、旧 H0 解析回归与仅出口压力改变的审计。
- `independent_validation.json`：原生日志、保存快照、最终场和 checkpoint 的独立复核。
- `all_balance_metrics.csv`：分流为无量纲 Qout/Qin；J、极差、std、最大等分偏差均使用比例而非百分数。
- `final_3D_ports.csv`：流量 m³/s，压力 Pa，3D−0D 差为比例。
- `wss_comparison.csv`：原始壁面面片 WSS 的面积统计，列名标明 Pa 和 μm²；包含人工延伸段。
- `data/`：旧/新 WSS 的真实壁面 VTP、原始面片 SI 数组及原网格面片/四面体编号。
- `figures/`：原始 WSS 同尺度比较、设计网格与实际分流比较，PNG/PDF；不新增平滑。
- `evidence/wss_core/`：字节不变的现有 WSS 核心；原路径和 SHA256 在 provenance 中，相关旧验证证据随附。
- `logs/`：本次命令输出，包含正式冻结前发现并修复的序列化测试失败，保留过程证据。

纯 0D 设计及完整 trace 在相邻 `../parameterized_0d_v1/best_feasible_balance/`。脚本位于 `../../scripts/`，科学方法说明位于 `../../docs/BEST_FEASIBLE_BALANCE_DESIGN.md`。

唯一 CFD 位于仓库的 `ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-best-feasible-balance-v1/`。完整中间时步留在本地和服务器 `run/1-procs/`；Git 不重复保存这些大文件，保留独立验收后的 frozen_flow、checkpoint、原生日志与各快照 SHA256。已有 dispatch 文件阻止再次启动该科学算例。
